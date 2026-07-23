# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import shutil
import subprocess
from bisect import bisect_left
from collections import Counter, OrderedDict, abc, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from functools import cache
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Optional, Sequence, cast

import biotite.sequence.align as align
import numpy as np
import pandas as pd
import pyarrow
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from numpy.typing import NDArray
from pyarrow import csv
from rdkit import Chem, DataStructs, RDConfig, rdBase
from rdkit.Chem import ChemicalFeatures, rdShapeAlign, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from tqdm import tqdm

from plinder.core.scores.entries import (
    EntryView,
    LigandView,
    SystemView,
    load_entry_views,
)
from plinder.core.scores.metrics import SCORE_NAMES
from plinder.core.structure.smallmols_similarity import mol2morgan_fp
from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger
from plinder.data import databases
from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig

LOG = setup_logger(__name__)

ECFP4_RADIUS = 2
ECFP4_NBITS = 1024
ECFP4_PARQUET_METADATA = {
    b"plinder.fingerprint": b"ECFP4",
    b"plinder.fingerprint.generator": b"RDKit Morgan",
    b"plinder.fingerprint.radius": b"2",
    b"plinder.fingerprint.nbits": b"1024",
    b"plinder.fingerprint.include_chirality": b"false",
}

_NON_CANONICAL_AA = str.maketrans({"J": "L", "U": "C", "O": "K"})


def _protein_similarity_lookups() -> tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    """Build byte-indexed validity and BLOSUM similarity lookup tables."""
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    alphabet = list(matrix.get_alphabet1())
    scores = matrix.score_matrix()
    valid = np.zeros(256, dtype=bool)
    similar = np.zeros((256, 256), dtype=bool)
    for i, aa_i in enumerate(alphabet):
        code_i = ord(aa_i)
        valid[code_i] = True
        for j, aa_j in enumerate(alphabet):
            pair_score = max(0, scores[i, j])
            normalization = max(scores[i, i], scores[j, j])
            similar[code_i, ord(aa_j)] = pair_score / normalization > 0.2
    return valid, similar


_VALID_PROTEIN_AA, _SIMILAR_PROTEIN_AA = _protein_similarity_lookups()


def write_ecfp4_fingerprint_table(ligands: pd.DataFrame, output_path: Path) -> None:
    """Write the fixed ECFP4/1024 fingerprint table with explicit metadata."""
    table = pa.Table.from_pandas(ligands, preserve_index=False)
    metadata = {**(table.schema.metadata or {}), **ECFP4_PARQUET_METADATA}
    pq.write_table(table.replace_schema_metadata(metadata), output_path)


SORT_ORDER = [
    ("similarity", "descending"),
    ("query_system", "ascending"),
    ("query_ligand_id", "ascending"),
    ("target_system", "ascending"),
    ("target_ligand_id", "ascending"),
]
INFO_COLUMNS = (
    "query_system",
    "query_ligand_id",
    "target_system",
    "target_ligand_id",
    "protein_mapping",
    "protein_mapper",
)
_ChainInstanceMapping = str
_ChainPairType = tuple[_ChainInstanceMapping, _ChainInstanceMapping]
_SimilarityScoreDictType = dict[str, float]
_Ligand3DCandidateType = dict[str, str | float | None]
_PocketDataType = tuple[
    dict[str, dict[int, int]],
    dict[str, dict[int, Counter[str]]],
    int,
    int,
    int,
]
_ProteinScoreResultType = tuple[
    dict[str, list[_ChainPairType]],
    _SimilarityScoreDictType,
    dict[_ChainPairType, pd.DataFrame],
    str,
]


def load_ligands_from_index(*, annotation: pd.DataFrame) -> pd.DataFrame:
    """Extract post-ingest ligand-similarity inputs from annotation rows."""
    columns = {
        "entry_pdb_id": "pdb_id",
        "system_id": "system_id",
        "ligand_rdkit_canonical_smiles": "ligand_rdkit_canonical_smiles",
        "ligand_unique_ccd_code": "ligand_ccd_code",
        "ligand_id": "ligand_id",
        "ligand_asym_id": "ligand_asym_id",
    }
    if annotation.empty:
        return pd.DataFrame(columns=list(columns.values()))
    eligible = annotation["system_type"].eq("holo") & annotation[
        "ligand_is_proper"
    ].fillna(False)
    ligands = annotation.loc[eligible, list(columns)].rename(columns=columns)
    return (
        ligands.dropna(subset=["ligand_id"])
        .drop_duplicates(subset=["ligand_id"])
        .reset_index(drop=True)
        .sort_values("ligand_id")
    )


def load_ligands_from_annotation_table(*, data_dir: Path) -> pd.DataFrame:
    """Load proper holo ligand instances from the collated annotation table."""
    columns = [
        "entry_pdb_id",
        "system_id",
        "system_type",
        "ligand_is_proper",
        "ligand_rdkit_canonical_smiles",
        "ligand_unique_ccd_code",
        "ligand_id",
        "ligand_asym_id",
    ]
    annotation = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
        columns=columns,
    )
    return load_ligands_from_index(annotation=annotation)


def annotate_cofactor_similarity(
    ligands: pd.DataFrame,
    *,
    cofactor_smiles: dict[str, str],
    threshold: float = 90.0,
    radius: int = 2,
    nbits: int = 1024,
) -> pd.DataFrame:
    """Annotate unique ligand SMILES against CCD cofactor structures."""
    cofactor_fingerprints: list[DataStructs.ExplicitBitVect] = []
    cofactor_codes: list[str] = []
    for code, smiles in sorted(cofactor_smiles.items()):
        try:
            cofactor_fingerprints.append(
                mol2morgan_fp(smiles, radius=radius, nbits=nbits)
            )
            cofactor_codes.append(code)
        except (TypeError, ValueError):
            LOG.warning(f"skipping invalid cofactor SMILES for {code}")
    if not cofactor_fingerprints:
        raise ValueError("no valid CCD cofactor structures are available")

    result = ligands.copy()
    maximum_similarities: list[float] = []
    closest_cofactors: list[str] = []
    total = len(result)
    started = perf_counter()
    for index, binary_fingerprint in enumerate(result["fingerprint"], start=1):
        fingerprint = DataStructs.CreateFromBinaryText(binary_fingerprint)
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprint, cofactor_fingerprints
        )
        best_index = int(np.argmax(similarities))
        maximum_similarities.append(float(similarities[best_index] * 100.0))
        closest_cofactors.append(cofactor_codes[best_index])
        if index % 5_000 == 0 or index == total:
            elapsed = perf_counter() - started
            rate = index / elapsed
            LOG.info(
                "cofactor similarity progress: "
                f"processed={index}/{total} rate={rate:.1f}/s "
                f"eta_seconds={(total - index) / rate:.1f}"
            )
    result["ligand_max_cofactor_similarity"] = np.asarray(
        maximum_similarities, dtype=np.float32
    )
    result["ligand_most_similar_cofactor"] = closest_cofactors
    result["ligand_is_cofactor_like"] = (
        result["ligand_max_cofactor_similarity"] >= threshold
    )
    return result


def compute_ligand_fingerprints(
    *,
    data_dir: Path,
    cofactor_similarity_threshold: float = 90.0,
) -> None:
    """Fingerprint the unique canonical SMILES in a completed ligand ingest."""
    ligands = load_ligands_from_annotation_table(data_dir=data_dir)
    for column in ligands.columns:
        LOG.info(
            f"compute_ligand_fingerprints: unique {column}="
            f"{ligands[column].nunique()}"
        )

    smiles_column = "ligand_rdkit_canonical_smiles"
    ligands = ligands.dropna(subset=[smiles_column])
    ligands = ligands[ligands[smiles_column].ne("")].copy()
    ligands_unique = (
        ligands[[smiles_column]]
        .drop_duplicates()
        .sort_values(smiles_column)
        .reset_index(drop=True)
    )
    ligands_unique.insert(
        0,
        "ligand_smiles_id",
        np.arange(len(ligands_unique), dtype=np.int32),
    )
    total = len(ligands_unique)
    started = perf_counter()
    binary_fingerprints: list[bytes] = []
    for index, smiles in enumerate(ligands_unique[smiles_column], start=1):
        fingerprint = mol2morgan_fp(
            smiles,
            radius=ECFP4_RADIUS,
            nbits=ECFP4_NBITS,
        )
        binary_fingerprints.append(DataStructs.BitVectToBinaryText(fingerprint))
        if index % 5_000 == 0 or index == total:
            elapsed = perf_counter() - started
            rate = index / elapsed
            LOG.info(
                "ligand fingerprint progress: "
                f"processed={index}/{total} rate={rate:.1f}/s "
                f"eta_seconds={(total - index) / rate:.1f}"
            )
    ligands_unique["fingerprint"] = binary_fingerprints

    from plinder.data.utils.annotations.ligand_utils import parse_cofactors

    component_path = data_dir / "dbs" / "components" / "components.parquet"
    if not component_path.is_file():
        raise FileNotFoundError(f"missing CCD component table: {component_path}")
    components = pd.read_parquet(
        component_path,
        columns=["binder_id", "canonical_smiles"],
    ).dropna(subset=["canonical_smiles"])
    cofactor_codes = parse_cofactors(data_dir)
    cofactor_smiles = dict(
        components.loc[
            components["binder_id"].isin(cofactor_codes),
            ["binder_id", "canonical_smiles"],
        ].itertuples(index=False, name=None)
    )
    missing_cofactors = cofactor_codes.difference(cofactor_smiles)
    if missing_cofactors:
        LOG.warning(
            f"{len(missing_cofactors)} cofactor codes have no CCD SMILES structure"
        )
    ligands_unique = annotate_cofactor_similarity(
        ligands_unique,
        cofactor_smiles=cofactor_smiles,
        threshold=cofactor_similarity_threshold,
        radius=ECFP4_RADIUS,
        nbits=ECFP4_NBITS,
    )

    output_dir = data_dir / "fingerprints"
    output_dir.mkdir(exist_ok=True, parents=True)
    write_ecfp4_fingerprint_table(
        ligands_unique,
        output_dir / "ligands_per_smiles.parquet",
    )

    # Fingerprints define the node universe, so any old score shards are stale.
    for path in (data_dir / "ligand_scores").glob("*.parquet"):
        path.unlink()


def ligand_scores(
    *,
    ligand_ids: list[int],
    data_dir: Path,
    output_path: Path,
    number_id_col: str = "ligand_smiles_id",
    minimum_similarity: float = 30.0,
) -> None:
    """Write all BulkTanimoto edges above ``minimum_similarity``."""
    fingerprint_path = data_dir / "fingerprints" / "ligands_per_smiles.parquet"
    fingerprint_metadata = pq.read_schema(fingerprint_path).metadata or {}
    if any(
        fingerprint_metadata.get(key) != value
        for key, value in ECFP4_PARQUET_METADATA.items()
    ):
        raise ValueError("ligand fingerprint metadata is not ECFP4/1024")
    all_ligands = pd.read_parquet(fingerprint_path)
    fingerprints = [
        DataStructs.CreateFromBinaryText(value) for value in all_ligands["fingerprint"]
    ]
    if len(fingerprints) != len(all_ligands):
        raise ValueError("ligands don't match fingerprints")
    node_ids = all_ligands[number_id_col].astype(int).tolist()
    if node_ids != list(range(len(all_ligands))):
        raise ValueError("ligand IDs must match contiguous fingerprint row indices")

    minimum_fraction = minimum_similarity / 100.0
    rows: list[dict[str, int | float]] = []
    for ligand_id in ligand_ids:
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[ligand_id], fingerprints
        )
        rows.extend(
            {
                "query_ligand_id": ligand_id,
                "target_ligand_id": target_id,
                "tanimoto_similarity_ecfp4_1024": similarity * 100.0,
            }
            for target_id, similarity in enumerate(similarities)
            if similarity >= minimum_fraction
        )
    table = pa.Table.from_pylist(
        rows,
        schema=schemas.TANIMOTO_SCORE_SCHEMA.with_metadata(ECFP4_PARQUET_METADATA),
    )
    pq.write_table(table, output_path)


def build_ligand_similarity_annotations(
    *,
    unique_ligands: pd.DataFrame,
    ligand_occurrences: pd.DataFrame,
    edges: pd.DataFrame,
    cluster_threshold: float = 90.0,
) -> pd.DataFrame:
    """Build deterministic Tanimoto components and PDB-frequency annotations."""
    node_ids = unique_ligands["ligand_smiles_id"].astype(int).tolist()
    parent = {node_id: node_id for node_id in node_ids}

    def find(node_id: int) -> int:
        while parent[node_id] != node_id:
            parent[node_id] = parent[parent[node_id]]
            node_id = parent[node_id]
        return node_id

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            parent[right_root] = left_root
        else:
            parent[left_root] = right_root

    selected_edges = edges[edges["tanimoto_similarity_ecfp4_1024"] >= cluster_threshold]
    for left, right in selected_edges[
        ["query_ligand_id", "target_ligand_id"]
    ].itertuples(index=False, name=None):
        left_id = int(left)
        right_id = int(right)
        if left_id not in parent or right_id not in parent:
            raise ValueError(
                f"ligand similarity edge references unknown node {left_id}, {right_id}"
            )
        union(left_id, right_id)

    components: dict[int, list[int]] = {}
    for node_id in node_ids:
        components.setdefault(find(node_id), []).append(node_id)
    ordered_components = sorted(
        components.values(), key=lambda members: (-len(members), members)
    )
    cluster_by_node = {
        node_id: f"c{cluster_index}"
        for cluster_index, members in enumerate(ordered_components)
        for node_id in members
    }
    pdb_ids_by_node = {
        int(node_id): set(group["pdb_id"].astype(str))
        for node_id, group in ligand_occurrences.groupby("ligand_smiles_id")
    }
    num_pdb_ids_by_node: dict[int, int] = {}
    for members in ordered_components:
        component_pdb_ids: set[str] = set()
        for node_id in members:
            component_pdb_ids.update(pdb_ids_by_node.get(node_id, set()))
        for node_id in members:
            num_pdb_ids_by_node[node_id] = len(component_pdb_ids)

    threshold_label = f"{cluster_threshold:g}".replace(".", "p")
    cluster_column = f"ligand_tanimoto_ecfp4_1024_{threshold_label}_cluster"
    count_column = f"{cluster_column}_num_pdb_ids"
    annotations = unique_ligands.drop(columns=["fingerprint"]).copy()
    annotations[cluster_column] = annotations["ligand_smiles_id"].map(cluster_by_node)
    annotations[count_column] = (
        annotations["ligand_smiles_id"].map(num_pdb_ids_by_node).astype(np.int32)
    )
    return annotations


def annotate_ligand_similarity(
    *, data_dir: Path, cluster_threshold: float = 90.0
) -> Path:
    """Collate score shards into ligand-level annotations for the final index."""
    fingerprint_dir = data_dir / "fingerprints"
    unique_ligands = pd.read_parquet(fingerprint_dir / "ligands_per_smiles.parquet")
    ligand_occurrences = load_ligands_from_annotation_table(data_dir=data_dir).merge(
        unique_ligands[["ligand_rdkit_canonical_smiles", "ligand_smiles_id"]],
        on="ligand_rdkit_canonical_smiles",
        how="inner",
        validate="many_to_one",
    )
    score_paths = sorted((data_dir / "ligand_scores").glob("*.parquet"))
    if len(unique_ligands) and not score_paths:
        raise FileNotFoundError("no BulkTanimoto score shards were generated")
    frames = [
        pd.read_parquet(
            path,
            filters=[("tanimoto_similarity_ecfp4_1024", ">=", cluster_threshold)],
        )
        for path in score_paths
    ]
    edges = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(
            columns=[
                "query_ligand_id",
                "target_ligand_id",
                "tanimoto_similarity_ecfp4_1024",
            ]
        )
    )
    expected_query_ids = set(unique_ligands["ligand_smiles_id"].astype(int))
    observed_query_ids = set(edges["query_ligand_id"].astype(int))
    if observed_query_ids != expected_query_ids:
        missing = sorted(expected_query_ids.difference(observed_query_ids))
        extra = sorted(observed_query_ids.difference(expected_query_ids))
        raise ValueError(
            "BulkTanimoto score shards do not cover the fingerprint set: "
            f"missing={missing[:10]}, extra={extra[:10]}"
        )
    annotations = build_ligand_similarity_annotations(
        unique_ligands=unique_ligands,
        ligand_occurrences=ligand_occurrences,
        edges=edges,
        cluster_threshold=cluster_threshold,
    )
    output_path = fingerprint_dir / "ligand_similarity_annotations.parquet"
    annotations.to_parquet(output_path, index=False)
    return output_path


PHARMACOPHORE_FEATURES = frozenset(
    {
        "Donor",
        "Acceptor",
        "NegIonizable",
        "PosIonizable",
        "ZnBinder",
        "Aromatic",
        "Hydrophobe",
        "LumpedHydrophobe",
    }
)
_PharmacophoreFeatures = tuple[Any, ...]
# RDKit ShapeAlign supplies default radii for the usual organic elements but
# raises for coordination metals such as HEM iron. Other elements receive a
# custom periodic-table radius without changing the molecule or canonical SDF.
SHAPE_DEFAULT_ATOMIC_NUMBERS = frozenset({1, 6, 7, 8, 9, 15, 16, 17, 35, 53})
LIGAND_3D_SCORE_ABILITY_CACHE_SIZE = 4096
_LIGAND_3D_SCORE_ABILITY_CACHE: OrderedDict[str, bool] = OrderedDict()


@cache
def _get_feature_scoring_context() -> tuple[Any, dict[str, FeatMaps.FeatMapParams]]:
    """Load RDKit's pharmacophore definitions only when shape scoring is used."""
    build_feature_factory = cast(Any, ChemicalFeatures).BuildFeatureFactory
    factory = build_feature_factory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    parameters = {
        family: FeatMaps.FeatMapParams() for family in factory.GetFeatureFamilies()
    }
    return factory, parameters


def align_molecules(
    reference: Chem.Mol,
    mobile: Chem.Mol,
    max_preiters: int = 100,
    max_postiters: int = 100,
) -> tuple[float, float]:
    """Shape-align ``mobile`` onto ``reference`` and return shape/color scores."""
    reference_options = _shape_input_options(reference)
    mobile_options = _shape_input_options(mobile)
    shape, color = rdShapeAlign.AlignMol(
        reference,
        mobile,
        reference_options,
        mobile_options,
        -1,
        -1,
        1.0,
        max_preiters,
        max_postiters,
    )
    return float(shape), float(color)


def _shape_input_options(molecule: Chem.Mol) -> Any:
    """Supply VdW radii for atoms outside ShapeAlign's default organic set."""
    options = rdShapeAlign.ShapeInputOptions()
    periodic_table = Chem.GetPeriodicTable()
    options.atomRadii = [
        (atom.GetIdx(), float(periodic_table.GetRvdw(atom.GetAtomicNum())))
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() not in SHAPE_DEFAULT_ATOMIC_NUMBERS
        and atom.GetAtomicNum() > 0
    ]
    return options


def get_feature_map_score(
    mol_1: Chem.Mol,
    mol_2: Chem.Mol,
    score_mode: int = FeatMaps.FeatMapScoreMode.All,
    features_1: Sequence[Any] | None = None,
    features_2: Sequence[Any] | None = None,
) -> float:
    """Calculate the normalized pharmacophore feature overlap."""
    factory, parameters = _get_feature_scoring_context()
    molecules_and_features = ((mol_1, features_1), (mol_2, features_2))
    feature_lists = [
        list(features)
        if features is not None
        else [
            feature
            for feature in factory.GetFeaturesForMol(molecule)
            if feature.GetFamily() in PHARMACOPHORE_FEATURES
        ]
        for molecule, features in molecules_and_features
    ]
    denominator = min(len(features) for features in feature_lists)
    if denominator == 0:
        return 0.0
    feature_map_type = cast(Any, FeatMaps.FeatMap)
    feature_map = feature_map_type(
        feats=feature_lists[0],
        weights=[1] * len(feature_lists[0]),
        params=parameters,
    )
    feature_map.scoreMode = score_mode
    return float(feature_map.ScoreFeats(feature_lists[1]) / denominator)


def get_sucos_score(
    mol_1: Chem.Mol,
    mol_2: Chem.Mol,
    score_mode: int = FeatMaps.FeatMapScoreMode.All,
    features_1: Sequence[Any] | None = None,
    features_2: Sequence[Any] | None = None,
) -> float:
    """Calculate SuCOS from feature overlap and shape protrusion distance."""
    feature_map_score = float(
        np.clip(
            get_feature_map_score(
                mol_1,
                mol_2,
                score_mode,
                features_1,
                features_2,
            ),
            0,
            1,
        )
    )
    protrude_distance = float(
        np.clip(
            rdShapeHelpers.ShapeProtrudeDist(mol_1, mol_2, allowReordering=False),
            0,
            1,
        )
    )
    return 0.5 * feature_map_score + 0.5 * (1 - protrude_distance)


def _pharmacophore_features(molecule: Chem.Mol) -> _PharmacophoreFeatures:
    """Detect the pharmacophore features of an unmodified reference ligand."""
    factory, _ = _get_feature_scoring_context()
    return tuple(
        feature
        for feature in factory.GetFeaturesForMol(molecule)
        if feature.GetFamily() in PHARMACOPHORE_FEATURES
    )


def _prepare_unsanitized_sdf_molecule(molecule: Chem.Mol, *, label: str) -> Chem.Mol:
    """Prepare a permissively parsed molecule without rejecting bad valence."""
    try:
        molecule.UpdatePropertyCache(strict=False)
    except Exception as exc:
        LOG.debug("property-cache update failed for %s: %s", label, exc)
        return molecule
    sanitize_ops = (
        Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
    )
    try:
        failed_operation = Chem.SanitizeMol(
            molecule,
            sanitizeOps=sanitize_ops,
            catchErrors=True,
        )
    except Exception as exc:
        LOG.debug("partial sanitization failed for %s: %s", label, exc)
    else:
        if failed_operation != Chem.SanitizeFlags.SANITIZE_NONE:
            LOG.debug(
                "partial sanitization stopped at operation %s for %s",
                failed_operation,
                label,
            )
    return molecule


def load_sdf_molecule(sdf_file: Path) -> Chem.Mol | None:
    """Load one SDF, falling back to permissive parsing for bad valence."""
    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromMolFile(str(sdf_file))
    except Exception:
        molecule = None
    if molecule is None:
        try:
            with rdBase.BlockLogs():
                molecule = Chem.MolFromMolFile(
                    str(sdf_file), sanitize=False, strictParsing=False
                )
        except Exception as exc:
            LOG.warning("failed loading ligand SDF %s: %s", sdf_file, exc)
            return None
        if molecule is not None:
            molecule = _prepare_unsanitized_sdf_molecule(molecule, label=str(sdf_file))
    if molecule is None or molecule.GetNumConformers() == 0:
        return None
    return molecule


def load_sdf_molecule_block(sdf: bytes, *, label: str) -> Chem.Mol | None:
    """Load one packed SDF record, with a permissive bad-valence fallback."""
    try:
        block = sdf.decode("utf-8")
    except (AttributeError, UnicodeDecodeError) as exc:
        LOG.warning("failed decoding packed ligand SDF %s: %s", label, exc)
        return None
    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromMolBlock(block)
    except Exception:
        molecule = None
    if molecule is None:
        try:
            with rdBase.BlockLogs():
                molecule = Chem.MolFromMolBlock(
                    block, sanitize=False, strictParsing=False
                )
        except Exception as exc:
            LOG.warning("failed loading packed ligand SDF %s: %s", label, exc)
            return None
        if molecule is not None:
            molecule = _prepare_unsanitized_sdf_molecule(molecule, label=label)
    if molecule is None or molecule.GetNumConformers() == 0:
        return None
    return molecule


def prepare_ligand_for_3d_scoring(molecule: Chem.Mol) -> Chem.Mol | None:
    """Clone a ligand after checking that it contains an explicit heavy atom."""
    if not any(atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms()):
        return None
    return Chem.Mol(molecule)


def canonical_ligand_sdf_path(data_dir: Path, *, pdb_id: str, asym_id: str) -> Path:
    """Return the canonical ASU SDF path for one ligand occurrence."""
    return (
        data_dir
        / "raw_entries"
        / pdb_id[-3:-1]
        / pdb_id
        / "ligand_files"
        / f"{asym_id}.sdf"
    )


def is_ligand_3d_score_able(sdf_file: Path) -> bool:
    """Test whether one canonical SDF supports shape, color, and SuCOS scoring."""
    molecule = load_sdf_molecule(sdf_file)
    if molecule is None:
        return False
    prepared = prepare_ligand_for_3d_scoring(molecule)
    if prepared is None:
        return False
    positions = prepared.GetConformer().GetPositions()
    if not np.isfinite(positions).all():
        LOG.warning(f"ligand SDF has non-finite coordinates ({sdf_file})")
        return False

    # Capability depends on the sanitized molecular graph, while each ASU copy
    # differs only in its coordinates. Loading every SDF above still validates
    # its conformer; reusing a successful graph-level probe avoids repeatedly
    # self-aligning large repeated ligands (e.g. symmetry-related oligos).
    try:
        cache_key = Chem.MolToSmiles(
            prepared,
            canonical=True,
            isomericSmiles=True,
        )
    except Exception:
        cache_key = None
    if cache_key is not None and cache_key in _LIGAND_3D_SCORE_ABILITY_CACHE:
        _LIGAND_3D_SCORE_ABILITY_CACHE.move_to_end(cache_key)
        return _LIGAND_3D_SCORE_ABILITY_CACHE[cache_key]

    reference = Chem.Mol(prepared)
    mobile = Chem.Mol(prepared)
    try:
        shape, color = align_molecules(reference, mobile)
        sucos = get_sucos_score(reference, mobile)
    except Exception as exc:
        LOG.warning(f"ligand SDF is not 3D-scoreable ({sdf_file}): {exc}")
        return False
    score_able = all(np.isfinite(value) for value in (shape, color, sucos))
    if score_able and cache_key is not None:
        _LIGAND_3D_SCORE_ABILITY_CACHE[cache_key] = True
        _LIGAND_3D_SCORE_ABILITY_CACHE.move_to_end(cache_key)
        while len(_LIGAND_3D_SCORE_ABILITY_CACHE) > LIGAND_3D_SCORE_ABILITY_CACHE_SIZE:
            _LIGAND_3D_SCORE_ABILITY_CACHE.popitem(last=False)
    return score_able


def annotate_ligand_3d_score_ability(
    ligands: pd.DataFrame, *, data_dir: Path
) -> pd.DataFrame:
    """Annotate each ligand occurrence using its canonical ASU conformation."""
    required = {"pdb_id", "ligand_asym_id"}
    missing = required.difference(ligands.columns)
    if missing:
        raise ValueError(f"missing ligand SDF identifiers: {sorted(missing)}")
    result = ligands.copy()
    scoreability_by_sdf: dict[Path, bool] = {}
    abilities = []
    for pdb_id, asym_id in result[["pdb_id", "ligand_asym_id"]].itertuples(
        index=False, name=None
    ):
        sdf_file = canonical_ligand_sdf_path(
            data_dir,
            pdb_id=str(pdb_id),
            asym_id=str(asym_id),
        )
        if sdf_file not in scoreability_by_sdf:
            scoreability_by_sdf[sdf_file] = is_ligand_3d_score_able(sdf_file)
        abilities.append(scoreability_by_sdf[sdf_file])
    result["ligand_is_3d_score_able"] = pd.array(abilities, dtype="boolean")
    return result


def get_sequence_similarity(seq_str1: str, seq_str2: str) -> tuple[float, float]:
    """
    Calculate the protein similarity score of an alignment.

    If the alignment contains more than two protein sequences,
    all pairwise scores are counted.

    Parameters
    ----------
    seq_str1 : str
        First protein sequence
    seq_str2 : str
        Second protein sequence

    Returns
    -------
    tuple[float, float]
        Sequence identity and sequence similarity score.
    """
    # Keep the established scoring semantics: remove gaps independently before
    # comparing positions, and map J/U/O to their ProteinSequence equivalents.
    # A byte lookup avoids constructing Biotite sequence/alignment objects for
    # every search hit, which dominates mapping time for dense query chains.
    s1 = seq_str1.translate(_NON_CANONICAL_AA).replace("-", "").encode("ascii")
    s2 = seq_str2.translate(_NON_CANONICAL_AA).replace("-", "").encode("ascii")
    if not s1 or len(s1) != len(s2):
        raise ValueError("aligned sequences must have equal non-zero lengths")
    codes1 = np.frombuffer(s1, dtype=np.uint8)
    codes2 = np.frombuffer(s2, dtype=np.uint8)
    if not np.all(_VALID_PROTEIN_AA[codes1]) or not np.all(_VALID_PROTEIN_AA[codes2]):
        raise ValueError("aligned sequence contains an unsupported residue")
    identity = float(np.mean(codes1 == codes2))
    similarity = float(np.mean(_SIMILAR_PROTEIN_AA[codes1, codes2]))
    return identity, similarity


def get_sequence_similarity_helper(seq_str1: str, seq_str2: str) -> float:
    try:
        return get_sequence_similarity(seq_str1, seq_str2)[1]
    except Exception:
        return 0


def _alignment_search_command(
    *,
    aln_type: str,
    query_db: Path,
    target_db: Path,
    search_db: Path,
    tmp_dir: Path,
    alignment_config: FoldseekConfig | MMSeqsConfig,
    threads: int,
    expand_exact_clusters: bool = False,
) -> list[str]:
    """Build a representative-target search from the configured filters."""
    command = [
        aln_type,
        "search",
        str(query_db),
        str(target_db),
        str(search_db),
        str(tmp_dir),
        "-a",
        "-e",
        f"{alignment_config.evalue}",
        "-s",
        f"{alignment_config.sensitivity}",
        "--max-seqs",
        f"{alignment_config.max_seqs}",
        "-c",
        f"{alignment_config.coverage}",
        "--cov-mode",
        "2",
        "--min-seq-id",
        f"{alignment_config.min_seq_id}",
    ]
    if aln_type == "mmseqs":
        if not isinstance(alignment_config, MMSeqsConfig):
            raise TypeError("mmseqs search requires MMSeqsConfig")
        if expand_exact_clusters:
            raise ValueError("MMseqs cluster expansion is a separate workflow")
    elif aln_type == "foldseek":
        if not isinstance(alignment_config, FoldseekConfig):
            raise TypeError("foldseek search requires FoldseekConfig")
        command.extend(["--sort-by-structure-bits", "0"])
        if expand_exact_clusters:
            command.extend(["--cluster-search", "1"])
    else:
        raise ValueError(f"unsupported alignment type: {aln_type}")
    command.extend(["--threads", str(threads)])
    return command


def run_alignment(
    aln_type: str,
    query_db: Path,
    target_db: Path,
    search_target_db: Path,
    search_db: Path,
    aln_file: Path,
    alignment_config: FoldseekConfig | MMSeqsConfig,
    cluster_alignment_db: Path | None = None,
    tmp_dir: Path = Path.cwd() / "tmp",
    remove_tmp: bool = True,
    threads: int = 1,
) -> None:
    def remove_search_results() -> None:
        for filename in search_db.parent.glob(f"{search_db.name}*"):
            if filename.is_dir():
                shutil.rmtree(filename)
            elif filename.is_file():
                filename.unlink()

    # A failed search may leave a partial result DB. Always start this
    # scratch-scoped operation from a clean prefix so Slurm retries work.
    remove_search_results()

    representative_result = (
        search_db.parent / f"{search_db.name}_representatives"
        if aln_type == "mmseqs"
        else search_db
    )
    search_commands = _alignment_search_command(
        aln_type=aln_type,
        query_db=query_db,
        target_db=search_target_db,
        search_db=representative_result,
        tmp_dir=tmp_dir,
        alignment_config=alignment_config,
        threads=threads,
        expand_exact_clusters=aln_type == "foldseek",
    )
    format_output = (
        "query,target,qlen,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,"
        "qcov,tcov,qaln,taln"
    )
    if aln_type == "foldseek":
        format_output += ",lddt"
    subprocess.check_call(search_commands, stdout=subprocess.DEVNULL)
    if aln_type == "mmseqs":
        if cluster_alignment_db is None:
            raise ValueError("MMseqs clustered search requires cluster alignments")
        expanded_result = search_db.parent / f"{search_db.name}_expanded"
        subprocess.check_call(
            [
                "mmseqs",
                "expandaln",
                str(query_db),
                str(target_db),
                str(representative_result),
                str(cluster_alignment_db),
                str(expanded_result),
                "--expansion-mode",
                "0",
                "--threads",
                str(threads),
            ],
            stdout=subprocess.DEVNULL,
        )
        # Re-align expanded members so scores and E-values are calculated for
        # the full target DB rather than copied from its representatives.
        subprocess.check_call(
            [
                "mmseqs",
                "align",
                str(query_db),
                str(target_db),
                str(expanded_result),
                str(search_db),
                "-a",
                "-e",
                f"{alignment_config.evalue}",
                "-c",
                f"{alignment_config.coverage}",
                "--cov-mode",
                "2",
                "--min-seq-id",
                f"{alignment_config.min_seq_id}",
                "--threads",
                str(threads),
            ],
            stdout=subprocess.DEVNULL,
        )

    # Foldseek's cluster-search DB contains the expanded member records needed
    # by convertalis, so a portable search bundle need not also ship full_db.
    conversion_target_db = search_target_db if aln_type == "foldseek" else target_db
    convert_commands = [
        aln_type,
        "convertalis",
        str(query_db),
        str(conversion_target_db),
        str(search_db),
        str(aln_file.with_suffix(".tsv")),
        "--format-mode",
        "4",
        "--format-output",
        format_output,
        "--threads",
        str(threads),
    ]
    subprocess.check_call(convert_commands, stdout=subprocess.DEVNULL)

    _stream_alignment_tsv_to_dataset(
        aln_file.with_suffix(".tsv"),
        aln_file.with_suffix(".parquet"),
        include_target_pdb_id=search_db != "pred",
    )

    # Cleanup
    if remove_tmp and tmp_dir.exists():
        shutil.rmtree(tmp_dir)
        aln_file.with_suffix(".tsv").unlink()
    remove_search_results()


def _pdb_id_from_alignment_identifier(identifier: Any) -> str:
    return (
        str(identifier)
        .replace("_xyz-enrich.cif.gz", "")
        .replace("_xyz-enrich.cif", "")
        .replace("_xyz-enrich", "")
        .replace(".cif.gz", "")
        .replace(".cif", "")
        .replace("pdb_0000", "")[:4]
    )


def _pdb_ids_from_alignment_identifiers(identifiers: Any) -> Any:
    """Extract four-character PDB IDs without a Python loop over Arrow rows."""
    without_foldseek_prefix = pc.replace_substring(
        identifiers,
        pattern="pdb_0000",
        replacement="",
    )
    return pc.utf8_slice_codeunits(without_foldseek_prefix, start=0, stop=4)


def _raw_alignment_schema(aln_type: str) -> pa.Schema:
    fields = [
        pa.field("query", pa.string()),
        pa.field("target", pa.string()),
        pa.field("qlen", pa.int64()),
        pa.field("fident", pa.float64()),
        pa.field("alnlen", pa.int64()),
        pa.field("qstart", pa.int64()),
        pa.field("qend", pa.int64()),
        pa.field("tstart", pa.int64()),
        pa.field("tend", pa.int64()),
        pa.field("evalue", pa.float64()),
        pa.field("bits", pa.int64()),
        pa.field("qcov", pa.float64()),
        pa.field("tcov", pa.float64()),
        pa.field("qaln", pa.string()),
        pa.field("taln", pa.string()),
    ]
    if aln_type == "foldseek":
        fields.append(pa.field("lddt", pa.float64()))
    elif aln_type != "mmseqs":
        raise ValueError(f"unknown alignment type: {aln_type}")
    fields.append(pa.field("target_pdb_id", pa.string()))
    return pa.schema(fields)


def _stream_alignment_tsv_to_dataset(
    tsv_path: Path,
    dataset_path: Path,
    *,
    include_target_pdb_id: bool,
) -> None:
    """Convert an alignment TSV without materializing the full search batch."""
    if dataset_path.exists():
        shutil.rmtree(dataset_path)
    reader = csv.open_csv(
        tsv_path,
        parse_options=csv.ParseOptions(delimiter="\t"),
        read_options=csv.ReadOptions(block_size=16 * 1024 * 1024),
    )
    for batch_index, batch in enumerate(reader):
        table = pyarrow.Table.from_batches([batch])
        table = table.append_column(
            "query_pdb_id",
            _pdb_ids_from_alignment_identifiers(table["query"]),
        )
        if include_target_pdb_id:
            table = table.append_column(
                "target_pdb_id",
                _pdb_ids_from_alignment_identifiers(table["target"]),
            )
        pq.write_to_dataset(
            table,
            dataset_path,
            partition_cols=["query_pdb_id"],
            basename_template=f"batch-{batch_index}-{{i}}.parquet",
        )


def combine_scores(
    q_t_scores: _SimilarityScoreDictType,
    q_t_mappings: dict[str, list[_ChainPairType]],
    protein_chain_mapper: str = "",
) -> dict[str, str | float]:
    def make_mapping(mapping: list[_ChainPairType]) -> str:
        return ";".join(f"{a}:{b}" if len(b) else a for a, b in mapping)

    q_t_scores_combined: dict[str, str | float] = {}
    for s in SCORE_NAMES:
        for suffix in ["_max", "_weighted_max", "_weighted_sum", ""]:
            if (
                f"{s}_foldseek{suffix}" not in q_t_scores
                and f"{s}_mmseqs{suffix}" not in q_t_scores
            ):
                continue
            s1, s2 = (
                float(q_t_scores.get(f"{s}_foldseek{suffix}", 0)),
                float(q_t_scores.get(f"{s}_mmseqs{suffix}", 0)),
            )
            source = "foldseek"
            if np.abs(s1 - s2) < 1e-6:
                q_t_scores_combined[f"{s}{suffix}"] = s1
                q_t_scores_combined[f"{s}{suffix}_source"] = "both"
            else:
                index = np.argmax([s1, s2])
                q_t_scores_combined[f"{s}{suffix}"] = [s1, s2][index]
                source = ["foldseek", "mmseqs"][index]
                q_t_scores_combined[f"{s}{suffix}_source"] = source
            source_name = "mmseqs" if source == "mmseqs" else "foldseek"
            n = f"{s}_{source_name}{suffix}"
            if n in q_t_mappings:
                q_t_scores_combined[f"{s}{suffix}_mapping"] = make_mapping(
                    q_t_mappings[n]
                )
    if protein_chain_mapper != "":
        q_t_scores_combined["protein_mapping"] = make_mapping(
            q_t_mappings.get(f"{protein_chain_mapper}_weighted_sum", [])
        )
        q_t_scores_combined["protein_mapper"] = (
            "foldseek" if "lddt" in protein_chain_mapper else "mmseqs"
        )
    return q_t_scores_combined


@dataclass
class Scorer:
    entries: dict[str, EntryView]
    source_to_full_db_file: dict[str, Path]
    db_dir: Path
    scores_dir: Path
    ligand_sdf_resolver: Callable[[LigandView], Path | None] | None = field(
        default=None, repr=False
    )
    foldseek_config: FoldseekConfig = field(default_factory=FoldseekConfig)
    mmseqs_config: MMSeqsConfig = field(default_factory=MMSeqsConfig)
    protein_chain_mappers: list[str] = field(
        default_factory=lambda: [
            "protein_lddt_qcov_foldseek",
            "protein_fident_qcov_mmseqs",
        ]
    )
    ligand_chain_mappers: list[str] = field(
        default_factory=lambda: [
            "pocket_fident_qcov_foldseek",
            "pocket_fident_qcov_mmseqs",
        ]
    )
    # default minimum threshold for every metric BEFORE scaling to 100
    minimum_threshold: float = 0.0
    # optional custom minimum threshold for each metric BEFORE scaling to 100
    minimum_thresholds: dict[str, float] = field(default_factory=dict)
    max_query_protein_chains: int = 30
    max_query_proper_ligand_chains: int = 30
    shape_score_threads: int = 1
    _ligand_mol_cache: dict[tuple[str, str], Chem.Mol | None] = field(
        default_factory=dict, init=False, repr=False
    )
    _ligand_shape_score_cache: dict[
        tuple[tuple[str, str], tuple[str, str]], _SimilarityScoreDictType
    ] = field(default_factory=dict, init=False, repr=False)
    _ligand_reference_feature_cache: dict[
        tuple[str, str], _PharmacophoreFeatures | None
    ] = field(default_factory=dict, init=False, repr=False)
    _ligand_sdf_block_cache: dict[tuple[str, str], bytes] = field(
        default_factory=dict, init=False, repr=False
    )
    _loaded_ligand_archive_entries: set[str] = field(
        default_factory=set, init=False, repr=False
    )
    _ligand_pocket_data_cache: dict[str, _PocketDataType] = field(
        default_factory=dict, init=False, repr=False
    )
    _protein_score_cache: dict[
        tuple[str, str, tuple[str, ...], tuple[str, ...]], _ProteinScoreResultType
    ] = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.shape_score_threads < 1:
            raise ValueError("shape_score_threads must be positive")
        if self.max_query_protein_chains < 1 or self.max_query_proper_ligand_chains < 1:
            raise ValueError("scoring system chain limits must be positive")
        self.db_dir.mkdir(exist_ok=True, parents=True)
        self.scores_dir.mkdir(exist_ok=True, parents=True)

    def system_is_query_scoreable(self, system: SystemView) -> bool:
        """Return whether a holo system is small enough to score as a query."""
        proper_ligand_count = sum(
            ligand.is_proper for ligand in system.ligands.values()
        )
        return (
            system.system_type == "holo"
            and 0 < len(system.protein_chains_asym_id) <= self.max_query_protein_chains
            and 0 < proper_ligand_count <= self.max_query_proper_ligand_chains
        )

    @staticmethod
    def system_is_target_scoreable(system: SystemView) -> bool:
        """Return whether a system may contribute scores as a target.

        Query limits exist to bound the work initiated by an entry.  Applying
        them to targets as well would remove large systems from the score graph
        entirely, so every protein-containing holo system remains available as
        a target. Only proper ligands are emitted later by
        :meth:`get_scores_holo`.
        """
        return system.system_type == "holo" and bool(system.protein_chains_asym_id)

    def resolve_ligand_sdf(self, data_dir: Path, ligand: LigandView) -> Path | None:
        """Resolve the canonical ASU ligand SDF for a ligand annotation."""
        candidates = [
            canonical_ligand_sdf_path(
                data_dir,
                pdb_id=ligand.pdb_id,
                asym_id=ligand.asym_id,
            ),
            data_dir
            / "ligand_archives"
            / ligand.pdb_id
            / "ligand_files"
            / f"{ligand.asym_id}.sdf",
        ]
        return next((path for path in candidates if path.is_file()), None)

    def _get_ligand_mol(self, data_dir: Path, ligand: LigandView) -> Chem.Mol | None:
        cache_key = (ligand.pdb_id, ligand.asym_id)
        if cache_key not in self._ligand_mol_cache:
            packed = self._ligand_sdf_block_cache.pop(cache_key, None)
            if packed is not None:
                molecule = load_sdf_molecule_block(
                    packed, label=f"{ligand.pdb_id}/{ligand.asym_id}"
                )
                self._ligand_mol_cache[cache_key] = (
                    prepare_ligand_for_3d_scoring(molecule)
                    if molecule is not None
                    else None
                )
                return self._ligand_mol_cache[cache_key]
            if self.ligand_sdf_resolver is None:
                sdf_file = self.resolve_ligand_sdf(data_dir, ligand)
            else:
                sdf_file = self.ligand_sdf_resolver(ligand)
            if sdf_file is None:
                LOG.warning(
                    "no ligand SDF found for "
                    f"{ligand.id} (canonical key={cache_key})"
                )
                self._ligand_mol_cache[cache_key] = None
            else:
                molecule = load_sdf_molecule(sdf_file)
                self._ligand_mol_cache[cache_key] = (
                    prepare_ligand_for_3d_scoring(molecule)
                    if molecule is not None
                    else None
                )
        return self._ligand_mol_cache[cache_key]

    def _preload_packed_ligand_sdfs(
        self,
        data_dir: Path,
        ligands: Sequence[LigandView],
        *,
        requested_entries_only: bool = False,
    ) -> None:
        """Bulk-read packed canonical SDFs once per requested PDB entry."""
        if self.ligand_sdf_resolver is not None:
            return
        requested_shards: set[str] = set()
        for ligand in ligands:
            cache_key = (ligand.pdb_id, ligand.asym_id)
            if (
                cache_key not in self._ligand_mol_cache
                and ligand.pdb_id not in self._loaded_ligand_archive_entries
            ):
                requested_shards.add(ligand.pdb_id[1:3])

        archive_root = data_dir / "ligand_archives"
        for shard in sorted(requested_shards):
            archive = archive_root / f"{shard}.parquet"
            if not archive.is_file():
                continue
            requested_pdb_ids = sorted(
                {
                    ligand.pdb_id
                    for ligand in ligands
                    if ligand.pdb_id[1:3] == shard
                    and ligand.pdb_id not in self._loaded_ligand_archive_entries
                }
            )
            # A shard is deliberately small (roughly one hundred PDB entries)
            # and normally one Parquet row group. Reading it once avoids
            # repeating the same NFS page read for later queries in this job.
            table = pq.read_table(
                archive,
                columns=["pdb_id", "ligand_asym_id", "sdf"],
                filters=(
                    [("pdb_id", "in", requested_pdb_ids)]
                    if requested_entries_only
                    else None
                ),
                use_threads=True,
            )
            for pdb_id, asym_id, sdf in zip(
                table["pdb_id"].to_pylist(),
                table["ligand_asym_id"].to_pylist(),
                table["sdf"].to_pylist(),
            ):
                self._ligand_sdf_block_cache[(str(pdb_id), str(asym_id))] = sdf
            self._loaded_ligand_archive_entries.update(
                str(pdb_id) for pdb_id in table["pdb_id"].to_pylist()
            )

    def _get_ligand_reference_features(
        self, data_dir: Path, ligand: LigandView
    ) -> _PharmacophoreFeatures | None:
        """Cache pharmacophore features for an unmodified reference ligand."""
        cache_key = (ligand.pdb_id, ligand.asym_id)
        if cache_key not in self._ligand_reference_feature_cache:
            molecule = self._get_ligand_mol(data_dir, ligand)
            self._ligand_reference_feature_cache[cache_key] = (
                _pharmacophore_features(molecule) if molecule is not None else None
            )
        return self._ligand_reference_feature_cache[cache_key]

    def get_ligand_pair_shape_scores(
        self,
        data_dir: Path,
        query_ligand: LigandView,
        target_ligand: LigandView,
        pocket_qcov: float,
    ) -> _SimilarityScoreDictType:
        """Shape-score one ligand pair after the pocket-coverage gate."""
        if not np.isfinite(pocket_qcov) or pocket_qcov <= 0:
            return {}
        if not query_ligand.is_3d_score_able or not target_ligand.is_3d_score_able:
            return {}

        query_key = (query_ligand.pdb_id, query_ligand.asym_id)
        target_key = (target_ligand.pdb_id, target_ligand.asym_id)
        cache_key = (query_key, target_key)
        if cache_key not in self._ligand_shape_score_cache:
            query_mol = self._get_ligand_mol(data_dir, query_ligand)
            target_mol = self._get_ligand_mol(data_dir, target_ligand)
            query_features = self._get_ligand_reference_features(data_dir, query_ligand)
            base_scores: _SimilarityScoreDictType = {}
            if (
                query_mol is not None
                and target_mol is not None
                and query_features is not None
            ):
                # ShapeAlign mutates the mobile conformer. Clone the cached
                # canonical ASU molecules once for this directed ligand pair,
                # then reuse its immutable scores across every system in which
                # those ligand instances occur.
                query_aligned = Chem.Mol(query_mol)
                target_aligned = Chem.Mol(target_mol)
                try:
                    shape, color = align_molecules(query_aligned, target_aligned)
                except Exception as exc:
                    LOG.warning(
                        "shape alignment failed for "
                        f"{query_ligand.id} to {target_ligand.id}: {exc}"
                    )
                else:
                    if all(np.isfinite(value) for value in (shape, color)):
                        base_scores.update(
                            {
                                "shape": float(np.clip(shape, 0, 1)),
                                "color": float(np.clip(color, 0, 1)),
                            }
                        )
                        try:
                            sucos_shape = get_sucos_score(
                                query_aligned,
                                target_aligned,
                                FeatMaps.FeatMapScoreMode.All,
                                query_features,
                                None,
                            )
                        except Exception as exc:
                            LOG.warning(
                                "SuCOS calculation failed for "
                                f"{query_ligand.id} to {target_ligand.id}: {exc}"
                            )
                        else:
                            if np.isfinite(sucos_shape):
                                base_scores["sucos_shape"] = sucos_shape
                            else:
                                LOG.warning(
                                    "SuCOS calculation returned a non-finite "
                                    f"score for {query_ligand.id} to "
                                    f"{target_ligand.id}"
                                )
                    else:
                        LOG.warning(
                            "shape alignment returned non-finite scores for "
                            f"{query_ligand.id} to {target_ligand.id}"
                        )
            self._ligand_shape_score_cache[cache_key] = base_scores

        scores = dict(self._ligand_shape_score_cache[cache_key])
        if "sucos_shape" in scores:
            scores["sucos_shape_pocket_qcov"] = scores["sucos_shape"] * pocket_qcov
        return scores

    def score_canonical_ligand_pairs(
        self,
        data_dir: Path,
        pairs: pd.DataFrame,
    ) -> pd.DataFrame:
        """Calculate base 3D metrics once for each canonical directed pair."""
        pair_columns = [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
        if pairs.empty:
            return pd.DataFrame(columns=schemas.LIGAND_3D_SCORE_SCHEMA.names)
        missing = sorted(set(pair_columns).difference(pairs.columns))
        if missing:
            raise ValueError(f"ligand 3D pair table is missing columns {missing}")
        pair_frame = pairs[pair_columns].drop_duplicates(ignore_index=True)
        ligand_views: dict[tuple[str, str], LigandView] = {}

        def ligand_view(pdb_id: str, asym_id: str) -> LigandView:
            key = (pdb_id, asym_id)
            if key not in ligand_views:
                ligand_views[key] = LigandView(
                    id=f"{pdb_id}__canonical__{asym_id}",
                    pdb_id=pdb_id,
                    system_id="",
                    instance_chain=f"0.{asym_id}",
                    asym_id=asym_id,
                    is_proper=True,
                    is_3d_score_able=True,
                    protein_chains_asym_id=[],
                    num_pocket_residues=0,
                    num_interactions=0,
                    num_unique_interactions=0,
                )
            return ligand_views[key]

        work: list[tuple[LigandView, LigandView]] = [
            (
                ligand_view(str(row.query_entry), str(row.query_ligand_asym_id)),
                ligand_view(str(row.target_entry), str(row.target_ligand_asym_id)),
            )
            for row in pair_frame.itertuples(index=False)
        ]
        self._preload_packed_ligand_sdfs(
            data_dir,
            list(ligand_views.values()),
            requested_entries_only=True,
        )

        def load_ligand(ligand: LigandView) -> Chem.Mol | None:
            return self._get_ligand_mol(data_dir, ligand)

        def score_pair(pair: tuple[LigandView, LigandView]) -> dict[str, float]:
            return self.get_ligand_pair_shape_scores(
                data_dir, pair[0], pair[1], pocket_qcov=1.0
            )

        with ThreadPoolExecutor(max_workers=self.shape_score_threads) as executor:
            list(executor.map(load_ligand, ligand_views.values()))
            for query_ligand in {
                (query.pdb_id, query.asym_id): query for query, _ in work
            }.values():
                self._get_ligand_reference_features(data_dir, query_ligand)
            scores = list(executor.map(score_pair, work))
        records = []
        for pair, score in zip(pair_frame.itertuples(index=False), scores):
            records.append(
                {
                    **{column: getattr(pair, column) for column in pair_columns},
                    "shape": score.get("shape"),
                    "color": score.get("color"),
                    "sucos_shape": score.get("sucos_shape"),
                }
            )
        return pd.DataFrame.from_records(
            records, columns=schemas.LIGAND_3D_SCORE_SCHEMA.names
        )

    def make_dbs(self) -> None:
        databases.make_sub_dbs(self.db_dir, self.source_to_full_db_file, self.entries)

    def get_config(
        self, search_db: str, aln_type: str
    ) -> FoldseekConfig | MMSeqsConfig:
        config: FoldseekConfig | MMSeqsConfig
        if aln_type == "foldseek":
            config = replace(self.foldseek_config)
        else:
            config = replace(self.mmseqs_config)
        if search_db in ["apo", "pred"]:
            config.coverage = 0.9
            config.min_seq_id = 0.9
        return config

    def run_alignments(
        self,
        entry_ids: list[str],
        search_db: str,
        output_folder: Path,
        overwrite: bool = False,
        threads: int = 1,
        alignment_types: Sequence[str] | None = None,
    ) -> None:
        output_folder.mkdir(exist_ok=True)
        failures: list[str] = []
        selected_alignment_types = list(alignment_types or ["mmseqs", "foldseek"])
        unsupported = sorted(set(selected_alignment_types) - {"foldseek", "mmseqs"})
        if unsupported:
            raise ValueError(f"unsupported alignment types: {unsupported}")
        for aln_type in selected_alignment_types:
            sub_db = output_folder / search_db / aln_type
            sub_db.mkdir(exist_ok=True, parents=True)
            db_ids = databases.get_db_ids(
                self.entries, "holo", aln_type, entry_ids=entry_ids
            )
            missing_query_ids = databases.make_sub_db(
                db_ids,
                self.source_to_full_db_file[f"holo_{aln_type}"],
                sub_db,
                aln_type,
            )
            if not db_ids or len(missing_query_ids) == len(db_ids):
                LOG.warning(
                    f"no {aln_type} query chains are available for "
                    f"{len(entry_ids)} entries"
                )
                continue
            tmp_dir = sub_db / f"tmp_{search_db}_{aln_type}"
            tmp_dir.mkdir(exist_ok=True, parents=True)
            aln_file = sub_db / f"aln_{search_db}.tsv"
            LOG.info("run_alignments calling run_alignment:")
            LOG.info(f"    sub_db={sub_db}")
            LOG.info(f"    aln_file={aln_file.with_suffix('.tsv')}")
            LOG.info(f"    aln_type={aln_type}")
            try:
                (
                    target_db,
                    search_target_db,
                    cluster_alignment_db,
                ) = databases.exact_search_database_paths(
                    self.db_dir,
                    search_db,
                    aln_type,
                )
                run_alignment(
                    aln_type=aln_type,
                    query_db=sub_db / sub_db.name,
                    target_db=target_db,
                    search_target_db=search_target_db,
                    cluster_alignment_db=cluster_alignment_db,
                    search_db=sub_db / "search",
                    aln_file=aln_file.with_suffix(".tsv"),
                    tmp_dir=tmp_dir / output_folder.stem,
                    alignment_config=self.get_config(search_db, aln_type),
                    threads=threads,
                )
            except Exception as e:
                scratch = (
                    output_folder / "scratch" / "scores" / "run_alignment_failures"
                )
                scratch.mkdir(exist_ok=True, parents=True)
                (scratch / f"{search_db}_{aln_type}.txt").write_text(f"{repr(e)}: {e}")
                LOG.error(f"scoring: Error for {search_db}_{aln_type}: {e}")
                failures.append(f"{search_db}_{aln_type}: {e}")
                continue
            aln_dir = self.db_dir / f"{search_db}_{aln_type}" / "aln"
            aln_dir.mkdir(exist_ok=True, parents=True)
            for pdb_id in tqdm(entry_ids):
                pdb_id_file = (
                    aln_file.with_suffix(".parquet") / f"query_pdb_id={pdb_id}"
                )
                target = aln_dir / f"{pdb_id}.parquet"
                local_output = (
                    output_folder
                    / "completed"
                    / search_db
                    / aln_type
                    / f"{pdb_id}.parquet"
                )
                local_output.parent.mkdir(exist_ok=True, parents=True)
                if pdb_id_file.exists():
                    pdb_id_df = pd.read_parquet(pdb_id_file)
                    pdb_id_df.to_parquet(local_output, index=False)
                else:
                    # Short chains can legitimately have no hit after the
                    # E-value filter.  A typed empty file is their durable
                    # searched/no-hit completion marker.
                    pq.write_table(
                        pa.Table.from_pylist(
                            [], schema=_raw_alignment_schema(aln_type)
                        ),
                        local_output,
                    )
                install_path = target.with_suffix(target.suffix + ".tmp")
                shutil.copyfile(local_output, install_path)
                install_path.replace(target)
                local_output.unlink(missing_ok=True)
        if failures:
            raise RuntimeError("alignment searches failed: " + "; ".join(failures))

    def map_alignment_files(
        self,
        data_dir: Path,
        pdb_id: str,
        search_db: str,
        overwrite: bool = True,
        scratch_dir: Path | None = None,
        mapped_db_dir: Path | None = None,
    ) -> list[Path]:
        """Map raw backend alignments without calculating system scores."""
        mapped_files: list[Path] = []
        failures: list[str] = []
        for aln_type in ["foldseek", "mmseqs"]:
            raw_file = (
                self.db_dir / f"{search_db}_{aln_type}" / "aln" / f"{pdb_id}.parquet"
            )
            if not raw_file.exists():
                LOG.info(f"map_alignment_files: raw_file={raw_file} does not exist")
                continue

            entries_to_load = {pdb_id}
            if search_db != "pred":
                entries_to_load |= set(
                    pd.read_parquet(raw_file, columns=["target_pdb_id"])[
                        "target_pdb_id"
                    ]
                )
            entries_to_load = entries_to_load.difference(self.entries)
            if entries_to_load:
                LOG.info(
                    f"loading {len(entries_to_load)} entries to map {pdb_id} "
                    f"against {search_db} with {aln_type}"
                )
                self.entries.update(
                    load_entry_views(pdb_ids=entries_to_load, data_dir=data_dir)
                )

            mapped_file = (
                (mapped_db_dir or self.db_dir)
                / f"{search_db}_{aln_type}"
                / "mapped_aln"
                / f"{pdb_id}.parquet"
            )
            mapped_file.parent.mkdir(exist_ok=True, parents=True)
            if not overwrite and mapped_file.exists():
                try:
                    mapped_columns = set(pq.read_schema(mapped_file).names)
                except (OSError, ValueError):
                    mapped_columns = set()
                if mapped_columns and schemas.mapped_alignment_schema_is_current(
                    mapped_columns, alignment_type=aln_type
                ):
                    mapped_files.append(mapped_file)
                    continue
                LOG.info(
                    "map_alignment_files: replacing mapped file with stale schema: "
                    f"{mapped_file}"
                )
            failure_file = (
                self.db_dir
                / "scratch"
                / "scores"
                / "map_alignment_df_failures"
                / f"{search_db}_{aln_type}_{pdb_id}.txt"
            )
            try:
                LOG.info(f"mapping {pdb_id} against {search_db} with {aln_type}")
                mapped = self.map_alignment_df(raw_file, aln_type, search_db)
                temporary_root = scratch_dir or mapped_file.parent
                temporary_root.mkdir(exist_ok=True, parents=True)
                temporary = temporary_root / (
                    f"{search_db}-{aln_type}-{pdb_id}.mapped.parquet"
                )
                mapped.to_parquet(temporary, index=True)
                install_path = mapped_file.with_suffix(mapped_file.suffix + ".tmp")
                shutil.copyfile(temporary, install_path)
                install_path.replace(mapped_file)
                temporary.unlink(missing_ok=True)
                mapped_files.append(mapped_file)
                failure_file.unlink(missing_ok=True)
            except Exception as exc:
                failure_file.parent.mkdir(exist_ok=True, parents=True)
                failure_file.write_text(f"{exc!r}: {exc}")
                failures.append(f"{search_db}_{aln_type}_{pdb_id}: {exc}")
        if failures:
            raise RuntimeError("alignment mapping failed: " + "; ".join(failures))
        return mapped_files

    def get_score_df(
        self,
        data_dir: Path,
        pdb_id: str,
        search_db: str,
        overwrite: bool = True,
        map_alignments: bool = True,
        scratch_dir: Path | None = None,
        source_to_aln_file: dict[str, Path] | None = None,
        defer_ligand_3d: bool = False,
    ) -> Path:
        """
        Convert aligmnent results to mapped alignment results. Then
        aggregate the mapped alignment results into the scores dataset.
        """
        if defer_ligand_3d and search_db != "holo":
            raise ValueError("ligand 3D scoring can only be deferred for holo scores")
        score_df_path = self.db_dir / f"search_db={search_db}" / f"{pdb_id}.parquet"
        candidate_path = (
            self.scores_dir
            / "ligand_3d_candidates"
            / f"search_db={search_db}"
            / f"shard={pdb_id[1:3]}"
            / f"{pdb_id}.parquet"
        )
        score_mode = b"deferred" if defer_ligand_3d else b"complete"
        cached_score_is_current = False
        if not overwrite and score_df_path.is_file():
            try:
                cached_schema = pq.read_schema(score_df_path)
                cached_columns = set(cached_schema.names)
                pq.read_metadata(score_df_path)
                cached_score_is_current = (
                    set(schemas.PROTEIN_SIMILARITY_SCHEMA.names).issubset(
                        cached_columns
                    )
                    and (cached_schema.metadata or {}).get(b"plinder.ligand_3d")
                    == score_mode
                )
                if defer_ligand_3d:
                    candidate_schema = pq.read_schema(candidate_path)
                    pq.read_metadata(candidate_path)
                    cached_score_is_current = cached_score_is_current and set(
                        schemas.LIGAND_3D_CANDIDATE_SCHEMA.names
                    ).issubset(candidate_schema.names)
            except (OSError, ValueError):
                cached_score_is_current = False
        if overwrite or not cached_score_is_current:
            LOG.info(f"get_score_df: aggregating scores for {pdb_id} to {search_db}")
        else:
            LOG.info(f"get_score_df: skipping existing {score_df_path}")
            return score_df_path
        if map_alignments:
            self.map_alignment_files(
                data_dir,
                pdb_id,
                search_db,
                overwrite=overwrite,
            )
        else:
            entries_to_load = {pdb_id}
            for aln_type in ["foldseek", "mmseqs"]:
                source = f"{search_db}_{aln_type}"
                mapped_file = (
                    source_to_aln_file[source]
                    if source_to_aln_file is not None and source in source_to_aln_file
                    else self.db_dir / source / "mapped_aln" / f"{pdb_id}.parquet"
                )
                if mapped_file.is_file():
                    target_entries = pd.read_parquet(
                        mapped_file,
                        columns=["target_entry"],
                        filters=[("query_entry", "==", pdb_id)],
                    )
                    if "target_entry" in target_entries.columns:
                        target_values = target_entries["target_entry"]
                    else:
                        target_values = pd.Series(
                            target_entries.index.get_level_values("target_entry")
                        )
                    entries_to_load.update(target_values.dropna().astype(str))
            entries_to_load.difference_update(self.entries)
            if entries_to_load:
                LOG.info(
                    f"loading {len(entries_to_load)} entries to score {pdb_id} "
                    f"against {search_db}"
                )
                self.entries.update(
                    load_entry_views(pdb_ids=entries_to_load, data_dir=data_dir)
                )

        try:
            score_df_path.parent.mkdir(exist_ok=True, parents=True)
            LOG.info(f"aggregating scores for {pdb_id} to {search_db}")
            # A canonical pair cannot recur across different query entries.
            # Bound this cache to one query while retaining the molecule cache
            # across a batch, where target ligands are commonly reused.
            self._ligand_shape_score_cache.clear()
            self._protein_score_cache.clear()
            ligand_3d_candidates: list[_Ligand3DCandidateType] | None = (
                [] if defer_ligand_3d else None
            )
            df = self.aggregate_scores(
                pdb_id,
                search_db=search_db,
                data_dir=None if defer_ligand_3d else data_dir,
                source_to_aln_file=source_to_aln_file,
                ligand_3d_candidates=ligand_3d_candidates,
            )
            if df is None or df.empty:
                df = pd.DataFrame(
                    {
                        name: pd.Series(dtype="object")
                        for name in schemas.PROTEIN_SIMILARITY_SCHEMA.names
                    }
                )
            temporary_root = scratch_dir or score_df_path.parent
            temporary_root.mkdir(exist_ok=True, parents=True)
            temporary = temporary_root / f"{search_db}-{pdb_id}.scores.parquet"
            score_schema = schemas.PROTEIN_SIMILARITY_SCHEMA.with_metadata(
                {b"plinder.ligand_3d": score_mode}
            )
            df.to_parquet(
                temporary,
                index=False,
                schema=score_schema,
            )
            candidate_temporary = (
                temporary_root / f"{search_db}-{pdb_id}.ligand-3d-candidates.parquet"
            )
            if ligand_3d_candidates is not None:
                candidate_table = pa.Table.from_pylist(
                    ligand_3d_candidates,
                    schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA,
                )
                pq.write_table(
                    candidate_table,
                    candidate_temporary,
                    compression="zstd",
                )
                candidate_path.parent.mkdir(exist_ok=True, parents=True)
                candidate_install = candidate_path.with_suffix(
                    candidate_path.suffix + ".tmp"
                )
                shutil.copyfile(candidate_temporary, candidate_install)
                candidate_install.replace(candidate_path)
            install_path = score_df_path.with_suffix(score_df_path.suffix + ".tmp")
            shutil.copyfile(temporary, install_path)
            install_path.replace(score_df_path)
            temporary.unlink(missing_ok=True)
            candidate_temporary.unlink(missing_ok=True)
            if ligand_3d_candidates is None:
                candidate_path.unlink(missing_ok=True)
        except Exception as e:
            scratch = data_dir / "scratch" / "scores" / "aggregate_scores_failures"
            scratch.mkdir(exist_ok=True, parents=True)
            (scratch / f"{search_db}_{pdb_id}.txt").write_text(f"{repr(e)}: {e}")
            LOG.error(
                f"scoring: Error in aggregate_scores: {pdb_id} searching against {search_db}: {repr(e)}"
            )
            raise
        return score_df_path

    def repair_score_df_targets(
        self,
        data_dir: Path,
        pdb_id: str,
        *,
        affected_target_entries: set[str],
        search_db: str = "holo",
        scratch_dir: Path | None = None,
    ) -> Path:
        """Replace only score rows whose target entry was reannotated."""
        if search_db != "holo":
            raise ValueError(
                "targeted score repair currently supports holo scores only"
            )
        if pdb_id in affected_target_entries:
            raise ValueError("an affected query entry requires a complete rescore")
        score_path = self.db_dir / f"search_db={search_db}" / f"{pdb_id}.parquet"
        candidate_path = (
            self.scores_dir
            / "ligand_3d_candidates"
            / f"search_db={search_db}"
            / f"shard={pdb_id[1:3]}"
            / f"{pdb_id}.parquet"
        )
        if not score_path.is_file() or not candidate_path.is_file():
            raise FileNotFoundError(
                f"targeted score repair requires existing score and candidate files: "
                f"score={score_path.is_file()} candidates={candidate_path.is_file()}"
            )

        requested_entries = {pdb_id, *affected_target_entries}
        requested_entries.difference_update(self.entries)
        if requested_entries:
            self.entries.update(
                load_entry_views(pdb_ids=requested_entries, data_dir=data_dir)
            )
        if pdb_id not in self.entries:
            raise KeyError(f"query entry is absent from the current index: {pdb_id}")
        target_system_ids = {
            system_id
            for target_pdb_id in affected_target_entries
            for system_id in (
                self.entries[target_pdb_id].systems
                if target_pdb_id in self.entries
                else {}
            )
        }
        source_to_aln_file = {
            f"{search_db}_{alignment_type}": data_dir
            / "alignments"
            / f"search_db={search_db}"
            / f"alignment_type={alignment_type}"
            / f"shard={pdb_id[1:3]}.parquet"
            for alignment_type in ["foldseek", "mmseqs"]
        }
        ligand_3d_candidates: list[_Ligand3DCandidateType] = []
        repaired = (
            self.aggregate_scores(
                pdb_id,
                search_db=search_db,
                data_dir=None,
                source_to_aln_file=source_to_aln_file,
                target_system_ids=target_system_ids,
                ligand_3d_candidates=ligand_3d_candidates,
            )
            if target_system_ids
            else None
        )

        def unaffected_target(values: pd.Series) -> pd.Series:
            return ~values.astype(str).str.split("__").str[0].isin(
                affected_target_entries
            )

        scores = pd.read_parquet(score_path)
        scores = scores[unaffected_target(scores["target_system"])]
        if repaired is not None and not repaired.empty:
            scores = pd.concat([scores, repaired], ignore_index=True)
        score_keys = [
            "query_system",
            "query_ligand_id",
            "target_system",
            "target_ligand_id",
            "metric",
        ]
        if scores.duplicated(score_keys).any():
            raise ValueError(
                f"targeted score repair produced duplicate rows for {pdb_id}"
            )
        scores = scores.sort_values(
            [column for column, _ in SORT_ORDER],
            ascending=[order == "ascending" for _, order in SORT_ORDER],
            ignore_index=True,
        )

        candidates = pd.read_parquet(candidate_path)
        candidates = candidates[
            ~candidates["target_entry"].astype(str).isin(affected_target_entries)
        ]
        if ligand_3d_candidates:
            candidates = pd.concat(
                [candidates, pd.DataFrame(ligand_3d_candidates)], ignore_index=True
            )
        candidate_keys = [
            "query_system",
            "query_ligand_id",
            "target_system",
            "target_ligand_id",
        ]
        if candidates.duplicated(candidate_keys).any():
            raise ValueError(
                f"targeted score repair produced duplicate candidates for {pdb_id}"
            )

        temporary_root = scratch_dir or score_path.parent
        temporary_root.mkdir(exist_ok=True, parents=True)
        score_temporary = temporary_root / f"{pdb_id}.repair-scores.parquet"
        candidate_temporary = temporary_root / f"{pdb_id}.repair-candidates.parquet"
        score_schema = schemas.PROTEIN_SIMILARITY_SCHEMA.with_metadata(
            {b"plinder.ligand_3d": b"deferred"}
        )
        scores.to_parquet(
            score_temporary,
            index=False,
            schema=score_schema,
        )
        candidates.to_parquet(
            candidate_temporary,
            index=False,
            schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA,
        )
        score_install = score_path.with_suffix(score_path.suffix + ".tmp")
        candidate_install = candidate_path.with_suffix(candidate_path.suffix + ".tmp")
        shutil.copyfile(score_temporary, score_install)
        shutil.copyfile(candidate_temporary, candidate_install)
        candidate_install.replace(candidate_path)
        score_install.replace(score_path)
        score_temporary.unlink(missing_ok=True)
        candidate_temporary.unlink(missing_ok=True)
        return score_path

    def load_alignments(
        self,
        source_to_aln_file: dict[str, Path],
        search_db: str = "holo",
        query_entry_ids: set[str] | None = None,
        target_entry_ids: set[str] | None = None,
    ) -> pd.DataFrame:
        data = []
        index_columns = [
            "query_entry",
            "target_entry",
            "query_chain_mapped",
            "target_chain_mapped",
            "source",
        ]
        for source, aln_file in source_to_aln_file.items():
            sdb, aln_type = source.split("_")
            if sdb != search_db:
                continue
            if not aln_file.exists():
                continue
            filters = []
            if query_entry_ids is not None:
                filters.append(("query_entry", "in", query_entry_ids))
            if target_entry_ids is not None:
                filters.append(("target_entry", "in", target_entry_ids))
            aln_df = pd.read_parquet(aln_file, filters=filters or None)
            if aln_df.empty:
                continue
            # Per-entry mapped files written by pandas restore these fields as
            # a MultiIndex. Release shards written by DuckDB expose the same
            # fields as regular columns. Normalize both layouts here so shard
            # reads can still use Parquet predicate pushdown.
            if set(index_columns).issubset(aln_df.columns):
                aln_df["source"] = aln_type
                aln_df = aln_df.set_index(index_columns)
            elif list(aln_df.index.names) != index_columns:
                raise ValueError(
                    f"unexpected mapped alignment schema in {aln_file}: "
                    f"index={aln_df.index.names}, columns={list(aln_df.columns)}"
                )
            aln_df["fident_qcov"] = aln_df["fident"] * aln_df["qcov"]
            aln_df["seqsim_qcov"] = aln_df["seqsim"] * aln_df["qcov"]
            if aln_type == "foldseek":
                aln_df["lddt_qcov"] = aln_df["lddt"] * aln_df["qcov"]
            # Legacy mapped files used four position-keyed nested columns.
            # Compact release shards already expose aligned residue arrays and
            # identity bytes, so they require no Python dictionary expansion.
            if "qrnum" in aln_df.columns:
                for column in ["qrnum", "trnum"]:
                    aln_df[column] = aln_df[column].apply(
                        lambda x: {int(i): int(r) for i, r in x}
                    )
                for column in ["qaa", "taa"]:
                    aln_df[column] = [
                        {
                            position: str(residue)
                            for position, residue in zip(qrnum, residues)
                        }
                        for qrnum, residues in zip(aln_df["qrnum"], aln_df[column])
                    ]
            data.append(aln_df)
        if len(data):
            df = pd.concat(data)
            return df.sort_index()
        else:
            return pd.DataFrame()

    def map_alignment_df(
        self, df_file: Path, aln_type: str, search_db: str
    ) -> pd.DataFrame:
        df = pd.read_parquet(df_file)
        if df.empty:
            for column in [
                "query_entry",
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
            ]:
                df[column] = pd.Series(dtype="string")
            for column in ["seqsim", "seqsim_qcov", "fident_qcov"]:
                df[column] = pd.Series(dtype="float64")
            if aln_type == "foldseek":
                df["lddt_qcov"] = pd.Series(dtype="float64")
            for column in [
                "query_pocket_residue_numbers",
                "target_pocket_residue_numbers",
                "pocket_residue_identity",
            ]:
                df[column] = pd.Series(dtype="object")
            df.drop(
                columns=["qaln", "taln", "evalue", "bits", "tcov"],
                inplace=True,
                errors="ignore",
            )
            df["source"] = pd.Series(dtype="string")
            return df.set_index(
                [
                    "query_entry",
                    "target_entry",
                    "query_chain_mapped",
                    "target_chain_mapped",
                    "source",
                ]
            )
        if aln_type == "foldseek":
            query_replacements = [
                "_xyz-enrich.cif.gz",
                "_xyz-enrich.cif",
                "_xyz-enrich",
                "pdb_0000",
                ".cif.gz",
                ".cif",
            ]
            for value in query_replacements:
                df["query"] = df["query"].str.replace(value, "", regex=False)
            if search_db == "pred":
                for value in ["-F1-model_v4.cif", "-F1-model_v4", "AF-"]:
                    df["target"] = df["target"].str.replace(value, "", regex=False)
            else:
                for value in query_replacements:
                    df["target"] = df["target"].str.replace(value, "", regex=False)
        query_parts = df["query"].str.split("_", n=1, expand=True)
        if aln_type == "foldseek":
            query_parts[1] = query_parts[1].str.replace(r"^MODEL_\d+_", "", regex=True)
        df["query_entry"] = query_parts[0]
        df["query_chain_mapped"] = [
            self.entries[entry].author_to_asym.get(author)
            for entry, author in zip(query_parts[0], query_parts[1])
        ]
        df = df.dropna(subset=["query_chain_mapped"]).reset_index(drop=True)
        if search_db == "pred":
            df["target_chain_mapped"] = "A"
            df["target_entry"] = df["target"].str.split("_", n=1).str[0]
        else:
            target_parts = df["target"].str.split("_", n=1, expand=True)
            if aln_type == "foldseek":
                target_parts[1] = target_parts[1].str.replace(
                    r"^MODEL_\d+_", "", regex=True
                )
            df["target_entry"] = target_parts[0]
            df["target_chain_mapped"] = [
                (
                    self.entries[entry].author_to_asym.get(author)
                    if entry in self.entries
                    else None
                )
                for entry, author in zip(target_parts[0], target_parts[1])
            ]
            df = df.dropna(subset=["target_chain_mapped"]).reset_index(drop=True)
        df["qaln"] = df["qaln"].str.upper()
        df["taln"] = df["taln"].str.upper()
        df["seqsim"] = [
            get_sequence_similarity_helper(qaln, taln)
            for qaln, taln in zip(df["qaln"], df["taln"])
        ]
        df["seqsim_qcov"] = df["seqsim"] * df["qcov"]
        df["fident_qcov"] = df["fident"] * df["qcov"]
        if aln_type == "foldseek":
            df["lddt_qcov"] = df["lddt"] * df["qcov"]
        pocket_mappings = [
            self._map_alignment_pocket_positions(
                query_entry=row.query_entry,
                target_entry=row.target_entry,
                query_chain=str(row.query_chain_mapped),
                target_chain=str(row.target_chain_mapped),
                qstart=int(row.qstart),
                tstart=int(row.tstart),
                qaln=str(row.qaln),
                taln=str(row.taln),
                aln_type=aln_type,
                search_db=search_db,
            )
            for row in df.itertuples(index=False)
        ]
        for index, column in enumerate(
            (
                "query_pocket_residue_numbers",
                "target_pocket_residue_numbers",
                "pocket_residue_identity",
            )
        ):
            df[column] = [mapping[index] for mapping in pocket_mappings]
        # Only ligand-level pocket score reconstruction is required downstream:
        # retain aligned residue numbers and equality flags, not alignment
        # positions or amino-acid letters. Raw search statistics remain in the
        # private search output and are omitted from the release representation.
        df.drop(
            columns=["qaln", "taln", "evalue", "bits", "tcov"],
            inplace=True,
            errors="ignore",
        )
        df["source"] = aln_type
        df.set_index(
            [
                "query_entry",
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
                "source",
            ],
            inplace=True,
        )
        return df

    def _map_alignment_pocket_positions(
        self,
        *,
        query_entry: str,
        target_entry: str,
        query_chain: str,
        target_chain: str,
        qstart: int,
        tstart: int,
        qaln: str,
        taln: str,
        aln_type: str,
        search_db: str,
    ) -> tuple[list[int], list[int], bytes]:
        """Map sparse query-pocket positions through one pairwise alignment."""
        # mmseqs operates on the SEQRES FASTA, so 1-based position
        #     equals the residue NUMBER (label_seq_id = Chain.residues key);
        # foldseek operates on the 3D structure, so position == 0-based
        #     resolved-residue INDEX;
        query_numbers: list[int] = []
        target_numbers: list[int] = []
        residue_identity = bytearray()
        q_i, t_i = qstart - 1, tstart - 1
        q_entry = self.entries[query_entry]
        q_i2n = q_entry.pocket_index_to_number_per_chain.get(query_chain, {})
        t_i2n: dict[int, int] = {}
        if search_db != "pred":
            t_entry = self.entries[target_entry]
            t_i2n = t_entry.pocket_index_to_number_per_chain.get(target_chain, {})
        alignment_length = min(len(qaln), len(taln))
        if not q_i2n or alignment_length == 0:
            return query_numbers, target_numbers, bytes(residue_identity)

        qaln = qaln[:alignment_length]
        taln = taln[:alignment_length]
        query_residue_positions = [
            position for position, residue in enumerate(qaln) if residue != "-"
        ]
        target_residue_positions = [
            position for position, residue in enumerate(taln) if residue != "-"
        ]

        # Locate only sparse query-pocket positions in the alignment instead of
        # walking every aligned character in Python.  Offsets are relative to
        # the first aligned residue reported by the search backend.
        if aln_type == "mmseqs":
            query_candidates = [
                (number - (q_i + 1), number) for number in set(q_i2n.values())
            ]
            target_pocket_numbers = set(t_i2n.values())
        else:
            query_candidates = [
                (index - q_i, number) for index, number in q_i2n.items()
            ]
            target_pocket_numbers = set()

        for query_offset, query_number in sorted(query_candidates):
            if not 0 <= query_offset < len(query_residue_positions):
                continue
            alignment_position = query_residue_positions[query_offset]
            if taln[alignment_position] == "-":
                continue
            target_index = t_i + bisect_left(
                target_residue_positions, alignment_position
            )
            if aln_type == "mmseqs":
                target_number = (
                    target_index + 1
                    if target_index + 1 in target_pocket_numbers
                    else None
                )
            else:
                target_number = t_i2n.get(target_index) if search_db != "pred" else None
            query_numbers.append(query_number)
            target_numbers.append(target_number if target_number is not None else -1)
            residue_identity.append(
                qaln[alignment_position] == taln[alignment_position]
            )
        return query_numbers, target_numbers, bytes(residue_identity)

    def map_row(self, parts: pd.Series, aln_type: str, search_db: str) -> pd.Series:
        """Map pocket positions for callers operating on one pandas row."""
        mapped = self._map_alignment_pocket_positions(
            query_entry=str(parts["query_entry"]),
            target_entry=str(parts["target_entry"]),
            query_chain=str(parts["query_chain_mapped"]),
            target_chain=str(parts["target_chain_mapped"]),
            qstart=int(parts["qstart"]),
            tstart=int(parts["tstart"]),
            qaln=str(parts["qaln"]),
            taln=str(parts["taln"]),
            aln_type=aln_type,
            search_db=search_db,
        )
        for column, values in zip(
            (
                "query_pocket_residue_numbers",
                "target_pocket_residue_numbers",
                "pocket_residue_identity",
            ),
            mapped,
        ):
            parts[column] = values
        return parts

    def get_protein_scores_pair(self, aln: pd.DataFrame) -> dict[str, float]:
        scores = {}
        for source in aln.index.unique():
            source_rows = aln.loc[source]
            if isinstance(source_rows, pd.DataFrame):
                rank_column = "lddt_qcov" if source == "foldseek" else "fident_qcov"
                data = source_rows.sort_values(
                    [rank_column, "qcov", "fident"],
                    ascending=False,
                    kind="stable",
                ).iloc[0]
            else:
                data = source_rows
            for score in ["qcov", "fident", "seqsim", "fident_qcov", "seqsim_qcov"]:
                scores[f"protein_{score}_{source}"] = data[score]
            if source == "foldseek":
                for score in ["lddt", "lddt_qcov"]:
                    scores[f"protein_{score}_{source}"] = data[score]
        return scores

    def get_protein_scores(
        self,
        query_target_entry_alignments: pd.DataFrame,
        query_system: SystemView,
        target_protein_chains: list[str],
        query_system_length: int,
        query_protein_chains: list[str] | None = None,
    ) -> tuple[
        dict[str, list[_ChainPairType]],
        _SimilarityScoreDictType,
        dict[_ChainPairType, pd.DataFrame],
        str,
    ]:
        scores: _SimilarityScoreDictType = defaultdict(float)
        mappings: dict[str, list[_ChainPairType]] = defaultdict(list)
        max_chain_lengths: dict[str, float] = defaultdict(float)
        protein_chain_mapper = ""
        if query_protein_chains is None:
            query_protein_chains = query_system.protein_chains_asym_id
        query_protein_chains = self.get_protein_receptor_chains(
            query_system.pdb_id, query_protein_chains
        )
        query_system_length = self.get_protein_chain_length(
            query_system.pdb_id, query_protein_chains
        )
        s_matrix = np.zeros(
            (
                len(query_protein_chains),
                len(target_protein_chains),
            )
        )
        for i, q_instance_chain in enumerate(query_protein_chains):
            q_chain = q_instance_chain.split(".")[1]
            q_chain_length = self.entries[query_system.pdb_id].chains[q_chain].length
            for j, t_instance_chain in enumerate(target_protein_chains):
                t_chain = t_instance_chain.split(".")[1]
                try:
                    aln = query_target_entry_alignments.loc[(q_chain, t_chain)]
                except KeyError:
                    LOG.debug("get_protein_scores indexing failed, aln is None")
                    aln = None
                if aln is not None:
                    pair = (q_instance_chain, t_instance_chain)
                    pair_scores = self.get_protein_scores_pair(aln)
                    if protein_chain_mapper == "":
                        for chain_mapper in self.protein_chain_mappers:
                            if chain_mapper in pair_scores:
                                protein_chain_mapper = chain_mapper
                                break
                    s_matrix[i, j] = pair_scores.get(protein_chain_mapper, 0)
                    for score in pair_scores:
                        if (
                            pair_scores[score] * q_chain_length
                            > scores[f"{score}_weighted_max"]
                        ):
                            scores[f"{score}_weighted_max"] = (
                                pair_scores[score] * q_chain_length
                            )
                            max_chain_lengths[f"{score}_weighted_max"] = q_chain_length
                            mappings[f"{score}_weighted_max"] = [pair]
                        if pair_scores[score] > scores[f"{score}_max"]:
                            scores[f"{score}_max"] = pair_scores[score]
                            mappings[f"{score}_max"] = [pair]
        for score in scores:
            if score.endswith("_weighted_max") and max_chain_lengths[score] > 0:
                scores[score] /= max_chain_lengths[score]

        # do chain mapping
        # Find the highest scores along with their indices in the matrix
        alns = {}
        while np.any(s_matrix):
            # changed score -> score_tag conflicting variable name
            score_val = float(np.amax(s_matrix))
            if score_val == 0:
                break
            q_idx, t_idx = np.unravel_index(np.argmax(s_matrix), s_matrix.shape)
            q_instance_chain, t_instance_chain = (
                query_protein_chains[q_idx],
                target_protein_chains[t_idx],
            )
            q_chain, t_chain = (
                q_instance_chain.split(".")[1],
                t_instance_chain.split(".")[1],
            )
            try:
                aln = query_target_entry_alignments.loc[(q_chain, t_chain)]
            except KeyError:
                break
            pair = (q_instance_chain, t_instance_chain)
            pair_scores = self.get_protein_scores_pair(aln)
            for score in pair_scores:
                scores[f"{score}_weighted_sum"] += (
                    pair_scores[score]
                    * self.entries[query_system.pdb_id].chains[q_chain].length
                )
            alns[pair] = aln

            mappings[f"{protein_chain_mapper}_weighted_sum"].append(pair)
            # Zero out the entire row and column for the max element
            s_matrix[q_idx, :] = 0
            s_matrix[:, t_idx] = 0
        for score in scores:
            if score.endswith("_weighted_sum") and query_system_length > 0:
                scores[score] /= query_system_length
        return mappings, scores, alns, protein_chain_mapper

    def get_pocket_pli_scores(
        self,
        alns: dict[_ChainPairType, pd.DataFrame],
        query_system: SystemView,
        target_system: SystemView | None = None,
    ) -> tuple[
        _SimilarityScoreDictType,
        _SimilarityScoreDictType,
    ]:
        (
            query_pocket,
            query_interactions,
            pocket_length,
            pli_length,
            unique_length,
        ) = self._protein_only_pocket_data(
            query_system.pdb_id,
            query_system.protein_chains_asym_id,
            query_system.pocket_residue_number_to_index,
            query_system.interactions_counter,
        )
        target_pocket = None
        target_interactions = None
        if target_system is not None:
            (
                target_pocket,
                target_interactions,
                _,
                _,
                _,
            ) = self._protein_only_pocket_data(
                target_system.pdb_id,
                target_system.protein_chains_asym_id,
                target_system.pocket_residue_number_to_index,
                target_system.interactions_counter,
            )
        return self._get_pocket_pli_scores(
            alns=alns,
            query_pocket=query_pocket,
            query_interactions=query_interactions,
            pocket_length=pocket_length,
            pli_length=pli_length,
            pli_unique_length=unique_length,
            target_pocket=target_pocket,
            target_interactions=target_interactions,
        )

    def _protein_only_pocket_data(
        self,
        pdb_id: str,
        receptor_chains: list[str],
        pocket: dict[str, dict[int, int]],
        interactions: dict[str, dict[int, Counter[str]]],
    ) -> tuple[
        dict[str, dict[int, int]],
        dict[str, dict[int, Counter[str]]],
        int,
        int,
        int,
    ]:
        """Restrict pocket and interaction data to polypeptide receptors."""
        allowed = set(self.get_protein_receptor_chains(pdb_id, receptor_chains))
        protein_pocket = {
            chain: residues for chain, residues in pocket.items() if chain in allowed
        }
        protein_interactions = {
            chain: residues
            for chain, residues in interactions.items()
            if chain in allowed
        }
        pocket_length = sum(len(residues) for residues in protein_pocket.values())
        interaction_length = sum(
            sum(counter.values())
            for residues in protein_interactions.values()
            for counter in residues.values()
        )
        unique_interaction_length = sum(
            len(counter)
            for residues in protein_interactions.values()
            for counter in residues.values()
        )
        return (
            protein_pocket,
            protein_interactions,
            pocket_length,
            interaction_length,
            unique_interaction_length,
        )

    def _get_pocket_pli_scores(
        self,
        *,
        alns: dict[_ChainPairType, pd.DataFrame],
        query_pocket: dict[str, dict[int, int]],
        query_interactions: dict[str, dict[int, Counter[str]]],
        pocket_length: int,
        pli_length: int,
        pli_unique_length: int,
        target_pocket: dict[str, dict[int, int]] | None,
        target_interactions: dict[str, dict[int, Counter[str]]] | None,
    ) -> tuple[
        _SimilarityScoreDictType,
        _SimilarityScoreDictType,
    ]:
        pocket_scores: _SimilarityScoreDictType = defaultdict(float)
        pli_scores: _SimilarityScoreDictType = defaultdict(float)
        has_target = target_pocket is not None
        target_pocket = target_pocket or {}
        target_interactions = target_interactions or {}
        # Compact maps contain only aligned query-pocket residue numbers,
        # target-pocket residue numbers (or -1), and amino-acid identity flags.
        # This is sufficient for exact reconstruction of every pocket metric.
        for q_instance_chain, t_instance_chain in alns:
            aln = alns[(q_instance_chain, t_instance_chain)]
            q_chain_pocket = query_pocket.get(q_instance_chain, {})
            q_chain_interactions = query_interactions.get(q_instance_chain, {})
            t_chain_pocket = target_pocket.get(t_instance_chain, {})
            t_chain_interactions = target_interactions.get(t_instance_chain, {})
            for source, aln_source in aln.iterrows():
                pocket_positions: abc.Iterable[tuple[int, int, bool]]
                if "query_pocket_residue_numbers" in aln_source.index:
                    compact_values = (
                        aln_source["query_pocket_residue_numbers"],
                        aln_source["target_pocket_residue_numbers"],
                        aln_source["pocket_residue_identity"],
                    )
                    # Pandas represents null list/binary Parquet cells as
                    # scalar NaN values. Such an alignment has no mapped
                    # pocket positions and contributes zero pocket coverage.
                    if not all(
                        isinstance(value, abc.Iterable) for value in compact_values
                    ):
                        continue
                    pocket_positions = zip(
                        *compact_values,
                    )
                else:
                    pocket_positions = (
                        (
                            query_number,
                            aln_source["trnum"].get(position, -1),
                            aln_source["qaa"][position] == aln_source["taa"][position],
                        )
                        for position, query_number in aln_source["qrnum"].items()
                    )
                for q_n, t_n, residues_are_identical in pocket_positions:
                    if q_n not in q_chain_pocket:
                        continue
                    if residues_are_identical:
                        pocket_scores[f"pocket_fident_{source}"] += 1
                    if has_target and t_n >= 0 and t_n in t_chain_pocket:
                        pocket_scores[f"pocket_qcov_{source}"] += 1
                        if residues_are_identical:
                            pocket_scores[f"pocket_fident_qcov_{source}"] += 1
                        if q_n in q_chain_interactions and t_n in t_chain_interactions:
                            pli_scores[f"pli_qcov_{source}"] += sum(
                                (
                                    q_chain_interactions[q_n]
                                    & t_chain_interactions[t_n]
                                ).values()
                            )
                            pli_scores[f"pli_unique_qcov_{source}"] += len(
                                set(q_chain_interactions[q_n].values())
                                & set(t_chain_interactions[t_n].values())
                            )
        if pocket_length:
            for score in pocket_scores:
                pocket_scores[score] /= pocket_length
        else:
            pocket_scores.clear()
        for score in list(pli_scores):
            denominator = pli_unique_length if "unique" in score else pli_length
            if denominator:
                pli_scores[score] /= denominator
            else:
                del pli_scores[score]
        return pocket_scores, pli_scores

    def get_ligand_pair_pocket_pli_scores(
        self,
        alns: dict[_ChainPairType, pd.DataFrame],
        query_ligand: LigandView,
        target_ligand: LigandView,
    ) -> tuple[
        _SimilarityScoreDictType,
        _SimilarityScoreDictType,
    ]:
        """Calculate directed pocket and PLI coverage for one ligand pair."""
        (
            query_pocket,
            query_interactions,
            pocket_length,
            pli_length,
            unique_length,
        ) = self._ligand_protein_only_pocket_data(query_ligand)
        (
            target_pocket,
            target_interactions,
            _,
            _,
            _,
        ) = self._ligand_protein_only_pocket_data(target_ligand)
        return self._get_pocket_pli_scores(
            alns=alns,
            query_pocket=query_pocket,
            query_interactions=query_interactions,
            pocket_length=pocket_length,
            pli_length=pli_length,
            pli_unique_length=unique_length,
            target_pocket=target_pocket,
            target_interactions=target_interactions,
        )

    def get_ligand_pocket_scores(
        self,
        alns: dict[_ChainPairType, pd.DataFrame],
        query_ligand: LigandView,
    ) -> _SimilarityScoreDictType:
        """Calculate directed pocket identity for a ligand against apo/pred."""
        (
            query_pocket,
            query_interactions,
            pocket_length,
            pli_length,
            unique_length,
        ) = self._ligand_protein_only_pocket_data(query_ligand)
        return self._get_pocket_pli_scores(
            alns=alns,
            query_pocket=query_pocket,
            query_interactions=query_interactions,
            pocket_length=pocket_length,
            pli_length=pli_length,
            pli_unique_length=unique_length,
            target_pocket=None,
            target_interactions=None,
        )[0]

    def _ligand_protein_only_pocket_data(self, ligand: LigandView) -> _PocketDataType:
        """Cache immutable protein-only pocket inputs per ligand occurrence."""
        if ligand.id not in self._ligand_pocket_data_cache:
            self._ligand_pocket_data_cache[ligand.id] = self._protein_only_pocket_data(
                ligand.pdb_id,
                ligand.protein_chains_asym_id,
                ligand.pocket_residue_number_to_index,
                ligand.interactions_counter,
            )
        return self._ligand_pocket_data_cache[ligand.id]

    def get_protein_receptor_chains(
        self, pdb_id: str, receptor_chains: list[str]
    ) -> list[str]:
        """Return only polypeptide instance chains from a receptor collection."""
        entry = self.entries.get(pdb_id)
        if entry is None:
            # LigandView instances constructed outside an EntryView are expected
            # to already carry the protein-only scoring representation.
            return sorted(set(receptor_chains))
        result = []
        for instance_chain in receptor_chains:
            asym_id = instance_chain.split(".", maxsplit=1)[-1]
            chain = entry.chains.get(asym_id)
            if chain is not None and chain.is_polypeptide:
                result.append(instance_chain)
        return sorted(set(result))

    def get_protein_chain_length(self, pdb_id: str, protein_chains: list[str]) -> int:
        """Return the total SEQRES length for an instance-chain collection."""
        return sum(
            self.entries[pdb_id].chains[instance_chain.split(".", 1)[1]].length
            for instance_chain in self.get_protein_receptor_chains(
                pdb_id, protein_chains
            )
        )

    def get_scores(
        self,
        search_db: str,
        query_system: SystemView,
        query_entry_alignments: pd.DataFrame,
        data_dir: Path | None = None,
        ligand_3d_candidates: list[_Ligand3DCandidateType] | None = None,
        query_ligand_ids: set[str] | None = None,
        target_system_ids: set[str] | None = None,
        target_ligand_ids: set[str] | None = None,
    ) -> abc.Generator[dict[str, str | float | None], None, None]:
        if search_db == "holo":
            return self.get_scores_holo(
                query_system,
                query_entry_alignments,
                data_dir=data_dir,
                ligand_3d_candidates=ligand_3d_candidates,
                query_ligand_ids=query_ligand_ids,
                target_system_ids=target_system_ids,
                target_ligand_ids=target_ligand_ids,
            )
        elif search_db == "apo" or search_db == "pred":
            return self.get_scores_apo_pred(
                query_system,
                query_entry_alignments,
                query_ligand_ids=query_ligand_ids,
                target_system_ids=target_system_ids,
            )
        else:
            raise ValueError(f"Invalid search_db: {search_db}")

    def get_scores_holo(
        self,
        query_system: SystemView,
        query_entry_alignments: pd.DataFrame,
        data_dir: Path | None = None,
        ligand_3d_candidates: list[_Ligand3DCandidateType] | None = None,
        query_ligand_ids: set[str] | None = None,
        target_system_ids: set[str] | None = None,
        target_ligand_ids: set[str] | None = None,
    ) -> abc.Generator[dict[str, str | float | None], None, None]:
        score_started = perf_counter()
        deferred_rows: list[
            tuple[
                dict[str, str | float | None],
                LigandView,
                LigandView,
                float,
            ]
        ] = []
        query_ligands = [
            ligand
            for ligand in query_system.ligands.values()
            if ligand.is_proper
            and (query_ligand_ids is None or ligand.id in query_ligand_ids)
        ]
        if not query_ligands:
            return
        query_instance_chains = sorted(
            {
                chain.split(".", 1)[1]
                for ligand in query_ligands
                for chain in self.get_protein_receptor_chains(
                    ligand.pdb_id, ligand.protein_chains_asym_id
                )
            }
        )
        if not query_instance_chains:
            return
        for target_entry in query_entry_alignments.index.get_level_values(
            "target_entry"
        ).unique():
            if target_entry not in self.entries:
                # No data for this target entry
                continue
            query_target_entry_alignments = query_entry_alignments.loc[target_entry]
            query_chain_mapped_values = set(
                query_target_entry_alignments.index.get_level_values(
                    "query_chain_mapped"
                )
            )
            if all(
                chain not in query_chain_mapped_values
                for chain in query_instance_chains
            ):
                # No alignments for this query system
                continue
            all_target_chains: set[str] = set()
            for q_chain in query_instance_chains:
                if q_chain in query_chain_mapped_values:
                    all_target_chains.update(
                        query_target_entry_alignments.loc[
                            q_chain
                        ].index.get_level_values("target_chain_mapped")
                    )
            for target_system_id in self.entries[target_entry].systems:
                target_system = self.entries[target_entry].systems[target_system_id]
                if not self.system_is_target_scoreable(target_system) or (
                    target_system_ids is not None
                    and target_system_id not in target_system_ids
                ):
                    continue
                target_system_protein_chains = self.get_protein_receptor_chains(
                    target_system.pdb_id, target_system.protein_chains_asym_id
                )
                if (
                    target_system_id == query_system.id
                    or not target_system_protein_chains
                    or all(
                        target_instance_chain.split(".")[1] not in all_target_chains
                        for target_instance_chain in target_system_protein_chains
                    )
                ):
                    # Same as query system or No alignments for this target system
                    continue
                target_ligands = [
                    ligand
                    for ligand in target_system.ligands.values()
                    if ligand.is_proper
                    and (target_ligand_ids is None or ligand.id in target_ligand_ids)
                ]
                for query_ligand in query_ligands:
                    # Keep every receptor chain in the directed query
                    # denominator. Chains without a hit remain zero-coverage
                    # rows in get_protein_scores(); dropping them here would
                    # inflate weighted similarities for partial alignments.
                    query_protein_chains = self.get_protein_receptor_chains(
                        query_ligand.pdb_id,
                        query_ligand.protein_chains_asym_id,
                    )
                    if not query_protein_chains:
                        continue
                    query_protein_length = self.get_protein_chain_length(
                        query_system.pdb_id, query_protein_chains
                    )
                    for target_ligand in target_ligands:
                        target_protein_chains = self.get_protein_receptor_chains(
                            target_ligand.pdb_id,
                            target_ligand.protein_chains_asym_id,
                        )
                        if not target_protein_chains:
                            continue
                        q_t_scores: dict[str, float] = {}

                        # Protein scores and mappings are restricted to the
                        # receptor chains belonging to this ligand pair. Many
                        # systems reuse the same receptor-chain sets, so avoid
                        # recomputing their assignment and aggregate scores.
                        protein_cache_key = (
                            query_system.pdb_id,
                            str(target_entry),
                            tuple(query_protein_chains),
                            tuple(target_protein_chains),
                        )
                        if protein_cache_key not in self._protein_score_cache:
                            self._protein_score_cache[
                                protein_cache_key
                            ] = self.get_protein_scores(
                                query_target_entry_alignments,
                                query_system,
                                target_protein_chains,
                                query_protein_length,
                                query_protein_chains=query_protein_chains,
                            )
                        (
                            q_t_mappings,
                            protein_scores,
                            alns,
                            protein_chain_mapper,
                        ) = self._protein_score_cache[protein_cache_key]
                        if not protein_scores:
                            continue
                        q_t_scores.update(protein_scores)

                        (
                            pocket_scores,
                            pli_scores,
                        ) = self.get_ligand_pair_pocket_pli_scores(
                            alns, query_ligand, target_ligand
                        )
                        q_t_scores.update(pocket_scores)
                        q_t_scores.update(pli_scores)
                        combined: dict[str, str | float | None] = {
                            **combine_scores(
                                q_t_scores, q_t_mappings, protein_chain_mapper
                            )
                        }

                        pocket_qcov_value = combined.get("pocket_qcov", 0.0)
                        pocket_qcov = (
                            float(pocket_qcov_value)
                            if isinstance(pocket_qcov_value, (int, float))
                            else 0.0
                        )
                        defer_shape = (
                            data_dir is not None
                            and pocket_qcov > 0
                            and self.shape_score_threads > 1
                        )
                        if data_dir is not None and pocket_qcov > 0 and not defer_shape:
                            combined.update(
                                self.get_ligand_pair_shape_scores(
                                    data_dir,
                                    query_ligand,
                                    target_ligand,
                                    pocket_qcov,
                                )
                            )
                        combined.update(
                            {
                                "query_system": query_system.id,
                                "query_ligand_id": query_ligand.id,
                                "target_system": target_system.id,
                                "target_ligand_id": target_ligand.id,
                            }
                        )
                        if (
                            ligand_3d_candidates is not None
                            and pocket_qcov > 0
                            and query_ligand.is_3d_score_able
                            and target_ligand.is_3d_score_able
                        ):
                            ligand_3d_candidates.append(
                                {
                                    "query_system": query_system.id,
                                    "query_ligand_id": query_ligand.id,
                                    "query_entry": query_ligand.pdb_id,
                                    "query_ligand_asym_id": query_ligand.asym_id,
                                    "target_system": target_system.id,
                                    "target_ligand_id": target_ligand.id,
                                    "target_entry": target_ligand.pdb_id,
                                    "target_ligand_asym_id": target_ligand.asym_id,
                                    "protein_mapping": combined.get("protein_mapping"),
                                    "protein_mapper": combined.get("protein_mapper"),
                                    "pocket_qcov": pocket_qcov,
                                }
                            )
                        if defer_shape:
                            deferred_rows.append(
                                (
                                    combined,
                                    query_ligand,
                                    target_ligand,
                                    pocket_qcov,
                                )
                            )
                        else:
                            yield combined

        if not deferred_rows:
            return
        assert data_dir is not None
        shape_data_dir = data_dir

        canonical_pairs: dict[
            tuple[tuple[str, str], tuple[str, str]], tuple[LigandView, LigandView]
        ] = {}
        canonical_ligands: dict[tuple[str, str], LigandView] = {}
        canonical_query_ligands: dict[tuple[str, str], LigandView] = {}
        for _, query_ligand, target_ligand, _ in deferred_rows:
            query_key = (query_ligand.pdb_id, query_ligand.asym_id)
            target_key = (target_ligand.pdb_id, target_ligand.asym_id)
            canonical_pairs.setdefault(
                (query_key, target_key), (query_ligand, target_ligand)
            )
            canonical_query_ligands.setdefault(query_key, query_ligand)
            canonical_ligands.setdefault(query_key, query_ligand)
            canonical_ligands.setdefault(target_key, target_ligand)

        core_seconds = perf_counter() - score_started

        def load_ligand(ligand: LigandView) -> Chem.Mol | None:
            return self._get_ligand_mol(shape_data_dir, ligand)

        def score_pair(pair: tuple[LigandView, LigandView]) -> None:
            self.get_ligand_pair_shape_scores(shape_data_dir, pair[0], pair[1], 1.0)

        # Preload each canonical molecule once, then calculate each unique
        # directed pair once. RDKit receives independent molecule clones in
        # get_ligand_pair_shape_scores(), so comparisons are thread-isolated.
        preload_started = perf_counter()
        self._preload_packed_ligand_sdfs(
            shape_data_dir, list(canonical_ligands.values())
        )
        ligands_to_load = [
            ligand
            for key, ligand in canonical_ligands.items()
            if key not in self._ligand_mol_cache
        ]
        with ThreadPoolExecutor(max_workers=self.shape_score_threads) as executor:
            list(executor.map(load_ligand, ligands_to_load))
            preload_seconds = perf_counter() - preload_started
            shape_started = perf_counter()
            # The reference conformer is not modified by ShapeAlign, so detect
            # its pharmacophore once instead of repeating that work per pair.
            for key, ligand in canonical_query_ligands.items():
                if key not in self._ligand_reference_feature_cache:
                    self._get_ligand_reference_features(shape_data_dir, ligand)
            pairs_to_score = [
                pair
                for key, pair in canonical_pairs.items()
                if key not in self._ligand_shape_score_cache
            ]
            list(executor.map(score_pair, pairs_to_score))
        shape_seconds = perf_counter() - shape_started
        LOG.info(
            "ligand scoring profile for %s: rows=%d canonical_pairs=%d new_pairs=%d "
            "core=%.2fs sdf_preload=%.2fs shape=%.2fs threads=%d",
            query_system.id,
            len(deferred_rows),
            len(canonical_pairs),
            len(pairs_to_score),
            core_seconds,
            preload_seconds,
            shape_seconds,
            self.shape_score_threads,
        )

        for combined, query_ligand, target_ligand, pocket_qcov in deferred_rows:
            combined.update(
                self.get_ligand_pair_shape_scores(
                    shape_data_dir, query_ligand, target_ligand, pocket_qcov
                )
            )
            yield combined

    def _get_suffixes(self, s: str) -> list[str]:
        suffixes = ["_weighted_sum", "_weighted_max", "_max"]
        if not s.startswith("protein"):
            suffixes = [""]
        return suffixes

    def get_column_mapr(self) -> dict[str, list[str]]:
        mapr = {"info": [x for x in INFO_COLUMNS]}
        for s in SCORE_NAMES:
            suffixes = self._get_suffixes(s)
            for suffix in suffixes:
                mapr[f"{s}{suffix}"] = [
                    f"{s}{suffix}",
                    f"{s}{suffix}_source",
                    f"{s}{suffix}_mapping",
                ]
        return mapr

    def aggregate_scores(
        self,
        pdb_id: str,
        search_db: str = "holo",
        data_dir: Path | None = None,
        source_to_aln_file: dict[str, Path] | None = None,
        query_system_ids: set[str] | None = None,
        query_ligand_ids: set[str] | None = None,
        target_system_ids: set[str] | None = None,
        target_ligand_ids: set[str] | None = None,
        ligand_3d_candidates: list[_Ligand3DCandidateType] | None = None,
    ) -> Optional[pd.DataFrame]:
        if source_to_aln_file is None:
            source_to_aln_file = {
                f"{search_db}_{aln_type}": self.db_dir
                / f"{search_db}_{aln_type}"
                / "mapped_aln"
                / f"{pdb_id}.parquet"
                for aln_type in ["foldseek", "mmseqs"]
            }
        target_entry_ids = None
        if target_system_ids is not None:
            target_entry_ids = {
                system_id.split("__", maxsplit=1)[0] for system_id in target_system_ids
            }
        alignments = self.load_alignments(
            search_db=search_db,
            source_to_aln_file=source_to_aln_file,
            query_entry_ids={pdb_id},
            target_entry_ids=target_entry_ids,
        )
        if alignments.empty:
            return None
        query_entry_alignments = alignments.loc[pdb_id]
        column_mapr = self.get_column_mapr()
        pdb_vals = []
        for system in self.entries[pdb_id].systems.values():
            if not self.system_is_query_scoreable(system) or (
                query_system_ids is not None and system.id not in query_system_ids
            ):
                continue
            for score_dict in self.get_scores(
                search_db,
                system,
                query_entry_alignments,
                data_dir=data_dir,
                ligand_3d_candidates=ligand_3d_candidates,
                query_ligand_ids=query_ligand_ids,
                target_system_ids=target_system_ids,
                target_ligand_ids=target_ligand_ids,
            ):
                # Keep nullable identifiers (notably target_ligand_id for
                # apo/pred) present so every search database shares a schema.
                info = {column: score_dict.get(column) for column in INFO_COLUMNS}
                grouped: dict[str, dict[str, str | float]] = {}
                # iterate over groups of columns instead of raw list of columns
                for metric, columns in column_mapr.items():
                    for column in columns:
                        score = score_dict.get(column)
                        if score is None:
                            continue
                        if metric == "info":
                            continue
                        else:
                            grouped.setdefault(metric, {})
                            abbr_column = column.replace(metric, "", 1).lstrip("_")
                            if not abbr_column:
                                abbr_column = "similarity"
                                if np.isnan(score):
                                    continue
                            grouped[metric][abbr_column] = score
                for metric, values in grouped.items():
                    if len(values) and values.get("similarity") is not None:
                        pdb_vals.append({**info, **values, **{"metric": metric}})
        if not len(pdb_vals):
            return None
        df = pd.DataFrame(pdb_vals)
        all_thresholds = [
            (m, self.minimum_thresholds.get(m, self.minimum_threshold))
            for m in df["metric"].unique()
        ]
        query = " or ".join(
            [f"(metric=='{m}' and similarity>={t})" for (m, t) in all_thresholds]
        )
        df = df.query(query).copy()
        df["protein_mapper"] = df["protein_mapper"].astype("category")
        for col in ["source", "metric"]:
            # PROTEIN_SIMILARITY_SCHEMA declares these dictionaries ordered.
            df[col] = df[col].astype(pd.CategoricalDtype(ordered=True))
        df["similarity"] = (df["similarity"] * 100).apply(round).astype(np.int8)
        columns = [col for (col, _) in SORT_ORDER]
        ascending = [direction == "ascending" for (_, direction) in SORT_ORDER]
        return df.sort_values(by=columns, ascending=ascending)

    def get_scores_apo_pred(
        self,
        query_system: SystemView,
        query_entry_alignments: pd.DataFrame,
        query_ligand_ids: set[str] | None = None,
        target_system_ids: set[str] | None = None,
    ) -> abc.Generator[dict[str, str | float | None], None, None]:
        query_ligands = [
            ligand
            for ligand in query_system.ligands.values()
            if ligand.is_proper
            and (query_ligand_ids is None or ligand.id in query_ligand_ids)
        ]
        for target_entry in query_entry_alignments.index.get_level_values(
            "target_entry"
        ).unique():
            query_target_entry_alignments = query_entry_alignments.loc[target_entry]
            query_chain_mapped_values = set(
                query_target_entry_alignments.index.get_level_values(
                    "query_chain_mapped"
                )
            )
            for query_ligand in query_ligands:
                query_protein_chains = self.get_protein_receptor_chains(
                    query_ligand.pdb_id,
                    query_ligand.protein_chains_asym_id,
                )
                if not query_protein_chains:
                    continue
                mapped_query_protein_chains = [
                    chain
                    for chain in query_protein_chains
                    if chain.split(".", 1)[1] in query_chain_mapped_values
                ]
                if not mapped_query_protein_chains:
                    continue
                query_protein_length = self.get_protein_chain_length(
                    query_system.pdb_id, query_protein_chains
                )
                target_chains: set[str] = set()
                for query_chain in mapped_query_protein_chains:
                    q_chain = query_chain.split(".", 1)[1]
                    target_chains.update(
                        query_target_entry_alignments.loc[
                            q_chain
                        ].index.get_level_values("target_chain_mapped")
                    )
                for t_chain in target_chains:
                    target_system_id = f"{target_entry}_{t_chain}"
                    if (
                        target_system_ids is not None
                        and target_system_id not in target_system_ids
                    ):
                        continue
                    if (
                        target_entry in self.entries
                        and t_chain in self.entries[target_entry].chains
                    ):
                        continue
                    q_t_scores: dict[str, float] = {}
                    # Protein score calculation
                    (
                        q_t_mappings,
                        protein_scores,
                        alns,
                        protein_chain_mapper,
                    ) = self.get_protein_scores(
                        query_target_entry_alignments,
                        query_system,
                        [f"0.{t_chain}"],
                        query_protein_length,
                        query_protein_chains=query_protein_chains,
                    )

                    if len(alns) == 0 or not len(protein_scores):
                        continue
                    q_t_scores.update(protein_scores)
                    # Pocket score calculation
                    q_t_scores.update(self.get_ligand_pocket_scores(alns, query_ligand))
                    q_t_scores_combined: dict[str, str | float | None] = {
                        **combine_scores(q_t_scores, q_t_mappings, protein_chain_mapper)
                    }
                    q_t_scores_combined["target_system"] = target_system_id
                    q_t_scores_combined["target_ligand_id"] = None
                    q_t_scores_combined["query_system"] = query_system.id
                    q_t_scores_combined["query_ligand_id"] = query_ligand.id
                    yield q_t_scores_combined
