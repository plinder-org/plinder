# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import shutil
import subprocess
from collections import Counter, abc, defaultdict
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
from typing import Any, Callable, Optional, cast

import biotite.sequence as seq
import biotite.sequence.align as align
import numpy as np
import pandas as pd
import pyarrow
import pyarrow as pa
import pyarrow.parquet as pq
from pyarrow import csv
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import ChemicalFeatures, rdMolAlign, rdShapeAlign, rdShapeHelpers
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
    ligands = annotation.loc[
        annotation["system_type"].eq("holo"), list(columns)
    ].rename(columns=columns)
    return (
        ligands.dropna(subset=["ligand_id"])
        .drop_duplicates(subset=["ligand_id"])
        .reset_index(drop=True)
        .sort_values("ligand_id")
    )


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
    for binary_fingerprint in result["fingerprint"]:
        fingerprint = DataStructs.CreateFromBinaryText(binary_fingerprint)
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprint, cofactor_fingerprints
        )
        best_index = int(np.argmax(similarities))
        maximum_similarities.append(float(similarities[best_index] * 100.0))
        closest_cofactors.append(cofactor_codes[best_index])
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
    ligands = (
        pd.read_parquet(data_dir / "ligands").drop_duplicates().reset_index(drop=True)
    )
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
    fingerprints = [
        mol2morgan_fp(smiles, radius=ECFP4_RADIUS, nbits=ECFP4_NBITS)
        for smiles in ligands_unique[smiles_column]
    ]
    ligands_unique["fingerprint"] = [
        DataStructs.BitVectToBinaryText(fingerprint) for fingerprint in fingerprints
    ]

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
    ligand_occurrences = pd.read_parquet(
        data_dir / "ligands",
        columns=["pdb_id", "ligand_rdkit_canonical_smiles"],
    ).merge(
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
    crippen_o3a = rdMolAlign.GetCrippenO3A(mobile, reference, maxIters=max_preiters)
    crippen_o3a.Align()
    shape, color = rdShapeAlign.AlignMol(
        reference,
        mobile,
        max_preiters=max_preiters,
        max_postiters=max_postiters,
    )
    return float(shape), float(color)


def get_feature_map_score(
    mol_1: Chem.Mol,
    mol_2: Chem.Mol,
    score_mode: int = FeatMaps.FeatMapScoreMode.All,
) -> float:
    """Calculate the normalized pharmacophore feature overlap."""
    factory, parameters = _get_feature_scoring_context()
    feature_lists = []
    for molecule in (mol_1, mol_2):
        feature_lists.append(
            [
                feature
                for feature in factory.GetFeaturesForMol(molecule)
                if feature.GetFamily() in PHARMACOPHORE_FEATURES
            ]
        )
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
) -> float:
    """Calculate SuCOS from feature overlap and shape protrusion distance."""
    feature_map_score = float(
        np.clip(get_feature_map_score(mol_1, mol_2, score_mode), 0, 1)
    )
    protrude_distance = float(
        np.clip(
            rdShapeHelpers.ShapeProtrudeDist(mol_1, mol_2, allowReordering=False),
            0,
            1,
        )
    )
    return 0.5 * feature_map_score + 0.5 * (1 - protrude_distance)


def load_sdf_molecule(sdf_file: Path) -> Chem.Mol | None:
    """Load and sanitize one SDF molecule, returning ``None`` on failure."""
    try:
        molecule = Chem.MolFromMolFile(str(sdf_file))
        if molecule is None:
            molecule = Chem.MolFromMolFile(
                str(sdf_file), sanitize=False, strictParsing=False
            )
            if molecule is not None:
                Chem.SanitizeMol(molecule)
        if molecule is None or molecule.GetNumConformers() == 0:
            return None
        return molecule
    except Exception as exc:
        LOG.warning(f"failed loading ligand SDF {sdf_file}: {exc}")
        return None


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
    reference = Chem.Mol(molecule)
    mobile = Chem.Mol(molecule)
    try:
        shape, color = align_molecules(reference, mobile)
        sucos = get_sucos_score(reference, mobile)
    except Exception as exc:
        LOG.warning(f"ligand SDF is not 3D-scoreable ({sdf_file}): {exc}")
        return False
    return all(np.isfinite(value) for value in (shape, color, sucos))


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
    # biotite's ProteinSequence handles B/Z/X natively; J/U/O still need mapping.
    non_canonical_aa = str.maketrans({"J": "L", "U": "C", "O": "K"})
    s1 = seq_str1.translate(non_canonical_aa).replace("-", "")
    s2 = seq_str2.translate(non_canonical_aa).replace("-", "")
    seq1 = seq.ProteinSequence(s1)
    seq2 = seq.ProteinSequence(s2)
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    trace = align.Alignment.trace_from_strings([s1, s2])
    ali = align.Alignment([seq1, seq2], trace)

    # "Similar" position = BLOSUM62 pair score > 20% of max self-score.
    # TODO: verify the definition for the normalization?
    codes = align.alignment.get_codes(ali)
    score_matrix = matrix.score_matrix()
    mask = (codes[0] != -1) & (codes[1] != -1)
    ci, cj = codes[0, mask], codes[1, mask]
    pair_scores = np.maximum(0, score_matrix[ci, cj])
    norm = np.maximum(score_matrix[ci, ci], score_matrix[cj, cj])
    similarity_score = float(np.mean(pair_scores / norm > 0.2))
    return (align.get_sequence_identity(ali), similarity_score)


def get_sequence_similarity_helper(seq_str1: str, seq_str2: str) -> float:
    try:
        return get_sequence_similarity(seq_str1, seq_str2)[1]
    except Exception:
        return 0


def run_alignment(
    aln_type: str,
    query_db: Path,
    target_db: Path,
    search_db: Path,
    aln_file: Path,
    alignment_config: FoldseekConfig | MMSeqsConfig,
    tmp_dir: Path = Path.cwd() / "tmp",
    remove_tmp: bool = True,
) -> None:
    if search_db.with_suffix(".dbtype").exists():
        search_db.with_suffix(".dbtype").unlink()

    search_commands = [
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
    if aln_type == "foldseek":
        search_commands += ["--sort-by-structure-bits", "0"]
    convert_commands = [
        aln_type,
        "convertalis",
        str(query_db),
        str(target_db),
        str(search_db),
        str(aln_file.with_suffix(".tsv")),
        "--format-mode",
        "4",
        "--format-output",
        "query,target,qlen,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln",
    ]
    if aln_type == "foldseek":
        convert_commands[-1] += ",lddt"
    subprocess.check_call(search_commands, stdout=subprocess.DEVNULL)
    subprocess.check_call(convert_commands, stdout=subprocess.DEVNULL)

    table = csv.read_csv(
        aln_file.with_suffix(".tsv"), parse_options=csv.ParseOptions(delimiter="\t")
    )
    table = table.append_column(
        "query_pdb_id",
        pyarrow.array(
            [
                str(x)
                .replace("_xyz-enrich.cif.gz", "")
                .replace("_xyz-enrich.cif", "")
                .replace("_xyz-enrich", "")
                .replace(".cif.gz", "")
                .replace(".cif", "")
                .replace("pdb_0000", "")[:4]
                for x in table["query"]
            ]
        ),
    )
    if search_db != "pred":
        table = table.append_column(
            "target_pdb_id",
            pyarrow.array(
                [
                    str(x)
                    .replace("_xyz-enrich.cif.gz", "")
                    .replace("_xyz-enrich.cif", "")
                    .replace("_xyz-enrich", "")
                    .replace(".cif.gz", "")
                    .replace(".cif", "")
                    .replace("pdb_0000", "")[:4]
                    for x in table["target"]
                ]
            ),
        )
    if aln_file.with_suffix(".parquet").exists():
        shutil.rmtree(aln_file.with_suffix(".parquet"))
    pq.write_to_dataset(
        table, aln_file.with_suffix(".parquet"), partition_cols=["query_pdb_id"]
    )

    # Cleanup
    if remove_tmp and tmp_dir.exists():
        shutil.rmtree(tmp_dir)
        aln_file.with_suffix(".tsv").unlink()
    for filename in search_db.parent.glob(f"{search_db.name}*"):
        if filename.is_dir():
            shutil.rmtree(filename)
        elif filename.is_file():
            filename.unlink()


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
    minimum_threshold: int = 0
    # optional custom minimum threshold for each metric BEFORE scaling to 100
    minimum_thresholds: dict[str, int] = field(
        default_factory=lambda: {"pli_qcov": 0, "pocket_qcov": 0, "pli_unique_qcov": 0}
    )
    _ligand_mol_cache: dict[tuple[str, str], Chem.Mol | None] = field(
        default_factory=dict, init=False, repr=False
    )

    def __post_init__(self) -> None:
        self.db_dir.mkdir(exist_ok=True, parents=True)
        self.scores_dir.mkdir(exist_ok=True, parents=True)

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
                self._ligand_mol_cache[cache_key] = load_sdf_molecule(sdf_file)
        return self._ligand_mol_cache[cache_key]

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

        query_mol = self._get_ligand_mol(data_dir, query_ligand)
        target_mol = self._get_ligand_mol(data_dir, target_ligand)
        if query_mol is None or target_mol is None:
            return {}

        # Both Crippen O3A and AlignMol mutate the mobile conformer. Always
        # clone cached molecules so one comparison cannot affect another.
        query_aligned = Chem.Mol(query_mol)
        target_aligned = Chem.Mol(target_mol)
        try:
            shape, color = align_molecules(query_aligned, target_aligned)
        except Exception as exc:
            LOG.warning(
                "shape alignment failed for "
                f"{query_ligand.id} to {target_ligand.id}: {exc}"
            )
            return {}
        if not all(np.isfinite(value) for value in (shape, color)):
            LOG.warning(
                "shape alignment returned non-finite scores for "
                f"{query_ligand.id} to {target_ligand.id}"
            )
            return {}

        scores: _SimilarityScoreDictType = {
            "shape": float(np.clip(shape, 0, 1)),
            "color": float(np.clip(color, 0, 1)),
        }
        try:
            sucos_shape = get_sucos_score(query_aligned, target_aligned)
        except Exception as exc:
            LOG.warning(
                "SuCOS calculation failed for "
                f"{query_ligand.id} to {target_ligand.id}: {exc}"
            )
            return scores
        if not np.isfinite(sucos_shape):
            LOG.warning(
                "SuCOS calculation returned a non-finite score for "
                f"{query_ligand.id} to {target_ligand.id}"
            )
            return scores
        scores["sucos_shape"] = sucos_shape
        scores["sucos_shape_pocket_qcov"] = sucos_shape * pocket_qcov
        return scores

    def make_dbs(self) -> None:
        databases.make_sub_dbs(self.db_dir, self.source_to_full_db_file, self.entries)

    @staticmethod
    def get_config(search_db: str, aln_type: str) -> FoldseekConfig | MMSeqsConfig:
        config: FoldseekConfig | MMSeqsConfig
        if aln_type == "foldseek":
            config = FoldseekConfig()
        else:
            config = MMSeqsConfig()
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
    ) -> None:
        output_folder.mkdir(exist_ok=True)
        for aln_type in ["mmseqs", "foldseek"]:
            sub_db = output_folder / search_db / aln_type
            sub_db.mkdir(exist_ok=True, parents=True)
            db_ids = databases.get_db_ids(
                self.entries, "holo", aln_type, entry_ids=entry_ids
            )
            databases.make_sub_db(
                db_ids,
                self.source_to_full_db_file[f"holo_{aln_type}"],
                sub_db,
                aln_type,
            )
            tmp_dir = sub_db / f"tmp_{search_db}_{aln_type}"
            tmp_dir.mkdir(exist_ok=True, parents=True)
            aln_file = sub_db / f"aln_{search_db}.tsv"
            LOG.info("run_alignments calling run_alignment:")
            LOG.info(f"    sub_db={sub_db}")
            LOG.info(f"    aln_file={aln_file.with_suffix('.tsv')}")
            LOG.info(f"    aln_type={aln_type}")
            try:
                run_alignment(
                    aln_type=aln_type,
                    query_db=sub_db / sub_db.name,
                    target_db=self.db_dir
                    / f"{search_db}_{aln_type}"
                    / f"{search_db}_{aln_type}",
                    search_db=sub_db / "search",
                    aln_file=aln_file.with_suffix(".tsv"),
                    tmp_dir=tmp_dir / output_folder.stem,
                    alignment_config=self.get_config(search_db, aln_type),
                )
            except Exception as e:
                scratch = (
                    output_folder / "scratch" / "scores" / "run_alignment_failures"
                )
                scratch.mkdir(exist_ok=True, parents=True)
                (scratch / f"{search_db}_{aln_type}.txt").write_text(f"{repr(e)}: {e}")
                LOG.error(f"scoring: Error for {search_db}_{aln_type}: {e}")
                continue
            aln_dir = self.db_dir / f"{search_db}_{aln_type}" / "aln"
            aln_dir.mkdir(exist_ok=True, parents=True)
            for pdb_id in tqdm(entry_ids):
                pdb_id_file = (
                    aln_file.with_suffix(".parquet") / f"query_pdb_id={pdb_id}"
                )
                if not pdb_id_file.exists():
                    scratch = (
                        output_folder / "scratch" / "scores" / "run_alignments_failures"
                    )
                    scratch.mkdir(exist_ok=True, parents=True)
                    (scratch / f"{search_db}_{aln_type}_{pdb_id}.txt").write_text("")
                    continue
                pdb_id_df = pd.read_parquet(pdb_id_file)
                if not pdb_id_df.empty:
                    pdb_id_df.to_parquet(aln_dir / f"{pdb_id}.parquet")
                else:
                    scratch = (
                        output_folder / "scratch" / "scores" / "run_alignments_empty"
                    )
                    scratch.mkdir(exist_ok=True, parents=True)
                    (scratch / f"{search_db}_{aln_type}_{pdb_id}.txt").write_text("")

    def get_score_df(
        self, data_dir: Path, pdb_id: str, search_db: str, overwrite: bool = True
    ) -> Path:
        """
        Convert aligmnent results to mapped alignment results. Then
        aggregate the mapped alignment results into the scores dataset.
        """
        score_df_path = self.db_dir / f"search_db={search_db}" / f"{pdb_id}.parquet"
        if overwrite or not score_df_path.exists():
            LOG.info(f"get_score_df: aggregating scores for {pdb_id} to {search_db}")
        else:
            LOG.info(f"get_score_df: skipping existing {score_df_path}")
            return score_df_path
        for aln_type in ["foldseek", "mmseqs"]:
            pdb_id_file = (
                self.db_dir / f"{search_db}_{aln_type}" / "aln" / f"{pdb_id}.parquet"
            )
            if not pdb_id_file.exists():
                LOG.info(f"get_score_df: pdb_id_file={pdb_id_file} does not exist")
                continue

            # self.entries = {}
            entries_to_load = {pdb_id}
            if search_db != "pred" and pdb_id_file.exists():
                entries_to_load |= set(
                    pd.read_parquet(pdb_id_file, columns=["target_pdb_id"])[
                        "target_pdb_id"
                    ]
                )
            entries_to_load = entries_to_load.difference(self.entries.keys())
            LOG.info(f"entries_to_load pdb_id={pdb_id} {len(entries_to_load)}")
            LOG.info(
                f"loading {len(entries_to_load)} (additional) entries for {pdb_id}"
            )
            if entries_to_load:
                self.entries.update(
                    load_entry_views(pdb_ids=entries_to_load, data_dir=data_dir)
                )
            pdb_file = (
                self.db_dir
                / f"{search_db}_{aln_type}"
                / "mapped_aln"
                / f"{pdb_id}.parquet"
            )
            pdb_file.parent.mkdir(exist_ok=True, parents=True)
            if overwrite or not pdb_file.exists():
                try:
                    LOG.info(
                        f"mapping aligment df for {pdb_id} to {search_db} for {aln_type}"
                    )
                    self.map_alignment_df(pdb_id_file, aln_type, search_db).to_parquet(
                        pdb_file, index=True
                    )
                except Exception as e:
                    scratch = (
                        Path(*pdb_id_file.parts[:3])
                        / "scratch"
                        / "scores"
                        / "map_alignment_df_failures"
                    )
                    scratch.mkdir(exist_ok=True, parents=True)
                    (scratch / f"{search_db}_{aln_type}_{pdb_id}.txt").write_text(
                        f"{repr(e)}: {e}"
                    )
                    LOG.error(
                        f"scoring: Error in map_alignment_df: {pdb_id} searching against {search_db} with {aln_type}: {repr(e)}"
                    )
                    continue
            else:
                LOG.info(f"skipping creating {pdb_file} because it already exists")

        try:
            score_df_path.parent.mkdir(exist_ok=True, parents=True)
            LOG.info(f"aggregating scores for {pdb_id} to {search_db}")
            df = self.aggregate_scores(pdb_id, search_db=search_db, data_dir=data_dir)
            if df is not None and not df.empty:
                df.to_parquet(score_df_path, index=False)
        except Exception as e:
            scratch = (
                Path(*score_df_path.parts[:3])
                / "scratch"
                / "scores"
                / "aggregate_scores_failures"
            )
            scratch.mkdir(exist_ok=True, parents=True)
            (scratch / f"{search_db}_{pdb_id}.txt").write_text(f"{repr(e)}: {e}")
            LOG.error(
                f"scoring: Error in aggregate_scores: {pdb_id} searching against {search_db}: {repr(e)}"
            )
        return score_df_path

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
            aln_df["qrnum"] = aln_df["qrnum"].apply(
                lambda x: dict([(int(i), int(r)) for i, r in x])
            )
            aln_df["trnum"] = aln_df["trnum"].apply(
                lambda x: dict([(int(i), int(r)) for i, r in x])
            )
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
        if aln_type == "foldseek":
            df["query"] = df["query"].replace(
                {
                    "_xyz-enrich.cif.gz": "",
                    "_xyz-enrich.cif": "",
                    "_xyz-enrich": "",
                    "pdb_0000": "",
                    ".cif.gz": "",
                    ".cif": "",
                },
                regex=True,
            )
            if search_db == "pred":
                df["target"] = df["target"].replace(
                    {"-F1-model_v4.cif": "", "-F1-model_v4": "", "AF-": ""}, regex=True
                )
            else:
                df["target"] = df["target"].replace(
                    {
                        "_xyz-enrich.cif.gz": "",
                        "_xyz-enrich.cif": "",
                        "_xyz-enrich": "",
                        "pdb_0000": "",
                        ".cif.gz": "",
                        ".cif": "",
                    },
                    regex=True,
                )
        df["query_chain_mapped"] = (
            df["query"]
            .str.split("_", expand=True)
            .apply(lambda x: self.entries[x[0]].author_to_asym.get(x[1], None), axis=1)
        )
        df = df.dropna(subset=["query_chain_mapped"]).reset_index(drop=True)
        df["query_entry"] = df["query"].str.split("_", expand=True)[0]
        if search_db == "pred":
            df["target_chain_mapped"] = "A"
        else:
            df["target_chain_mapped"] = (
                df["target"]
                .str.split("_", expand=True)
                .apply(
                    lambda x: (
                        self.entries[x[0]].author_to_asym.get(x[1], None)
                        if x[0] in self.entries
                        else None
                    ),
                    axis=1,
                )
            )
            df = df.dropna(subset=["target_chain_mapped"]).reset_index(drop=True)
        df["target_entry"] = df["target"].str.split("_", expand=True)[0]
        df["qaln"] = df["qaln"].str.upper()
        df["taln"] = df["taln"].str.upper()
        df["seqsim"] = df[["qaln", "taln"]].apply(
            lambda x: get_sequence_similarity_helper(x["qaln"], x["taln"]), axis=1
        )
        df["seqsim_qcov"] = df["seqsim"] * df["qcov"]
        df["fident_qcov"] = df["fident"] * df["qcov"]
        if aln_type == "foldseek":
            df["lddt_qcov"] = df["lddt"] * df["qcov"]
        df = df.apply(
            lambda x: self.map_row(x, aln_type=aln_type, search_db=search_db), axis=1
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

    def map_row(self, parts: pd.Series, aln_type: str, search_db: str) -> pd.Series:
        # mmseqs operates on the SEQRES FASTA, so 1-based position
        #     equals the residue NUMBER (label_seq_id = Chain.residues key);
        # foldseek operates on the 3D structure, so position == 0-based
        #     resolved-residue INDEX;
        parts["qrnum"] = []
        parts["trnum"] = []
        q_i, t_i = int(parts["qstart"]) - 1, int(parts["tstart"]) - 1
        q_i2n: dict[int, int] = {}
        t_i2n: dict[int, int] = {}
        if aln_type == "foldseek":
            q_entry = self.entries[parts["query_entry"]]
            q_i2n = q_entry.pocket_index_to_number_per_chain.get(
                parts["query_chain_mapped"], {}
            )
            if search_db != "pred":
                t_entry = self.entries[parts["target_entry"]]
                t_i2n = t_entry.pocket_index_to_number_per_chain.get(
                    parts["target_chain_mapped"], {}
                )
        for x, (q_a, t_a) in enumerate(zip(parts["qaln"], parts["taln"])):
            if q_a != "-" and t_a != "-":
                if aln_type == "mmseqs":
                    parts["qrnum"].append((x, q_i + 1))
                    if search_db != "pred":
                        parts["trnum"].append((x, t_i + 1))
                else:
                    q_n = q_i2n.get(q_i)
                    if q_n is not None:
                        parts["qrnum"].append((x, q_n))
                    if search_db != "pred":
                        t_n = t_i2n.get(t_i)
                        if t_n is not None:
                            parts["trnum"].append((x, t_n))
            if q_a != "-":
                q_i += 1
            if t_a != "-":
                t_i += 1
        return parts

    def get_protein_scores_pair(self, aln: pd.DataFrame) -> dict[str, float]:
        scores = {}
        for source, data in aln.iterrows():
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
                    aln = query_target_entry_alignments.loc[q_chain].loc[t_chain]
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
        # qrnum/trnum entries from map_row are uniformly (aln_position,
        # residue_number) pairs — see map_row's docstring. Source-agnostic.
        for q_instance_chain, t_instance_chain in alns:
            aln = alns[(q_instance_chain, t_instance_chain)]
            q_chain_pocket = query_pocket.get(q_instance_chain, {})
            q_chain_interactions = query_interactions.get(q_instance_chain, {})
            t_chain_pocket = target_pocket.get(t_instance_chain, {})
            t_chain_interactions = target_interactions.get(t_instance_chain, {})
            for source, aln_source in aln.iterrows():
                for i, q_n in aln_source["qrnum"].items():
                    if q_n not in q_chain_pocket:
                        continue
                    q_a = aln_source["qaln"][i]
                    t_a = aln_source["taln"][i]
                    t_n = aln_source["trnum"].get(i)
                    if q_a == t_a:
                        pocket_scores[f"pocket_fident_{source}"] += 1
                    if has_target and t_n is not None and t_n in t_chain_pocket:
                        pocket_scores[f"pocket_qcov_{source}"] += 1
                        if q_a == t_a:
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
        ) = self._protein_only_pocket_data(
            query_ligand.pdb_id,
            query_ligand.protein_chains_asym_id,
            query_ligand.pocket_residue_number_to_index,
            query_ligand.interactions_counter,
        )
        target_pocket, target_interactions, _, _, _ = self._protein_only_pocket_data(
            target_ligand.pdb_id,
            target_ligand.protein_chains_asym_id,
            target_ligand.pocket_residue_number_to_index,
            target_ligand.interactions_counter,
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
        ) = self._protein_only_pocket_data(
            query_ligand.pdb_id,
            query_ligand.protein_chains_asym_id,
            query_ligand.pocket_residue_number_to_index,
            query_ligand.interactions_counter,
        )
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
        query_ligand_ids: set[str] | None = None,
        target_system_ids: set[str] | None = None,
        target_ligand_ids: set[str] | None = None,
    ) -> abc.Generator[dict[str, str | float | None], None, None]:
        if search_db == "holo":
            return self.get_scores_holo(
                query_system,
                query_entry_alignments,
                data_dir=data_dir,
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
        query_ligand_ids: set[str] | None = None,
        target_system_ids: set[str] | None = None,
        target_ligand_ids: set[str] | None = None,
    ) -> abc.Generator[dict[str, str | float | None], None, None]:
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
                if target_system.system_type != "holo" or (
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
                        # receptor chains belonging to this ligand pair.
                        (
                            q_t_mappings,
                            protein_scores,
                            alns,
                            protein_chain_mapper,
                        ) = self.get_protein_scores(
                            query_target_entry_alignments,
                            query_system,
                            target_protein_chains,
                            query_protein_length,
                            query_protein_chains=query_protein_chains,
                        )
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
                        if data_dir is not None and pocket_qcov > 0:
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
        column_mapr = self.get_column_mapr()
        pdb_vals = []
        for system in self.entries[pdb_id].systems.values():
            if system.system_type != "holo" or (
                query_system_ids is not None and system.id not in query_system_ids
            ):
                continue
            for score_dict in self.get_scores(
                search_db,
                system,
                alignments.loc[pdb_id],
                data_dir=data_dir,
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
