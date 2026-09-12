# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""CCD-anchored ligand relationship databases (MMP + ECFP4/Tanimoto).

Builds the ligand relationship databases over the fixed Chemical Component
Dictionary universe rather than one ingest run's ligands, so the artifacts are
stable across releases and dataset ligands join by CCD identity first and fall
back to chemistry only when they have no CCD identity.

See ``make_ligand_dbs.md`` for the design and the wire-in checklist.

Notes
-----
The component table deliberately reuses the column names of the existing
unique-SMILES table (``ligand_smiles_id`` / ``ligand_rdkit_canonical_smiles``)
so :mod:`plinder.data.annotations.mmpdb_utils` and the ECFP4 writer in
:mod:`plinder.data.annotations.get_similarity_scores` can be reused verbatim on
this universe.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sqlite3
from collections import Counter, deque
from collections.abc import Collection, Iterable, Mapping, Sequence
from difflib import SequenceMatcher
from functools import cache, lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, NamedTuple, TypedDict, cast

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

if TYPE_CHECKING:
    from collections.abc import Callable

    import biotite.structure as struc
    from rdkit.Chem.rdchem import Mol

import pandas as pd

from plinder.core.structure.smallmols_similarity import parity_similarity
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

FRAGMENT_OPTIONS_METADATA = b"mmpdb_fragment_options"
FRAGMENT_COLUMNS = [
    "ligand_smiles_id",
    "normalized_smiles",
    "num_cuts",
    "enumeration_label",
    "variable_num_heavies",
    "variable_symmetry_class",
    "variable_smiles",
    "attachment_order",
    "constant_num_heavies",
    "constant_symmetry_class",
    "constant_smiles",
    "constant_with_H_smiles",
]
HYDROGEN = ("[*][H]", "1", "0", "N", 0)  # mmpdb's hydrogen variable fragment


# TODO: three ad-hoc bundled-CCD readers now exist - this one,
# ``protein_utils._ccd_parent_components`` and ``ccd_template``'s leaving-atom
# lookup.
# Hoist a shared accessor into ``cif_utils`` rather than adding a fourth.
@lru_cache(maxsize=1)
def _ccd_chem_comp_table() -> pd.DataFrame:
    """Return biotite's bundled CCD ``chem_comp`` category as a DataFrame.

    Only the columns this module consumes, all as strings: the CCD stores
    everything as text and nothing here does arithmetic on it.
    """
    from biotite.structure.info.ccd import get_ccd

    category = get_ccd()["chem_comp"]
    wanted = ["id", "type", "pdbx_release_status"]
    columns = {
        name: category[name].as_array(str) for name in wanted if name in category
    }
    return pd.DataFrame(columns)


def _ccd_smiles(comp_id: str) -> str | None:
    """Canonical SMILES for one CCD component, or None when underivable."""
    from plinder.data.annotations.ligand_utils import _get_ccd_smiles

    try:
        return _get_ccd_smiles(comp_id)
    except Exception as exc:  # noqa: BLE001 - one bad component must not abort
        LOG.warning(f"CCD SMILES failed for {comp_id}: {exc}")
        return None


def ccd_component_table(
    *,
    exclude_non_druglike: bool = True,
    exclude_polymer_linking: bool = False,
    limit: int | None = None,
) -> pd.DataFrame:
    """Build the active-CCD node universe with canonical SMILES.

    Parameters
    ----------
    exclude_non_druglike : bool, default=True
        Apply plinder's own :func:`~plinder.data.annotations.ligand_utils.is_excluded_mol`
        filter, the same rule that marks a dataset ligand an artifact. Keeps the
        CCD universe consistent with native ligand treatment rather than
        inventing a second notion of "real ligand". Of the dataset's four
        exclusion sources (these chemistry rules, the artifact bad-list, the
        cofactor list and the dummy codes) only this one applies here: listed
        artifacts and cofactors stay nodes, since a ligand may still be related
        to them; dummies and ions drop out for lacking a molecule.
    exclude_polymer_linking : bool, default=False
        Drop ``*-linking`` components (amino acids, nucleotides, saccharides).
        plinder treats these as polymer/receptor rather than ligand, so enable
        this for a ligand-only universe; the default keeps the broad chemical
        reference set.
    limit : int or None
        Process only the first *limit* released components (testing aid).

    Returns
    -------
    pd.DataFrame
        Columns ``ccd_id``, ``ligand_rdkit_canonical_smiles`` (stereo included,
        the key the dataset deduplicates ligands on and the one chemistry
        matching uses), ``ligand_smiles_id``, ``num_heavy_atoms`` and
        ``ccd_type``. The
        ``ligand_*`` names intentionally match the unique-SMILES table so the
        MMP and ECFP4 builders can be reused.

    Notes
    -----
    Only ``pdbx_release_status == "REL"`` components are kept, so obsolete
    entries never enter the databases. Placeholder codes are dropped via
    plinder's own :func:`~plinder.data.annotations.ligand_utils.lig_has_dummies`.
    """
    from rdkit import Chem

    from plinder.data.annotations.ligand_utils import is_excluded_mol, lig_has_dummies

    table = _ccd_chem_comp_table()
    released = table[table["pdbx_release_status"] == "REL"]
    # plinder's own dummy/unknown list, not a local copy: it covers UNK, UPL,
    # DN and N, which carry real SMILES and would otherwise enter the universe
    # as if they were genuine chemistry.
    comp_ids = [
        comp_id for comp_id in released["id"].tolist() if not lig_has_dummies(comp_id)
    ]
    if limit is not None:
        comp_ids = comp_ids[:limit]
    LOG.info(f"ccd_component_table: deriving SMILES for {len(comp_ids)} components")

    type_by_id = (
        dict(zip(released["id"], released["type"])) if "type" in released else {}
    )

    rows: list[dict[str, object]] = []
    skipped_no_smiles = 0
    skipped_non_druglike = 0
    skipped_polymer = 0
    for comp_id in comp_ids:
        ccd_type = str(type_by_id.get(comp_id, "")).lower()
        if exclude_polymer_linking and "linking" in ccd_type:
            skipped_polymer += 1
            continue
        smiles = _ccd_smiles(comp_id)
        if not smiles:
            skipped_no_smiles += 1
            continue
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            skipped_no_smiles += 1
            continue
        num_heavy = mol.GetNumHeavyAtoms()
        # Same rule plinder uses to call a dataset ligand an artifact, so the
        # reference universe and the annotation path agree on "real ligand".
        if exclude_non_druglike and is_excluded_mol(smiles):
            skipped_non_druglike += 1
            continue
        rows.append(
            {
                "ccd_id": comp_id,
                "ligand_rdkit_canonical_smiles": smiles,
                "num_heavy_atoms": int(num_heavy),
                # CCD ``type`` is mixed case in the bundle ("NON-POLYMER" vs
                # "non-polymer"); normalise so consumers can group on it.
                "ccd_type": ccd_type,
            }
        )
    LOG.info(
        f"ccd_component_table: kept {len(rows)}; skipped {skipped_no_smiles} "
        f"without SMILES, {skipped_non_druglike} non-druglike, "
        f"{skipped_polymer} polymer-linking"
    )
    components = pd.DataFrame(rows).sort_values("ccd_id", ignore_index=True)
    # Stable integer node id: assigned over the sorted CCD ids so the universe
    # is reproducible for a given CCD bundle. Must be 0-based and contiguous -
    # ``get_similarity_scores.ligand_scores`` indexes the fingerprint list
    # directly by node id and rejects any other numbering.
    components["ligand_smiles_id"] = range(len(components))
    return components


def ccd_universe_signature(components: pd.DataFrame) -> str:
    """Digest of the component universe, for the build manifest."""
    digest = hashlib.sha256()
    for ccd_id, smiles in components[
        ["ccd_id", "ligand_rdkit_canonical_smiles"]
    ].itertuples(index=False, name=None):
        digest.update(str(ccd_id).encode())
        digest.update(b"\0")
        digest.update(str(smiles).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def ccd_dbs_dir(data_dir: Path) -> Path:
    """Directory holding the CCD-anchored database artifacts."""
    return Path(data_dir) / "ccd_dbs"


def ccd_fingerprint_path(data_dir: Path) -> Path:
    """Path of the CCD ECFP4 table, in the layout the Tanimoto scorer expects.

    :func:`plinder.data.annotations.get_similarity_scores.ligand_scores`
    resolves ``<data_dir>/fingerprints/ligands_per_smiles.parquet`` internally,
    so mirroring that layout under ``ccd_dbs`` lets the existing scorer run on
    the CCD universe unchanged - no new parameter, no edit to a shared file.
    """
    return ccd_dbs_dir(data_dir) / "fingerprints" / "ligands_per_smiles.parquet"


def build_ccd_tanimoto_scores(
    components: pd.DataFrame,
    data_dir: Path,
    *,
    minimum_similarity: float = 30.0,
    batch_size: int = 5_000,
) -> Path:
    """Compute ECFP4 Tanimoto edges over the CCD universe.

    Parameters
    ----------
    components : pd.DataFrame
        Output of :func:`ccd_component_table`; supplies the node ids.
    data_dir : Path
        Release directory; scores land in ``<data_dir>/ccd_dbs/ligand_scores``.
    minimum_similarity : float
        Percentage cutoff below which edges are dropped.
    batch_size : int
        Node ids scored per output shard, mirroring the ingest scatter so a
        50k-node universe never materialises an all-pairs matrix at once.

    Returns
    -------
    Path
        Directory holding the per-batch score shards.

    Notes
    -----
    Delegates to the same ``ligand_scores`` routine the ingest pipeline uses,
    pointed at the CCD ``fingerprints`` layout, so the edge semantics and the
    ECFP4 metadata check are identical to the dataset-level scores.
    """
    from plinder.data.annotations import get_similarity_scores

    ccd_dir = ccd_dbs_dir(data_dir)
    scores_dir = ccd_dir / "ligand_scores"
    scores_dir.mkdir(parents=True, exist_ok=True)
    # Shards are numbered per batch, so a rebuild over a smaller universe would
    # otherwise leave higher-numbered shards behind and any consumer globbing
    # this directory would read them as live edges referencing dead node ids.
    for stale in scores_dir.glob("ccd_scores_*.parquet"):
        stale.unlink()
    node_ids = [int(value) for value in components["ligand_smiles_id"]]
    for start in range(0, len(node_ids), batch_size):
        batch = node_ids[start : start + batch_size]
        output_path = scores_dir / f"ccd_scores_{start // batch_size:04d}.parquet"
        get_similarity_scores.ligand_scores(
            ligand_ids=batch,
            data_dir=ccd_dir,
            output_path=output_path,
            minimum_similarity=minimum_similarity,
        )
    LOG.info(
        f"build_ccd_tanimoto_scores: scored {len(node_ids)} nodes into {scores_dir}"
    )
    return scores_dir


def build_ccd_ecfp_db(components: pd.DataFrame, output_path: Path) -> Path:
    """Write the ECFP4/1024 fingerprint table over the CCD universe.

    Mirrors ``fingerprints/ligands_per_smiles.parquet`` (same columns, same
    parquet metadata) so the existing similarity/clustering machinery can run
    against this node universe unchanged.
    """
    from rdkit import DataStructs

    from plinder.core.structure.smallmols_similarity import mol2morgan_fp
    from plinder.data.annotations.get_similarity_scores import (
        ECFP4_NBITS,
        ECFP4_RADIUS,
        write_ecfp4_fingerprint_table,
    )

    table = components.copy()
    binary_fingerprints: list[bytes] = []
    for smiles in table["ligand_rdkit_canonical_smiles"]:
        fingerprint = mol2morgan_fp(smiles, radius=ECFP4_RADIUS, nbits=ECFP4_NBITS)
        binary_fingerprints.append(DataStructs.BitVectToBinaryText(fingerprint))
    table["fingerprint"] = binary_fingerprints

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(".parquet.tmp")
    temporary_path.unlink(missing_ok=True)
    write_ecfp4_fingerprint_table(table, temporary_path)
    temporary_path.replace(output_path)
    LOG.info(f"build_ccd_ecfp_db: wrote {len(table)} fingerprints to {output_path}")
    return output_path


def build_ccd_mmp_db(
    components: pd.DataFrame,
    output_path: Path,
    *,
    scratch_dir: Path,
    threads: int = 4,
    fragments_path: Path | None = None,
) -> Path:
    """Write matched molecular pairs over the CCD universe.

    Reuses the existing mmpdb flow (``smi_split -> fragment ->
    fragdb_partition -> index``) unchanged; only the input universe differs.

    Parameters
    ----------
    fragments_path : Path or None
        When given, persist every fragmentation of the universe there, keyed by
        constant part; :func:`query_ccd_mmp_pairs` pairs molecules outside the
        universe against it without a rebuild.

    Notes
    -----
    The universe is CCD components only; composite ligands are deliberately not
    added. Pairing joins fragmentations on their shared constant part, and
    residues of one class leave the same constant when their side chain is cut -
    every amino acid leaves the backbone - so composites of that class all pair
    with one another: recall without precision. Fragmentation also enumerates
    cuts over acyclic single bonds, of which a composite has one per linkage
    plus its side chains, so cost grows combinatorially with residue count while
    those shared constants make the pair count quadratic in the number of
    composites. Relate composites by their residue graph instead
    (:func:`composite_similarity`), and use this database only to describe a
    residue substitution, never to decide whether two composites are related.
    """
    from tempfile import TemporaryDirectory

    from plinder.data.annotations.mmpdb_utils import (
        _generate_pair_files,
        _write_pair_parquet,
    )

    executable = shutil.which("mmpdb")
    if executable is None:
        raise RuntimeError(
            "mmpdb is required to build the CCD MMP database; install the "
            "PLINDER data dependencies"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scratch_dir = Path(scratch_dir)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(".parquet.tmp")
    temporary_path.unlink(missing_ok=True)
    with TemporaryDirectory(prefix="plinder-ccd-mmp-", dir=scratch_dir) as work:
        work_dir = Path(work)
        pair_files = _generate_pair_files(
            ligands=components,
            work_dir=work_dir,
            threads=threads,
            executable=executable,
        )
        _write_pair_parquet(
            pair_files=pair_files,
            ligands=components,
            output_path=temporary_path,
        )
        if fragments_path is not None:
            _export_fragmentations(work_dir, Path(fragments_path))
    temporary_path.replace(output_path)
    LOG.info(f"build_ccd_mmp_db: wrote matched pairs to {output_path}")
    return output_path


def _export_fragmentations(work_dir: Path, fragments_path: Path) -> Path:
    """Persist the universe's fragmentations as the constant-keyed pair dictionary.

    ``mmpdb index`` forms pairs by joining fragmentations on their constant part;
    the ``ligands.NNNN.fragdb`` shards hold exactly those rows. Sorted by
    constant so parquet row-group statistics let a lookup skip most of the file,
    with the fragment options in the parquet metadata so queries fragment the
    same way.
    """
    frames: list[pd.DataFrame] = []
    options: dict[str, Any] | None = None
    record_columns = {
        "ligand_smiles_id": "r.title AS ligand_smiles_id",
        "normalized_smiles": "r.normalized_smiles",
    }
    select = ", ".join(
        record_columns.get(name, f"f.{name}") for name in FRAGMENT_COLUMNS
    )
    for database in sorted(work_dir.glob("ligands.[0-9][0-9][0-9][0-9].fragdb")):
        connection = sqlite3.connect(database)
        try:
            frames.append(
                pd.read_sql_query(
                    f"SELECT {select} FROM fragmentation f "
                    "JOIN record r ON r.id = f.record_id",
                    connection,
                )
            )
            names = [row[1] for row in connection.execute("PRAGMA table_info(options)")]
            values = connection.execute("SELECT * FROM options").fetchone()
        finally:
            connection.close()
        options = {
            name: value
            for name, value in zip(names, values)
            if name not in {"id", "version"}
        }
    if options is None:
        raise RuntimeError("mmpdb fragmentation left no fragment databases behind")
    fragments = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=FRAGMENT_COLUMNS)
    )
    fragments["ligand_smiles_id"] = fragments["ligand_smiles_id"].astype(int)
    fragments = fragments.sort_values(
        ["constant_smiles", "num_cuts"], ignore_index=True
    )
    table = pa.Table.from_pandas(
        fragments, preserve_index=False
    ).replace_schema_metadata({FRAGMENT_OPTIONS_METADATA: json.dumps(options).encode()})
    fragments_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = fragments_path.with_suffix(".parquet.tmp")
    pq.write_table(table, temporary_path, compression="zstd", row_group_size=50_000)
    temporary_path.replace(fragments_path)
    LOG.info(f"_export_fragmentations: wrote {len(fragments)} rows to {fragments_path}")
    return fragments_path


def _fragment_dictionary(
    fragments_path: Path, *, constants: Collection[str], molecules: Collection[str]
) -> tuple[pd.DataFrame, dict[str, list[int]]]:
    """Fragmentation rows a query can pair with, and the universe's molecule map.

    Rows share a constant with a query fragmentation, or are one-cut rows whose
    constant plus hydrogen is a query molecule. The map (normalised SMILES ->
    node ids) serves the opposite hydrogen case.
    """
    frames = []
    if constants:
        frames.append(
            pq.read_table(
                fragments_path,
                filters=[("constant_smiles", "in", sorted(set(constants)))],
            ).to_pandas()
        )
    if molecules:
        frames.append(
            pq.read_table(
                fragments_path,
                filters=[
                    ("num_cuts", "==", 1),
                    ("constant_with_H_smiles", "in", sorted(set(molecules))),
                ],
            ).to_pandas()
        )
    rows = (
        pd.concat(frames, ignore_index=True).drop_duplicates()
        if frames
        else pq.read_schema(fragments_path).empty_table().to_pandas()
    )
    nodes = (
        pq.read_table(fragments_path, columns=["ligand_smiles_id", "normalized_smiles"])
        .to_pandas()
        .drop_duplicates()
    )
    molecule_ids = {
        smiles: [int(value) for value in group]
        for smiles, group in nodes.groupby("normalized_smiles")["ligand_smiles_id"]
    }
    return rows, molecule_ids


def query_ccd_mmp_pairs(
    queries: Mapping[str, str],
    data_dir: Path,
    *,
    max_variable_heavies: int | None = 10,
    max_heavies_transf: int | None = None,
) -> pd.DataFrame:
    """Matched molecular pairs between molecules outside the universe and the CCD.

    Runs mmpdb's pairing rule against the persisted fragmentation dictionary: a
    query and a component are a pair when a fragmentation of each leaves the
    same constant part (cut count and symmetry class included) with different
    variable parts; a one-cut constant whose hydrogen-capped form is itself a
    molecule of the other side pairs against the hydrogen fragment, as in
    ``mmpdb index``; two enumerated-chirality constants never pair. SMIRKS come
    from ``mmpdblib.index_algorithm.cansmirks``, so transformations are written
    exactly as in ``ccd_mmp_pairs.parquet``. Defaults are mmpdb's index defaults
    (``--max-variable-heavies 10``).

    Parameters
    ----------
    queries : mapping
        ``{query_id: smiles}`` for the molecules to place (novel or external).
    data_dir : Path
        Release directory holding ``ccd_dbs``.
    max_variable_heavies, max_heavies_transf : int or None
        mmpdb index filters on the variable part and on the change in its size.

    Returns
    -------
    pd.DataFrame
        One row per (query, component, shared core, transformation) with
        ``query_id``, ``query_smiles``, ``ligand_smiles_id``, ``ccd_id``,
        ``ligand_smiles``, ``transformation`` (query >> component),
        ``shared_core_smiles``, ``num_cuts``, ``shared_core_num_heavy_atoms``.

    Notes
    -----
    Components with no fragmentation of their own (very small molecules) cannot
    be hydrogen partners, since the molecule map is read from the dictionary.

    TODO (discuss): the pair relation is symmetric; the directed quantity plinder
    uses is core coverage per side (``shared_core_fraction`` in
    ``ccd_mmp_pairs.parquet``), which this output lacks. Also whether mmpdb's
    index defaults (10-atom variable cap, no coverage floor) and the downstream
    5-atom core floor are the cutoffs wanted for the CCD reference.
    """
    from mmpdblib import (
        fragment_algorithm,
        fragment_records,
        fragment_types,
        index_algorithm,
    )

    fragments_path = ccd_dbs_dir(data_dir) / "ccd_mmp_fragments.parquet"
    if not fragments_path.is_file():
        raise FileNotFoundError(f"missing CCD fragmentation table: {fragments_path}")
    options = json.loads(
        pq.read_schema(fragments_path).metadata[FRAGMENT_OPTIONS_METADATA]
    )
    fragment_filter = fragment_types.get_fragment_filter(
        fragment_types.FragmentOptions(**options)
    )
    records = {}
    for query_id, smiles in queries.items():
        record = fragment_records.make_fragment_record_from_smiles(
            smiles, fragment_filter
        )
        if getattr(record, "errmsg", None):
            LOG.warning(
                f"query_ccd_mmp_pairs: cannot fragment {query_id}: {record.errmsg}"
            )
            continue
        records[query_id] = record
    constants = {
        fragmentation.constant_smiles
        for record in records.values()
        for fragmentation in record.fragmentations
    }
    molecules = {record.normalized_smiles for record in records.values()}
    dictionary, molecule_ids = _fragment_dictionary(
        fragments_path, constants=constants, molecules=molecules
    )
    groups = {
        key: group
        for key, group in dictionary.groupby(
            ["num_cuts", "constant_smiles", "constant_symmetry_class"], sort=False
        )
    }
    capped = dictionary[dictionary["num_cuts"] == 1].groupby("constant_with_H_smiles")
    relabel_cache = index_algorithm.RelabelCache()
    no_enumeration = fragment_algorithm.EnumerationLabel.NO_ENUMERATION
    rows: list[dict[str, object]] = []

    def pair(
        query_id: str,
        query_smiles: str,
        query_molecule: str,
        constant: tuple[int, str, str, int],
        query_side: tuple[str, str, str, str, int],
        target_id: int,
        target_molecule: str,
        target_side: tuple[str, str, str, str, int],
    ) -> None:
        """Append one pair if mmpdb's index filters admit it."""
        num_cuts, constant_smiles, constant_symmetry_class, constant_heavies = constant
        (
            query_smiles_,
            query_symmetry,
            query_order,
            query_label,
            query_heavies,
        ) = query_side
        (
            target_smiles,
            target_symmetry,
            target_order,
            target_label,
            target_heavies,
        ) = target_side
        if max_variable_heavies is not None and (
            query_heavies > max_variable_heavies
            or target_heavies > max_variable_heavies
        ):
            return
        if target_smiles == query_smiles_ and target_order == query_order:
            return
        if target_molecule == query_molecule:  # the query is that component
            return
        if query_label != no_enumeration and target_label != no_enumeration:
            return
        if (
            max_heavies_transf is not None
            and abs(target_heavies - query_heavies) > max_heavies_transf
        ):
            return
        smirks, shared_core = index_algorithm.cansmirks(
            num_cuts,
            query_smiles_,
            query_symmetry,
            query_order,
            constant_smiles,
            constant_symmetry_class,
            target_smiles,
            target_symmetry,
            target_order,
            relabel_cache,
        )
        rows.append(
            {
                "query_id": query_id,
                "query_smiles": query_smiles,
                "ligand_smiles_id": int(target_id),
                "transformation": smirks,
                "shared_core_smiles": shared_core,
                "num_cuts": int(num_cuts),
                "shared_core_num_heavy_atoms": int(constant_heavies),
            }
        )

    for query_id, record in records.items():
        for query in record.fragmentations:
            constant = (
                query.num_cuts,
                query.constant_smiles,
                query.constant_symmetry_class,
                query.constant_num_heavies,
            )
            query_side = (
                query.variable_smiles,
                query.variable_symmetry_class,
                query.attachment_order,
                query.enumeration_label,
                query.variable_num_heavies,
            )
            group = groups.get(constant[:3])
            if group is not None:
                for target in group.itertuples(index=False):
                    pair(
                        query_id,
                        record.input_smiles,
                        record.normalized_smiles,
                        constant,
                        query_side,
                        target.ligand_smiles_id,
                        target.normalized_smiles,
                        (
                            target.variable_smiles,
                            target.variable_symmetry_class,
                            target.attachment_order,
                            target.enumeration_label,
                            target.variable_num_heavies,
                        ),
                    )
            if query.num_cuts == 1:
                # the query minus this substituent is itself a component
                for target_id in molecule_ids.get(query.constant_with_H_smiles, []):
                    pair(
                        query_id,
                        record.input_smiles,
                        record.normalized_smiles,
                        constant,
                        query_side,
                        target_id,
                        query.constant_with_H_smiles,
                        HYDROGEN,
                    )
        if record.normalized_smiles in capped.groups:
            # the query is a component minus one substituent
            for target in capped.get_group(record.normalized_smiles).itertuples(
                index=False
            ):
                pair(
                    query_id,
                    record.input_smiles,
                    record.normalized_smiles,
                    (
                        1,
                        target.constant_smiles,
                        target.constant_symmetry_class,
                        target.constant_num_heavies,
                    ),
                    HYDROGEN,
                    target.ligand_smiles_id,
                    target.normalized_smiles,
                    (
                        target.variable_smiles,
                        target.variable_symmetry_class,
                        target.attachment_order,
                        target.enumeration_label,
                        target.variable_num_heavies,
                    ),
                )
    pairs = pd.DataFrame(
        rows,
        columns=[
            "query_id",
            "query_smiles",
            "ligand_smiles_id",
            "transformation",
            "shared_core_smiles",
            "num_cuts",
            "shared_core_num_heavy_atoms",
        ],
    ).drop_duplicates()
    components = load_ccd_components(data_dir)[
        ["ligand_smiles_id", "ccd_id", "ligand_rdkit_canonical_smiles"]
    ].rename(columns={"ligand_rdkit_canonical_smiles": "ligand_smiles"})
    return pairs.merge(components, on="ligand_smiles_id", how="left").sort_values(
        ["query_id", "ligand_smiles_id", "shared_core_smiles", "transformation"],
        ignore_index=True,
    )


class CcdMatch(TypedDict):
    """One ligand resolved against the CCD universe."""

    match_kind: str
    matched_ccd_ids: list[str]
    component_ccd_ids: list[str]


class CcdIndex(NamedTuple):
    """Lookups over the universe that :func:`match_ligand_to_ccd` needs."""

    ccd_ids: frozenset[str]
    smiles_to_ccd_ids: dict[
        str, list[str]
    ]  # canonical SMILES, stereo included -> codes
    released_ids: frozenset[str]


def ccd_index(components: pd.DataFrame) -> CcdIndex:
    """Build the matcher's lookups from :func:`ccd_component_table` output."""
    # both sides through canonical_smiles: re-canonicalising is not idempotent
    # for a few dative/cage SMILES, so the stored string alone would not match
    by_smiles: dict[str, list[str]] = {}
    for ccd_id, smiles in zip(
        components["ccd_id"].astype(str),
        components["ligand_rdkit_canonical_smiles"].astype(str),
    ):
        by_smiles.setdefault(canonical_smiles(smiles) or smiles, []).append(ccd_id)
    released = _ccd_chem_comp_table().query("pdbx_release_status == 'REL'")["id"]
    return CcdIndex(
        frozenset(components["ccd_id"].astype(str)),
        by_smiles,
        frozenset(released.astype(str)),
    )


def canonical_smiles(smiles: str | None) -> str | None:
    """RDKit canonical SMILES with stereo, the form the CCD table stores."""
    from rdkit import Chem

    from plinder.core.utils.sanitize import mol_from_smiles

    if not smiles:
        return None
    mol = mol_from_smiles(smiles)
    return None if mol is None else str(Chem.MolToSmiles(mol))


def match_ligand_to_ccd(
    ccd_code: str, ligand_smiles: str | None, index: CcdIndex
) -> CcdMatch:
    """Resolve one ligand against the CCD universe by code, then by chemistry.

    ``ligand_smiles`` must come through :func:`canonical_smiles`, as the index
    keys do. ``match_kind`` is one of

    ``exact``
        The ligand is a universe component: its single code is one, or its
        whole-molecule canonical SMILES equals a component's. The second test
        is what lets one molecule match however it was deposited, e.g. lactose
        as ``LAT`` or as the linked sugars ``GAL-BGC``. Stereo is part of the
        identity, as in the dataset's own ligand key: the ligand's SMILES
        carries the stereo perceived from its coordinates, the component's the
        stereo of its ideal coordinates, so cellobiose and maltose never match
        lactose.
    ``excluded``
        A single code that is a released CCD component the universe filters
        out (ions, artifacts, placeholders): known, not novel.
    ``novel``
        No component is this molecule. Relating it to the universe is the
        fingerprint's job (:func:`query_ccd_tanimoto`, :func:`query_ccd_mmp_pairs`).

    ``component_ccd_ids`` lists a composite's codes that are universe
    components, as annotation only: a missing component says nothing about
    what the whole molecule is.
    """
    codes = [code for code in str(ccd_code or "").split("-") if code]
    matched: list[str] = []
    if len(codes) == 1 and codes[0] in index.ccd_ids:
        matched = [codes[0]]
    elif ligand_smiles:
        matched = sorted(index.smiles_to_ccd_ids.get(ligand_smiles, []))
    if matched:
        kind = "exact"
    elif len(codes) == 1 and codes[0] in index.released_ids:
        kind = "excluded"
    else:
        kind = "novel"
    return {
        "match_kind": kind,
        "matched_ccd_ids": matched,
        "component_ccd_ids": [code for code in codes if code in index.ccd_ids],
    }


def load_ccd_components(data_dir: Path) -> pd.DataFrame:
    """Load the persisted component universe the artifacts were built from.

    Always prefer this over rebuilding with :func:`ccd_component_table` when
    joining against existing artifacts: node ids are positions in *this* table,
    so a table built with different filters would silently mis-map every id.
    """
    path = ccd_dbs_dir(data_dir) / "ccd_components.parquet"
    if not path.is_file():
        raise FileNotFoundError(f"missing CCD component table: {path}")
    return pd.read_parquet(path)


def build_ligand_ccd_match(
    ligands: pd.DataFrame,
    components: pd.DataFrame,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Map dataset ligands onto the CCD universe.

    Parameters
    ----------
    ligands : pd.DataFrame
        ``ligand_id``, ``ligand_ccd_code`` and ``ligand_smiles`` (the
        whole-molecule SMILES; composites are one connected molecule), plus
    components : pd.DataFrame
        Output of :func:`ccd_component_table`.
    output_path : Path or None
        When given, write the join table as parquet.

    Returns
    -------
    pd.DataFrame
        ``ligand_id``, ``ligand_ccd_code``, ``match_kind``, ``matched_ccd_ids``,
        ``ccd_node_ids`` (nodes of the matched components), ``component_ccd_ids``.

    Notes
    -----
    Emitted as a sidecar join table rather than annotation columns, so building
    these databases never contends with the index schema.
    """
    index = ccd_index(components)
    node_by_ccd = dict(
        zip(components["ccd_id"].astype(str), components["ligand_smiles_id"])
    )
    canonical_of = {
        smiles: canonical_smiles(smiles)
        for smiles in ligands["ligand_smiles"].dropna().unique()
    }
    rows: list[dict[str, object]] = []
    for ligand_id, ccd_code, smiles in zip(
        ligands["ligand_id"], ligands["ligand_ccd_code"], ligands["ligand_smiles"]
    ):
        match = match_ligand_to_ccd(str(ccd_code), canonical_of.get(smiles), index)
        rows.append(
            {
                "ligand_id": ligand_id,
                "ligand_ccd_code": ccd_code,
                "match_kind": match["match_kind"],
                "matched_ccd_ids": match["matched_ccd_ids"],
                "ccd_node_ids": [int(node_by_ccd[c]) for c in match["matched_ccd_ids"]],
                "component_ccd_ids": match["component_ccd_ids"],
            }
        )
    matches = pd.DataFrame(rows)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = output_path.with_suffix(".parquet.tmp")
        temporary_path.unlink(missing_ok=True)
        matches.to_parquet(temporary_path, index=False)
        temporary_path.replace(output_path)
    return matches


def make_ligand_ccd_match(*, data_dir: Path) -> Path:
    """Write ``ccd_dbs/ligand_ccd_match.parquet`` for the collated annotation table."""
    annotation_path = data_dir / "index" / "annotation_table.parquet"
    if not annotation_path.is_file():
        raise FileNotFoundError(
            f"collate the annotation table first: {annotation_path}"
        )
    columns = [
        "ligand_id",
        "ligand_ccd_code",
        "ligand_smiles",
    ]
    available = set(pq.read_schema(annotation_path).names)
    ligands = pd.read_parquet(
        annotation_path, columns=[c for c in columns if c in available]
    ).drop_duplicates()
    output_path = ccd_dbs_dir(data_dir) / "ligand_ccd_match.parquet"
    build_ligand_ccd_match(
        ligands, load_ccd_components(data_dir), output_path=output_path
    )
    return output_path


def query_ccd_tanimoto(
    queries: Mapping[str, str], data_dir: Path, *, minimum_similarity: float = 30.0
) -> pd.DataFrame:
    """ECFP4/1024 Tanimoto neighbours of molecules outside the universe.

    The same fingerprint, table and percentage cutoff as
    :func:`~plinder.data.annotations.get_similarity_scores.ligand_scores`, so a
    novel or external ligand (a composite included, through its whole-molecule
    SMILES) is related to the CCD exactly as dataset ligands are to each other.

    Returns
    -------
    pd.DataFrame
        ``query_id``, ``query_smiles``, ``ligand_smiles_id``,
        ``tanimoto_similarity_ecfp4_1024`` (percent), ``ccd_id``; best hit first.
    """
    from rdkit import DataStructs

    from plinder.core.structure.smallmols_similarity import mol2morgan_fp
    from plinder.data.annotations.get_similarity_scores import (
        ECFP4_NBITS,
        ECFP4_PARQUET_METADATA,
        ECFP4_RADIUS,
    )

    fingerprint_path = ccd_fingerprint_path(data_dir)
    metadata = pq.read_schema(fingerprint_path).metadata or {}
    if any(metadata.get(key) != value for key, value in ECFP4_PARQUET_METADATA.items()):
        raise ValueError("CCD fingerprint metadata is not ECFP4/1024")
    table = pd.read_parquet(
        fingerprint_path, columns=["ligand_smiles_id", "fingerprint"]
    )
    fingerprints = [DataStructs.CreateFromBinaryText(v) for v in table["fingerprint"]]
    node_ids = table["ligand_smiles_id"].astype(int).tolist()
    minimum_fraction = minimum_similarity / 100.0
    rows: list[dict[str, object]] = []
    for query_id, smiles in queries.items():
        try:
            fingerprint = mol2morgan_fp(smiles, radius=ECFP4_RADIUS, nbits=ECFP4_NBITS)
        except Exception as exc:  # noqa: BLE001 - one bad query must not abort
            LOG.warning(f"query_ccd_tanimoto: cannot fingerprint {query_id}: {exc}")
            continue
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprint, fingerprints)
        rows.extend(
            {
                "query_id": query_id,
                "query_smiles": smiles,
                "ligand_smiles_id": node_id,
                "tanimoto_similarity_ecfp4_1024": similarity * 100.0,
            }
            for node_id, similarity in zip(node_ids, similarities)
            if similarity >= minimum_fraction
        )
    hits = pd.DataFrame(
        rows,
        columns=[
            "query_id",
            "query_smiles",
            "ligand_smiles_id",
            "tanimoto_similarity_ecfp4_1024",
        ],
    )
    components = load_ccd_components(data_dir)[["ligand_smiles_id", "ccd_id"]]
    return hits.merge(components, on="ligand_smiles_id", how="left").sort_values(
        ["query_id", "tanimoto_similarity_ecfp4_1024", "ligand_smiles_id"],
        ascending=[True, False, True],
        ignore_index=True,
    )


def ccd_parity_path(data_dir: Path) -> Path:
    """Path of the CCD-vs-CCD PARITY-like score table."""
    return ccd_dbs_dir(data_dir) / "ccd_parity_scores.parquet"


@cache
def _ccd_named_mol(ccd_id: str) -> Mol | None:
    """RDKit molecule of a CCD component carrying its atom names, or None."""
    from plinder.data.annotations.ligand_utils import _get_ccd_mol

    return _get_ccd_mol(ccd_id)


def _ccd_heavy_atoms(ccd_id: str) -> tuple[tuple[str, ...], frozenset[frozenset[str]]]:
    """Heavy-atom names of a CCD component and its bonds as name pairs."""
    from plinder.core.structure.ccd_template import ccd_component_template

    template = ccd_component_template(ccd_id)
    if template is None:
        return (), frozenset()
    return template.heavy, frozenset(frozenset((a, b)) for a, b, _ in template.bonds)


def _atom_names(mol: Mol) -> list[str]:
    """PDB atom names when the molecule carries them, else atom indices."""
    names = []
    for index in range(mol.GetNumAtoms()):
        info = mol.GetAtomWithIdx(index).GetPDBResidueInfo()
        names.append(info.GetName().strip() if info else str(index))
    return names


def _parity_rows(
    pairs: Sequence[tuple[Mol | None, Mol | None, bool]], *, threads: int
) -> pd.DataFrame:
    """PARITY-like score (percent) and largest-fragment atom names per molecule pair.

    Each pair is ``(first, second, keep_alternatives)``; the scored fragment is
    always the best of every equally good MCES, and with ``keep_alternatives``
    the others are stored after it.
    """
    from concurrent.futures import ThreadPoolExecutor

    from plinder.core.structure.smallmols_similarity import rascal_parity_match

    def row(
        pair: tuple[Mol | None, Mol | None, bool],
    ) -> tuple[float, list[list[str]], list[list[str]]]:
        first, second, keep_alternatives = pair
        if first is None or second is None:
            return float("nan"), [], []
        match = rascal_parity_match(first, second, all_best=keep_alternatives)
        names_1, names_2 = _atom_names(first), _atom_names(second)
        fragments = (match.atoms, *(match.alternatives if keep_alternatives else ()))
        return (
            match.score(first, second) * 100.0,
            [[names_1[i] for i in fragment] for fragment in fragments],
            [[names_2[j] for j in fragment.values()] for fragment in fragments],
        )

    with ThreadPoolExecutor(threads) as pool:
        rows = list(pool.map(row, pairs))
    return pd.DataFrame(
        rows, columns=["parity_similarity", "fragment_atoms_1", "fragment_atoms_2"]
    )


class CcdParityTable:
    """Mono-to-mono PARITY-like matches keyed by CCD code pair, for composite scoring.

    Rows carry ``ccd_id_1``, ``ccd_id_2``, ``parity_similarity`` (percent) and
    the largest-fragment atom names ``fragment_atoms_1`` -> ``fragment_atoms_2``,
    one list per equally good fragment, the scored one first.
    """

    def __init__(self, rows: pd.DataFrame) -> None:
        self._rows: dict[tuple[str, str], tuple[float, list[dict[str, str]]]] = {
            (first, second): (
                score / 100.0,
                [dict(zip(one, other)) for one, other in zip(names_1, names_2)],
            )
            for first, second, score, names_1, names_2 in rows[
                [
                    "ccd_id_1",
                    "ccd_id_2",
                    "parity_similarity",
                    "fragment_atoms_1",
                    "fragment_atoms_2",
                ]
            ].itertuples(index=False, name=None)
            if score == score
        }

    @classmethod
    def load(cls, data_dir: Path) -> CcdParityTable:
        """The precompiled table of :func:`build_ccd_parity_scores`."""
        codes = load_ccd_components(data_dir).set_index("ligand_smiles_id")["ccd_id"]
        table = pd.read_parquet(ccd_parity_path(data_dir))
        table["ccd_id_1"] = table["ligand_smiles_id_1"].map(codes)
        table["ccd_id_2"] = table["ligand_smiles_id_2"].map(codes)
        return cls(table)

    @classmethod
    def compute(cls, codes: Iterable[str], *, threads: int = 1) -> CcdParityTable:
        """Match every pair of *codes* now, for small residue sets and tests."""
        unique = sorted(set(codes))
        mols = {code: _ccd_named_mol(code) for code in unique}
        pairs = [(a, b) for k, a in enumerate(unique) for b in unique[k + 1 :]]
        rows = _parity_rows(
            [(mols[a], mols[b], True) for a, b in pairs], threads=threads
        )
        rows.insert(0, "ccd_id_2", [b for _, b in pairs])
        rows.insert(0, "ccd_id_1", [a for a, _ in pairs])
        return cls(rows)

    def _lookup(self, first: str, second: str) -> tuple[float, list[dict[str, str]]]:
        if (first, second) in self._rows:
            return self._rows[(first, second)]
        score, mappings = self._rows.get((second, first), (0.0, []))
        return score, [{q: p for p, q in mapping.items()} for mapping in mappings]

    def similarity(self, first: str, second: str) -> float:
        """Node kernel: 1 for identical codes, the table score, 0 when absent."""
        return 1.0 if first == second else self._lookup(first, second)[0]

    def mappings(self, first: str, second: str) -> list[dict[str, str]]:
        """Equally good largest fragments of *first* onto *second*, scored one first."""
        if first == second:
            return [{name: name for name in _ccd_heavy_atoms(first)[0]}]
        return self._lookup(first, second)[1]

    def mapping(self, first: str, second: str) -> dict[str, str]:
        """The scored fragment of *first* onto *second*; empty when absent."""
        return next(iter(self.mappings(first, second)), {})


def build_ccd_parity_scores(
    components: pd.DataFrame, data_dir: Path, *, threads: int = 4
) -> Path:
    """Score every ECFP4 edge of the CCD universe with :func:`rascal_parity_score`.

    The Tanimoto shards written by :func:`build_ccd_tanimoto_scores` are the
    prefilter: each unordered pair above the ECFP4 cutoff is scored once, so the
    graph-based score never runs over the all-pairs matrix. Molecules are built
    from the CCD atoms so the stored fragment carries CCD atom names, which is
    what :func:`composite_parity` joins on. A pair with an unbuildable molecule
    scores NaN rather than being dropped.

    Returns
    -------
    Path
        Parquet with ``ligand_smiles_id_1 < ligand_smiles_id_2``,
        ``tanimoto_similarity_ecfp4_1024``, ``parity_similarity`` (percent) and
        the largest-fragment atom names ``fragment_atoms_1`` -> ``fragment_atoms_2``.
    """
    scores_dir = ccd_dbs_dir(data_dir) / "ligand_scores"
    edges = pd.concat(
        [
            pd.read_parquet(shard)
            for shard in sorted(scores_dir.glob("ccd_scores_*.parquet"))
        ],
        ignore_index=True,
    )
    edges = edges[edges["query_ligand_id"] < edges["target_ligand_id"]].rename(
        columns={
            "query_ligand_id": "ligand_smiles_id_1",
            "target_ligand_id": "ligand_smiles_id_2",
        }
    )
    codes = components.set_index("ligand_smiles_id")["ccd_id"]
    # only residues carry linkages, so only their pairs need alternative fragments
    linking = components.set_index("ligand_smiles_id")["ccd_type"].str.contains(
        "LINKING", case=False
    )
    edges = edges.sort_values(
        ["ligand_smiles_id_1", "ligand_smiles_id_2"], ignore_index=True
    )
    pairs = [
        (
            _ccd_named_mol(codes[i]),
            _ccd_named_mol(codes[j]),
            bool(linking[i] and linking[j]),
        )
        for i, j in edges[["ligand_smiles_id_1", "ligand_smiles_id_2"]].itertuples(
            index=False, name=None
        )
    ]
    table = pd.concat([edges, _parity_rows(pairs, threads=threads)], axis=1)
    output_path = ccd_parity_path(data_dir)
    temporary_path = output_path.with_suffix(".parquet.tmp")
    table.to_parquet(temporary_path, index=False)
    temporary_path.replace(output_path)
    LOG.info(f"build_ccd_parity_scores: scored {len(table)} edges into {output_path}")
    return output_path


def query_ccd_parity(
    queries: Mapping[str, str],
    data_dir: Path,
    *,
    minimum_similarity: float = 30.0,
    threads: int = 4,
) -> pd.DataFrame:
    """PARITY-like scores of molecules outside the universe against their ECFP4 neighbours.

    Extends the precompiled CCD-vs-CCD table to ligands that only exist at
    dataset build time (novel molecules, composites through their whole-molecule
    SMILES): the same Tanimoto prefilter as :func:`query_ccd_tanimoto`, then the
    same score as :func:`build_ccd_parity_scores` on each hit.

    Returns
    -------
    pd.DataFrame
        :func:`query_ccd_tanimoto` columns plus ``parity_similarity`` (percent);
        best PARITY-like hit first.
    """
    hits = query_ccd_tanimoto(queries, data_dir, minimum_similarity=minimum_similarity)
    from plinder.core.utils.sanitize import mol_from_smiles

    smiles = load_ccd_components(data_dir).set_index("ligand_smiles_id")[
        "ligand_rdkit_canonical_smiles"
    ]
    mols = {value: mol_from_smiles(value) for value in set(hits["query_smiles"])}
    pairs = [
        (mols[query], mol_from_smiles(smiles[node]), False)
        for query, node in hits[["query_smiles", "ligand_smiles_id"]].itertuples(
            index=False, name=None
        )
    ]
    hits["parity_similarity"] = _parity_rows(pairs, threads=threads)[
        "parity_similarity"
    ]
    return hits.sort_values(
        ["query_id", "parity_similarity", "ligand_smiles_id"],
        ascending=[True, False, True],
        ignore_index=True,
    )


def make_ccd_ligand_dbs(
    *,
    data_dir: Path,
    scratch_dir: Path,
    threads: int = 4,
    force_update: bool = False,
    limit: int | None = None,
    build_tanimoto_scores: bool = True,
    minimum_similarity: float = 30.0,
) -> Path:
    """Build CCD-anchored MMP, ECFP4 and PARITY-like databases with cached settings.

    Parameters
    ----------
    data_dir : Path
        Release directory; artifacts land in ``<data_dir>/ccd_dbs``.
    scratch_dir : Path
        Working area for the mmpdb fragment/index run.
    threads : int
        Parallelism for the mmpdb fragmentation and indexing steps.
    force_update : bool
        Rebuild even when the manifest matches.
    limit : int or None
        Restrict to the first *limit* released components (testing aid).
    build_tanimoto_scores : bool
        Build the Tanimoto edges and their PARITY-like scores.
    minimum_similarity : float
        ECFP4 Tanimoto percentage cutoff for pairs retained in both score tables.

    Returns
    -------
    Path
        The ``ccd_dbs`` directory.

    Notes
    -----
    The manifest records the component universe, mmpdb version, and score
    settings. Matching builds reuse their files; changed settings rebuild them.
    """
    import json

    from plinder.data.annotations.mmpdb_utils import _mmpdb_version

    output_dir = ccd_dbs_dir(data_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    components_path = output_dir / "ccd_components.parquet"
    # Written into a ``fingerprints/`` subdirectory so the existing Tanimoto
    # scorer (which resolves ``<data_dir>/fingerprints/ligands_per_smiles.parquet``
    # internally) runs against this universe with ``data_dir=ccd_dbs`` and no
    # change to its signature.
    ecfp_path = ccd_fingerprint_path(data_dir)
    mmp_path = output_dir / "ccd_mmp_pairs.parquet"
    fragments_path = output_dir / "ccd_mmp_fragments.parquet"
    manifest_path = output_dir / "ccd_dbs.manifest.json"

    components = ccd_component_table(limit=limit)
    if components.empty:
        raise ValueError("the CCD component universe is empty")
    manifest = {
        "ccd_universe_signature": ccd_universe_signature(components),
        "num_components": int(len(components)),
        "mmpdb_version": _mmpdb_version(),
        "minimum_similarity": minimum_similarity,
        "build_tanimoto_scores": build_tanimoto_scores,
    }

    artifacts: tuple[Path, ...] = (
        components_path,
        ecfp_path,
        mmp_path,
        fragments_path,
        manifest_path,
    )
    if build_tanimoto_scores:
        artifacts += (ccd_parity_path(data_dir),)
    if not force_update and all(path.is_file() for path in artifacts):
        try:
            cached = json.loads(manifest_path.read_text())
        except (OSError, ValueError):
            cached = None
        if cached == manifest:
            LOG.info("make_ccd_ligand_dbs: artifacts already match the build settings")
            return output_dir

    temporary_components = components_path.with_suffix(".parquet.tmp")
    temporary_components.unlink(missing_ok=True)
    components.to_parquet(temporary_components, index=False)
    temporary_components.replace(components_path)

    build_ccd_ecfp_db(components, ecfp_path)
    build_ccd_mmp_db(
        components,
        mmp_path,
        scratch_dir=Path(scratch_dir),
        threads=threads,
        fragments_path=fragments_path,
    )
    if build_tanimoto_scores:
        build_ccd_tanimoto_scores(
            components,
            data_dir,
            minimum_similarity=minimum_similarity,
        )
        build_ccd_parity_scores(components, data_dir, threads=threads)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    LOG.info(f"make_ccd_ligand_dbs: built {len(components)} components -> {output_dir}")
    return output_dir


def ccd_code_sequence(ccd_code: str) -> list[str]:
    """Ordered component codes of a (possibly composite) ligand code.

    No CCD id contains ``-`` (ids are at most five characters), so splitting on
    the hyphen recovers the components unambiguously.
    """
    return [code for code in str(ccd_code).split("-") if code]


# Linear sequence identity is the weak point of composite comparison. It is
# correct for genuinely linear composites (peptides, oligonucleotides) but wrong
# in two ways for branched ones, chiefly glycans:
#   1. A glycan is a tree, not a chain, so linear identity cannot express branch
#      topology - two different topologies over the same residues score 1.0.
#   2. ``ccd_code`` is the chain's residue order, i.e. one arbitrary traversal.
#      The same glycan deposited with a different residue order linearises
#      differently and scores < 1.0 against itself.
# Prefer :func:`composite_similarity`, which follows the bond graph instead: a
# branched composite and the linear reading of its own code string are identical
# to this function and distinct to that one. This one stays for records that
# carry only a code string and no structure.
#
# TODO: neither comparison is chirality-aware, so stereoisomeric residues are
# interchangeable here (GLC/GAL/MAN all score 1.0 against each other). That is
# consistent with the rest of the stack - plinder's ECFP4 fingerprints set
# ``include_chirality: false`` and MMP fragmentation is likewise achiral - so
# making this one component stereo-aware would introduce an inconsistency rather
# than remove one. Revisit only as a stack-wide change: chirality-aware ECFP4
# does separate stereoisomeric hexoses, which the achiral form scores identical,
# so the signal is available if the whole similarity stack moves together. Note that MCES (RascalMCES) and distance-geometry bounds matrices
# were both measured stereo-blind, so neither is a shortcut to it.
def component_sequence_identity(query: Sequence[str], target: Sequence[str]) -> float:
    """Fractional identity between two component-code sequences.

    Parameters
    ----------
    query, target : sequence of str
        Ordered CCD component codes, e.g. ``["SER", "TPO", "GLY"]`` for a
        peptide, ``["DA", "DT", "DG"]`` for DNA, ``["NAG", "BMA"]`` for a glycan.

    Returns
    -------
    float
        Identity in ``[0, 1]``: matched positions over the longer sequence.

    Notes
    -----
    Applies to **any covalently-linked biopolymer composite** - peptide, RNA/DNA
    or oligosaccharide - because in every case the discriminating signal is the
    component sequence, not the chemistry:

    - mmpdb cuts rotatable single bonds. The peptide backbone amide and the
      nucleic-acid phosphodiester are both such bonds, so every peptide and
      every oligonucleotide shatters into the same repeating residue fragments.
    - ECFP4's radius-2 window is dominated by the repeating backbone motif
      (``-NH-CO-CH-`` or sugar-phosphate). Nucleic acids are *worse* than
      peptides here: with only four canonical nucleotides the fingerprint space
      is far more degenerate than with twenty amino acids.
    - Worst of all for sugars: plinder builds ECFP4 with
      ``include_chirality=false``, so hexose stereoisomers (glucose / galactose
      / mannose) have **identical connectivity and therefore identical
      fingerprints**. Nothing but the component code separates them.

    Alignment uses :class:`difflib.SequenceMatcher` over tokens, since residue
    codes are multi-character and the protein-sequence aligners do not apply.

    Caveat: this models a **linear** sequence. Branched glycans are trees
    flattened into residue order, so identity here is a weak proxy for them -
    exact-set or topology-aware comparison is more appropriate for branched
    composites.
    """
    query = list(query)
    target = list(target)
    if not query or not target:
        return 0.0
    matcher = SequenceMatcher(a=query, b=target, autojunk=False)
    matched = sum(block.size for block in matcher.get_matching_blocks())
    return matched / max(len(query), len(target))


def match_composite_by_sequence(
    ccd_code: str,
    reference_codes: Mapping[str, str],
    *,
    min_identity: float = 0.8,
) -> list[tuple[str, float]]:
    """Rank reference ligands whose component-code sequence matches the query.

    Works for any biopolymer composite (peptide, RNA/DNA, oligosaccharide); see
    :func:`component_sequence_identity` for why sequence beats chemistry here.

    Parameters
    ----------
    ccd_code : str
        Query ligand code (composite, hyphen-joined).
    reference_codes : mapping
        ``{reference_id: ccd_code}`` to compare against.
    min_identity : float
        Minimum fractional identity to report.

    Returns
    -------
    list of (reference_id, identity)
        Sorted by descending identity; empty when nothing clears the threshold.

    Notes
    -----
    The sequence arm of the matching strategy: name first, chemistry never. Use
    it instead of fingerprint or MMP similarity for any linked-residue
    composite.
    """
    query_sequence = ccd_code_sequence(ccd_code)
    if not query_sequence:
        return []
    matches: list[tuple[str, float]] = []
    for reference_id, reference_code in reference_codes.items():
        identity = component_sequence_identity(
            query_sequence, ccd_code_sequence(reference_code)
        )
        if identity >= min_identity:
            matches.append((reference_id, identity))
    matches.sort(key=lambda item: (-item[1], item[0]))
    return matches


class ResidueGraph(NamedTuple):
    """Rooted residue-level graph of a composite ligand.

    ``labels[i]`` is residue *i*'s CCD code, ``adjacency[i]`` the residues it is
    covalently bonded to, and ``roots`` the residues that only accept linkages.
    This is the connectivity ``ccd_code`` cannot express: the code string is one
    arbitrary traversal, so any comparison built on it is order-dependent and
    blind to branching.

    The rooting matters as much as the edges. Undirected, a linear and a
    branched trisaccharide are the same three-node path, so a comparison free to
    pick its own roots scores them identical. ``roots`` fixes the traversal
    chemically instead - see :func:`residue_graph_from_atoms`.

    ``donors[i]`` holds the residues donating a linkage into residue *i*, i.e.
    the edges directed by chemistry. It orients traversal of a macrocycle, which
    has no terminus to root at, and may be left empty for a graph built without
    structural information.

    ``atoms[i]`` is the heavy-atom names present in residue *i* and ``links``
    each inter-residue bond as ``(i, j, atom_in_i, atom_in_j)``; both empty for
    a graph built without structural information.
    """

    labels: tuple[str, ...]
    adjacency: tuple[frozenset[int], ...]
    roots: tuple[int, ...]
    donors: tuple[frozenset[int], ...] = ()
    atoms: tuple[frozenset[str], ...] = ()
    links: tuple[tuple[int, int, str, str], ...] = ()

    def donors_of(self, residue: int) -> frozenset[int]:
        """Residues donating a linkage into *residue*, empty when undirected."""
        return self.donors[residue] if self.donors else frozenset()

    @property
    def is_connected(self) -> bool:
        """Whether every residue is reachable from the first one."""
        if not self.labels:
            return True
        seen = {0}
        stack = [0]
        while stack:
            for neighbour in self.adjacency[stack.pop()]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        return len(seen) == len(self.labels)

    @property
    def is_tree(self) -> bool:
        edges = sum(len(neighbours) for neighbours in self.adjacency) // 2
        return edges == len(self.labels) - 1 and self.is_connected


def residue_graph_from_atoms(atoms: struc.AtomArray) -> ResidueGraph:
    """Build the rooted residue graph of a ligand from its bond graph.

    Parameters
    ----------
    atoms : struc.AtomArray
        Ligand atoms carrying ``bonds``, as produced by ``build_biounit``.

    Returns
    -------
    ResidueGraph
        Residues as nodes, inter-residue covalent bonds as edges, rooted at the
        residues that only accept linkages.

    Raises
    ------
    ValueError
        If the atoms carry no bond graph, which is the whole input here.

    Notes
    -----
    Linkages are directed by the chemistry of the bond rather than by residue
    numbering, which is not comparable across entries. In both families plinder
    sees, one side donates a carbon and the other accepts through a heteroatom:
    a glycosidic bond runs anomeric ``C1`` into the next ring's oxygen, a
    peptide bond runs carbonyl ``C`` into the next residue's nitrogen. Directing
    every linkage donor -> acceptor therefore roots the graph at the reducing
    end of a glycan and the C-terminus of a peptide, both canonical and both
    derived from the structure alone.

    Cyclic residue graphs (cyclic peptides, macrocycles) have no such terminus;
    every residue is then reported as a root, and comparison falls back to the
    best alignment over root pairs.
    """
    import biotite.structure as structure

    if atoms.bonds is None:
        raise ValueError("residue_graph_from_atoms requires bonded atoms")

    starts = structure.get_residue_starts(atoms, add_exclusive_stop=True)
    residue_of_atom: dict[int, int] = {}
    labels: list[str] = []
    for index, (start, stop) in enumerate(zip(starts[:-1], starts[1:])):
        labels.append(str(atoms.res_name[start]))
        for atom_index in range(int(start), int(stop)):
            residue_of_atom[atom_index] = index

    elements = [str(value).upper() for value in atoms.element]
    names = [str(value) for value in atoms.atom_name]
    present: list[set[str]] = [set() for _ in labels]
    for atom_index, residue in residue_of_atom.items():
        if elements[atom_index] not in ("H", "D"):
            present[residue].add(names[atom_index])
    adjacency: list[set[int]] = [set() for _ in labels]
    donors: list[set[int]] = [set() for _ in labels]
    donates: set[int] = set()
    links: list[tuple[int, int, str, str]] = []
    for first, second, _ in atoms.bonds.as_array():
        left = residue_of_atom.get(int(first))
        right = residue_of_atom.get(int(second))
        if left is None or right is None or left == right:
            continue
        adjacency[left].add(right)
        adjacency[right].add(left)
        links.append((left, right, names[int(first)], names[int(second)]))
        left_element, right_element = elements[int(first)], elements[int(second)]
        if left_element == "C" and right_element != "C":
            donates.add(left)
            donors[right].add(left)
        elif right_element == "C" and left_element != "C":
            donates.add(right)
            donors[left].add(right)

    roots = tuple(index for index in range(len(labels)) if index not in donates)
    graph = ResidueGraph(
        tuple(labels),
        tuple(frozenset(entry) for entry in adjacency),
        roots or tuple(range(len(labels))),
        tuple(frozenset(entry) for entry in donors),
        tuple(frozenset(entry) for entry in present),
        tuple(links),
    )
    if len(labels) > 1 and not graph.is_connected:
        # plinder builds ligands as connected components of the covalent graph,
        # so a disconnected result means unrelated residues were passed in
        # together and any similarity computed from it would be meaningless.
        LOG.warning(
            "residue graph is disconnected: %d residues (%s) with no single "
            "covalent component; these are separate ligands, not a composite",
            len(labels),
            "-".join(labels),
        )
    return graph


def ccd_code_residue_graph(ccd_code: str) -> ResidueGraph:
    """Residue graph implied by a hyphenated ``ccd_code``, assuming linearity.

    Fallback for records that carry only the code string. It reproduces the
    linear assumption of :func:`component_sequence_identity`, so prefer
    :func:`residue_graph_from_atoms` whenever the structure is available.
    """
    labels = ccd_code_sequence(ccd_code)
    adjacency = [
        frozenset(
            neighbour
            for neighbour in (index - 1, index + 1)
            if 0 <= neighbour < len(labels)
        )
        for index in range(len(labels))
    ]
    return ResidueGraph(tuple(labels), tuple(adjacency), (0,) if labels else ())


def _exact_label_similarity(first: str, second: str) -> float:
    return 1.0 if first == second else 0.0


def _rooted_children(graph: ResidueGraph, root: int) -> tuple[tuple[int, ...], ...]:
    """Spanning tree from *root*, as a child list per residue.

    Descends against the direction of donation, so from a residue the traversal
    reaches the residues that donate into it and never the one it donates to.
    On a tree this reproduces the graph exactly; on a macrocycle, where every
    residue both donates and accepts, it unrolls the ring into a path whose
    direction is set by the linkage chemistry rather than by residue numbering,
    which is what makes the score invariant to rotating the residue order.
    Edges whose direction is undetermined are traversed either way.
    """
    children: list[tuple[int, ...]] = [() for _ in graph.labels]
    seen = {root}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        descend = []
        for neighbour in sorted(graph.adjacency[node]):
            if neighbour in seen or node in graph.donors_of(neighbour):
                continue
            seen.add(neighbour)
            descend.append(neighbour)
            queue.append(neighbour)
        children[node] = tuple(descend)
    return tuple(children)


def _align_rooted(
    query: ResidueGraph,
    target: ResidueGraph,
    children_query: tuple[tuple[int, ...], ...],
    children_target: tuple[tuple[int, ...], ...],
    root_query: int,
    root_target: int,
    similarity: Callable[[str, str], float],
) -> tuple[float, list[tuple[int, int]]]:
    """Best alignment of two spanning trees: summed node score and aligned pairs.

    Each residue match is weighted by how well the two residues' linkage counts
    agree. The spanning tree necessarily drops any ring-closing edge, so without
    this a macrocycle and the open chain of the same residues would align
    perfectly; degree is where that edge survives.
    """
    from scipy.optimize import linear_sum_assignment

    memo: dict[tuple[int, int], tuple[float, list[tuple[int, int]]]] = {}

    def linkage_agreement(a: int, b: int) -> float:
        # Smoothed, so a free residue (degree 0) is merely a poor match for a
        # linked one rather than an impossible one. An unsmoothed ratio zeroes
        # that pair outright, which would score every mono-residue ligand 0.0
        # against every composite containing it - worse than the whole-molecule
        # fingerprint this is meant to improve on.
        degree_query = len(query.adjacency[a])
        degree_target = len(target.adjacency[b])
        return (1 + min(degree_query, degree_target)) / (
            1 + max(degree_query, degree_target)
        )

    def score(a: int, b: int) -> tuple[float, list[tuple[int, int]]]:
        cached = memo.get((a, b))
        if cached is not None:
            return cached
        total = similarity(query.labels[a], target.labels[b]) * linkage_agreement(a, b)
        pairs = [(a, b)]
        kids_query, kids_target = children_query[a], children_target[b]
        if kids_query and kids_target:
            scored = [[score(ca, cb) for cb in kids_target] for ca in kids_query]
            matrix = np.array([[value for value, _ in row] for row in scored])
            rows, columns = linear_sum_assignment(-matrix)
            for row, column in zip(rows, columns):
                total += float(matrix[row, column])
                pairs += scored[row][column][1]
        memo[(a, b)] = (total, pairs)
        return total, pairs

    return score(root_query, root_target)


def _best_alignment(
    query: ResidueGraph, target: ResidueGraph, similarity: Callable[[str, str], float]
) -> tuple[float, list[tuple[int, int]]]:
    """Best root pairing: summed node score and the one-to-one residue pairs.

    A single residue has no topology to orient, so against a composite it is
    tried at every residue rather than only at the roots.
    """
    if not query.labels or not target.labels:
        return 0.0, []
    roots_a = query.roots if len(target.labels) > 1 else range(len(query.labels))
    roots_b = target.roots if len(query.labels) > 1 else range(len(target.labels))
    return max(
        (
            _align_rooted(
                query,
                target,
                _rooted_children(query, root_a),
                _rooted_children(target, root_b),
                root_a,
                root_b,
                similarity,
            )
            for root_a in roots_a
            for root_b in roots_b
        ),
        key=lambda item: item[0],
    )


def composite_similarity(
    query: ResidueGraph,
    target: ResidueGraph,
    *,
    node_similarity: Callable[[str, str], float] | None = None,
) -> float:
    """Similarity of two composite ligands following their connectivity.

    A rooted tree alignment: residues match residues by a node score and
    subtrees are aligned recursively, each level resolved by an optimal (1:1)
    assignment between the two child sets. Unlike sequence identity this follows
    the actual bond graph, so branch topology is respected and the arbitrary
    residue order in ``ccd_code`` is irrelevant.

    This is *not* Feature Trees (Rarey & Dixon 1998), though it shares the
    tree-alignment shape and, with a graded ``node_similarity``, the fuzzy node
    comparison. FTrees' defining step is a split search in which one node may
    match a contiguous *set* of nodes; matching here is strictly one residue to
    one residue.

    That 1:1 restriction is a real limitation, not a free simplification. The
    CCD assigns single codes to molecules that are chemically composites
    (``LAT`` for the same molecule as ``GAL-BGC``, ``SUC`` for ``GLC-FRU``,
    ``GSH`` for a tripeptide), so the same chemistry does appear at two
    granularities, and this function scores such pairs 0.0 - it can only see
    residue labels, which share nothing across a regrouping.

    This function therefore measures topology agreement *given* a decomposition,
    and is only meaningful between ligands decomposed comparably. For a measure
    that does not assume that, use :func:`ligand_features` with
    :func:`ligand_similarity`, where molecular and residue descriptions share
    one feature space and the decomposition is not privileged.

    Parameters
    ----------
    query, target : ResidueGraph
        Rooted residue graphs, e.g. from :func:`residue_graph_from_atoms`.
    node_similarity : callable, optional
        ``(code_a, code_b) -> [0, 1]``. Defaults to exact CCD-code identity;
        pass :func:`ccd_node_similarity` to let chemically similar residues
        match partially, so ``A-TYR-B`` and ``A-PHE-B`` score as near-misses
        rather than as a plain mismatch. Note that a graded kernel also raises
        the floor for unrelated composites - every amino acid shares backbone
        chemistry - so thresholds tuned for exact matching are too permissive
        under it.

    Returns
    -------
    float
        Score in ``[0, 1]``, normalised by the larger residue count, so a
        subgraph match is penalised by the residues it leaves unmatched. Each
        residue match is additionally weighted by how well the two residues'
        linkage counts agree, which is what separates a macrocycle from the
        open chain of the same residues.

    Notes
    -----
    Only the roots recorded on each graph are tried, which is what keeps a
    linear and a branched composite distinguishable; a single residue, having
    no topology, is tried at every residue of the other graph. Residue graphs run to a
    handful of nodes, so the assignment at each level is negligible.

    Comparison runs over the spanning tree from each root (see
    :func:`_rooted_children`), so a macrocycle is handled like any other graph
    rather than recursing around its ring.
    """
    best, _ = _best_alignment(query, target, node_similarity or _exact_label_similarity)
    return best / max(len(query.labels), len(target.labels))


def _present_atoms(graph: ResidueGraph, residue: int) -> frozenset[str]:
    """Heavy atoms of a residue: the structure's when recorded, else the component's."""
    if graph.atoms:
        return graph.atoms[residue]
    return frozenset(_ccd_heavy_atoms(graph.labels[residue])[0])


def _present_bonds(graph: ResidueGraph, residue: int) -> int:
    """Component bonds whose both atoms are present in the residue."""
    present = _present_atoms(graph, residue)
    return sum(bond <= present for bond in _ccd_heavy_atoms(graph.labels[residue])[1])


def ligand_parity(
    query: ResidueGraph,
    target: ResidueGraph,
    table: CcdParityTable,
    query_mol: Mol,
    target_mol: Mol,
) -> float:
    """PARITY-like similarity of two ligands, assembled or whole-molecule as fits.

    Two composites are assembled from the mono table by :func:`composite_parity`.
    When either ligand is a single residue the whole molecules go through
    :func:`rascal_parity_score`: a mono may span several residues of the other
    side, which one-to-one residue alignment cannot express, and one molecule
    of residue size keeps the engine in its fast range.
    """
    from plinder.core.structure.smallmols_similarity import rascal_parity_score

    if min(len(query.labels), len(target.labels)) <= 1:
        return rascal_parity_score(query_mol, target_mol, stereo=False)
    return composite_parity(query, target, table)


def _rewired(mapping: dict[str, str], atom: str, image: str) -> dict[str, str]:
    """*mapping* with ``atom -> image``, the displaced images swapped to stay one-to-one."""
    rewired = dict(mapping)
    holder = next((p for p, q in mapping.items() if q == image), None)
    if holder is not None:
        if atom in mapping:
            rewired[holder] = mapping[atom]
        else:
            del rewired[holder]
    rewired[atom] = image
    return rewired


def _fragment_bonds(first: str, second: str, mapping: dict[str, str]) -> int:
    """Bonds of *first* inside the fragment whose images are bonds of *second*."""
    bonds_second = _ccd_heavy_atoms(second)[1]
    return sum(
        frozenset(mapping[p] for p in bond) in bonds_second
        for bond in _ccd_heavy_atoms(first)[1]
        if bond <= mapping.keys()
    )


def composite_parity(
    query: ResidueGraph, target: ResidueGraph, table: CcdParityTable
) -> float:
    """PARITY-like similarity of two composites assembled from mono-to-mono matches.

    Residues are aligned one-to-one with the table as node kernel; each aligned
    pair contributes its largest matched fragment, restricted to the atoms the
    structure actually holds (a linked residue has lost its leaving atom); a
    linkage counts, and joins the two fragments, only when both mappings carry
    its atoms onto the partner linkage. The assembled match is a common
    subgraph of the two molecules, so the score sits at or just below the
    whole-molecule :func:`rascal_parity_score`. No stereo term, as in PARITY.

    Residues match one-to-one, so a single residue that spans several residues
    of the other side (chitotriose against NAG-NAG-NAG) is credited for one of
    them; :func:`ligand_parity` routes such pairs to the whole-molecule engine.

    Where the table holds several equally good fragments for a residue pair,
    the one carrying the most linkage atoms onto the partner's linkage is taken,
    then the biggest, then the first in the table's fixed order, so the choice
    is deterministic and no atom is used twice on either side. A fragment that
    misses a linkage atom is also tried rewired onto it (a beta1-4 against a
    beta1-6 chain keeps the chain connected at the cost of the two displaced
    hydroxyl bonds, as the whole-molecule engine does).
    """
    _, pairs = _best_alignment(query, target, table.similarity)
    partner = dict(pairs)
    target_links = {(i, j): (p, q) for i, j, p, q in target.links}
    target_links |= {(j, i): (q, p) for i, j, p, q in target.links}
    # link atoms each residue should carry onto its partner's linkage
    wanted: dict[int, dict[str, str]] = {a: {} for a, _ in pairs}
    for i, j, p, q in query.links:
        link = target_links.get((partner.get(i, -1), partner.get(j, -1)))
        if link:
            wanted[i][p], wanted[j][q] = link
    fragments: dict[int, dict[str, str]] = {}
    bonds = 0
    for a, b in pairs:
        present_a, present_b = _present_atoms(query, a), _present_atoms(target, b)
        candidates = [
            {p: q for p, q in mapping.items() if p in present_a and q in present_b}
            for mapping in table.mappings(query.labels[a], target.labels[b])
        ] or [{}]
        candidates += [
            _rewired(mapping, p, r)
            for mapping in list(candidates)
            for p, r in wanted[a].items()
            if mapping.get(p) != r and p in present_a and r in present_b
        ]
        scored = [
            (
                sum(mapping.get(p) == r for p, r in wanted[a].items()),
                len(mapping),
                _fragment_bonds(query.labels[a], target.labels[b], mapping),
            )
            for mapping in candidates
        ]
        best = max(range(len(candidates)), key=lambda k: scored[k])
        fragments[a] = candidates[best]
        bonds += scored[best][2]
    # a matched linkage adds its bond and merges the two fragments it joins
    joined = {a: a for a in fragments}

    def root(a: int) -> int:
        while joined[a] != a:
            a = joined[a]
        return a

    for i, j, p, q in query.links:
        if i not in partner or j not in partner:
            continue
        link = target_links.get((partner[i], partner[j]))
        if link and fragments[i].get(p) == link[0] and fragments[j].get(q) == link[1]:
            bonds += 1
            joined[root(i)] = root(j)
    merged: Counter[int] = Counter()
    for a, mapping in fragments.items():
        merged[root(a)] += len(mapping)
    graphs = (query, target)
    atom_total = sum(
        len(_present_atoms(g, i)) for g in graphs for i in range(len(g.labels))
    )
    bond_total = sum(_present_bonds(g, i) for g in graphs for i in range(len(g.labels)))
    bond_total += len(query.links) + len(target.links)
    return parity_similarity(
        max(merged.values(), default=0), bonds, atom_total, bond_total
    )


@lru_cache(maxsize=8192)
def _ccd_component_fingerprint(ccd_id: str) -> Any | None:
    """ECFP4/1024 of one CCD component, cached; None when underivable.

    Uses plinder's own ``mol2morgan_fp`` at the native ECFP4 radius and size, so
    residue-level similarity here is the same measure the ligand similarity
    stack uses elsewhere.
    """
    from plinder.core.structure.smallmols_similarity import mol2morgan_fp
    from plinder.data.annotations.get_similarity_scores import (
        ECFP4_NBITS,
        ECFP4_RADIUS,
    )

    smiles = _ccd_smiles(ccd_id)
    if not smiles:
        return None
    try:
        fingerprint = mol2morgan_fp(smiles, radius=ECFP4_RADIUS, nbits=ECFP4_NBITS)
    except Exception as exc:  # noqa: BLE001 - one bad component must not abort
        LOG.warning(f"CCD fingerprint failed for {ccd_id}: {exc}")
        return None
    return cast("object | None", fingerprint)


def ccd_node_similarity(first: str, second: str) -> float:
    """Graded similarity between two CCD components, for use as a node kernel.

    Pass to :func:`composite_similarity` to let chemically similar residues
    match partially instead of all-or-nothing, so ``A-TYR-B`` and ``A-PHE-B``
    score as near-misses rather than as unrelated composites.

    Parameters
    ----------
    first, second : str
        CCD component codes.

    Returns
    -------
    float
        ``1.0`` for identical codes, otherwise the ECFP4/1024 Tanimoto of the
        two components. ``0.0`` when either has no derivable chemistry, which
        keeps an unresolvable residue from inventing similarity.

    Notes
    -----
    Achiral, like the rest of plinder's fingerprint stack, so stereoisomeric
    residues score 1.0 against each other - see the TODO above
    :func:`component_sequence_identity`.
    """
    from rdkit import DataStructs

    if first == second:
        return 1.0
    left, right = _ccd_component_fingerprint(first), _ccd_component_fingerprint(second)
    if left is None or right is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(left, right))


def residue_graph_features(graph: ResidueGraph, *, iterations: int = 2) -> Counter[str]:
    """Weisfeiler-Lehman features of a residue graph.

    Each residue's label is refined by the labels of the residues donating into
    it and those it donates to, repeatedly; every label seen at every round is a
    feature. Features depend only on labels and linkage, never on residue
    numbering, so relabelling a ring or reordering a ``ccd_code`` cannot change
    them, while a different cyclic order or branch topology does.

    Parameters
    ----------
    graph : ResidueGraph
        Residue graph, from :func:`residue_graph_from_atoms`.
    iterations : int
        Refinement rounds. Two reaches each residue's neighbours-of-neighbours,
        which is the residue-level analogue of ECFP4's radius 2.

    Returns
    -------
    Counter of str
        Feature strings, ``r:``-prefixed to share a space with molecular bits,
        counted by how many residues carry each. Multiplicity is kept because
        residue composition is exact and unhashed: a decapeptide genuinely has
        two alanines where a pentapeptide has one, and a presence-only feature
        set scores those two molecules identical.
    """
    labels = dict(enumerate(graph.labels))
    features = Counter(f"r:{label}" for label in labels.values())
    for _ in range(iterations):
        refined: dict[int, str] = {}
        for node, label in labels.items():
            donors = sorted(labels[other] for other in graph.donors_of(node))
            accepts = sorted(
                labels[other]
                for other in graph.adjacency[node]
                if node in graph.donors_of(other)
            )
            undirected = sorted(
                labels[other]
                for other in graph.adjacency[node]
                if other not in graph.donors_of(node)
                and node not in graph.donors_of(other)
            )
            refined[
                node
            ] = f"{label}<{','.join(donors)}>{','.join(accepts)}~{','.join(undirected)}"
        labels = refined
        features.update(f"r:{label}" for label in labels.values())
    return features


def ligand_features(smiles: str | None, graph: ResidueGraph | None) -> Counter[str]:
    """The feature set describing one ligand: molecular bits plus residue graph.

    Both descriptions live in one space, so similarity is a single Tanimoto over
    whatever each ligand has. Nothing branches on how a ligand was decomposed:
    the molecular bits are identical for the same chemistry however it was
    carved into residues, and the residue features add the topology that
    fingerprints cannot see.

    Parameters
    ----------
    smiles : str or None
        Whole-molecule SMILES (plinder's ``ligand_rdkit_canonical_smiles``).
    graph : ResidueGraph or None
        Residue graph, when the structure is available.

    Returns
    -------
    Counter of str
        ``m:``-prefixed molecular bits and ``r:``-prefixed residue features.

    Notes
    -----
    Molecular bits are kept presence-only, matching plinder's ECFP4 tables and
    :func:`~plinder.data.annotations.get_similarity_scores.ligand_scores`, so
    this measure stays consistent with every other similarity in the codebase.
    Counting them was measured and rejected: a larger molecule's multiplicities
    inflate the denominator, sharpening the penalty on any size difference, and
    since folded bits outnumber residue features it also dilutes the residue
    signal that separates a resequenced macrocycle. Residue features are counted
    instead, where multiplicity is exact rather than hashed.
    """
    from plinder.core.structure.smallmols_similarity import mol2morgan_fp
    from plinder.data.annotations.get_similarity_scores import (
        ECFP4_NBITS,
        ECFP4_RADIUS,
    )

    features: Counter[str] = Counter()
    if smiles:
        try:
            fingerprint = mol2morgan_fp(smiles, radius=ECFP4_RADIUS, nbits=ECFP4_NBITS)
            features.update({f"m:{bit}": 1 for bit in fingerprint.GetOnBits()})
        except Exception as exc:  # noqa: BLE001 - a bad SMILES must not abort
            LOG.warning(f"molecular features failed for {smiles!r}: {exc}")
    if graph is not None and graph.labels:
        features.update(residue_graph_features(graph))
    return features


def ligand_similarity(query: Counter[str], target: Counter[str]) -> float:
    """Similarity between two :func:`ligand_features` sets.

    ``sum(min) / sum(max)`` over the pooled features - the count-aware form of
    Tanimoto, which it equals exactly when every count is 1. Counts only enter
    through the residue features, so a molecule and a tandem repeat of it are no
    longer identical, while any pair whose residue multiplicities already agree
    scores exactly as it did presence-only.
    """
    if not query or not target:
        return 0.0
    shared = set(query) | set(target)
    lower = sum(min(query.get(key, 0), target.get(key, 0)) for key in shared)
    upper = sum(max(query.get(key, 0), target.get(key, 0)) for key in shared)
    return lower / upper if upper else 0.0
