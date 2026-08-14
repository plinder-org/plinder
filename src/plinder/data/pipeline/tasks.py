# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import hashlib
import json
import os
import time
from concurrent.futures import (
    ALL_COMPLETED,
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)
from pathlib import Path
from shutil import copyfile, rmtree
from string import ascii_lowercase, digits
from textwrap import dedent
from typing import Any, Sequence

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from omegaconf import DictConfig
from tqdm import tqdm

from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger
from plinder.data import clusters, databases
from plinder.data.annotations import get_similarity_scores, mmpdb_utils
from plinder.data.pipeline import collate, io, utils
from plinder.data.pipeline.ingest import (
    balance_entries,
    completed_entry_metrics,
    completed_interface_metrics,
    discover_entries,
    ingest_pdb_batch,
    normalize_ingest_mode,
    normalize_pdb_id,
)

LOG = setup_logger(__name__)
ALIGNMENT_CHAIN_LOOKUP_RELATIVE = Path("index/alignment_chain_lookup.parquet")
ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE = Path(
    "index/alignment_chain_lookup.manifest.json"
)
INTERFACE_HALF_REPRESENTATIVES_RELATIVE = Path(
    "index/interface_half_representatives.parquet"
)
INTERFACE_REPRESENTATIVES_RELATIVE = Path("index/interface_representatives.parquet")
INTERFACE_MEMBERSHIP_RELATIVE = Path("index/interface_membership.parquet")
INTERFACE_REPRESENTATIVES_MANIFEST_RELATIVE = Path(
    "index/interface_representatives.manifest.json"
)
LIGAND_POCKET_REPRESENTATIVES_RELATIVE = Path(
    "index/ligand_pocket_representatives.parquet"
)
LIGAND_POCKET_MEMBERSHIP_RELATIVE = Path("index/ligand_pocket_membership.parquet")
LIGAND_POCKET_REPRESENTATIVES_MANIFEST_RELATIVE = Path(
    "index/ligand_pocket_representatives.manifest.json"
)
STAGES = [
    "download_rcsb_files",
    "download_alternative_datasets",
    "make_entries",
    "collate_entries",
    "make_dbs",
    "make_canonical_ligand_archives",
    "finalize_ligand_archives",
    "compute_ligand_fingerprints",
    "make_ligand_scores",
    "annotate_ligand_similarity",
    "make_ligand_mmp_pairs",
    "make_sub_dbs",
    "run_batch_searches",
    "map_batch_alignments",
    "collate_alignments",
    "finalize_alignments",
    "plan_interface_scores",
    "make_interface_scores",
    "finalize_interface_scores",
    "make_batch_scores",
    "collate_ligand_3d_candidates",
    "plan_ligand_3d_scores",
    "make_ligand_3d_scores",
    "collate_ligand_3d_scores",
    "merge_ligand_3d_scores",
    "finalize_scores",
    "export_sucos_shape_pocket_qcov",
    "finalize_sucos_export",
    "collate_partitions",
    "make_linked_apo_structures",
    "plan_clusters",
    "make_symmetric_edge_fragments",
    "make_symmetric_edge_shards",
    "make_component_reductions",
    "merge_component_reductions",
    "make_set_covers",
    "make_directed_set_covers",
    "summarize_clusters",
    "finalize_index",
]


def _file_content_signature(path: Path) -> dict[str, int | str]:
    """Return a stable signature for a createdb file manifest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {"size": path.stat().st_size, "sha256": digest.hexdigest()}


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text())
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def scatter_download_rcsb_files(
    *,
    data_dir: Path,
    batch_size: int,
    two_char_codes: list[str],
) -> list[list[str]]:
    """
    Split the task of rsyncing to RCSB along
    the middle two character codes in PDB IDs.
    Check both CIFs and validation reports.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    batch_size : int
        how many codes to put in a chunk
    two_char_codes : list[str], default=[]
        only consider particular codes

    Returns
    -------
    codes : list[list[str]]
        list of lists of chunks of two character codes
    """
    if len(two_char_codes):
        codes = sorted(two_char_codes)
    else:
        codes = io.list_rcsb(kind="cif")
        cif_codes = io.get_missing_two_char_codes(
            kind="cif",
            data_dir=data_dir / "ingest",
            two_char_codes=codes,
        )
        codes = io.list_rcsb(kind="val")
        val_codes = io.get_missing_two_char_codes(
            kind="val",
            data_dir=data_dir / "reports",
            two_char_codes=codes,
        )
        codes = sorted(set(cif_codes).union(val_codes))
    LOG.info(f"scatter_download_rcsb_files: found {len(codes)} two character codes")
    return [codes[pos : pos + batch_size] for pos in range(0, len(codes), batch_size)]


def download_rcsb_files(
    *,
    data_dir: Path,
    two_char_codes: list[str],
) -> None:
    """
    Download both CIF files and validation reports
    for a list of two character codes.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    two_char_codes : list[str]
        the list of two character codes to download
    """
    for code in two_char_codes:
        LOG.info(f"downloading cifs for two_char_code={code}")
        io.rsync_rcsb(
            kind="cif",
            data_dir=data_dir / "ingest",
            two_char_code=code,
        )
        LOG.info(f"downloading reports for two_char_code={code}")
        io.rsync_rcsb(
            kind="val",
            data_dir=data_dir / "reports",
            two_char_code=code,
        )


def download_alternative_datasets(
    *,
    data_dir: Path,
    threads: int,
    force_update: bool,
) -> None:
    """
    Download the alternative datasets that don't fit neatly
    into the rcsb rsync setup used downstream in entry and
    annotation generation.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    force_update : bool
        if True, force re-download

    """
    kws = dict(data_dir=data_dir, force_update=force_update)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures: list[Future[Any]] = [
            executor.submit(io.download_cofactors, **kws),
            executor.submit(io.download_seqres_data, **kws),
            executor.submit(io.download_components_cif, **kws),
            executor.submit(io.download_affinity_data, **kws),
        ]
        wait(futures, return_when=ALL_COMPLETED)
        for future in futures:
            exc = future.exception()
            if exc is not None:
                raise exc


def make_dbs(
    *,
    data_dir: Path,
    sub_databases: list[str],
    cpu: int,
    cif_root: Path | None = None,
    seqres_path: Path | None = None,
    scratch_dir: Path | None = None,
    create: bool = True,
    index: bool = True,
    force_update: bool = False,
    build_dir: Path | None = None,
) -> None:
    """
    Make the foldseek and mmseqs dbs

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    """
    input_dirs = {}
    if "apo" in sub_databases or "holo" in sub_databases:
        input_dirs["foldseek"] = cif_root or data_dir / "ingest"
        input_dirs["mmseqs"] = seqres_path or io.download_seqres_data(data_dir=data_dir)
        if create:
            for source in input_dirs.values():
                if not source.exists():
                    raise FileNotFoundError(f"missing database source: {source}")
    if "pred" in sub_databases:
        input_dirs["pred_mmseqs"] = io.download_uniprot_fasta_data(data_dir=data_dir)
        input_dirs["pred_foldseek"] = io.download_alphafold_cif_files(data_dir=data_dir)
    for db, source in input_dirs.items():
        output_dir = data_dir / "dbs" / db
        working_output_dir = build_dir / db if build_dir is not None else output_dir
        LOG.info(f"make_dbs: making {db} in {working_output_dir}")
        database_type = db.split("_")[-1]
        tmp_dir = (scratch_dir or data_dir / "scratch" / "databases") / db
        tmp_dir.mkdir(exist_ok=True, parents=True)
        database_path = output_dir / database_type
        complete = databases.created_database_is_complete(database_path, database_type)
        source_signature = (
            _file_content_signature(source)
            if database_type == "foldseek" and source.is_file()
            else None
        )
        input_marker = output_dir / f"{database_type}.createdb-input.json"
        source_is_current = (
            source_signature is None or _read_json(input_marker) == source_signature
        )
        rebuild = force_update or not complete or not source_is_current
        if create and rebuild:
            working_marker = working_output_dir / f"{database_type}.createdb-input.json"
            working_marker.unlink(missing_ok=True)
            databases.create_db(
                source,
                working_output_dir,
                database_type,
                threads=cpu,
            )
            if source_signature is not None:
                _write_json_atomic(working_marker, source_signature)
        elif create:
            LOG.info(f"make_dbs: reusing completed {database_path}")
        if index:
            index_output_dir = working_output_dir if rebuild else output_dir
            databases.create_db_index(
                index_output_dir,
                database_type,
                tmp_dir=tmp_dir,
                threads=cpu,
            )
        if build_dir is not None and create and rebuild:
            working_database = working_output_dir / database_type
            if not databases.created_database_is_complete(
                working_database, database_type
            ):
                raise RuntimeError(f"createdb output is incomplete: {working_database}")
            databases.install_database_directory(working_output_dir, output_dir)


def scatter_make_entries(
    *,
    data_dir: Path,
    cif_root: Path,
    validation_root: Path,
    batch_size: int,
    two_char_codes: list[str],
    pdb_ids: list[str],
    force_update: bool,
    discovery_threads: int = 8,
    interface_min_residues: int = 7,
    interface_annotate_prodigy: bool = True,
    ingest_mode: str = "all",
) -> list[list[str]]:
    """Discover and size-balance source entries for V3 annotation."""
    ingest_mode = normalize_ingest_mode(ingest_mode)
    selected_pdb_ids = [normalize_pdb_id(pdb_id) for pdb_id in pdb_ids]
    selected_codes = _selected_context_codes(two_char_codes, selected_pdb_ids)
    entries = discover_entries(
        cif_root,
        validation_root,
        check_validation=False,
        threads=discovery_threads,
        two_char_codes=selected_codes,
        pdb_ids=selected_pdb_ids,
    )
    if not force_update:
        entries = [
            entry
            for entry in entries
            if (
                completed_interface_metrics(
                    data_dir,
                    entry.pdb_id,
                    expected_interface_min_residues=interface_min_residues,
                    expected_annotate_prodigy=interface_annotate_prodigy,
                )
                if ingest_mode == "interfaces"
                else completed_entry_metrics(
                    data_dir,
                    entry.pdb_id,
                    expected_interface_min_residues=(
                        interface_min_residues if ingest_mode == "all" else None
                    ),
                    expected_annotate_prodigy=(
                        interface_annotate_prodigy if ingest_mode == "all" else None
                    ),
                    expected_ingest_mode=ingest_mode,
                )
            )
            is None
        ]
    LOG.info(f"scatter_make_entries: found {len(entries)} PDBs in {cif_root}")
    return [
        [entry.pdb_id for entry in batch]
        for batch in balance_entries(entries, batch_size=batch_size)
    ]


def make_entries(
    *,
    data_dir: Path,
    pdb_ids: list[str],
    cif_root: Path,
    validation_root: Path,
    force_update: bool,
    annotation_cfg: DictConfig,
    entry_cfg: DictConfig,
    interface_cfg: DictConfig,
    cpu: int = 1,
    ingest_mode: str = "all",
) -> list[str]:
    """Run the same resumable V3 batch implementation used by Slurm."""
    del cpu  # Entry annotation is intentionally sequential within each worker.
    normalized_ids = [normalize_pdb_id(pdb_id) for pdb_id in pdb_ids]
    hash_id = utils.hash_contents(normalized_ids)
    metrics_path, _ = ingest_pdb_batch(
        pdb_ids=normalized_ids,
        output_root=data_dir,
        cif_root=cif_root,
        validation_root=validation_root,
        force=force_update,
        job_id=f"metaflow-{hash_id}",
        annotation_cfg=annotation_cfg,
        entry_cfg=entry_cfg,
        interface_cfg=interface_cfg,
        mode=normalize_ingest_mode(ingest_mode),
    )
    payload = json.loads(metrics_path.read_text())
    failed = [
        str(entry["pdb_id"])
        for entry in payload["entries"]
        if entry["status"] == "failed"
    ]
    LOG.info(f"would rerun {len(failed)} entries")
    return failed


def _selected_context_codes(two_char_codes: list[str], pdb_ids: list[str]) -> list[str]:
    """Apply context precedence consistently across entry-derived stages."""
    if pdb_ids:
        return sorted({normalize_pdb_id(pdb_id)[1:3] for pdb_id in pdb_ids})
    return sorted(str(code).lower() for code in two_char_codes)


def scatter_collate_entries(
    *,
    data_dir: Path,
    batch_size: int,
) -> list[list[str]]:
    """Plan V3 collation and scatter deterministic two-character shards."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    plan = collate.plan_collation(data_dir)
    codes = [str(code) for code in plan["codes"]]
    return [codes[pos : pos + batch_size] for pos in range(0, len(codes), batch_size)]


def collate_entries(
    *,
    data_dir: Path,
    two_char_codes: list[str],
    cpu: int,
    memory_limit: str,
) -> None:
    """Collate one batch of V3 entry shards through the shared implementation."""
    for code in two_char_codes:
        collate.collate_shard(
            data_dir,
            code,
            threads=cpu,
            memory_limit=memory_limit,
            scratch_dir=data_dir / "scratch" / "collation",
        )


def finalize_entry_collation(
    *,
    data_dir: Path,
    cpu: int,
    memory_limit: str,
) -> dict[str, Any]:
    """Validate and fail-closed install the sharded V3 annotation index."""
    return collate.finalize_collation(
        data_dir,
        threads=cpu,
        memory_limit=memory_limit,
        scratch_dir=data_dir / "scratch" / "collation-finalize",
    )


def scatter_make_canonical_ligand_archives(
    *,
    data_dir: Path,
    two_char_codes: list[str],
    pdb_ids: list[str],
    batch_size: int,
) -> list[list[str]]:
    """
    Scatter two-character codes for canonical ligand archive generation.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    batch_size : int
        how many codes to put in a chunk
    two_char_codes : list[str], default=[]
        only consider particular codes
    pdb_ids : list[str], default=[]
        if set, derive codes from these IDs and ignore ``two_char_codes``

    Returns
    -------
    chunks : list[list[str]]
        batches of two character codes
    """
    entry_dir = data_dir / "raw_entries"
    selected_codes = _selected_context_codes(two_char_codes, pdb_ids)
    if selected_codes:
        codes = selected_codes
    else:
        codes = sorted(os.listdir(entry_dir.as_posix()))
    LOG.info(
        "scatter_make_canonical_ligand_archives: "
        f"found {len(codes)} two character codes"
    )
    return [codes[pos : pos + batch_size] for pos in range(0, len(codes), batch_size)]


def make_canonical_ligand_archives(
    *,
    data_dir: Path,
    two_char_codes: list[str],
    scratch_dir: Path | None = None,
) -> None:
    """Pack canonical ASU ligand SDFs into one queryable Parquet per shard."""
    schema = pa.schema(
        [
            pa.field("pdb_id", pa.string(), nullable=False),
            pa.field("ligand_asym_id", pa.string(), nullable=False),
            pa.field("sdf", pa.binary(), nullable=False),
        ]
    )
    for code in two_char_codes:
        entry_dir = data_dir / "raw_entries" / code
        archive = data_dir / "ligand_archives" / f"{code}.parquet"
        archive.parent.mkdir(exist_ok=True, parents=True)
        records = []
        for entry_parquet in sorted(entry_dir.glob("*.parquet")):
            ligand_dir = entry_dir / entry_parquet.stem / "ligand_files"
            for ligand_file in sorted(ligand_dir.glob("*.sdf")):
                records.append(
                    {
                        "pdb_id": entry_parquet.stem,
                        "ligand_asym_id": ligand_file.stem,
                        "sdf": ligand_file.read_bytes(),
                    }
                )
        table = pa.Table.from_pylist(records, schema=schema)
        temporary_root = scratch_dir or archive.parent
        temporary_root.mkdir(exist_ok=True, parents=True)
        temporary = temporary_root / f"{code}.parquet"
        if temporary == archive:
            temporary = archive.with_suffix(".tmp.parquet")
        pq.write_table(
            table,
            temporary,
            compression="zstd",
            compression_level=6,
            use_dictionary=["pdb_id", "ligand_asym_id"],
            row_group_size=2_048,
            write_statistics=True,
        )
        install = archive.with_suffix(".tmp.parquet")
        if temporary != install:
            copyfile(temporary, install)
            temporary.unlink(missing_ok=True)
        install.replace(archive)


def make_sub_dbs(
    *,
    data_dir: Path,
    sub_databases: list[str],
    cpu: int = 1,
    scratch_dir: Path | None = None,
) -> None:
    """
    Get the list of all pdb IDs to load all the entries
    for full sub-database generation context. Explicitly
    don't support two_char_codes forwarding to get_local_contents
    to avoid an issue where sub dbs are created without full
    context.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    """
    entries = None
    identifiers_by_database = None
    if set(sub_databases) <= {"holo", "apo"}:
        identifiers_by_database = {}
        for search_db in sub_databases:
            chains = (
                _protein_scoring_chains(data_dir)
                if search_db == "holo"
                else _apo_scoring_chains(data_dir)
            )
            if search_db == "holo":
                chain_auth_ids = pd.read_parquet(
                    data_dir / "index" / "entry_chains.parquet",
                    columns=["entry_pdb_id", "chain_asym_id", "chain_auth_id"],
                )
                chains = chains.merge(
                    chain_auth_ids,
                    on=["entry_pdb_id", "chain_asym_id"],
                    how="left",
                    validate="one_to_one",
                )
            chains = chains[chains["chain_auth_id"].notna()]
            identifiers_by_database[f"{search_db}_foldseek"] = {
                f"pdb_0000{row.entry_pdb_id}_xyz-enrich_{row.chain_auth_id}"
                for row in chains.itertuples(index=False)
            }
            identifiers_by_database[f"{search_db}_mmseqs"] = {
                f"{row.entry_pdb_id}_{row.chain_auth_id}"
                for row in chains.itertuples(index=False)
            }
    else:
        from plinder.core.scores.entries import entry_views_from_df

        entries = entry_views_from_df(
            pd.read_parquet(data_dir / "index" / "annotation_table.parquet"),
            entry_chains=pd.read_parquet(data_dir / "index" / "entry_chains.parquet"),
        )
    db_dir = data_dir / "dbs" / "subdbs"
    db_dir.mkdir(exist_ok=True)
    LOG.info("making sub-databases for scoring")
    db_sources = utils.get_db_sources(data_dir=data_dir, sub_databases=sub_databases)
    databases.make_sub_dbs(
        db_dir,
        db_sources,
        entries,
        identifiers_by_database=identifiers_by_database,
        tmp_dir=scratch_dir,
        threads=cpu,
    )
    if set(sub_databases).intersection({"holo", "apo"}):
        make_alignment_chain_lookup(
            data_dir=data_dir,
            scratch_dir=scratch_dir,
            threads=cpu,
        )


def _interface_representative_source_signature(
    data_dir: Path,
) -> dict[str, int | str]:
    source = data_dir / "index" / "interface_annotation_table.parquet"
    stat = source.stat()
    return {
        "name": source.name,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _interface_representative_output_signature(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {"name": path.name, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _completed_interface_representatives(
    data_dir: Path,
) -> dict[str, Any] | None:
    manifest = _read_json(data_dir / INTERFACE_REPRESENTATIVES_MANIFEST_RELATIVE)
    expected = {
        "half_interfaces": (
            INTERFACE_HALF_REPRESENTATIVES_RELATIVE,
            schemas.INTERFACE_HALF_REPRESENTATIVE_SCHEMA,
        ),
        "interfaces": (
            INTERFACE_REPRESENTATIVES_RELATIVE,
            schemas.INTERFACE_REPRESENTATIVE_SCHEMA,
        ),
        "membership": (
            INTERFACE_MEMBERSHIP_RELATIVE,
            schemas.INTERFACE_MEMBERSHIP_SCHEMA,
        ),
    }
    try:
        if manifest is None or manifest.get(
            "source"
        ) != _interface_representative_source_signature(data_dir):
            return None
        for key, (relative, schema) in expected.items():
            path = data_dir / relative
            if not pq.read_schema(path).equals(schema) or manifest.get(
                "outputs", {}
            ).get(key) != _interface_representative_output_signature(path):
                return None
    except (OSError, TypeError, ValueError):
        return None
    return manifest


def make_interface_representatives(
    *,
    data_dir: Path,
    scratch_dir: Path | None,
    threads: int,
    force_update: bool = False,
) -> dict[str, Any]:
    """Normalize exact assembly-copy interfaces before mapping and scoring."""
    current = _completed_interface_representatives(data_dir)
    if not force_update and current is not None:
        LOG.info(
            "make_interface_representatives: reusing %d representatives for "
            "%d interfaces",
            current["representative_interface_count"],
            current["interface_count"],
        )
        return current

    import duckdb

    source_signature = _interface_representative_source_signature(data_dir)
    started = time.monotonic()
    LOG.info("make_interface_representatives: normalizing interface masks")
    source = data_dir / "index" / "interface_annotation_table.parquet"
    working_root = (scratch_dir or data_dir / "scratch") / "interface-representatives"
    if working_root.exists():
        rmtree(working_root)
    working_root.mkdir(parents=True)
    temporary_paths = {
        "half_interfaces": working_root / INTERFACE_HALF_REPRESENTATIVES_RELATIVE.name,
        "interfaces": working_root / INTERFACE_REPRESENTATIVES_RELATIVE.name,
        "membership": working_root / INTERFACE_MEMBERSHIP_RELATIVE.name,
    }
    connection = duckdb.connect()
    connection.sql(f"SET threads={max(1, threads)}")
    connection.sql("SET memory_limit='64GB'")
    connection.sql(f"SET temp_directory='{working_root.as_posix()}'")
    connection.sql("SET preserve_insertion_order=false")
    connection.sql(
        dedent(
            f"""
            CREATE TEMP TABLE interface_sides AS
            SELECT
                entry_pdb_id,
                system_id,
                1::TINYINT AS side,
                interface_chain_1 AS instance_chain_id,
                split_part(interface_chain_1, '.', 2) AS chain_asym_id,
                interface_chain_1_residue_numbers AS residue_numbers,
                interface_chain_1_residue_indices AS residue_indices
            FROM read_parquet('{source.as_posix()}')
            UNION ALL
            SELECT
                entry_pdb_id,
                system_id,
                2::TINYINT AS side,
                interface_chain_2 AS instance_chain_id,
                split_part(interface_chain_2, '.', 2) AS chain_asym_id,
                interface_chain_2_residue_numbers AS residue_numbers,
                interface_chain_2_residue_indices AS residue_indices
            FROM read_parquet('{source.as_posix()}');

            CREATE TEMP TABLE half_representatives AS
            SELECT
                min(system_id || '::side=' || side::VARCHAR)::VARCHAR
                    AS half_interface_id,
                entry_pdb_id,
                arg_min(
                    instance_chain_id,
                    system_id || '::side=' || side::VARCHAR
                )::VARCHAR AS instance_chain_id,
                chain_asym_id,
                residue_numbers,
                residue_indices
            FROM interface_sides
            GROUP BY
                entry_pdb_id,
                chain_asym_id,
                residue_numbers,
                residue_indices;

            CREATE TEMP TABLE side_membership AS
            SELECT
                sides.entry_pdb_id,
                sides.system_id,
                sides.side,
                representatives.half_interface_id
            FROM interface_sides AS sides
            INNER JOIN half_representatives AS representatives
              ON sides.entry_pdb_id = representatives.entry_pdb_id
             AND sides.chain_asym_id = representatives.chain_asym_id
             AND sides.residue_numbers = representatives.residue_numbers
             AND sides.residue_indices = representatives.residue_indices;

            CREATE TEMP TABLE interface_members AS
            SELECT
                entry_pdb_id,
                system_id,
                max(half_interface_id) FILTER (WHERE side = 1)
                    AS side_1_half_interface_id,
                max(half_interface_id) FILTER (WHERE side = 2)
                    AS side_2_half_interface_id
            FROM side_membership
            GROUP BY entry_pdb_id, system_id;

            CREATE TEMP TABLE representative_groups AS
            SELECT
                entry_pdb_id,
                least(side_1_half_interface_id, side_2_half_interface_id)
                    AS half_interface_1_id,
                greatest(side_1_half_interface_id, side_2_half_interface_id)
                    AS half_interface_2_id,
                min(system_id)::VARCHAR AS representative_system_id
            FROM interface_members
            GROUP BY entry_pdb_id, half_interface_1_id, half_interface_2_id;

            COPY (
                SELECT * FROM half_representatives
                ORDER BY half_interface_id
            ) TO '{temporary_paths["half_interfaces"].as_posix()}' (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000
            );

            COPY (
                SELECT
                    representative_system_id,
                    entry_pdb_id,
                    half_interface_1_id,
                    half_interface_2_id
                FROM representative_groups
                ORDER BY representative_system_id
            ) TO '{temporary_paths["interfaces"].as_posix()}' (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000
            );

            COPY (
                SELECT
                    members.system_id,
                    groups.representative_system_id,
                    members.side_1_half_interface_id,
                    members.side_2_half_interface_id
                FROM interface_members AS members
                INNER JOIN representative_groups AS groups
                  ON members.entry_pdb_id = groups.entry_pdb_id
                 AND least(
                        members.side_1_half_interface_id,
                        members.side_2_half_interface_id
                     ) = groups.half_interface_1_id
                 AND greatest(
                        members.side_1_half_interface_id,
                        members.side_2_half_interface_id
                     ) = groups.half_interface_2_id
                ORDER BY members.system_id
            ) TO '{temporary_paths["membership"].as_posix()}' (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000
            );
            """
        )
    )
    connection.close()
    LOG.info(
        "make_interface_representatives: DuckDB normalization complete "
        "elapsed_seconds=%.1f",
        time.monotonic() - started,
    )

    expected = {
        "half_interfaces": (
            INTERFACE_HALF_REPRESENTATIVES_RELATIVE,
            schemas.INTERFACE_HALF_REPRESENTATIVE_SCHEMA,
        ),
        "interfaces": (
            INTERFACE_REPRESENTATIVES_RELATIVE,
            schemas.INTERFACE_REPRESENTATIVE_SCHEMA,
        ),
        "membership": (
            INTERFACE_MEMBERSHIP_RELATIVE,
            schemas.INTERFACE_MEMBERSHIP_SCHEMA,
        ),
    }
    for key, (_, schema) in expected.items():
        observed = pq.read_schema(temporary_paths[key])
        if not observed.equals(schema):
            rmtree(working_root)
            raise ValueError(
                f"interface representative {key} has unexpected schema: {observed}"
            )
    membership_rows = pq.ParquetFile(temporary_paths["membership"]).metadata.num_rows
    source_rows = pq.ParquetFile(source).metadata.num_rows
    if membership_rows != source_rows:
        rmtree(working_root)
        raise ValueError(
            "interface membership is incomplete: "
            f"rows={membership_rows}/{source_rows}"
        )
    if _interface_representative_source_signature(data_dir) != source_signature:
        rmtree(working_root)
        raise RuntimeError("interface annotation changed while making representatives")

    outputs: dict[str, dict[str, int | str]] = {}
    for key, (relative, _) in expected.items():
        target = data_dir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        install = target.with_suffix(target.suffix + ".tmp")
        copyfile(temporary_paths[key], install)
        install.replace(target)
        outputs[key] = _interface_representative_output_signature(target)
    rmtree(working_root)
    payload = {
        "status": "complete",
        "source": source_signature,
        "outputs": outputs,
        "interface_count": source_rows,
        "representative_interface_count": pq.ParquetFile(
            data_dir / INTERFACE_REPRESENTATIVES_RELATIVE
        ).metadata.num_rows,
        "representative_half_interface_count": pq.ParquetFile(
            data_dir / INTERFACE_HALF_REPRESENTATIVES_RELATIVE
        ).metadata.num_rows,
    }
    _write_json_atomic(data_dir / INTERFACE_REPRESENTATIVES_MANIFEST_RELATIVE, payload)
    LOG.info(
        "make_interface_representatives: complete interfaces=%d "
        "representatives=%d half_representatives=%d reduction=%.2fx "
        "elapsed_seconds=%.1f",
        payload["interface_count"],
        payload["representative_interface_count"],
        payload["representative_half_interface_count"],
        (
            payload["interface_count"] / payload["representative_interface_count"]
            if payload["representative_interface_count"]
            else 1.0
        ),
        time.monotonic() - started,
    )
    return payload


def _ligand_pocket_representative_source_signatures(
    data_dir: Path,
) -> dict[str, dict[str, int | str]]:
    signatures: dict[str, dict[str, int | str]] = {}
    for name in ["annotation_table.parquet", "entry_chains.parquet"]:
        path = data_dir / "index" / name
        stat = path.stat()
        signatures[name] = {
            "name": name,
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    return signatures


def _completed_ligand_pocket_representatives(
    data_dir: Path,
) -> dict[str, Any] | None:
    manifest = _read_json(data_dir / LIGAND_POCKET_REPRESENTATIVES_MANIFEST_RELATIVE)
    expected = {
        "representatives": (
            LIGAND_POCKET_REPRESENTATIVES_RELATIVE,
            schemas.LIGAND_POCKET_REPRESENTATIVE_SCHEMA,
        ),
        "membership": (
            LIGAND_POCKET_MEMBERSHIP_RELATIVE,
            schemas.LIGAND_POCKET_MEMBERSHIP_SCHEMA,
        ),
    }
    try:
        if manifest is None or manifest.get(
            "sources"
        ) != _ligand_pocket_representative_source_signatures(data_dir):
            return None
        for key, (relative, schema) in expected.items():
            path = data_dir / relative
            if not pq.read_schema(path).equals(schema) or manifest.get(
                "outputs", {}
            ).get(key) != _interface_representative_output_signature(path):
                return None
    except (OSError, TypeError, ValueError):
        return None
    return manifest


def make_ligand_pocket_representatives(
    *,
    data_dir: Path,
    scratch_dir: Path | None,
    threads: int,
    force_update: bool = False,
) -> dict[str, Any]:
    """Normalize exact protein-pocket copies while retaining ligand membership."""
    current = _completed_ligand_pocket_representatives(data_dir)
    if not force_update and current is not None:
        LOG.info(
            "make_ligand_pocket_representatives: reusing %d representatives "
            "for %d scoreable ligands",
            current["representative_ligand_count"],
            current["ligand_count"],
        )
        return current

    import duckdb

    sources = _ligand_pocket_representative_source_signatures(data_dir)
    annotation = data_dir / "index" / "annotation_table.parquet"
    chains = data_dir / "index" / "entry_chains.parquet"
    working_root = (scratch_dir or data_dir / "scratch") / "ligand-representatives"
    if working_root.exists():
        rmtree(working_root)
    working_root.mkdir(parents=True)
    temporary_paths = {
        "representatives": working_root / LIGAND_POCKET_REPRESENTATIVES_RELATIVE.name,
        "membership": working_root / LIGAND_POCKET_MEMBERSHIP_RELATIVE.name,
    }
    started = time.monotonic()
    LOG.info("make_ligand_pocket_representatives: normalizing ligand pockets")
    connection = duckdb.connect()
    connection.sql(f"SET threads={max(1, threads)}")
    connection.sql("SET memory_limit='64GB'")
    connection.sql(f"SET temp_directory='{working_root.as_posix()}'")
    connection.sql("SET preserve_insertion_order=false")
    connection.sql(
        dedent(
            f"""
            CREATE TEMP TABLE protein_receptors AS
            SELECT entry_pdb_id, chain_asym_id
            FROM read_parquet('{chains.as_posix()}')
            WHERE chain_receptor_type = 'protein';

            CREATE TEMP TABLE eligible_ligands AS
            SELECT
                entry_pdb_id,
                system_id,
                ligand_id,
                ligand_asym_id,
                coalesce(ligand_is_3d_score_able, false)
                    AS ligand_is_3d_score_able,
                ligand_protein_chains_asym_id,
                ligand_neighboring_residues,
                ligand_interactions
            FROM read_parquet('{annotation.as_posix()}')
            WHERE system_type = 'holo' AND ligand_is_proper;

            CREATE TEMP TABLE receptor_sets AS
            SELECT
                ligands.entry_pdb_id,
                ligands.ligand_id,
                list_sort(list_distinct(list(receptors.chain_asym_id)))
                    AS receptor_chain_asym_ids
            FROM eligible_ligands AS ligands,
                 unnest(ligands.ligand_protein_chains_asym_id)
                    AS instances(instance_chain)
            INNER JOIN protein_receptors AS receptors
              ON ligands.entry_pdb_id = receptors.entry_pdb_id
             AND split_part(instance_chain, '.', 2) = receptors.chain_asym_id
            GROUP BY ligands.entry_pdb_id, ligands.ligand_id;

            CREATE TEMP TABLE canonical_receptor_sets AS
            SELECT
                *,
                min(ligand_id) OVER (
                    PARTITION BY entry_pdb_id, receptor_chain_asym_ids
                )::VARCHAR AS receptor_set_id
            FROM receptor_sets;

            CREATE TEMP TABLE normalized_pockets AS
            SELECT
                ligands.ligand_id,
                list_sort(list_distinct(list(
                    split_part(neighbor, '.', 2)
                ))) AS pocket_residues
            FROM eligible_ligands AS ligands
            INNER JOIN canonical_receptor_sets AS receptor_sets USING (ligand_id),
                 unnest(ligands.ligand_neighboring_residues)
                    AS residues(neighbor)
            WHERE list_contains(
                receptor_sets.receptor_chain_asym_ids,
                split_part(split_part(neighbor, '.', 2), '_', 1)
            )
            GROUP BY ligands.ligand_id;

            CREATE TEMP TABLE normalized_interactions AS
            SELECT
                ligands.ligand_id,
                list_sort(list(split_part(interaction, '.', 2))) AS interactions
            FROM eligible_ligands AS ligands
            INNER JOIN canonical_receptor_sets AS receptor_sets USING (ligand_id),
                 unnest(ligands.ligand_interactions)
                    AS interaction_rows(interaction)
            WHERE list_contains(
                receptor_sets.receptor_chain_asym_ids,
                split_part(split_part(interaction, '.', 2), '_', 1)
            )
            GROUP BY ligands.ligand_id;

            CREATE TEMP TABLE normalized_ligands AS
            SELECT
                ligands.entry_pdb_id,
                ligands.system_id,
                ligands.ligand_id,
                ligands.ligand_asym_id,
                ligands.ligand_is_3d_score_able,
                receptors.receptor_set_id,
                receptors.receptor_chain_asym_ids,
                coalesce(pockets.pocket_residues, []::VARCHAR[])
                    AS pocket_residues,
                coalesce(interactions.interactions, []::VARCHAR[])
                    AS interactions
            FROM eligible_ligands AS ligands
            INNER JOIN canonical_receptor_sets AS receptors USING (ligand_id)
            LEFT JOIN normalized_pockets AS pockets USING (ligand_id)
            LEFT JOIN normalized_interactions AS interactions USING (ligand_id);

            CREATE TEMP TABLE representative_groups AS
            SELECT
                min(ligand_id)::VARCHAR AS representative_ligand_id,
                arg_min(system_id, ligand_id)::VARCHAR AS representative_system_id,
                entry_pdb_id,
                ligand_asym_id,
                ligand_is_3d_score_able,
                receptor_set_id,
                receptor_chain_asym_ids,
                pocket_residues,
                interactions
            FROM normalized_ligands
            GROUP BY
                entry_pdb_id,
                ligand_asym_id,
                ligand_is_3d_score_able,
                receptor_set_id,
                receptor_chain_asym_ids,
                pocket_residues,
                interactions;

            COPY (
                SELECT * FROM representative_groups
                ORDER BY representative_ligand_id
            ) TO '{temporary_paths["representatives"].as_posix()}' (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000
            );

            COPY (
                SELECT
                    ligands.system_id,
                    ligands.ligand_id,
                    representatives.representative_system_id,
                    representatives.representative_ligand_id
                FROM normalized_ligands AS ligands
                INNER JOIN representative_groups AS representatives
                  ON ligands.entry_pdb_id = representatives.entry_pdb_id
                 AND ligands.ligand_asym_id = representatives.ligand_asym_id
                 AND ligands.ligand_is_3d_score_able
                        = representatives.ligand_is_3d_score_able
                 AND ligands.receptor_set_id = representatives.receptor_set_id
                 AND ligands.receptor_chain_asym_ids
                        = representatives.receptor_chain_asym_ids
                 AND ligands.pocket_residues = representatives.pocket_residues
                 AND ligands.interactions = representatives.interactions
                ORDER BY ligands.ligand_id
            ) TO '{temporary_paths["membership"].as_posix()}' (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000
            );
            """
        )
    )
    connection.close()
    LOG.info(
        "make_ligand_pocket_representatives: DuckDB normalization complete "
        "elapsed_seconds=%.1f",
        time.monotonic() - started,
    )

    expected = {
        "representatives": (
            LIGAND_POCKET_REPRESENTATIVES_RELATIVE,
            schemas.LIGAND_POCKET_REPRESENTATIVE_SCHEMA,
        ),
        "membership": (
            LIGAND_POCKET_MEMBERSHIP_RELATIVE,
            schemas.LIGAND_POCKET_MEMBERSHIP_SCHEMA,
        ),
    }
    for key, (_, schema) in expected.items():
        observed = pq.read_schema(temporary_paths[key])
        if not observed.equals(schema):
            rmtree(working_root)
            raise ValueError(
                f"ligand pocket representative {key} has unexpected schema: "
                f"{observed}"
            )
    ligand_count = pq.ParquetFile(temporary_paths["membership"]).metadata.num_rows
    representative_count = pq.ParquetFile(
        temporary_paths["representatives"]
    ).metadata.num_rows
    if _ligand_pocket_representative_source_signatures(data_dir) != sources:
        rmtree(working_root)
        raise RuntimeError("entry indexes changed while making ligand representatives")

    outputs: dict[str, dict[str, int | str]] = {}
    for key, (relative, _) in expected.items():
        target = data_dir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        install = target.with_suffix(target.suffix + ".tmp")
        copyfile(temporary_paths[key], install)
        install.replace(target)
        outputs[key] = _interface_representative_output_signature(target)
    rmtree(working_root)
    payload = {
        "status": "complete",
        "sources": sources,
        "outputs": outputs,
        "ligand_count": ligand_count,
        "representative_ligand_count": representative_count,
    }
    _write_json_atomic(
        data_dir / LIGAND_POCKET_REPRESENTATIVES_MANIFEST_RELATIVE,
        payload,
    )
    LOG.info(
        "make_ligand_pocket_representatives: complete ligands=%d "
        "representatives=%d reduction=%.2fx elapsed_seconds=%.1f",
        ligand_count,
        representative_count,
        ligand_count / representative_count if representative_count else 1.0,
        time.monotonic() - started,
    )
    return payload


def _completed_alignment_chain_lookup(
    data_dir: Path,
) -> dict[str, int | str] | None:
    lookup = data_dir / ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    manifest = data_dir / ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE
    try:
        if (
            _completed_ligand_pocket_representatives(data_dir) is None
            or _completed_interface_representatives(data_dir) is None
        ):
            return None
        stat = lookup.stat()
        columns = set(pq.read_schema(lookup).names)
        input_signatures = _alignment_chain_lookup_input_signatures(data_dir)
    except (OSError, TypeError, ValueError):
        return None
    expected_columns = {
        "entry_pdb_id",
        "chain_asym_id",
        "chain_auth_id",
        "selected_residue_numbers",
        "selected_residue_indices",
    }
    output: dict[str, int | str] = {
        "name": lookup.name,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if not expected_columns.issubset(columns):
        return None
    try:
        payload = json.loads(manifest.read_text())
    except (OSError, TypeError, ValueError):
        # Allow an in-progress V3 release to adopt a lookup created before the
        # input manifest existed, but only when both inputs are older. The next
        # make_alignment_chain_lookup() call records exact signatures.
        newest_input_mtime = max(
            int(signature["mtime_ns"]) for signature in input_signatures.values()
        )
        if stat.st_mtime_ns < newest_input_mtime:
            return None
        return output
    if payload.get("inputs") != input_signatures or payload.get("output") != output:
        return None
    return output


def _alignment_chain_lookup_input_signatures(
    data_dir: Path,
) -> dict[str, dict[str, int | str]]:
    signatures: dict[str, dict[str, int | str]] = {}
    for name in [
        LIGAND_POCKET_REPRESENTATIVES_RELATIVE.name,
        INTERFACE_HALF_REPRESENTATIVES_RELATIVE.name,
        "entry_chains.parquet",
    ]:
        path = data_dir / "index" / name
        stat = path.stat()
        signatures[name] = {
            "name": name,
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    return signatures


def _write_alignment_chain_lookup_manifest(data_dir: Path) -> None:
    lookup = data_dir / ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    stat = lookup.stat()
    manifest = data_dir / ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE
    temporary = manifest.with_suffix(manifest.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            {
                "inputs": _alignment_chain_lookup_input_signatures(data_dir),
                "output": {
                    "name": lookup.name,
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    temporary.replace(manifest)


def _refresh_representative_source_manifests(data_dir: Path) -> None:
    """Record index-only enrichments without rebuilding unchanged representatives."""
    ligand_manifest_path = data_dir / LIGAND_POCKET_REPRESENTATIVES_MANIFEST_RELATIVE
    ligand_manifest = _read_json(ligand_manifest_path)
    interface_manifest_path = data_dir / INTERFACE_REPRESENTATIVES_MANIFEST_RELATIVE
    interface_manifest = _read_json(interface_manifest_path)
    if ligand_manifest is None or interface_manifest is None:
        raise RuntimeError(
            "representative source manifests disappeared during finalization"
        )
    ligand_manifest["sources"] = _ligand_pocket_representative_source_signatures(
        data_dir
    )
    interface_manifest["source"] = _interface_representative_source_signature(data_dir)
    _write_json_atomic(ligand_manifest_path, ligand_manifest)
    _write_json_atomic(interface_manifest_path, interface_manifest)


def make_alignment_chain_lookup(
    *,
    data_dir: Path,
    scratch_dir: Path | None,
    threads: int,
    force_update: bool = False,
) -> Path:
    """Build the compact chain and selected-residue map used by every shard."""
    make_ligand_pocket_representatives(
        data_dir=data_dir,
        scratch_dir=scratch_dir,
        threads=threads,
        force_update=force_update,
    )
    make_interface_representatives(
        data_dir=data_dir,
        scratch_dir=scratch_dir,
        threads=threads,
        force_update=force_update,
    )
    lookup = data_dir / ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    if not force_update and _completed_alignment_chain_lookup(data_dir) is not None:
        manifest = data_dir / ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE
        if not manifest.is_file():
            _write_alignment_chain_lookup_manifest(data_dir)
        LOG.info(f"make_alignment_chain_lookup: reusing {lookup}")
        return lookup

    import duckdb

    input_signatures = _alignment_chain_lookup_input_signatures(data_dir)
    working_root = (scratch_dir or data_dir / "scratch") / "alignment-chain-lookup"
    working_root.mkdir(exist_ok=True, parents=True)
    temporary = working_root / lookup.name
    temporary.unlink(missing_ok=True)
    ligand_representatives = (
        data_dir / LIGAND_POCKET_REPRESENTATIVES_RELATIVE
    ).as_posix()
    interfaces = (data_dir / INTERFACE_HALF_REPRESENTATIVES_RELATIVE).as_posix()
    chains = (data_dir / "index" / "entry_chains.parquet").as_posix()
    con = duckdb.connect()
    con.sql(f"set threads={max(1, threads)};")
    con.sql("set memory_limit='64GB';")
    con.sql("set preserve_insertion_order=false;")
    con.sql(f"set temp_directory='{working_root.as_posix()}';")
    con.sql(
        dedent(
            f"""
            COPY (
                WITH protein_chains AS (
                    SELECT entry_pdb_id, chain_asym_id, chain_auth_id
                    FROM read_parquet('{chains}')
                    WHERE chain_receptor_type = 'protein'
                      AND chain_auth_id IS NOT NULL
                ),
                ligand_pocket_residues AS (
                    SELECT
                        representatives.entry_pdb_id,
                        split_part(neighbor, '_', 1) AS chain_asym_id,
                        CAST(split_part(neighbor, '_', 2) AS INTEGER)
                            AS residue_number,
                        CAST(split_part(neighbor, '_', 3) AS INTEGER)
                            AS residue_index
                    FROM read_parquet('{ligand_representatives}')
                        AS representatives,
                    UNNEST(representatives.pocket_residues)
                        AS residues(neighbor)
                ),
                interface_residues AS (
                    SELECT
                        entry_pdb_id,
                        chain_asym_id,
                        unnest(residue_numbers)
                            AS residue_number,
                        unnest(residue_indices)
                            AS residue_index
                    FROM read_parquet('{interfaces}')
                ),
                selected_residues AS (
                    SELECT * FROM ligand_pocket_residues
                    UNION ALL
                    SELECT * FROM interface_residues
                ),
                canonical_selected_residues AS (
                    SELECT
                        entry_pdb_id,
                        chain_asym_id,
                        residue_index,
                        min(residue_number) AS residue_number
                    FROM selected_residues
                    GROUP BY entry_pdb_id, chain_asym_id, residue_index
                ),
                selected_mapping AS (
                    SELECT
                        entry_pdb_id,
                        chain_asym_id,
                        list(residue_number ORDER BY residue_index)
                            AS selected_residue_numbers,
                        list(residue_index ORDER BY residue_index)
                            AS selected_residue_indices
                    FROM canonical_selected_residues
                    GROUP BY entry_pdb_id, chain_asym_id
                )
                SELECT
                    c.entry_pdb_id,
                    c.chain_asym_id,
                    c.chain_auth_id,
                    coalesce(p.selected_residue_numbers, []::INTEGER[])
                        AS selected_residue_numbers,
                    coalesce(p.selected_residue_indices, []::INTEGER[])
                        AS selected_residue_indices
                FROM protein_chains c
                LEFT JOIN selected_mapping p
                    USING (entry_pdb_id, chain_asym_id)
                ORDER BY c.entry_pdb_id, c.chain_asym_id
            ) TO '{temporary.as_posix()}'
            (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100_000);
            """
        )
    )
    lookup_is_unchanged = False
    if lookup.is_file():
        lookup_schema = pq.read_schema(lookup)
        temporary_schema = pq.read_schema(temporary)
        if lookup_schema.equals(temporary_schema):
            difference_count = con.sql(
                dedent(
                    f"""
                    SELECT count(*) FROM (
                        (SELECT * FROM read_parquet('{lookup.as_posix()}')
                         EXCEPT ALL
                         SELECT * FROM read_parquet('{temporary.as_posix()}'))
                        UNION ALL
                        (SELECT * FROM read_parquet('{temporary.as_posix()}')
                         EXCEPT ALL
                         SELECT * FROM read_parquet('{lookup.as_posix()}'))
                    )
                    """
                )
            ).fetchone()
            lookup_is_unchanged = difference_count == (0,)
    con.close()
    if _alignment_chain_lookup_input_signatures(data_dir) != input_signatures:
        temporary.unlink(missing_ok=True)
        raise RuntimeError("entry indexes changed while building chain lookup")
    if lookup_is_unchanged:
        temporary.unlink(missing_ok=True)
        _write_alignment_chain_lookup_manifest(data_dir)
        LOG.info(
            "make_alignment_chain_lookup: refreshed unchanged lookup input signatures"
        )
        return lookup
    lookup.parent.mkdir(exist_ok=True, parents=True)
    install_path = lookup.with_suffix(lookup.suffix + ".tmp")
    copyfile(temporary, install_path)
    install_path.replace(lookup)
    temporary.unlink(missing_ok=True)
    _write_alignment_chain_lookup_manifest(data_dir)
    return lookup


def compute_ligand_fingerprints(
    *,
    data_dir: Path,
    cofactor_similarity_threshold: float = 90.0,
) -> None:
    """Fingerprint unique ligand SMILES and annotate cofactor similarity."""
    LOG.info("compute_ligand_fingerprints: running")
    get_similarity_scores.compute_ligand_fingerprints(
        data_dir=data_dir,
        cofactor_similarity_threshold=cofactor_similarity_threshold,
    )


def scatter_make_ligand_scores(
    *,
    data_dir: Path,
    batch_size: int,
    number_id_col: str = "ligand_smiles_id",
) -> list[list[int]]:
    """Scatter the unique-SMILES fingerprint node IDs."""
    ligands = pd.read_parquet(
        data_dir / "fingerprints" / "ligands_per_smiles.parquet",
        columns=[number_id_col],
    )[number_id_col].to_list()
    LOG.info(f"scatter_make_ligand_scores: found {len(ligands)} ligands")
    chunks = [
        ligands[pos : pos + batch_size] for pos in range(0, len(ligands), batch_size)
    ]
    return chunks or [[]]


def make_ligand_scores(
    *,
    data_dir: Path,
    ligand_ids: list[int],
    minimum_similarity: float = 30.0,
    number_id_col: str = "ligand_smiles_id",
) -> None:
    if not ligand_ids:
        LOG.info("make_ligand_scores: no ligand nodes to score")
        return
    hashid = utils.hash_contents([str(i) for i in ligand_ids])
    output_path = data_dir / "ligand_scores" / f"{hashid}.parquet"
    output_path.parent.mkdir(exist_ok=True, parents=True)
    get_similarity_scores.ligand_scores(
        ligand_ids=ligand_ids,
        data_dir=data_dir,
        output_path=output_path,
        number_id_col=number_id_col,
        minimum_similarity=minimum_similarity,
    )


def annotate_ligand_similarity(*, data_dir: Path) -> None:
    """Write ligand identifiers and cofactor annotations."""
    get_similarity_scores.annotate_ligand_similarity(data_dir=data_dir)


def make_ligand_mmp_pairs(
    *,
    data_dir: Path,
    scratch_dir: Path,
    threads: int,
    force_update: bool = False,
) -> Path:
    """Write the unique-SMILES matched-molecular-pair release table."""
    return mmpdb_utils.make_ligand_mmp_pairs(
        data_dir=data_dir,
        scratch_dir=scratch_dir,
        threads=threads,
        force_update=force_update,
    )


def _interface_scoring_chain_keys(data_dir: Path) -> pd.DataFrame:
    """Return entry/asym keys used by a published protein interface."""
    interface_path = data_dir / "index" / "interface_annotation_table.parquet"
    if not interface_path.is_file():
        return pd.DataFrame(
            columns=["entry_pdb_id", "chain_asym_id", "chain_is_interface"]
        )
    interfaces = pd.read_parquet(
        interface_path,
        columns=["entry_pdb_id", "interface_chain_1", "interface_chain_2"],
    )
    interface_keys = pd.concat(
        [
            interfaces[["entry_pdb_id", column]].rename(
                columns={column: "instance_chain"}
            )
            for column in ["interface_chain_1", "interface_chain_2"]
        ],
        ignore_index=True,
    )
    interface_keys["chain_asym_id"] = (
        interface_keys.pop("instance_chain").astype(str).str.split(".", n=1).str[-1]
    )
    interface_keys = interface_keys.drop_duplicates()
    interface_keys["chain_is_interface"] = True
    return interface_keys


def _protein_scoring_chains(data_dir: Path) -> pd.DataFrame:
    """Return protein chains used by a ligand receptor or protein interface."""
    chain_path = data_dir / "index" / "entry_chains.parquet"
    chains = pd.read_parquet(
        chain_path,
        columns=[
            "entry_pdb_id",
            "chain_asym_id",
            "chain_receptor_type",
            "chain_is_holo",
        ],
    )
    chains = chains.merge(
        _interface_scoring_chain_keys(data_dir),
        on=["entry_pdb_id", "chain_asym_id"],
        how="left",
    )
    return chains.loc[
        chains["chain_receptor_type"].fillna("").astype(str).eq("protein")
        & (
            chains["chain_is_holo"].fillna(False).astype(bool)
            | chains["chain_is_interface"].eq(True)
        )
    ].copy()


def _apo_scoring_chains(data_dir: Path) -> pd.DataFrame:
    """Return reconstructable chains without a proper ligand receptor."""
    from plinder.data.linked_apo import ligand_holo_chain_keys

    chains = pd.read_parquet(
        data_dir / "index" / "entry_chains.parquet",
        columns=[
            "entry_pdb_id",
            "chain_asym_id",
            "chain_auth_id",
            "chain_entity_id",
            "chain_receptor_type",
            "chain_is_ligand_like",
        ],
    )
    holo_keys = ligand_holo_chain_keys(data_dir / "index" / "annotation_table.parquet")
    holo_keys["chain_is_holo"] = True
    chains = chains.merge(
        holo_keys,
        on=["entry_pdb_id", "chain_asym_id"],
        how="left",
        validate="one_to_one",
    )
    holo = chains["chain_is_holo"].eq(True)
    entity_ids = chains["chain_entity_id"].astype("string")
    usable_entity = entity_ids.notna() & entity_ids.str.strip().ne("")
    holo_entities = chains.loc[
        holo & usable_entity, ["entry_pdb_id", "chain_entity_id"]
    ].drop_duplicates()
    holo_entities["entity_is_holo"] = True
    chains = chains.merge(
        holo_entities,
        on=["entry_pdb_id", "chain_entity_id"],
        how="left",
        validate="many_to_one",
    )
    chains = chains.loc[
        chains["chain_receptor_type"].fillna("").astype(str).eq("protein")
        & ~holo
        & ~chains["chain_is_ligand_like"].fillna(False).astype(bool)
        & chains["entity_is_holo"].ne(True)
        & chains["chain_auth_id"].notna()
    ].copy()
    receptor_chains = pd.read_parquet(
        data_dir / "index" / "entry_biounit_chains.parquet",
        columns=["entry_pdb_id", "chain_asym_id", "chain_role"],
    )
    receptor_chains = receptor_chains.loc[
        receptor_chains["chain_role"].fillna("").astype(str).str.lower().eq("receptor"),
        ["entry_pdb_id", "chain_asym_id"],
    ].drop_duplicates()
    return chains.merge(
        receptor_chains,
        on=["entry_pdb_id", "chain_asym_id"],
        how="inner",
        validate="one_to_one",
    )


def scatter_protein_scoring(
    *,
    data_dir: Path,
    batch_size: int,
    two_char_codes: list[str],
    pdb_ids: list[str],
) -> list[list[str]]:
    """Split protein-containing V3 entries into score-generation batches.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    batch_size : int
        how many codes to put in a chunk
    two_char_codes : list[str], default=[]
        only consider particular codes
    pdb_ids : list[str], default=[]
        only consider particular pdb IDs

    Returns
    -------
    codes : list[list[str]]
        list of lists of chunks of two character codes
    """
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    chains = _protein_scoring_chains(data_dir)
    protein_entries = set(chains["entry_pdb_id"].astype(str))
    selected_pdb_ids = [normalize_pdb_id(pdb_id) for pdb_id in pdb_ids]
    if selected_pdb_ids:
        selected = sorted(protein_entries.intersection(selected_pdb_ids))
    else:
        selected_codes = {str(code).lower() for code in two_char_codes}
        selected = sorted(
            pdb_id
            for pdb_id in protein_entries
            if not selected_codes or pdb_id[1:3] in selected_codes
        )
    LOG.info(f"scatter_protein_scoring: found {len(selected)} protein PDB IDs")
    chunks = [
        selected[pos : pos + batch_size] for pos in range(0, len(selected), batch_size)
    ]
    return chunks or [[]]


def run_batch_searches(
    *,
    data_dir: Path,
    pdb_ids: list[str],
    scorer_cfg: DictConfig,
    foldseek_cfg: DictConfig,
    mmseqs_cfg: DictConfig,
    cpu: int,
    scratch_dir: Path | None = None,
    alignment_types: Sequence[str] | None = None,
    force_update: bool = False,
) -> None:
    selected_alignment_types = list(alignment_types or ["foldseek", "mmseqs"])
    for search_db in scorer_cfg.sub_databases:
        for alignment_type in selected_alignment_types:
            output_dir = (
                data_dir / "dbs" / "subdbs" / f"{search_db}_{alignment_type}" / "aln"
            )
            pending = [
                pdb_id
                for pdb_id in pdb_ids
                if force_update or not (output_dir / f"{pdb_id}.parquet").is_file()
            ]
            if not pending:
                LOG.info(
                    f"run_batch_searches: all {len(pdb_ids)} {search_db} "
                    f"{alignment_type} queries are complete"
                )
                continue
            LOG.info(
                f"run_batch_searches: searching {len(pending)} of {len(pdb_ids)} "
                f"{search_db} queries with {alignment_type}"
            )
            scorer, entry_ids, batch_db_dir = utils.get_scorer(
                data_dir=data_dir,
                pdb_ids=pending,
                scorer_cfg=scorer_cfg,
                load_entries=True,
                foldseek_cfg=foldseek_cfg,
                mmseqs_cfg=mmseqs_cfg,
                scratch_dir=scratch_dir,
            )
            try:
                scorer.run_alignments(
                    entry_ids=entry_ids,
                    output_folder=batch_db_dir,
                    overwrite=True,
                    search_db=search_db,
                    threads=cpu,
                    alignment_types=[alignment_type],
                )
            finally:
                rmtree(batch_db_dir)
            missing_outputs = [
                pdb_id
                for pdb_id in pending
                if not (output_dir / f"{pdb_id}.parquet").is_file()
            ]
            if missing_outputs:
                raise RuntimeError(
                    f"{search_db} {alignment_type} search produced no output for "
                    f"{len(missing_outputs)} query entries: "
                    f"{missing_outputs[:10]}"
                )
        if search_db == "holo" and alignment_types is None:
            missing = [
                pdb_id
                for pdb_id in pdb_ids
                if not any(
                    (
                        data_dir
                        / "dbs"
                        / "subdbs"
                        / f"{search_db}_{alignment_type}"
                        / "aln"
                        / f"{pdb_id}.parquet"
                    ).is_file()
                    for alignment_type in selected_alignment_types
                )
            ]
            if missing:
                raise RuntimeError(
                    "holo search did not produce results from either backend for "
                    f"{len(missing)} query entries: {missing[:10]}"
                )


def scatter_missing_alignment_mappings(
    *, data_dir: Path, batch_size: int, search_db: str = "holo"
) -> list[list[str]]:
    """Scatter query shards whose raw alignments are not mapped and published."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    raw_root = data_dir / "dbs" / "subdbs"
    shards = sorted(
        {
            path.stem[1:3]
            for alignment_type in ["foldseek", "mmseqs"]
            for path in (raw_root / f"{search_db}_{alignment_type}" / "aln").glob(
                "*.parquet"
            )
            if not path.name.endswith(".tmp.parquet")
        }
    )
    missing = [
        shard
        for shard in shards
        if not alignment_mapping_shard_is_current(
            data_dir=data_dir,
            search_db=search_db,
            shard=shard,
        )
    ]
    chunks = [
        missing[pos : pos + batch_size] for pos in range(0, len(missing), batch_size)
    ]
    return chunks or [[]]


def _alignment_input_signatures(
    *, data_dir: Path, shard: str, search_db: str = "holo"
) -> dict[str, list[dict[str, int | str]]]:
    signatures: dict[str, list[dict[str, int | str]]] = {}
    for alignment_type in ["foldseek", "mmseqs"]:
        source_dir = (
            data_dir / "dbs" / "subdbs" / f"{search_db}_{alignment_type}" / "aln"
        )
        sources = sorted(
            path
            for path in source_dir.glob("*.parquet")
            if path.stem[1:3] == shard and not path.name.endswith(".tmp.parquet")
        )
        signatures[alignment_type] = [
            {
                "name": path.name,
                "size": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
            }
            for path in sources
        ]
    return signatures


def _alignment_release_path(
    *, data_dir: Path, search_db: str, alignment_type: str, shard: str
) -> Path:
    return (
        data_dir
        / "alignments"
        / f"search_db={search_db}"
        / f"alignment_type={alignment_type}"
        / f"shard={shard}.parquet"
    )


def _alignment_mapping_manifest_path(
    *, data_dir: Path, shard: str, search_db: str = "holo"
) -> Path:
    root = data_dir / "alignments" / "manifests"
    if search_db != "holo":
        root = root / f"search_db={search_db}"
    return root / f"shard={shard}.json"


def _alignment_target_entry_ids(sources: Sequence[Path]) -> set[str]:
    """Read the distinct target IDs once for a complete query shard."""
    if not sources:
        return set()
    targets: set[str] = set()
    scanner = ds.dataset(sources, format="parquet").scanner(
        columns=["target_pdb_id"],
        batch_size=262_144,
        use_threads=False,
    )
    for batch in scanner.to_batches():
        unique = pc.unique(batch.column(0))
        targets.update(str(value) for value in unique.to_pylist() if value is not None)
    return targets


def alignment_mapping_shard_is_current(
    *, data_dir: Path, shard: str, search_db: str = "holo"
) -> bool:
    """Validate a shard manifest against raw inputs and published outputs."""
    manifest_path = _alignment_mapping_manifest_path(
        data_dir=data_dir,
        search_db=search_db,
        shard=shard,
    )
    try:
        payload = json.loads(manifest_path.read_text())
    except (OSError, TypeError, ValueError):
        return False
    inputs = _alignment_input_signatures(
        data_dir=data_dir,
        search_db=search_db,
        shard=shard,
    )
    lookup_signature = _completed_alignment_chain_lookup(data_dir)
    if lookup_signature is None:
        return False
    if (
        payload.get("shard") != shard
        or (search_db != "holo" and payload.get("search_db") != search_db)
        or payload.get("inputs") != inputs
        or payload.get("alignment_chain_lookup") != lookup_signature
    ):
        return False
    outputs = payload.get("outputs")
    if not isinstance(outputs, dict):
        return False
    for alignment_type, source_signatures in inputs.items():
        output = _alignment_release_path(
            data_dir=data_dir,
            search_db=search_db,
            alignment_type=alignment_type,
            shard=shard,
        )
        if not source_signatures:
            if outputs.get(alignment_type) is not None:
                return False
            continue
        if not output.is_file() or not utils._mapped_alignment_file_is_current(
            output, alignment_type=alignment_type
        ):
            return False
        stat = output.stat()
        if outputs.get(alignment_type) != {
            "name": output.name,
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }:
            return False
    return True


def map_batch_alignments(
    *,
    data_dir: Path,
    shards: list[str],
    scorer_cfg: DictConfig,
    force_update: bool,
    scratch_dir: Path | None = None,
    search_db: str = "holo",
) -> None:
    """Map raw backend hits directly into atomic query-shard release files."""
    if search_db not in scorer_cfg.sub_databases:
        raise ValueError(f"alignment database is not enabled: {search_db}")
    maximum_rows = int(getattr(scorer_cfg, "max_alignment_rows_per_query", 5_000_000))
    for shard in shards:
        if not force_update and alignment_mapping_shard_is_current(
            data_dir=data_dir,
            search_db=search_db,
            shard=shard,
        ):
            LOG.info(
                "map_batch_alignments: %s shard %s is complete",
                search_db,
                shard,
            )
            continue
        inputs = _alignment_input_signatures(
            data_dir=data_dir,
            search_db=search_db,
            shard=shard,
        )
        pdb_ids = sorted(
            {
                Path(str(signature["name"])).stem
                for signatures in inputs.values()
                for signature in signatures
            }
        )
        if not pdb_ids:
            continue
        rows_by_query: dict[str, dict[str, int]] = {pdb_id: {} for pdb_id in pdb_ids}
        for alignment_type, signatures in inputs.items():
            source_root = (
                data_dir / "dbs" / "subdbs" / f"{search_db}_{alignment_type}" / "aln"
            )
            for signature in signatures:
                source = source_root / str(signature["name"])
                pdb_id = source.stem
                rows_by_query[pdb_id][alignment_type] = int(
                    pq.ParquetFile(source).metadata.num_rows
                )
        skipped_queries = {
            pdb_id: {
                "reason": "raw_alignment_row_budget_exceeded",
                "maximum_rows": maximum_rows,
                "total_rows": sum(rows.values()),
                "rows_by_alignment_type": rows,
            }
            for pdb_id, rows in rows_by_query.items()
            if sum(rows.values()) > maximum_rows
        }
        if skipped_queries:
            LOG.warning(
                "map_batch_alignments: shard %s skips %s oversized queries: %s",
                shard,
                len(skipped_queries),
                sorted(skipped_queries),
            )
        mapped_pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in skipped_queries]
        lookup_signature = _completed_alignment_chain_lookup(data_dir)
        if lookup_signature is None:
            raise FileNotFoundError(
                "missing or stale alignment chain lookup; run make-sub-dbs first"
            )
        raw_root = data_dir / "dbs" / "subdbs"
        raw_sources = [
            raw_root / f"{search_db}_{alignment_type}" / "aln" / str(signature["name"])
            for alignment_type, signatures in inputs.items()
            for signature in signatures
            if Path(str(signature["name"])).stem not in skipped_queries
        ]
        mapping_entry_ids = set(mapped_pdb_ids)
        if search_db != "pred":
            mapping_entry_ids.update(_alignment_target_entry_ids(raw_sources))
        working_root = (scratch_dir or data_dir / "scratch" / "mapping") / (
            f"{search_db}-{shard}"
        )
        if working_root.exists():
            rmtree(working_root)
        mapped_db_dir = working_root / "mapped"
        scorer, entry_ids, _ = utils.get_scorer(
            data_dir=data_dir,
            pdb_ids=mapped_pdb_ids,
            scorer_cfg=scorer_cfg,
            load_entries=False,
            scratch_dir=working_root,
        )
        from plinder.core.scores.entries import load_alignment_entry_views

        scorer.entries.update(
            load_alignment_entry_views(
                lookup_path=data_dir / ALIGNMENT_CHAIN_LOOKUP_RELATIVE,
                pdb_ids=mapping_entry_ids,
            )
        )
        missing_mapping_entries = mapping_entry_ids.difference(scorer.entries)
        if missing_mapping_entries:
            raise ValueError(
                "alignment chain lookup is missing protein entries: "
                f"{sorted(missing_mapping_entries)[:10]}"
            )
        try:
            for pdb_id in tqdm(entry_ids):
                expected = sum(
                    any(
                        str(signature["name"]) == f"{pdb_id}.parquet"
                        for signature in inputs[alignment_type]
                    )
                    for alignment_type in ["foldseek", "mmseqs"]
                )
                mapped = scorer.map_alignment_files(
                    data_dir,
                    pdb_id,
                    search_db,
                    overwrite=True,
                    scratch_dir=working_root / "temporary",
                    mapped_db_dir=mapped_db_dir,
                )
                if len(mapped) != expected:
                    raise RuntimeError(
                        f"{search_db} alignment mapping for {pdb_id} produced "
                        f"{len(mapped)} of {expected} available backends"
                    )
            outputs: dict[str, dict[str, int | str] | None] = {}
            for alignment_type, source_signatures in inputs.items():
                target = _alignment_release_path(
                    data_dir=data_dir,
                    search_db=search_db,
                    alignment_type=alignment_type,
                    shard=shard,
                )
                if not source_signatures:
                    target.unlink(missing_ok=True)
                    outputs[alignment_type] = None
                    continue
                local_sources = sorted(
                    (
                        mapped_db_dir / f"{search_db}_{alignment_type}" / "mapped_aln"
                    ).glob("*.parquet")
                )
                _write_alignment_release_shard(
                    sources=local_sources,
                    target=target,
                    alignment_type=alignment_type,
                    temp_dir=working_root / "collate" / alignment_type,
                    threads=1,
                    memory_limit="7GB",
                )
                stat = target.stat()
                outputs[alignment_type] = {
                    "name": target.name,
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            manifest = _alignment_mapping_manifest_path(
                data_dir=data_dir,
                search_db=search_db,
                shard=shard,
            )
            manifest.parent.mkdir(exist_ok=True, parents=True)
            temporary_manifest = manifest.with_suffix(".tmp.json")
            temporary_manifest.write_text(
                json.dumps(
                    {
                        "search_db": search_db,
                        "shard": shard,
                        "alignment_chain_lookup": lookup_signature,
                        "inputs": inputs,
                        "outputs": outputs,
                        "skipped_queries": skipped_queries,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
            temporary_manifest.replace(manifest)
        finally:
            if working_root.exists():
                rmtree(working_root)


def scatter_missing_scores(
    *,
    data_dir: Path,
    batch_size: int,
    scorer_cfg: DictConfig,
    search_dbs: Sequence[str] = ("holo",),
) -> list[list[str]]:
    from plinder.data.pipeline.score import dropped_query_ids

    present = utils.get_pdb_ids_in_scoring_dataset(data_dir=data_dir)
    score_work = data_dir / "manifests" / "protein_scoring_work.parquet"
    if score_work.is_file():
        complete_holo: list[str] = []
        for pdb_id in present["holo"]:
            score_path = (
                data_dir / "dbs" / "subdbs" / "search_db=holo" / f"{pdb_id}.parquet"
            )
            candidate_path = (
                data_dir
                / "scores"
                / "ligand_3d_candidates"
                / "search_db=holo"
                / f"shard={pdb_id[1:3]}"
                / f"{pdb_id}.parquet"
            )
            try:
                metadata = pq.read_schema(score_path).metadata or {}
            except (OSError, ValueError):
                continue
            if (
                metadata.get(b"plinder.ligand_3d") in {b"deferred", b"complete"}
                and candidate_path.is_file()
            ):
                complete_holo.append(pdb_id)
        present["holo"] = complete_holo
    rerun: set[str] = set()
    dropped = dropped_query_ids(data_dir)
    for search_db in search_dbs:
        score_mode = b"deferred" if search_db == "holo" else b"complete"
        present[search_db] = [
            pdb_id
            for pdb_id in present.get(search_db, [])
            if get_similarity_scores.score_cache_is_current(
                data_dir
                / "dbs"
                / "subdbs"
                / f"search_db={search_db}"
                / f"{pdb_id}.parquet",
                ligand_3d_mode=score_mode,
                minimum_threshold=float(scorer_cfg.minimum_threshold),
                minimum_thresholds=dict(scorer_cfg.minimum_thresholds),
                holo_protein_scores_mode=(b"excluded" if search_db == "holo" else None),
            )
        ]
        manifest_root = data_dir / "alignments" / "manifests"
        if search_db != "holo":
            manifest_root = manifest_root / f"search_db={search_db}"
        mapped_queries: set[str] = set()
        for manifest_path in sorted(manifest_root.glob("shard=*.json")):
            shard = manifest_path.stem.removeprefix("shard=")
            if not alignment_mapping_shard_is_current(
                data_dir=data_dir,
                search_db=search_db,
                shard=shard,
            ):
                continue
            payload = json.loads(manifest_path.read_text())
            mapped_queries.update(
                Path(str(signature["name"])).stem
                for alignment_type in ["foldseek", "mmseqs"]
                for signature in payload["inputs"][alignment_type]
            )
        rerun.update(
            mapped_queries.difference(present.get(search_db, [])).difference(dropped)
        )
    run = sorted(rerun)
    if score_work.is_file():
        planned = pd.read_parquet(
            score_work, columns=["pdb_id", "score_batch_index", "estimated_work"]
        )
        planned = planned[planned["pdb_id"].isin(run)].sort_values(
            ["score_batch_index", "estimated_work"], ascending=[True, False]
        )
        chunks = [
            group["pdb_id"].astype(str).tolist()
            for _, group in planned.groupby("score_batch_index", sort=True)
        ]
        return chunks or [[]]
    chunks = [run[pos : pos + batch_size] for pos in range(0, len(run), batch_size)]
    return chunks or [[]]


def make_batch_scores(
    *,
    data_dir: Path,
    pdb_ids: list[str],
    scorer_cfg: DictConfig,
    force_update: bool,
    scratch_dir: Path | None = None,
    threads: int = 1,
    defer_ligand_3d: bool = True,
) -> None:
    scorer, entry_ids, _ = utils.get_scorer(
        data_dir=data_dir,
        pdb_ids=pdb_ids,
        scorer_cfg=scorer_cfg,
        load_entries=False,
    )
    if threads < 1:
        raise ValueError("threads must be positive")
    scorer.shape_score_threads = threads
    for search_db in scorer_cfg.sub_databases:
        if search_db != "holo":
            from plinder.core.scores.entries import load_entry_views

            missing_entries = set(entry_ids).difference(scorer.entries)
            if missing_entries:
                scorer.entries.update(
                    load_entry_views(pdb_ids=missing_entries, data_dir=data_dir)
                )
            entries_by_shard: dict[str, list[str]] = {}
            for pdb_id in entry_ids:
                entries_by_shard.setdefault(pdb_id[1:3], []).append(pdb_id)
            for shard, shard_entry_ids in entries_by_shard.items():
                source_to_aln_file = {
                    f"{search_db}_{alignment_type}": _alignment_release_path(
                        data_dir=data_dir,
                        search_db=search_db,
                        alignment_type=alignment_type,
                        shard=shard,
                    )
                    for alignment_type in ["foldseek", "mmseqs"]
                }
                alignments = scorer.load_alignments(
                    source_to_aln_file=source_to_aln_file,
                    search_db=search_db,
                    query_entry_ids=set(shard_entry_ids),
                )
                for pdb_id in tqdm(shard_entry_ids):
                    if alignments.empty:
                        query_alignments = alignments
                    else:
                        try:
                            query_alignments = alignments.loc[pdb_id]
                        except KeyError:
                            query_alignments = pd.DataFrame()
                    scorer.get_score_df(
                        data_dir,
                        pdb_id,
                        search_db=search_db,
                        overwrite=force_update,
                        map_alignments=False,
                        scratch_dir=scratch_dir,
                        source_to_aln_file=source_to_aln_file,
                        defer_ligand_3d=False,
                        query_entry_alignments=query_alignments,
                    )
            continue
        for pdb_id in tqdm(entry_ids):
            source_to_aln_file = {
                f"{search_db}_{alignment_type}": _alignment_release_path(
                    data_dir=data_dir,
                    search_db=search_db,
                    alignment_type=alignment_type,
                    shard=pdb_id[1:3],
                )
                for alignment_type in ["foldseek", "mmseqs"]
            }
            scorer.get_score_df(
                data_dir,
                pdb_id,
                search_db=search_db,
                overwrite=force_update,
                map_alignments=False,
                scratch_dir=scratch_dir,
                source_to_aln_file=source_to_aln_file,
                defer_ligand_3d=defer_ligand_3d and search_db == "holo",
            )


def repair_batch_scores(
    *,
    data_dir: Path,
    repairs: list[dict[str, Any]],
    scorer_cfg: DictConfig,
    scratch_dir: Path,
    threads: int = 1,
) -> None:
    """Repair complete affected queries and target-only rows in place."""
    if threads < 1:
        raise ValueError("threads must be positive")
    pdb_ids = [
        str(repair["pdb_id"])
        for repair in repairs
        if str(repair["repair_mode"]) != "drop"
    ]
    scorer = None
    if pdb_ids:
        scorer, _, _ = utils.get_scorer(
            data_dir=data_dir,
            pdb_ids=pdb_ids,
            scorer_cfg=scorer_cfg,
            load_entries=False,
        )
        scorer.shape_score_threads = threads
    started = time.perf_counter()
    for index, repair in enumerate(repairs, start=1):
        pdb_id = str(repair["pdb_id"])
        mode = str(repair["repair_mode"])
        query_started = time.perf_counter()
        LOG.info(
            "score repair progress: query=%d/%d pdb_id=%s mode=%s targets=%d",
            index,
            len(repairs),
            pdb_id,
            mode,
            len(repair.get("target_pdb_ids", [])),
        )
        if mode == "drop":
            (
                data_dir / "dbs" / "subdbs" / "search_db=holo" / f"{pdb_id}.parquet"
            ).unlink(missing_ok=True)
            (
                data_dir
                / "scores"
                / "ligand_3d_candidates"
                / "search_db=holo"
                / f"shard={pdb_id[1:3]}"
                / f"{pdb_id}.parquet"
            ).unlink(missing_ok=True)
        elif mode == "full":
            if scorer is None:
                raise RuntimeError("score repair unexpectedly lacks a scorer")
            scorer.get_score_df(
                data_dir,
                pdb_id,
                search_db="holo",
                overwrite=True,
                map_alignments=False,
                scratch_dir=scratch_dir,
                source_to_aln_file={
                    f"holo_{alignment_type}": _alignment_release_path(
                        data_dir=data_dir,
                        search_db="holo",
                        alignment_type=alignment_type,
                        shard=pdb_id[1:3],
                    )
                    for alignment_type in ["foldseek", "mmseqs"]
                },
                defer_ligand_3d=True,
            )
        elif mode in {"targets", "bounded"}:
            if scorer is None:
                raise RuntimeError("score repair unexpectedly lacks a scorer")
            scorer.repair_score_df_targets(
                data_dir,
                pdb_id,
                affected_target_entries=set(map(str, repair["target_pdb_ids"])),
                scratch_dir=scratch_dir,
                allow_missing=mode == "bounded",
                query_system_ids=(
                    set(map(str, repair["query_system_ids"]))
                    if mode == "bounded"
                    else None
                ),
                query_ligand_ids=(
                    set(map(str, repair["query_ligand_ids"]))
                    if mode == "bounded"
                    else None
                ),
                target_system_ids=(
                    set(map(str, repair["target_system_ids"]))
                    if mode == "bounded"
                    else None
                ),
                target_ligand_ids=(
                    set(map(str, repair["target_ligand_ids"]))
                    if mode == "bounded"
                    else None
                ),
            )
        else:
            raise ValueError(f"unknown score repair mode: {mode!r}")
        LOG.info(
            "score repair complete: query=%d/%d pdb_id=%s "
            "query_seconds=%.1f elapsed_seconds=%.1f",
            index,
            len(repairs),
            pdb_id,
            time.perf_counter() - query_started,
            time.perf_counter() - started,
        )


def _ligand_3d_candidate_shard_paths(data_dir: Path, shard: str) -> tuple[Path, Path]:
    root = data_dir / "scores" / "ligand_3d_candidate_shards"
    return root / f"shard={shard}.parquet", root / f"shard={shard}.json"


def _ligand_3d_pair_candidate_shard_path(data_dir: Path, shard: str) -> Path:
    return (
        data_dir
        / "scores"
        / "ligand_3d_pair_candidate_shards"
        / f"shard={shard}.parquet"
    )


def _ligand_3d_candidate_input_signatures(
    data_dir: Path, pdb_ids: Sequence[str]
) -> list[dict[str, int | str]]:
    signatures: list[dict[str, int | str]] = []
    for pdb_id in pdb_ids:
        path = (
            data_dir
            / "scores"
            / "ligand_3d_candidates"
            / "search_db=holo"
            / f"shard={pdb_id[1:3]}"
            / f"{pdb_id}.parquet"
        )
        stat = path.stat()
        schema = pq.read_schema(path)
        missing = sorted(
            set(schemas.LIGAND_3D_CANDIDATE_SCHEMA.names).difference(schema.names)
        )
        if missing:
            raise ValueError(f"candidate file {path} is missing columns {missing}")
        signatures.append(
            {
                "pdb_id": pdb_id,
                "path": str(path.resolve()),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "rows": pq.ParquetFile(path).metadata.num_rows,
            }
        )
    return signatures


def scatter_ligand_3d_candidate_shards(
    *, data_dir: Path, batch_size: int
) -> list[list[str]]:
    """Group protein-score query shards for candidate-file consolidation."""
    from plinder.data.pipeline.score import published_scoring_query_ids

    if batch_size < 1:
        raise ValueError("batch size must be positive")
    shards = sorted({pdb_id[1:3] for pdb_id in published_scoring_query_ids(data_dir)})
    return [
        shards[start : start + batch_size]
        for start in range(0, len(shards), batch_size)
    ] or [[]]


def collate_ligand_3d_candidates(
    *,
    data_dir: Path,
    shards: list[str],
    scratch_dir: Path,
    threads: int,
    replacement_query_ids: set[str] | None = None,
    source_query_ids: set[str] | None = None,
) -> list[Path]:
    """Consolidate or patch per-PDB positive-pocket candidate shards."""
    from plinder.data.pipeline.score import published_scoring_query_ids

    if threads < 1:
        raise ValueError("threads must be positive")
    import duckdb

    patch_existing = replacement_query_ids is not None
    if patch_existing:
        replacement_query_ids = set(map(str, replacement_query_ids or set()))
        active_queries = set(map(str, source_query_ids or set()))
        if active_queries.difference(replacement_query_ids):
            raise ValueError("candidate patch sources must be replacement queries")
    elif source_query_ids is not None:
        raise ValueError("source_query_ids requires replacement_query_ids")
    else:
        active_queries = published_scoring_query_ids(data_dir)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    outputs: list[Path] = []
    for shard in shards:
        if len(shard) != 2 or any(
            char not in ascii_lowercase + digits for char in shard
        ):
            raise ValueError(f"invalid ligand 3D candidate shard: {shard!r}")
        pdb_ids = sorted(pdb_id for pdb_id in active_queries if pdb_id[1:3] == shard)
        replaced_pdb_ids = sorted(
            pdb_id for pdb_id in replacement_query_ids or set() if pdb_id[1:3] == shard
        )
        if patch_existing and not replaced_pdb_ids:
            continue
        if not patch_existing and not pdb_ids:
            continue
        try:
            inputs = _ligand_3d_candidate_input_signatures(data_dir, pdb_ids)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"ligand 3D candidates are incomplete for shard {shard}: {exc}"
            ) from exc
        output, manifest = _ligand_3d_candidate_shard_paths(data_dir, shard)
        pair_output = _ligand_3d_pair_candidate_shard_path(data_dir, shard)
        if patch_existing and not output.is_file():
            raise FileNotFoundError(
                f"candidate patch requires existing shard: {output}"
            )
        output_is_current = False
        pair_output_is_current = False
        payload: dict[str, Any] = {}
        if not patch_existing and output.is_file() and manifest.is_file():
            try:
                payload = json.loads(manifest.read_text())
                stat = output.stat()
                current_output = {
                    "path": str(output.resolve()),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "rows": pq.ParquetFile(output).metadata.num_rows,
                }
                output_is_current = (
                    payload.get("shard") == shard
                    and payload.get("inputs") == inputs
                    and payload.get("output") == current_output
                )
                if output_is_current and pair_output.is_file():
                    pair_stat = pair_output.stat()
                    current_pair_output = {
                        "path": str(pair_output.resolve()),
                        "size": pair_stat.st_size,
                        "mtime_ns": pair_stat.st_mtime_ns,
                        "rows": pq.ParquetFile(pair_output).metadata.num_rows,
                    }
                    pair_schema = pq.read_schema(pair_output)
                    pair_output_is_current = (
                        set(schemas.LIGAND_3D_PAIR_CANDIDATE_SCHEMA.names).issubset(
                            pair_schema.names
                        )
                        and payload.get("pair_output") == current_pair_output
                    )
            except (OSError, TypeError, ValueError):
                output_is_current = False
                pair_output_is_current = False

        if output_is_current and pair_output_is_current:
            outputs.append(output)
            continue

        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
        temporary = scratch_dir / f"shard={shard}.parquet"
        if not output_is_current:
            source_paths = [Path(str(item["path"])) for item in inputs]
            paths_sql = ", ".join(f"'{path.as_posix()}'" for path in source_paths)
            temporary.unlink(missing_ok=True)
            if patch_existing:
                connection.register(
                    "replacement_queries",
                    pd.DataFrame({"query_entry": replaced_pdb_ids}),
                )
                replacement_sql = (
                    f"SELECT * FROM read_parquet([{paths_sql}])"
                    if source_paths
                    else f"SELECT * FROM read_parquet('{output.as_posix()}') "
                    "WHERE false"
                )
                source_sql = f"""
                    SELECT existing.*
                    FROM read_parquet('{output.as_posix()}') AS existing
                    ANTI JOIN replacement_queries USING (query_entry)
                    UNION ALL BY NAME
                    {replacement_sql}
                """
            else:
                source_sql = f"SELECT * FROM read_parquet([{paths_sql}])"
            connection.sql(
                f"""
                COPY (
                    SELECT * FROM ({source_sql})
                    ORDER BY
                        query_entry,
                        query_ligand_asym_id,
                        target_entry,
                        target_ligand_asym_id,
                        query_system,
                        target_system
                ) TO '{temporary.as_posix()}' (FORMAT PARQUET, COMPRESSION ZSTD)
                """
            )
            observed_rows = pq.ParquetFile(temporary).metadata.num_rows
            expected_rows = connection.sql(
                f"SELECT count(*) FROM ({source_sql})"
            ).fetchone()[0]
            if observed_rows != expected_rows:
                connection.close()
                raise ValueError(
                    f"ligand 3D candidate shard {shard} has {observed_rows} rows; "
                    f"expected {expected_rows}"
                )
            output.parent.mkdir(exist_ok=True, parents=True)
            install = output.with_suffix(output.suffix + ".tmp")
            copyfile(temporary, install)
            install.replace(output)
            temporary.unlink(missing_ok=True)

        pair_temporary = scratch_dir / f"pair-shard={shard}.parquet"
        pair_temporary.unlink(missing_ok=True)
        connection.sql(
            f"""
            COPY (
                SELECT
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id,
                    max(pocket_qcov)::DOUBLE AS pocket_qcov
                FROM read_parquet('{output.as_posix()}')
                GROUP BY
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id
                ORDER BY
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id
            ) TO '{pair_temporary.as_posix()}' (
                FORMAT PARQUET,
                COMPRESSION ZSTD
            )
            """
        )
        connection.close()
        pair_output.parent.mkdir(exist_ok=True, parents=True)
        pair_install = pair_output.with_suffix(pair_output.suffix + ".tmp")
        copyfile(pair_temporary, pair_install)
        pair_install.replace(pair_output)
        pair_temporary.unlink(missing_ok=True)

        output_stat = output.stat()
        observed_rows = pq.ParquetFile(output).metadata.num_rows
        pair_output_stat = pair_output.stat()
        payload = {
            "shard": shard,
            "inputs": (
                {
                    "base": "existing packed candidate shard",
                    "replaced_query_ids": replaced_pdb_ids,
                    "replacements": inputs,
                }
                if patch_existing
                else inputs
            ),
            "output": {
                "path": str(output.resolve()),
                "size": output_stat.st_size,
                "mtime_ns": output_stat.st_mtime_ns,
                "rows": observed_rows,
            },
            "pair_output": {
                "path": str(pair_output.resolve()),
                "size": pair_output_stat.st_size,
                "mtime_ns": pair_output_stat.st_mtime_ns,
                "rows": pq.ParquetFile(pair_output).metadata.num_rows,
            },
        }
        manifest.parent.mkdir(exist_ok=True, parents=True)
        temporary_manifest = manifest.with_suffix(".tmp.json")
        temporary_manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        temporary_manifest.replace(manifest)
        outputs.append(output)
    return outputs


def make_ligand_3d_scores(
    *,
    data_dir: Path,
    pairs: pd.DataFrame,
    batch_index: int,
    scorer_cfg: DictConfig,
    force_update: bool,
    scratch_dir: Path,
    threads: int,
    output_path: Path | None = None,
) -> Path:
    """Score one balanced batch of unique canonical ASU ligand pairs."""
    if batch_index < 0 or threads < 1:
        raise ValueError("batch index must be non-negative and threads positive")
    pair_columns = [
        "query_entry",
        "query_ligand_asym_id",
        "target_entry",
        "target_ligand_asym_id",
    ]
    missing = sorted(set(pair_columns).difference(pairs.columns))
    if missing:
        raise ValueError(f"ligand 3D work batch is missing columns {missing}")
    expected_pairs = pd.MultiIndex.from_frame(
        pairs[pair_columns].drop_duplicates(ignore_index=True)
    )
    if len(expected_pairs) != len(pairs):
        raise ValueError("ligand 3D work batch contains duplicate canonical pairs")
    output = output_path or (
        data_dir / "scores" / "ligand_3d_pairs" / f"{batch_index}.parquet"
    )
    if output.is_file() and not force_update:
        try:
            existing_schema = pq.read_schema(output)
            existing_pairs = pd.MultiIndex.from_frame(
                pd.read_parquet(output, columns=pair_columns)
            )
            if set(schemas.LIGAND_3D_SCORE_SCHEMA.names).issubset(
                existing_schema.names
            ) and existing_pairs.equals(expected_pairs):
                return output
        except (OSError, ValueError):
            pass

    scorer, _, _ = utils.get_scorer(
        data_dir=data_dir,
        pdb_ids=[],
        scorer_cfg=scorer_cfg,
        load_entries=False,
        scratch_dir=scratch_dir,
    )
    scorer.shape_score_threads = threads
    scores = scorer.score_canonical_ligand_pairs(data_dir, pairs)
    table = pa.Table.from_pandas(
        scores,
        schema=schemas.LIGAND_3D_SCORE_SCHEMA,
        preserve_index=False,
    )
    scratch_dir.mkdir(exist_ok=True, parents=True)
    temporary = scratch_dir / f"ligand-3d-{output.stem}.parquet"
    pq.write_table(
        table,
        temporary,
        compression="zstd",
    )
    output.parent.mkdir(exist_ok=True, parents=True)
    install = output.with_suffix(output.suffix + ".tmp")
    copyfile(temporary, install)
    install.replace(output)
    temporary.unlink(missing_ok=True)
    return output


def scatter_ligand_3d_score_batches(*, data_dir: Path) -> list[list[int]]:
    """Return the frozen canonical-pair batch indices."""
    from plinder.data.pipeline.score import _load_plan

    plan = _load_plan(data_dir)
    if "ligand_3d_batch_count" not in plan:
        raise ValueError("ligand 3D score batches have not been planned")
    return [[index] for index in range(int(plan["ligand_3d_batch_count"]))]


def scatter_ligand_3d_query_shards(
    *, data_dir: Path, batch_size: int
) -> list[list[str]]:
    """Group query shards that contain planned canonical ligand pairs."""
    from plinder.data.pipeline.score import published_scoring_query_ids

    if batch_size < 1:
        raise ValueError("batch size must be positive")
    shards = sorted({pdb_id[1:3] for pdb_id in published_scoring_query_ids(data_dir)})
    return [
        shards[start : start + batch_size]
        for start in range(0, len(shards), batch_size)
    ] or [[]]


def scatter_ligand_3d_merge(*, data_dir: Path, batch_size: int) -> list[list[str]]:
    """Group stable query shards for direct packed-score publication."""
    return scatter_ligand_3d_query_shards(data_dir=data_dir, batch_size=batch_size)


def collate_ligand_3d_scores(
    *,
    data_dir: Path,
    shards: list[str],
    scratch_dir: Path,
    threads: int,
) -> list[Path]:
    """Repartition balanced pair-score batches for fast query-PDB lookup."""
    from plinder.data.pipeline.score import (
        LIGAND_3D_WORK_RELATIVE,
        active_scoring_query_ids,
    )

    if threads < 1:
        raise ValueError("threads must be positive")
    import duckdb

    pair_columns = [
        "query_entry",
        "query_ligand_asym_id",
        "target_entry",
        "target_ligand_asym_id",
    ]
    work_path = data_dir / LIGAND_3D_WORK_RELATIVE
    pair_dir = data_dir / "scores" / "ligand_3d_pairs"
    output_dir = data_dir / "scores" / "ligand_3d_by_query"
    output_dir.mkdir(exist_ok=True, parents=True)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    active_queries = active_scoring_query_ids(data_dir)
    outputs: list[Path] = []
    for shard in shards:
        if len(shard) != 2:
            raise ValueError(f"invalid ligand 3D query shard: {shard!r}")
        shard_queries = sorted(
            pdb_id for pdb_id in active_queries if pdb_id[1:3] == shard
        )
        if not shard_queries:
            continue
        query_values = ", ".join(f"'{pdb_id}'" for pdb_id in shard_queries)
        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
        batch_rows = connection.sql(
            f"""
            SELECT DISTINCT ligand_3d_batch_index
            FROM read_parquet('{work_path.as_posix()}')
            WHERE substr(query_entry, 2, 2) = '{shard}'
              AND query_entry IN ({query_values})
            ORDER BY ligand_3d_batch_index
            """
        ).fetchall()
        batch_paths = [pair_dir / f"{int(row[0])}.parquet" for row in batch_rows]
        missing = [path for path in batch_paths if not path.is_file()]
        if missing:
            connection.close()
            raise FileNotFoundError(
                f"missing ligand 3D pair batches for shard {shard}: {missing[:10]}"
            )
        if not batch_paths:
            connection.close()
            continue
        paths_sql = ", ".join(f"'{path.as_posix()}'" for path in batch_paths)
        keys = ", ".join(pair_columns)
        counts = connection.sql(
            f"""
            WITH expected AS (
                SELECT {keys}
                FROM read_parquet('{work_path.as_posix()}')
                WHERE substr(query_entry, 2, 2) = '{shard}'
                  AND query_entry IN ({query_values})
            ), observed AS (
                SELECT scored.*
                FROM read_parquet([{paths_sql}]) AS scored
                INNER JOIN expected USING ({keys})
            )
            SELECT
                (SELECT count(*) FROM expected),
                (SELECT count(*) FROM observed),
                (SELECT count(*) - count(DISTINCT ({keys})) FROM observed)
            """
        ).fetchone()
        if counts is None:
            connection.close()
            raise RuntimeError(f"failed to validate ligand 3D query shard {shard}")
        if counts[0] != counts[1] or counts[2]:
            connection.close()
            raise ValueError(
                f"ligand 3D query shard {shard} has invalid coverage: "
                f"expected={counts[0]}, observed={counts[1]}, duplicates={counts[2]}"
            )
        temporary = scratch_dir / f"{shard}.parquet"
        temporary.unlink(missing_ok=True)
        connection.sql(
            f"""
            COPY (
                WITH expected AS (
                    SELECT {keys}
                    FROM read_parquet('{work_path.as_posix()}')
                    WHERE substr(query_entry, 2, 2) = '{shard}'
                      AND query_entry IN ({query_values})
                )
                SELECT scored.*
                FROM read_parquet([{paths_sql}]) AS scored
                INNER JOIN expected USING ({keys})
                ORDER BY {keys}
            ) TO '{temporary.as_posix()}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
        connection.close()
        output = output_dir / f"{shard}.parquet"
        install = output.with_suffix(output.suffix + ".tmp")
        copyfile(temporary, install)
        install.replace(output)
        temporary.unlink(missing_ok=True)
        outputs.append(output)
    return outputs


def _ligand_3d_query_shard_is_ready(*, data_dir: Path, shard: str) -> bool:
    """Return whether every canonical-pair batch needed by a query shard exists."""
    import duckdb

    from plinder.data.pipeline.score import (
        LIGAND_3D_WORK_RELATIVE,
        active_scoring_query_ids,
    )

    work_path = data_dir / LIGAND_3D_WORK_RELATIVE
    shard_queries = sorted(
        pdb_id for pdb_id in active_scoring_query_ids(data_dir) if pdb_id[1:3] == shard
    )
    if not shard_queries:
        return False
    query_values = ", ".join(f"'{pdb_id}'" for pdb_id in shard_queries)
    batch_rows = duckdb.sql(
        f"""
        SELECT DISTINCT ligand_3d_batch_index
        FROM read_parquet('{work_path.as_posix()}')
        WHERE substr(query_entry, 2, 2) = '{shard}'
          AND query_entry IN ({query_values})
        """
    ).fetchall()
    pair_dir = data_dir / "scores" / "ligand_3d_pairs"
    return bool(batch_rows) and all(
        (pair_dir / f"{int(row[0])}.parquet").is_file() for row in batch_rows
    )


def merge_ligand_3d_scores(
    *,
    data_dir: Path,
    shards: list[str],
    scorer_cfg: DictConfig,
    force_update: bool,
    scratch_dir: Path,
    threads: int = 1,
    reuse_cached_pairs: bool = False,
    replacement_query_ids: set[str] | None = None,
) -> list[Path]:
    """Publish complete score shards or patch selected query entries."""
    from plinder.data.pipeline.score import published_scoring_query_ids

    if threads < 1:
        raise ValueError("threads must be positive")
    import duckdb

    metric_columns = ["shape", "color", "sucos_shape"]
    all_ligand_metrics = [*metric_columns, "sucos_shape_pocket_qcov"]
    pair_keys = [
        "query_entry",
        "query_ligand_asym_id",
        "target_entry",
        "target_ligand_asym_id",
    ]
    candidate_keys = [
        "query_system",
        "query_ligand_id",
        "target_system",
        "target_ligand_id",
    ]
    outputs: list[Path] = []
    scratch_dir.mkdir(exist_ok=True, parents=True)
    output_dir = data_dir / "scores" / "search_db=holo"
    output_dir.mkdir(exist_ok=True, parents=True)
    published_queries = published_scoring_query_ids(data_dir)
    patch_existing = replacement_query_ids is not None
    replacement_query_ids = set(map(str, replacement_query_ids or set()))
    if patch_existing and not reuse_cached_pairs:
        raise ValueError("query replacement requires cached ligand-pair scores")
    active_queries = (
        published_queries.intersection(replacement_query_ids)
        if patch_existing
        else published_queries
    )
    thresholds = {
        metric: float(
            scorer_cfg.minimum_thresholds.get(metric, scorer_cfg.minimum_threshold)
        )
        for metric in all_ligand_metrics
    }

    for shard in shards:
        if len(shard) != 2 or any(
            char not in ascii_lowercase + digits for char in shard
        ):
            raise ValueError(f"invalid ligand 3D query shard: {shard!r}")
        output = output_dir / f"{shard}.parquet"
        replaced_pdb_ids = sorted(
            pdb_id for pdb_id in replacement_query_ids if pdb_id[1:3] == shard
        )
        if patch_existing and not replaced_pdb_ids:
            continue
        if patch_existing and not output.is_file():
            raise FileNotFoundError(
                f"score patch requires existing published shard: {output}"
            )
        if output.is_file() and not force_update and not patch_existing:
            schema = pq.read_schema(output)
            if set(schemas.PROTEIN_SIMILARITY_SCHEMA.names).issubset(schema.names):
                outputs.append(output)
                continue
        pair_path = data_dir / "scores" / "ligand_3d_by_query" / f"{shard}.parquet"
        if reuse_cached_pairs:
            if not pair_path.is_file():
                raise FileNotFoundError(
                    f"missing cached canonical ligand-pair scores for shard {shard}"
                )
        else:
            if not _ligand_3d_query_shard_is_ready(data_dir=data_dir, shard=shard):
                raise RuntimeError(
                    f"ligand 3D query shard {shard} is not ready; retry this merge "
                    "task after all canonical-pair batches are complete"
                )
            collate_ligand_3d_scores(
                data_dir=data_dir,
                shards=[shard],
                scratch_dir=scratch_dir / "pairs" / shard,
                threads=threads,
            )
        candidate_path = (
            data_dir
            / "scores"
            / "ligand_3d_candidate_shards"
            / f"shard={shard}.parquet"
        )
        pdb_ids = sorted(pdb_id for pdb_id in active_queries if pdb_id[1:3] == shard)
        base_paths = [
            data_dir / "dbs" / "subdbs" / "search_db=holo" / f"{pdb_id}.parquet"
            for pdb_id in pdb_ids
        ]
        missing_inputs = [
            path
            for path in [candidate_path, pair_path, *base_paths]
            if not path.is_file()
        ]
        if missing_inputs:
            raise FileNotFoundError(
                f"missing score inputs for query shard {shard}: {missing_inputs[:10]}"
            )
        base_paths_sql = ", ".join(f"'{path.as_posix()}'" for path in base_paths)
        pair_keys_sql = ", ".join(pair_keys)
        candidate_keys_sql = ", ".join(candidate_keys)
        shard_scratch = scratch_dir / shard
        shard_scratch.mkdir(exist_ok=True, parents=True)
        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{shard_scratch.as_posix()}'")
        if patch_existing:
            connection.register(
                "replacement_queries",
                pd.DataFrame({"query_entry": replaced_pdb_ids}),
            )
            candidate_sql = f"""
                SELECT candidates.*
                FROM read_parquet('{candidate_path.as_posix()}') AS candidates
                INNER JOIN replacement_queries USING (query_entry)
            """
        else:
            candidate_sql = f"SELECT * FROM read_parquet('{candidate_path.as_posix()}')"

        validation = connection.sql(
            f"""
            WITH candidates AS (
                {candidate_sql}
            ), pairs AS (
                SELECT *
                FROM read_parquet('{pair_path.as_posix()}')
            ), merged AS (
                SELECT candidates.*, pairs.query_entry AS matched_query_entry
                FROM candidates
                LEFT JOIN pairs USING ({pair_keys_sql})
            )
            SELECT
                (SELECT count(*) - count(DISTINCT ({candidate_keys_sql}))
                 FROM candidates) AS duplicate_candidates,
                (SELECT count(*) FROM merged WHERE matched_query_entry IS NULL)
                    AS missing_pairs,
                (SELECT count(*) - count(DISTINCT ({pair_keys_sql})) FROM pairs)
                    AS duplicate_pairs
            """
        ).fetchone()
        if validation is None:
            connection.close()
            raise RuntimeError(f"failed to validate ligand 3D query shard {shard}")
        if any(int(value) for value in validation):
            connection.close()
            raise ValueError(
                f"invalid ligand 3D query shard {shard}: "
                f"duplicate_candidates={validation[0]}, "
                f"missing_pairs={validation[1]}, duplicate_pairs={validation[2]}"
            )

        metric_selects = []
        for metric in all_ligand_metrics:
            value = (
                "sucos_shape * pocket_qcov"
                if metric == "sucos_shape_pocket_qcov"
                else metric
            )
            metric_selects.append(
                dedent(
                    f"""
                    SELECT
                        query_system,
                        query_ligand_id,
                        target_system,
                        target_ligand_id,
                        protein_mapping,
                        NULL::VARCHAR AS mapping,
                        protein_mapper,
                        NULL::VARCHAR AS source,
                        '{metric}'::VARCHAR AS metric,
                        round(({value}) * 100)::TINYINT AS similarity
                    FROM merged
                    WHERE {value} IS NOT NULL
                      AND {value} >= {thresholds[metric]:.17g}
                    """
                )
            )
        metric_union = "\nUNION ALL\n".join(metric_selects)
        excluded_metrics = ", ".join(f"'{metric}'" for metric in all_ligand_metrics)
        if patch_existing:
            replacement_base_sql = (
                f"""
                SELECT *
                FROM read_parquet(
                    [{base_paths_sql}],
                    hive_partitioning = false
                )
                WHERE metric NOT IN ({excluded_metrics})
                  AND NOT starts_with(cast(metric AS VARCHAR), 'protein_')
                """
                if base_paths
                else f"""
                SELECT *
                FROM read_parquet(
                    '{output.as_posix()}', hive_partitioning = false
                )
                WHERE false
                """
            )
            base_scores_sql = f"""
                SELECT existing.*
                FROM read_parquet(
                    '{output.as_posix()}', hive_partitioning = false
                ) AS existing
                ANTI JOIN replacement_queries
                  ON split_part(existing.query_system, '__', 1)
                   = replacement_queries.query_entry
                WHERE NOT starts_with(
                    cast(existing.metric AS VARCHAR), 'protein_'
                )
                UNION ALL BY NAME
                {replacement_base_sql}
            """
        else:
            base_scores_sql = f"""
                SELECT *
                FROM read_parquet(
                    [{base_paths_sql}],
                    hive_partitioning = false
                )
                WHERE metric NOT IN ({excluded_metrics})
                  AND NOT starts_with(cast(metric AS VARCHAR), 'protein_')
            """
        local_output = shard_scratch / f"{shard}.parquet"
        local_output.unlink(missing_ok=True)
        connection.sql(
            dedent(
                f"""
                COPY (
                    WITH candidates AS (
                        {candidate_sql}
                    ), pairs AS (
                        SELECT *
                        FROM read_parquet('{pair_path.as_posix()}')
                    ), merged AS (
                        SELECT candidates.*, pairs.shape, pairs.color,
                               pairs.sucos_shape
                        FROM candidates
                        INNER JOIN pairs USING ({pair_keys_sql})
                    ), ligand_scores AS (
                        {metric_union}
                    ), base_scores AS (
                        {base_scores_sql}
                    )
                    SELECT * FROM base_scores
                    UNION ALL
                    SELECT * FROM ligand_scores
                    ORDER BY
                        similarity DESC,
                        query_system,
                        query_ligand_id,
                        target_system,
                        target_ligand_id
                ) TO '{local_output.as_posix()}' (
                    FORMAT PARQUET,
                    COMPRESSION ZSTD,
                    ROW_GROUP_SIZE 500000
                )
                """
            )
        )
        connection.close()
        install = output.with_suffix(output.suffix + ".tmp")
        copyfile(local_output, install)
        install.replace(output)
        local_output.unlink(missing_ok=True)
        outputs.append(output)
    return outputs


def _write_alignment_release_shard(
    *,
    sources: list[Path],
    target: Path,
    alignment_type: str,
    temp_dir: Path,
    threads: int,
    memory_limit: str,
) -> None:
    """Sort and atomically install mapped parts as one query-shard Parquet."""
    import duckdb

    con = duckdb.connect()
    temp_dir.mkdir(exist_ok=True, parents=True)
    con.sql(f"set temp_directory='{temp_dir.as_posix()}';")
    con.sql(f"set threads={threads};")
    con.sql(f"set memory_limit='{memory_limit}';")
    target.parent.mkdir(exist_ok=True, parents=True)
    non_empty_sources = [
        path for path in sources if pq.ParquetFile(path).metadata.num_rows > 0
    ]
    temporary = temp_dir / target.name
    temporary.unlink(missing_ok=True)
    release_columns = [
        "query_entry",
        "target_entry",
        "query_chain_mapped",
        "target_chain_mapped",
        "source",
        "qcov",
        "fident",
        "seqsim",
        "query_selected_residue_numbers",
        "target_selected_residue_numbers",
        "selected_residue_identity",
    ]
    if alignment_type == "foldseek":
        release_columns.append("lddt")
    if non_empty_sources:
        source_sql = ", ".join(f"'{path.as_posix()}'" for path in non_empty_sources)
        select_columns = ", ".join(
            (
                f"CAST({column} AS INTEGER[]) AS {column}"
                if column
                in {
                    "query_selected_residue_numbers",
                    "target_selected_residue_numbers",
                }
                else f"CAST({column} AS BLOB) AS {column}"
                if column == "selected_residue_identity"
                else column
            )
            for column in release_columns
        )
        con.sql(
            dedent(
                f"""
                COPY (
                    SELECT {select_columns}
                    FROM read_parquet([{source_sql}], union_by_name = true)
                    ORDER BY query_entry, target_entry,
                             query_chain_mapped, target_chain_mapped, source
                ) TO '{temporary.as_posix()}'
                (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100_000);
                """
            )
        )
    else:
        pq.write_table(
            pa.Table.from_pylist(
                [],
                schema=schemas.mapped_alignment_schema(alignment_type=alignment_type),
            ),
            temporary,
            compression="zstd",
        )
    install_path = target.with_suffix(target.suffix + ".tmp")
    copyfile(temporary, install_path)
    install_path.replace(target)
    temporary.unlink(missing_ok=True)


def scatter_collate_alignments(*, data_dir: Path) -> list[list[str]]:
    """Return mapped alignment shards that require release collation."""
    mapped_files = (data_dir / "dbs" / "subdbs").glob("*_*/mapped_aln/*.parquet")
    shards = sorted(
        {
            path.stem[-3:-1]
            for path in mapped_files
            if not path.name.endswith(".tmp.parquet")
        }
    )
    # Preserve a join branch when there is no work, as required by Metaflow.
    return [[shard] for shard in shards] or [[]]


def collate_alignments(
    *,
    data_dir: Path,
    partition: list[str],
    scratch_dir: Path | None = None,
    threads: int = 1,
    memory_limit: str = "7GB",
) -> None:
    """Collate mapped Foldseek/MMseqs hits into query-addressable shards."""
    if not partition:
        LOG.info("collate_alignments: no mapped alignments found")
        return
    [shard] = partition

    scratch_root = scratch_dir or data_dir / "scratch" / "duckdb" / "alignments"
    temp_dir = scratch_root / shard
    for search_db in ["holo", "apo", "pred"]:
        for alignment_type in ["foldseek", "mmseqs"]:
            source_dir = (
                data_dir
                / "dbs"
                / "subdbs"
                / f"{search_db}_{alignment_type}"
                / "mapped_aln"
            )
            sources = sorted(
                path
                for path in source_dir.glob("*.parquet")
                if path.stem[-3:-1] == shard and not path.name.endswith(".tmp.parquet")
            )
            if not sources:
                continue
            target = _alignment_release_path(
                data_dir=data_dir,
                search_db=search_db,
                alignment_type=alignment_type,
                shard=shard,
            )
            _write_alignment_release_shard(
                sources=sources,
                target=target,
                alignment_type=alignment_type,
                temp_dir=temp_dir / f"{search_db}-{alignment_type}",
                threads=threads,
                memory_limit=memory_limit,
            )


def scatter_collate_partitions() -> list[list[str]]:
    partitions = [[i] for i in digits + ascii_lowercase] + [["apo"], ["pred"]]
    return partitions


def collate_partitions(
    *,
    data_dir: Path,
    partition: list[str],
    scratch_dir: Path | None = None,
    threads: int = 1,
    memory_limit: str = "7GB",
) -> None:
    """
    Collate the batch results from make_batch_scores into a partitioned
    sorted dataset using duckdb.

    Parameters
    ----------
    data_dir : Path
        plinder root dir
    partition : list[str]
        partitions to re-write
    """
    import duckdb

    part: str
    [part] = partition
    if threads < 1:
        raise ValueError("threads must be positive")
    con = duckdb.connect()
    temp_dir = (scratch_dir or data_dir / "scratch" / "duckdb") / part
    temp_dir.mkdir(exist_ok=True, parents=True)
    con.sql(f"set temp_directory='{temp_dir.as_posix()}';")
    con.sql(f"set threads={threads};")
    con.sql(f"set memory_limit='{memory_limit}';")

    search_db = "holo"
    src = f"*{part}.parquet"
    tgt = f"{part}.parquet"
    if part in ["apo", "pred"]:
        search_db = part
        src = "*.parquet"
    score_dir = data_dir / "scores" / f"search_db={search_db}"
    source_dir = data_dir / "dbs" / "subdbs" / f"search_db={search_db}"
    source = f"{source_dir}/{src}"
    target = score_dir / tgt

    if not list(source_dir.glob(src)):
        LOG.info(f"collate_partitions: no source files matching {source}")
        return
    score_dir.mkdir(exist_ok=True, parents=True)
    local_target = temp_dir / tgt
    local_target.unlink(missing_ok=True)
    con.sql(
        dedent(
            f"""
                COPY
                    (select * from '{source}')
                TO
                    '{local_target.as_posix()}'
                (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500_000);
            """
        )
    )
    con.close()
    install = target.with_suffix(target.suffix + ".tmp")
    copyfile(local_target, install)
    install.replace(target)
    local_target.unlink(missing_ok=True)


def make_linked_apo_structures(
    *,
    data_dir: Path,
    scratch_dir: Path | None = None,
    threads: int = 1,
    memory_limit: str = "7GB",
) -> Path:
    """Publish ranked deposited apo-chain links for each holo system."""
    from plinder.data.linked_apo import (
        build_apo_candidate_manifest,
        write_linked_apo_structure_table,
    )

    index_dir = data_dir / "index"
    inputs = {
        "protein_scores": data_dir / "scores/search_db=apo/apo.parquet",
        "annotation": index_dir / "annotation_table.parquet",
        "entry_chains": index_dir / "entry_chains.parquet",
        "biounit_chains": index_dir / "entry_biounit_chains.parquet",
        "entry_metadata": index_dir / "entry_metadata.parquet",
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    candidates = build_apo_candidate_manifest(
        inputs["entry_chains"],
        biounit_chains=inputs["biounit_chains"],
        entry_metadata=inputs["entry_metadata"],
        annotation=inputs["annotation"],
    )
    manifest = data_dir / "manifests/apo_candidates.parquet"
    manifest.parent.mkdir(exist_ok=True, parents=True)
    temporary_manifest = manifest.with_suffix(manifest.suffix + ".tmp")
    candidates.to_parquet(temporary_manifest, index=False)
    temporary_manifest.replace(manifest)

    output = index_dir / "linked_apo_structures.parquet"
    write_linked_apo_structure_table(
        inputs["protein_scores"],
        annotation=inputs["annotation"],
        candidates=manifest,
        output_path=output,
        scratch_dir=scratch_dir,
        threads=threads,
        memory_limit=memory_limit,
    )
    linked_rows = pq.ParquetFile(output).metadata.num_rows
    LOG.info(
        "make_linked_apo_structures: selected %d links from " "%d apo-chain candidates",
        linked_rows,
        len(candidates),
    )
    return output


def scatter_component_reduction_sources(
    *,
    data_dir: Path,
    metrics: list[str],
    batch_size: int,
    entity_type: clusters.ClusterEntity = "ligand",
) -> list[list[str]]:
    """Scatter each physical score shard once, independent of threshold count."""
    if batch_size < 1:
        raise ValueError("component reduction source batch size must be positive")
    sources: set[Path] = set()
    for metric in metrics:
        sources.update(
            clusters.component_score_sources(
                data_dir=data_dir,
                metric=metric,
                entity_type=entity_type,
            )
        )
    ordered = [str(path) for path in sorted(sources)]
    return [
        ordered[start : start + batch_size]
        for start in range(0, len(ordered), batch_size)
    ] or [[]]


def make_symmetric_edge_fragments(
    *,
    data_dir: Path,
    batches: list[dict[str, Any]],
    scratch_dir: Path,
    threads: int,
    force_update: bool,
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Map raw score batches into hash-partitioned canonical-pair fragments."""
    scratch_dir.mkdir(exist_ok=True, parents=True)
    for batch_index, batch in enumerate(batches):
        if not force_update and clusters.symmetric_edge_fragment_batch_is_complete(
            data_dir=data_dir,
            batch=batch,
            entity_type=entity_type,
        ):
            LOG.info(
                "symmetric fragment batch already complete: progress=%d/%d key=%s",
                batch_index + 1,
                len(batches),
                batch["key"],
            )
            continue
        sources = [Path(str(value["path"])) for value in batch["sources"]]
        local_sources: list[Path] = []
        batch_scratch = scratch_dir / f"batch-{batch['key']}"
        batch_scratch.mkdir(exist_ok=True, parents=True)
        try:
            for source_index, source in enumerate(sources):
                local_source = batch_scratch / f"{source_index:04d}-{source.name}"
                copyfile(source, local_source)
                local_sources.append(local_source)
            LOG.info(
                "symmetric fragment batch copied: progress=%d/%d key=%s " "sources=%d",
                batch_index + 1,
                len(batches),
                batch["key"],
                len(sources),
            )
            clusters.write_symmetric_edge_fragment_batch(
                data_dir=data_dir,
                batch=batch,
                scratch_dir=batch_scratch,
                threads=threads,
                force_update=force_update,
                read_paths=local_sources,
                entity_type=entity_type,
            )
        finally:
            for local_source in local_sources:
                local_source.unlink(missing_ok=True)


def make_symmetric_edge_shards(
    *,
    data_dir: Path,
    metric_buckets: list[tuple[str, int]],
    scratch_dir: Path,
    threads: int,
    force_update: bool,
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Merge reciprocal fragment values into compact minimum edge shards."""
    for index, (metric, bucket) in enumerate(metric_buckets, start=1):
        LOG.info(
            "symmetric edge shard start: progress=%d/%d metric=%s bucket=%d",
            index,
            len(metric_buckets),
            metric,
            bucket,
        )
        clusters.write_symmetric_edge_shard(
            data_dir=data_dir,
            metric=metric,
            bucket=bucket,
            scratch_dir=scratch_dir / f"{metric}-{bucket:03d}",
            threads=threads,
            force_update=force_update,
            entity_type=entity_type,
        )


def scatter_make_set_covers(
    *,
    data_dir: Path,
    metrics: list[str],
    thresholds: list[int],
    stop_on_cluster: int,
    skip_existing_clusters: bool,
    entity_type: clusters.ClusterEntity = "ligand",
) -> list[list[tuple[str, int]]]:
    """Scatter ligand-Tanimoto set-cover work after component reduction."""
    values = [
        [(metric, threshold)]
        for metric in metrics
        if entity_type == "ligand" and metric == "tanimoto_similarity_ecfp4_1024"
        for threshold in thresholds
    ]
    if stop_on_cluster:
        values = values[:stop_on_cluster]
    if not skip_existing_clusters:
        return values
    pending = []
    for item in values:
        metric, threshold = item[0]
        output = (
            clusters._sampling_root(data_dir, entity_type)
            / "set_cover"
            / f"metric={metric}"
            / f"threshold={threshold}.parquet"
        )
        if not clusters.set_cover_is_complete(output):
            pending.append(item)
    return pending or [[]]


def scatter_make_directed_set_covers(
    *,
    data_dir: Path,
    metrics: list[str],
    thresholds: list[int],
    stop_on_cluster: int,
    skip_existing: bool,
    entity_type: clusters.ClusterEntity = "ligand",
) -> list[list[tuple[str, int]]]:
    """Scatter directed centroid-cover work after connectivity publication."""
    values = [
        [(metric, threshold)]
        for metric in metrics
        if metric != "tanimoto_similarity_ecfp4_1024"
        for threshold in thresholds
    ]
    if stop_on_cluster:
        values = values[:stop_on_cluster]
    if not skip_existing:
        return values
    pending = []
    for item in values:
        metric, threshold = item[0]
        output = (
            clusters._sampling_root(data_dir, entity_type)
            / "directed_set_cover"
            / f"metric={metric}"
            / f"threshold={threshold}.parquet"
        )
        if not clusters.directed_set_cover_is_complete(
            output,
            entity_type=entity_type,
        ):
            pending.append(item)
    return pending or [[]]


_COMPONENT_REDUCTION_CONTEXT: dict[str, Any] | None = None


def _initialize_component_reduction_worker(context: dict[str, Any]) -> None:
    """Install immutable per-source state once in each metric worker."""
    global _COMPONENT_REDUCTION_CONTEXT
    _COMPONENT_REDUCTION_CONTEXT = context


def _reduce_component_metric(metric_index: int, metric: str) -> dict[str, Any]:
    """Reduce reciprocal and any-direction connectivity from one metric shard."""
    context = _COMPONENT_REDUCTION_CONTEXT
    if context is None:
        raise RuntimeError("component reduction worker context is not initialized")
    metric_started = time.time()
    chemical = metric == "tanimoto_similarity_ecfp4_1024"
    if chemical:
        manifests = [
            clusters.make_score_component_reduction(
                data_dir=context["data_dir"],
                metric=metric,
                thresholds=context["thresholds"],
                source_path=context["source"],
                read_path=context["local_source"],
                all_nodes=context["nodes"],
                eligible_systems=context["eligible_systems"],
                force_update=context["force_update"],
                entity_type=context["entity_type"],
            )
        ]
    else:
        manifests = [
            clusters.make_directed_cover_component_reduction(
                data_dir=context["data_dir"],
                metric=metric,
                thresholds=context["thresholds"],
                source_path=context["source"],
                read_path=context["local_source"],
                all_nodes=context["nodes"],
                eligible_systems=context["eligible_systems"],
                force_update=context["force_update"],
                entity_type=context["entity_type"],
            )
        ]
    return {
        "metric_index": metric_index,
        "metric": metric,
        "output_rows": sum(
            int(output["rows"])
            for manifest in manifests
            for output in manifest["outputs"]
        ),
        "elapsed_seconds": time.time() - metric_started,
    }


def make_component_reductions(
    *,
    data_dir: Path,
    source_paths: list[str],
    metrics: list[str],
    thresholds: list[int],
    scratch_dir: Path,
    force_update: bool,
    metric_workers: int = 1,
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Map reciprocal-minimum edge sources to exact component reductions."""
    if not source_paths:
        return
    if metric_workers < 1:
        raise ValueError("component reduction metric workers must be positive")
    chemical_metric = "tanimoto_similarity_ecfp4_1024"
    nonchemical_metrics = [metric for metric in metrics if metric != chemical_metric]
    generic_nodes: list[str] = []
    generic_systems: set[str] | None = None
    if nonchemical_metrics and entity_type != "interface":
        generic_nodes, generic_systems = clusters.component_node_universe(
            data_dir=data_dir,
            metric=nonchemical_metrics[0],
            entity_type=entity_type,
        )
    interface_nodes: dict[str, list[str]] = {}
    if entity_type == "interface":
        for metric in nonchemical_metrics:
            interface_nodes[metric], _ = clusters.component_node_universe(
                data_dir=data_dir,
                metric=metric,
                entity_type=entity_type,
            )
    chemical_nodes: list[str] = []
    if chemical_metric in metrics:
        chemical_nodes, _ = clusters.component_node_universe(
            data_dir=data_dir,
            metric=chemical_metric,
            entity_type=entity_type,
        )
    scratch_dir.mkdir(exist_ok=True, parents=True)

    for source_index, source_value in enumerate(source_paths):
        source = Path(source_value)
        source_metric = source.parent.name.removeprefix("metric=")
        if source_metric not in metrics:
            raise ValueError(
                f"component source metric is not selected: {source_metric}"
            )
        is_chemical = source_metric == chemical_metric
        source_metrics = [source_metric]
        nodes = (
            chemical_nodes
            if is_chemical
            else interface_nodes.get(source_metric, generic_nodes)
        )
        eligible_systems = None if is_chemical else generic_systems
        pending_metrics = [
            metric
            for metric in source_metrics
            if force_update
            or not (
                clusters.score_component_reduction_is_complete(
                    data_dir=data_dir,
                    metric=metric,
                    thresholds=thresholds,
                    source_path=source,
                    all_nodes=nodes,
                    eligible_systems=eligible_systems,
                    entity_type=entity_type,
                )
                if is_chemical
                else clusters.directed_cover_component_reduction_is_complete(
                    data_dir=data_dir,
                    metric=metric,
                    thresholds=thresholds,
                    source_path=source,
                    all_nodes=nodes,
                    eligible_systems=eligible_systems,
                    entity_type=entity_type,
                )
            )
        ]
        if not pending_metrics:
            LOG.info(
                "component reduction source already complete: source=%s metrics=%d",
                source.name,
                len(source_metrics),
            )
            continue
        source_started = time.time()
        LOG.info(
            "component reduction source start: source=%s size_bytes=%d "
            "pending_metrics=%d/%d nodes=%d",
            source.name,
            source.stat().st_size,
            len(pending_metrics),
            len(source_metrics),
            len(nodes),
        )
        local_source = scratch_dir / f"{source_index}-{source.name}"
        copyfile(source, local_source)
        try:
            LOG.info(
                "component reduction copied source to scratch: source=%s "
                "elapsed_seconds=%.1f",
                source.name,
                time.time() - source_started,
            )
            context = {
                "data_dir": data_dir,
                "thresholds": thresholds,
                "source": source,
                "local_source": local_source,
                "nodes": nodes,
                "eligible_systems": eligible_systems,
                "force_update": force_update,
                "entity_type": entity_type,
            }

            def log_metric_start(metric_index: int, metric: str) -> None:
                LOG.info(
                    "component reduction metric start: source=%s metric=%s "
                    "progress=%d/%d workers=%d",
                    source.name,
                    metric,
                    metric_index,
                    len(pending_metrics),
                    min(metric_workers, len(pending_metrics)),
                )

            def log_metric_complete(result: dict[str, Any]) -> None:
                LOG.info(
                    "component reduction metric complete: source=%s metric=%s "
                    "output_rows=%d elapsed_seconds=%.1f",
                    source.name,
                    result["metric"],
                    result["output_rows"],
                    result["elapsed_seconds"],
                )

            workers = min(metric_workers, len(pending_metrics))
            if workers == 1:
                for metric_index, metric in enumerate(pending_metrics, start=1):
                    log_metric_start(metric_index, metric)
                    _initialize_component_reduction_worker(context)
                    log_metric_complete(_reduce_component_metric(metric_index, metric))
            else:
                for metric_index, metric in enumerate(pending_metrics, start=1):
                    log_metric_start(metric_index, metric)
                with ProcessPoolExecutor(
                    max_workers=workers,
                    initializer=_initialize_component_reduction_worker,
                    initargs=(context,),
                ) as executor:
                    futures = [
                        executor.submit(_reduce_component_metric, metric_index, metric)
                        for metric_index, metric in enumerate(pending_metrics, start=1)
                    ]
                    wait(futures, return_when=ALL_COMPLETED)
                    for future in futures:
                        log_metric_complete(future.result())
        finally:
            local_source.unlink(missing_ok=True)
        LOG.info(
            "component reduction source complete: source=%s metrics=%d "
            "elapsed_seconds=%.1f",
            source.name,
            len(pending_metrics),
            time.time() - source_started,
        )


def merge_component_reductions(
    *,
    data_dir: Path,
    metrics: list[str],
    thresholds: list[int],
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Merge source reductions into internal connectivity labels."""
    started = time.time()
    for index, metric in enumerate(metrics, start=1):
        metric_started = time.time()
        LOG.info(
            "merge_component_reductions: metric=%s progress=%d/%d",
            metric,
            index,
            len(metrics),
        )
        if metric == "tanimoto_similarity_ecfp4_1024":
            clusters.merge_score_component_reductions(
                data_dir=data_dir,
                metric=metric,
                thresholds=thresholds,
                entity_type=entity_type,
            )
        else:
            clusters.merge_directed_cover_component_reductions(
                data_dir=data_dir,
                metric=metric,
                thresholds=thresholds,
                entity_type=entity_type,
            )
        elapsed = time.time() - started
        rate = index / elapsed
        LOG.info(
            "merge_component_reductions: metric=%s complete "
            "elapsed_seconds=%.1f overall_eta_seconds=%.1f",
            metric,
            time.time() - metric_started,
            (len(metrics) - index) / rate,
        )


def make_set_covers(
    *,
    data_dir: Path,
    metric_threshold: list[tuple[str, int]],
    skip_existing_clusters: bool,
    scratch_dir: Path | None = None,
    threads: int = 1,
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Compute one ligand-Tanimoto set cover after component reduction."""
    if not metric_threshold:
        LOG.info("make_set_covers: all set covers are cached")
        return
    [(metric, threshold)] = metric_threshold
    clusters.make_set_cover(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        skip_existing_clusters=skip_existing_clusters,
        scratch_dir=scratch_dir,
        threads=threads,
        entity_type=entity_type,
    )


def make_directed_set_covers(
    *,
    data_dir: Path,
    metric_threshold: list[tuple[str, int]],
    skip_existing: bool,
    scratch_dir: Path | None = None,
    threads: int = 1,
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Compute one directed cover for annotation and training-set sampling."""
    if not metric_threshold:
        LOG.info("make_directed_set_covers: all covers are cached")
        return
    [(metric, threshold)] = metric_threshold
    clusters.make_directed_set_cover(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        skip_existing=skip_existing,
        scratch_dir=scratch_dir,
        threads=threads,
        entity_type=entity_type,
    )


def summarize_clusters(
    *,
    data_dir: Path,
    metrics: list[str],
    thresholds: list[int],
    entity_type: clusters.ClusterEntity = "ligand",
) -> None:
    """Validate and summarize every cluster artifact before index publication."""
    from plinder.data.pipeline.score import summarize_clustering_artifacts

    summarize_clustering_artifacts(
        data_dir,
        metrics=metrics,
        thresholds=thresholds,
        entity_type=entity_type,
    )


def finalize_index(*, data_dir: Path) -> None:
    """Merge locally generated cluster IDs into the annotation indexes."""
    lookup_was_current = _completed_alignment_chain_lookup(data_dir) is not None
    utils.finalize_index(data_dir=data_dir)
    # finalize_index only adds ligand and cluster annotations; it preserves the
    # entry, chain, pocket, and interface fields from which the normalized
    # inputs were built. Refresh those source signatures so an otherwise valid
    # mapped release does not become stale merely because clusters were published.
    if lookup_was_current:
        _refresh_representative_source_manifests(data_dir)
        _write_alignment_chain_lookup_manifest(data_dir)
    collate.finalize_repair_marker(data_dir)
