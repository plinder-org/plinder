# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

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

from plinder.core.utils import gcs, schemas
from plinder.core.utils.log import setup_logger
from plinder.data import clusters, databases, splits
from plinder.data.pipeline import collate, io, utils
from plinder.data.pipeline.ingest import (
    balance_entries,
    completed_entry_metrics,
    discover_entries,
    ingest_pdb_batch,
    normalize_pdb_id,
)
from plinder.data.utils.annotations import get_similarity_scores

LOG = setup_logger(__name__)
ALIGNMENT_CHAIN_LOOKUP_RELATIVE = Path("index/alignment_chain_lookup.parquet")
ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE = Path(
    "index/alignment_chain_lookup.manifest.json"
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
    "make_sub_dbs",
    "run_batch_searches",
    "map_batch_alignments",
    "collate_alignments",
    "finalize_alignments",
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
    "make_component_reductions",
    "merge_component_reductions",
    "make_communities",
    "summarize_clusters",
    "finalize_index",
    "make_mmp_index",
    "make_splits",
    "make_links",
    "make_linked_structures",
    "score_linked_structures",
]


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
        if create and (force_update or not complete):
            databases.create_db(
                source,
                working_output_dir,
                database_type,
                threads=cpu,
            )
        elif create:
            LOG.info(f"make_dbs: reusing completed {database_path}")
        if index:
            index_output_dir = (
                working_output_dir if force_update or not complete else output_dir
            )
            databases.create_db_index(
                index_output_dir,
                database_type,
                tmp_dir=tmp_dir,
                threads=cpu,
            )
        if build_dir is not None and create and (force_update or not complete):
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
) -> list[list[str]]:
    """Discover and size-balance source entries for V3 annotation."""
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
            if completed_entry_metrics(data_dir, entry.pdb_id) is None
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
    cpu: int = 1,
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
    if set(sub_databases) == {"holo"}:
        chains = pd.read_parquet(
            data_dir / "index" / "entry_chains.parquet",
            columns=[
                "entry_pdb_id",
                "chain_auth_id",
                "chain_receptor_type",
                "chain_is_holo",
            ],
        )
        chains = chains[
            chains["chain_receptor_type"].fillna("").astype(str).eq("protein")
            & chains["chain_is_holo"].fillna(False).astype(bool)
            & chains["chain_auth_id"].notna()
        ]
        identifiers_by_database = {
            "holo_foldseek": {
                f"pdb_0000{row.entry_pdb_id}_xyz-enrich_{row.chain_auth_id}"
                for row in chains.itertuples(index=False)
            },
            "holo_mmseqs": {
                f"{row.entry_pdb_id}_{row.chain_auth_id}"
                for row in chains.itertuples(index=False)
            },
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
    if set(sub_databases) == {"holo"}:
        make_alignment_chain_lookup(
            data_dir=data_dir,
            scratch_dir=scratch_dir,
            threads=cpu,
        )


def _completed_alignment_chain_lookup(
    data_dir: Path,
) -> dict[str, int | str] | None:
    lookup = data_dir / ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    manifest = data_dir / ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE
    try:
        stat = lookup.stat()
        columns = set(pq.read_schema(lookup).names)
        input_signatures = _alignment_chain_lookup_input_signatures(data_dir)
    except (OSError, TypeError, ValueError):
        return None
    expected_columns = {
        "entry_pdb_id",
        "chain_asym_id",
        "chain_auth_id",
        "pocket_residue_numbers",
        "pocket_residue_indices",
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
    for name in ["annotation_table.parquet", "entry_chains.parquet"]:
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


def make_alignment_chain_lookup(
    *,
    data_dir: Path,
    scratch_dir: Path | None,
    threads: int,
    force_update: bool = False,
) -> Path:
    """Build the compact author-chain and pocket map used by every map shard."""
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
    annotation = (data_dir / "index" / "annotation_table.parquet").as_posix()
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
                pocket_residues AS (
                    SELECT
                        a.entry_pdb_id,
                        split_part(split_part(neighbor, '_', 1), '.', 2)
                            AS chain_asym_id,
                        CAST(split_part(neighbor, '_', 2) AS INTEGER)
                            AS residue_number,
                        CAST(split_part(neighbor, '_', 3) AS INTEGER)
                            AS residue_index
                    FROM read_parquet('{annotation}') AS a,
                    UNNEST(a.ligand_neighboring_residues) AS residues(neighbor)
                    WHERE a.ligand_is_proper
                ),
                canonical_pocket_residues AS (
                    SELECT
                        entry_pdb_id,
                        chain_asym_id,
                        residue_index,
                        min(residue_number) AS residue_number
                    FROM pocket_residues
                    GROUP BY entry_pdb_id, chain_asym_id, residue_index
                ),
                pocket_mapping AS (
                    SELECT
                        entry_pdb_id,
                        chain_asym_id,
                        list(residue_number ORDER BY residue_index)
                            AS pocket_residue_numbers,
                        list(residue_index ORDER BY residue_index)
                            AS pocket_residue_indices
                    FROM canonical_pocket_residues
                    GROUP BY entry_pdb_id, chain_asym_id
                )
                SELECT
                    c.entry_pdb_id,
                    c.chain_asym_id,
                    c.chain_auth_id,
                    coalesce(p.pocket_residue_numbers, []::INTEGER[])
                        AS pocket_residue_numbers,
                    coalesce(p.pocket_residue_indices, []::INTEGER[])
                        AS pocket_residue_indices
                FROM protein_chains c
                LEFT JOIN pocket_mapping p
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


def annotate_ligand_similarity(
    *, data_dir: Path, cluster_threshold: float = 90.0
) -> None:
    """Add 90%-component frequencies after all BulkTanimoto shards finish."""
    get_similarity_scores.annotate_ligand_similarity(
        data_dir=data_dir,
        cluster_threshold=cluster_threshold,
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
    chain_path = data_dir / "index" / "entry_chains.parquet"
    chains = pd.read_parquet(
        chain_path,
        columns=["entry_pdb_id", "chain_receptor_type", "chain_is_holo"],
    )
    protein_entries = set(
        chains.loc[
            chains["chain_receptor_type"].fillna("").astype(str).eq("protein")
            & chains["chain_is_holo"].fillna(False).astype(bool),
            "entry_pdb_id",
        ].astype(str)
    )
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
            target_database = (
                data_dir
                / "dbs"
                / "subdbs"
                / f"{search_db}_{alignment_type}"
                / f"{search_db}_{alignment_type}"
            )
            identifiers = databases.database_identifiers(target_database)
            if alignment_type == "foldseek":
                eligible_queries = {
                    identifier.replace("pdb_0000", "", 1)[:4]
                    for identifier in identifiers
                }
            else:
                eligible_queries = {
                    identifier.split("_", maxsplit=1)[0] for identifier in identifiers
                }
            output_dir = (
                data_dir / "dbs" / "subdbs" / f"{search_db}_{alignment_type}" / "aln"
            )
            pending = [
                pdb_id
                for pdb_id in pdb_ids
                if pdb_id in eligible_queries
                and (force_update or not (output_dir / f"{pdb_id}.parquet").is_file())
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
                    f"{len(missing_outputs)} eligible query entries: "
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
    *, data_dir: Path, batch_size: int
) -> list[list[str]]:
    """Scatter query shards whose raw alignments are not mapped and published."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    raw_root = data_dir / "dbs" / "subdbs"
    shards = sorted(
        {
            path.stem[1:3]
            for alignment_type in ["foldseek", "mmseqs"]
            for path in (raw_root / f"holo_{alignment_type}" / "aln").glob("*.parquet")
            if not path.name.endswith(".tmp.parquet")
        }
    )
    missing = [
        shard
        for shard in shards
        if not alignment_mapping_shard_is_current(data_dir=data_dir, shard=shard)
    ]
    chunks = [
        missing[pos : pos + batch_size] for pos in range(0, len(missing), batch_size)
    ]
    return chunks or [[]]


def _alignment_input_signatures(
    *, data_dir: Path, shard: str
) -> dict[str, list[dict[str, int | str]]]:
    signatures: dict[str, list[dict[str, int | str]]] = {}
    for alignment_type in ["foldseek", "mmseqs"]:
        source_dir = data_dir / "dbs" / "subdbs" / f"holo_{alignment_type}" / "aln"
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


def _alignment_mapping_manifest_path(*, data_dir: Path, shard: str) -> Path:
    return data_dir / "alignments" / "manifests" / f"shard={shard}.json"


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


def alignment_mapping_shard_is_current(*, data_dir: Path, shard: str) -> bool:
    """Validate a shard manifest against raw inputs and published outputs."""
    manifest_path = _alignment_mapping_manifest_path(data_dir=data_dir, shard=shard)
    try:
        payload = json.loads(manifest_path.read_text())
    except (OSError, TypeError, ValueError):
        return False
    inputs = _alignment_input_signatures(data_dir=data_dir, shard=shard)
    lookup_signature = _completed_alignment_chain_lookup(data_dir)
    if lookup_signature is None:
        return False
    if (
        payload.get("shard") != shard
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
            search_db="holo",
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
) -> None:
    """Map raw backend hits directly into atomic query-shard release files."""
    if list(scorer_cfg.sub_databases) != ["holo"]:
        raise ValueError("V3 sharded alignment mapping currently supports holo only")
    maximum_rows = int(getattr(scorer_cfg, "max_alignment_rows_per_query", 5_000_000))
    for shard in shards:
        if not force_update and alignment_mapping_shard_is_current(
            data_dir=data_dir, shard=shard
        ):
            LOG.info(f"map_batch_alignments: shard {shard} is complete")
            continue
        inputs = _alignment_input_signatures(data_dir=data_dir, shard=shard)
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
            source_root = data_dir / "dbs" / "subdbs" / f"holo_{alignment_type}" / "aln"
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
            raw_root / f"holo_{alignment_type}" / "aln" / str(signature["name"])
            for alignment_type, signatures in inputs.items()
            for signature in signatures
            if Path(str(signature["name"])).stem not in skipped_queries
        ]
        mapping_entry_ids = set(mapped_pdb_ids) | _alignment_target_entry_ids(
            raw_sources
        )
        working_root = (scratch_dir or data_dir / "scratch" / "mapping") / shard
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
                    "holo",
                    overwrite=True,
                    scratch_dir=working_root / "temporary",
                    mapped_db_dir=mapped_db_dir,
                )
                if len(mapped) != expected:
                    raise RuntimeError(
                        f"holo alignment mapping for {pdb_id} produced "
                        f"{len(mapped)} of {expected} available backends"
                    )
            outputs: dict[str, dict[str, int | str] | None] = {}
            for alignment_type, source_signatures in inputs.items():
                target = _alignment_release_path(
                    data_dir=data_dir,
                    search_db="holo",
                    alignment_type=alignment_type,
                    shard=shard,
                )
                if not source_signatures:
                    target.unlink(missing_ok=True)
                    outputs[alignment_type] = None
                    continue
                local_sources = sorted(
                    (mapped_db_dir / f"holo_{alignment_type}" / "mapped_aln").glob(
                        "*.parquet"
                    )
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
            manifest = _alignment_mapping_manifest_path(data_dir=data_dir, shard=shard)
            manifest.parent.mkdir(exist_ok=True, parents=True)
            temporary_manifest = manifest.with_suffix(".tmp.json")
            temporary_manifest.write_text(
                json.dumps(
                    {
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
    mapped_queries: set[str] = set()
    for manifest_path in sorted(
        (data_dir / "alignments" / "manifests").glob("shard=*.json")
    ):
        shard = manifest_path.stem.removeprefix("shard=")
        if not alignment_mapping_shard_is_current(data_dir=data_dir, shard=shard):
            continue
        payload = json.loads(manifest_path.read_text())
        mapped_queries.update(
            Path(str(signature["name"])).stem
            for alignment_type in ["foldseek", "mmseqs"]
            for signature in payload["inputs"][alignment_type]
        )
    rerun = mapped_queries.difference(present["holo"]).difference(
        dropped_query_ids(data_dir)
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
    for repair in repairs:
        pdb_id = str(repair["pdb_id"])
        mode = str(repair["repair_mode"])
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
        elif mode == "targets":
            if scorer is None:
                raise RuntimeError("score repair unexpectedly lacks a scorer")
            scorer.repair_score_df_targets(
                data_dir,
                pdb_id,
                affected_target_entries=set(map(str, repair["target_pdb_ids"])),
                scratch_dir=scratch_dir,
            )
        else:
            raise ValueError(f"unknown score repair mode: {mode!r}")


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
    from plinder.data.pipeline.score import active_scoring_query_ids

    if batch_size < 1:
        raise ValueError("batch size must be positive")
    shards = sorted({pdb_id[1:3] for pdb_id in active_scoring_query_ids(data_dir)})
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
) -> list[Path]:
    """Consolidate per-PDB positive-pocket candidates into query shards."""
    from plinder.data.pipeline.score import active_scoring_query_ids

    if threads < 1:
        raise ValueError("threads must be positive")
    import duckdb

    active_queries = active_scoring_query_ids(data_dir)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    outputs: list[Path] = []
    for shard in shards:
        if len(shard) != 2 or any(
            char not in ascii_lowercase + digits for char in shard
        ):
            raise ValueError(f"invalid ligand 3D candidate shard: {shard!r}")
        pdb_ids = sorted(pdb_id for pdb_id in active_queries if pdb_id[1:3] == shard)
        if not pdb_ids:
            continue
        try:
            inputs = _ligand_3d_candidate_input_signatures(data_dir, pdb_ids)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"ligand 3D candidates are incomplete for shard {shard}: {exc}"
            ) from exc
        output, manifest = _ligand_3d_candidate_shard_paths(data_dir, shard)
        pair_output = _ligand_3d_pair_candidate_shard_path(data_dir, shard)
        output_is_current = False
        pair_output_is_current = False
        payload: dict[str, Any] = {}
        if output.is_file() and manifest.is_file():
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
            connection.sql(
                f"""
                COPY (
                    SELECT *
                    FROM read_parquet([{paths_sql}])
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
            expected_rows = sum(int(item["rows"]) for item in inputs)
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
            "inputs": inputs,
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
    from plinder.data.pipeline.score import active_scoring_query_ids

    if batch_size < 1:
        raise ValueError("batch size must be positive")
    shards = sorted({pdb_id[1:3] for pdb_id in active_scoring_query_ids(data_dir)})
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
) -> list[Path]:
    """Publish complete V3 scores directly as immutable query shards."""
    from plinder.data.pipeline.score import active_scoring_query_ids

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
    active_queries = active_scoring_query_ids(data_dir)
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
        if output.is_file() and not force_update:
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

        validation = connection.sql(
            f"""
            WITH candidates AS (
                SELECT *
                FROM read_parquet('{candidate_path.as_posix()}')
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
        local_output = shard_scratch / f"{shard}.parquet"
        local_output.unlink(missing_ok=True)
        connection.sql(
            dedent(
                f"""
                COPY (
                    WITH candidates AS (
                        SELECT *
                        FROM read_parquet('{candidate_path.as_posix()}')
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
                        SELECT *
                        FROM read_parquet(
                            [{base_paths_sql}],
                            hive_partitioning = false
                        )
                        WHERE metric NOT IN ({excluded_metrics})
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
        "query_pocket_residue_numbers",
        "target_pocket_residue_numbers",
        "pocket_residue_identity",
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
                    "query_pocket_residue_numbers",
                    "target_pocket_residue_numbers",
                }
                else f"CAST({column} AS BLOB) AS {column}"
                if column == "pocket_residue_identity"
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


def scatter_component_reduction_sources(
    *,
    data_dir: Path,
    metrics: list[str],
    batch_size: int,
) -> list[list[str]]:
    """Scatter each physical score shard once, independent of threshold count."""
    if batch_size < 1:
        raise ValueError("component reduction source batch size must be positive")
    sources: set[Path] = set()
    for metric in metrics:
        sources.update(
            clusters.component_score_sources(data_dir=data_dir, metric=metric)
        )
    ordered = [str(path) for path in sorted(sources)]
    return [
        ordered[start : start + batch_size]
        for start in range(0, len(ordered), batch_size)
    ] or [[]]


def scatter_make_communities(
    *,
    data_dir: Path,
    metrics: list[str],
    thresholds: list[int],
    stop_on_cluster: int,
    skip_existing_clusters: bool,
) -> list[list[tuple[str, int]]]:
    """Scatter PLM work after exact weak component labels are available."""
    values = [[(metric, threshold)] for metric in metrics for threshold in thresholds]
    if stop_on_cluster:
        values = values[:stop_on_cluster]
    if not skip_existing_clusters:
        return values
    pending = []
    for item in values:
        metric, threshold = item[0]
        output = (
            data_dir
            / "ligand_clusters"
            / "cluster=communities"
            / "directed=False"
            / f"metric={metric}"
            / f"threshold={threshold}.parquet"
        )
        if not output.is_file():
            pending.append(item)
    return pending or [[]]


_COMPONENT_REDUCTION_CONTEXT: dict[str, Any] | None = None


def _initialize_component_reduction_worker(context: dict[str, Any]) -> None:
    """Install immutable per-source state once in each metric worker."""
    global _COMPONENT_REDUCTION_CONTEXT
    _COMPONENT_REDUCTION_CONTEXT = context


def _reduce_component_metric(metric_index: int, metric: str) -> dict[str, Any]:
    """Reduce one metric using state shared at process initialization."""
    context = _COMPONENT_REDUCTION_CONTEXT
    if context is None:
        raise RuntimeError("component reduction worker context is not initialized")
    metric_started = time.time()
    manifest = clusters.make_score_component_reduction(
        data_dir=context["data_dir"],
        metric=metric,
        thresholds=context["thresholds"],
        source_path=context["source"],
        read_path=context["local_source"],
        all_nodes=context["nodes"],
        eligible_systems=context["eligible_systems"],
        force_update=context["force_update"],
    )
    return {
        "metric_index": metric_index,
        "metric": metric,
        "output_rows": sum(int(output["rows"]) for output in manifest["outputs"]),
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
) -> None:
    """Map score sources to exact weak/strong reductions using local scratch."""
    if not source_paths:
        return
    if metric_workers < 1:
        raise ValueError("component reduction metric workers must be positive")
    chemical_metric = "tanimoto_similarity_ecfp4_1024"
    chemical_sources = set(
        clusters.component_score_sources(data_dir=data_dir, metric=chemical_metric)
        if chemical_metric in metrics
        else []
    )
    nonchemical_metrics = [metric for metric in metrics if metric != chemical_metric]
    generic_nodes: list[str] = []
    generic_systems: set[str] | None = None
    if nonchemical_metrics:
        generic_nodes, generic_systems = clusters.component_node_universe(
            data_dir=data_dir,
            metric=nonchemical_metrics[0],
        )
    chemical_nodes: list[str] = []
    if chemical_metric in metrics:
        chemical_nodes, _ = clusters.component_node_universe(
            data_dir=data_dir,
            metric=chemical_metric,
        )
    scratch_dir.mkdir(exist_ok=True, parents=True)

    for source_index, source_value in enumerate(source_paths):
        source = Path(source_value)
        is_chemical = source in chemical_sources
        source_metrics = [chemical_metric] if is_chemical else nonchemical_metrics
        nodes = chemical_nodes if is_chemical else generic_nodes
        eligible_systems = None if is_chemical else generic_systems
        pending_metrics = [
            metric
            for metric in source_metrics
            if force_update
            or not clusters.score_component_reduction_is_complete(
                data_dir=data_dir,
                metric=metric,
                thresholds=thresholds,
                source_path=source,
                all_nodes=nodes,
                eligible_systems=eligible_systems,
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
) -> None:
    """Merge every expected source reduction and publish ligand components."""
    started = time.time()
    for index, metric in enumerate(metrics, start=1):
        metric_started = time.time()
        LOG.info(
            "merge_component_reductions: metric=%s progress=%d/%d",
            metric,
            index,
            len(metrics),
        )
        clusters.merge_score_component_reductions(
            data_dir=data_dir,
            metric=metric,
            thresholds=thresholds,
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


def make_communities(
    *,
    data_dir: Path,
    metric_threshold: list[tuple[str, int]],
    skip_existing_clusters: bool,
    scratch_dir: Path | None = None,
    threads: int = 1,
) -> None:
    """Compute one community result after component publication."""
    if not metric_threshold:
        LOG.info("make_communities: all communities are cached")
        return
    [(metric, threshold)] = metric_threshold
    clusters.make_communities(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        skip_existing_clusters=skip_existing_clusters,
        scratch_dir=scratch_dir,
        threads=threads,
    )


def summarize_clusters(
    *, data_dir: Path, metrics: list[str], thresholds: list[int]
) -> None:
    """Validate and summarize every cluster artifact before index publication."""
    from plinder.data.pipeline.score import summarize_clustering_artifacts

    summarize_clustering_artifacts(
        data_dir,
        metrics=metrics,
        thresholds=thresholds,
    )


def finalize_index(*, data_dir: Path) -> None:
    """Merge locally generated cluster IDs into the annotation indexes."""
    lookup_was_current = _completed_alignment_chain_lookup(data_dir) is not None
    utils.finalize_index(data_dir=data_dir)
    # finalize_index only adds ligand and cluster annotations; it preserves the
    # entry, chain, and pocket fields from which the lookup was built. Refresh
    # that one input signature so an otherwise valid mapped release does not
    # become stale merely because cluster columns were published.
    if lookup_was_current:
        _write_alignment_chain_lookup_manifest(data_dir)
    utils.create_nonredundant_dataset(data_dir=data_dir)


def make_mmp_index(
    *,
    data_dir: Path,
) -> None:
    """
    Get the list of all pdb IDs to load all the entries
    for mmp indexing.
    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    """

    from plinder.data.utils.annotations.mmpdb_utils import (
        add_mmp_clusters_to_data,
        make_mmp_index_from_annotation_table,
    )

    LOG.info("making annotation table (and non-redundant) indexes")
    utils.create_nonredundant_dataset(data_dir=data_dir)

    LOG.info("making mmp index for all entries")
    annotation_index = data_dir / "index" / "annotation_table.parquet"
    columns = ["system_id", "ligand_rdkit_canonical_smiles", "ligand_unique_ccd_code"]
    annotation_df = pd.read_parquet(annotation_index, columns=columns)
    mmp_df_path = make_mmp_index_from_annotation_table(data_dir, annotation_df)
    load_mmp_df = pd.read_csv(mmp_df_path, compression="gzip", header=None, sep="\t")
    load_mmp_df.columns = ["SMILES1", "SMILES2", "id1", "id2", "V1>>V2", "CONSTANT"]
    mmp_data = add_mmp_clusters_to_data(
        load_mmp_df,
        annotation_df,
        cluster_folder=data_dir / "clusters",
    )
    mmp_data.to_parquet(data_dir / "mmp" / "plinder_mmp_series.parquet", index=False)


def scatter_make_splits(
    *,
    data_dir: Path,
    split_config_dir: str,
) -> list[list[tuple[DictConfig, str]]]:
    # defaults to empty string so skip it
    configs: list[list[tuple[DictConfig, str]]]
    if not len(split_config_dir):
        configs = [[]]
    # allow configs living in cloud buckets configured by PLINDER_BUCKET
    elif split_config_dir.startswith("gs:"):
        bucket_name = Path(split_config_dir).parts[1]
        configs = [
            [
                (
                    splits.get_config(
                        gcs.download_as_str(
                            gcs_path=cloud_path,
                            bucket_name=bucket_name,
                        )
                    ),
                    cloud_path,
                )
            ]
            for cloud_path in gcs.list_dir(
                gcs_path=str(split_config_dir),
                bucket_name=bucket_name,
            )
        ]
    else:
        # support relative local split_config_dir from data_dir
        # and absolute split_config_dir
        split_dir = Path(split_config_dir or "splits")
        if not split_dir.is_absolute():
            split_dir = data_dir / split_config_dir
        configs = [
            [
                (
                    splits.get_config(path.read_text()),
                    path.as_posix(),
                )
            ]
            for path in split_dir.rglob("*.yaml")
        ]
    for tup in configs:
        if len(tup[0]):
            LOG.info(f"scatter_make_splits: config={tup[0][1]}")
    return configs


def make_splits(
    *,
    data_dir: Path,
    cfg_and_path: list[tuple[DictConfig, str]],
) -> None:
    [(cfg, path)] = cfg_and_path
    splits.split(data_dir=data_dir, cfg=cfg, relpath=path)


def scatter_make_links(
    *,
    data_dir: Path,
    search_dbs: list[str],
) -> list[list[str]]:
    return [[obj] for obj in search_dbs]


def make_links(
    *,
    data_dir: Path,
    search_dbs: list[str],
    cpu: int = 8,
) -> None:
    from plinder.data.save_linked_structures import make_linked_structures_data_file

    save_dir = data_dir / "assignments"
    linked_structures = data_dir / "linked_staging"
    for search_db in search_dbs:
        output_file = linked_structures / f"{search_db}_links.parquet"
        make_linked_structures_data_file(
            data_dir=data_dir,
            search_db=search_db,
            superposed_folder=save_dir,
            output_file=output_file,
            num_processes=cpu,
        )


def make_linked_structures(
    *,
    data_dir: Path,
    search_dbs: list[str],
    cpu: int = 8,
    force_update: bool = False,
) -> None:
    import multiprocessing

    linked_structures = data_dir / "linked_staging"
    for search_db in search_dbs:
        if search_db == "holo":
            continue
        source_structures = linked_structures / "source" / search_db
        source_structures.mkdir(exist_ok=True, parents=True)
        df = pd.read_parquet(
            linked_structures / f"{search_db}_links.parquet", columns=["id"]
        )
        LOG.info(
            f"make_linked_structures: collecting {df['id'].nunique()} {search_db} linked structures"
        )
        func = None
        if search_db == "apo":
            func = utils.apo_file_from_link_id
        elif search_db == "pred":
            func = utils.pred_file_from_link_id
        if func is not None:
            args = [
                (data_dir, source_structures, link_id, force_update)
                for link_id in df["id"].unique()
            ]
            with multiprocessing.get_context("spawn").Pool(cpu) as p:
                p.starmap(func, args)
        utils.pack_source_structures(data_dir, search_db)


def scatter_score_linked_structures(
    *,
    data_dir: Path,
    search_dbs: list[str],
    batch_size: int,
) -> list[list[tuple[str, str]]]:
    items = []
    for search_db in search_dbs:
        links = pd.read_parquet(
            data_dir / "linked_staging" / f"{search_db}_links.parquet"
        )
        items.extend(
            [
                (search_db, system_id)
                for system_id in sorted(links["reference_system_id"])
            ]
        )
    return [items[pos : pos + batch_size] for pos in range(0, len(items), batch_size)]


def score_linked_structures(
    *,
    data_dir: Path,
    search_dbs: list[str],
    system_ids: list[tuple[str, str]],
    cpu: int = 8,
    force_update: bool = False,
) -> None:
    import multiprocessing

    from plinder.data.save_linked_structures import (
        system_save_and_score_representatives,
    )

    linked_structures = data_dir / "linked_staging"
    grouped = {
        search_db: [tup[1] for tup in system_ids if tup[0] == search_db]
        for search_db in search_dbs
    }
    dfs = []
    for search_db in search_dbs:
        df = pd.read_parquet(linked_structures / f"{search_db}_links.parquet")
        slc = df[df["reference_system_id"].isin(grouped[search_db])]
        if not slc.empty:
            dfs.append(slc.copy())
            dfs[-1]["kind"] = search_db
    if not len(dfs):
        LOG.info("no linked structures to make")
        return
    links = pd.concat(dfs).reset_index(drop=True)

    with multiprocessing.get_context("spawn").Pool(cpu) as p:
        p.starmap(
            system_save_and_score_representatives,
            [
                (system, group, data_dir, search_db, linked_structures, force_update)
                for (search_db, system), group in links.groupby(
                    [
                        "kind",
                        "reference_system_id",
                    ]
                )
            ],
        )
