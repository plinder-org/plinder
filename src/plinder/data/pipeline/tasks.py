# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import os
from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from pathlib import Path
from shutil import rmtree
from string import ascii_lowercase, digits
from textwrap import dedent
from typing import Any
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
from omegaconf import DictConfig
from tqdm import tqdm

from plinder.core.scores.metrics import is_ligand_level_metric
from plinder.core.utils import gcs
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
STAGES = [
    "download_rcsb_files",
    "download_alternative_datasets",
    "make_dbs",
    "make_entries",
    "collate_entries",
    "make_canonical_ligand_archives",
    "make_ligands",
    "compute_ligand_fingerprints",
    "make_ligand_scores",
    "annotate_ligand_similarity",
    "make_sub_dbs",
    "run_batch_searches",
    "make_batch_scores",
    "collate_alignments",
    "collate_partitions",
    "make_components_and_communities",
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


def make_dbs(*, data_dir: Path, sub_databases: list[str], cpu: int) -> None:
    """
    Make the foldseek and mmseqs dbs

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    """
    input_dirs = {}
    if "apo" in sub_databases or "holo" in sub_databases:
        input_dirs["foldseek"] = data_dir / "ingest"
        input_dirs["mmseqs"] = io.download_seqres_data(data_dir=data_dir)
    if "pred" in sub_databases:
        input_dirs["pred_mmseqs"] = io.download_uniprot_fasta_data(data_dir=data_dir)
        input_dirs["pred_foldseek"] = io.download_alphafold_cif_files(data_dir=data_dir)
    for db, source in input_dirs.items():
        output_dir = data_dir / "dbs" / db
        LOG.info(f"make_dbs: making {db} in {output_dir}")
        databases.make_db(source, output_dir, db.split("_")[-1], threads=cpu)


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
) -> None:
    """
    Archive canonical ASU ligand SDFs for a two-character code.
    """
    for code in two_char_codes:
        entry_dir = data_dir / "raw_entries" / code
        archive = data_dir / "ligand_archives" / f"{code}.zip"
        archive.parent.mkdir(exist_ok=True, parents=True)
        with ZipFile(archive.as_posix(), "w", compression=ZIP_DEFLATED) as zip_archive:
            for entry_parquet in sorted(entry_dir.glob("*.parquet")):
                ligand_dir = entry_dir / entry_parquet.stem / "ligand_files"
                for ligand_file in sorted(ligand_dir.glob("*.sdf")):
                    zip_archive.write(
                        ligand_file,
                        ligand_file.relative_to(entry_dir),
                    )


def make_sub_dbs(
    *,
    data_dir: Path,
    sub_databases: list[str],
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
    from plinder.core.scores.entries import entry_views_from_df

    entries = entry_views_from_df(
        pd.read_parquet(data_dir / "index" / "annotation_table.parquet"),
        entry_chains=pd.read_parquet(data_dir / "index" / "entry_chains.parquet"),
    )
    db_dir = data_dir / "dbs" / "subdbs"
    db_dir.mkdir(exist_ok=True)
    LOG.info("making sub-databases for scoring")
    db_sources = utils.get_db_sources(data_dir=data_dir, sub_databases=sub_databases)
    databases.make_sub_dbs(db_dir, db_sources, entries)


def scatter_make_ligands(
    *,
    data_dir: Path,
    batch_size: int,
    two_char_codes: list[str],
    pdb_ids: list[str],
) -> list[list[str]]:
    """ """
    pdb_ids = utils.get_local_contents(
        data_dir=data_dir / "ingest",
        two_char_codes=two_char_codes,
        pdb_ids=pdb_ids,
        as_four_char_ids=True,
    )
    LOG.info(f"scatter_make_ligands: found {len(pdb_ids)} PDBs from ingest")
    pdb_ids = [
        pdb_id
        for pdb_id in pdb_ids
        if utils.entry_exists(
            entry_dir=data_dir / "raw_entries",
            pdb_id=pdb_id,
        )
    ]
    LOG.info(f"scatter_make_ligands: found {len(pdb_ids)} PDBs from raw_entries")
    return [
        pdb_ids[pos : pos + batch_size] for pos in range(0, len(pdb_ids), batch_size)
    ]


def make_ligands(
    *,
    data_dir: Path,
    pdb_ids: list[str],
) -> None:
    """ """
    if not pdb_ids:
        LOG.info("make_ligands: no entries to process")
        return
    annotations = [
        pd.read_parquet(data_dir / "raw_entries" / pdb_id[-3:-1] / f"{pdb_id}.parquet")
        for pdb_id in pdb_ids
    ]
    annotation = pd.concat(annotations, ignore_index=True)
    hashed_contents = utils.hash_contents(pdb_ids)
    output_dir = data_dir / "ligands"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / f"{hashed_contents}.parquet"
    LOG.info("make_ligands: running save_ligand_batch")
    utils.save_ligand_batch(
        data_dir=data_dir,
        annotation=annotation,
        output_path=output_path,
    )


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
    utils.update_index_ligand_3d_score_ability(data_dir=data_dir)


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
    """
    Split all the PDB IDs in the dataset
    to be used in score generation.

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
    pdb_ids = utils.get_local_contents(
        data_dir=data_dir / "ingest",
        two_char_codes=two_char_codes,
        pdb_ids=pdb_ids,
        as_four_char_ids=True,
    )
    pdb_ids = [
        pdb_id
        for pdb_id in pdb_ids
        if utils.entry_exists(
            entry_dir=data_dir / "raw_entries",
            pdb_id=pdb_id,
        )
    ]
    LOG.info(f"scatter_make_scorers: found {len(pdb_ids)} pdb IDs")
    return [
        pdb_ids[pos : pos + batch_size] for pos in range(0, len(pdb_ids), batch_size)
    ]


def run_batch_searches(
    *,
    data_dir: Path,
    # TODO: : use large batches for run_batch_searches
    pdb_ids: list[str],
    scorer_cfg: DictConfig,
    cpu: int,
) -> None:
    scorer, entry_ids, batch_db_dir = utils.get_scorer(
        data_dir=data_dir,
        pdb_ids=pdb_ids,
        scorer_cfg=scorer_cfg,
        load_entries=True,
    )
    # TODO: convert get_similarity_scores.run_alignment to
    #       accept FoldseekConfig / MMSeqsConfig
    for search_db in scorer_cfg.sub_databases:
        LOG.info(f"make_scorers: run_alignments for {search_db}")
        scorer.run_alignments(
            entry_ids=entry_ids,
            output_folder=batch_db_dir,
            overwrite=True,
            search_db=search_db,
            threads=cpu,
        )
    rmtree(batch_db_dir)


def scatter_missing_scores(
    *,
    data_dir: Path,
    batch_size: int,
) -> list[list[str]]:
    present = utils.get_pdb_ids_in_scoring_dataset(data_dir=data_dir)
    alns = utils.get_alns(data_dir=data_dir, mapped=True)
    rerun = set()
    for search_db, pdb_ids in present.items():
        aln_foldseek = alns[search_db]["foldseek"]
        aln_mmseqs = alns[search_db]["mmseqs"]
        aln = set(aln_foldseek).union(aln_mmseqs)
        rerun |= aln.difference(pdb_ids)
    run = sorted(rerun)
    return [run[pos : pos + batch_size] for pos in range(0, len(run), batch_size)]


def make_batch_scores(
    *,
    data_dir: Path,
    pdb_ids: list[str],
    scorer_cfg: DictConfig,
    force_update: bool,
) -> None:
    scorer, entry_ids, _ = utils.get_scorer(
        data_dir=data_dir,
        pdb_ids=pdb_ids,
        scorer_cfg=scorer_cfg,
        load_entries=False,
    )
    for search_db in scorer_cfg.sub_databases:
        for pdb_id in tqdm(entry_ids):
            scorer.get_score_df(
                data_dir, pdb_id, search_db=search_db, overwrite=force_update
            )


def scatter_collate_alignments(*, data_dir: Path) -> list[list[str]]:
    """Return deterministic PDB two-character shards with mapped alignments."""
    mapped_files = (data_dir / "dbs" / "subdbs").glob("*_*/mapped_aln/*.parquet")
    shards = sorted({path.stem[-3:-1] for path in mapped_files})
    # Preserve a join branch when there is no work, as required by Metaflow.
    return [[shard] for shard in shards] or [[]]


def collate_alignments(*, data_dir: Path, partition: list[str]) -> None:
    """Collate mapped Foldseek/MMseqs hits into query-addressable shards."""
    if not partition:
        LOG.info("collate_alignments: no mapped alignments found")
        return
    [shard] = partition

    import duckdb

    con = duckdb.connect()
    temp_dir = data_dir / "scratch" / "duckdb" / "alignments" / shard
    temp_dir.mkdir(exist_ok=True, parents=True)
    con.sql(f"set temp_directory='{temp_dir.as_posix()}';")
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
                if path.stem[-3:-1] == shard
            )
            if not sources:
                continue
            target = (
                data_dir
                / "alignments"
                / f"search_db={search_db}"
                / f"alignment_type={alignment_type}"
                / f"shard={shard}.parquet"
            )
            target.parent.mkdir(exist_ok=True, parents=True)
            source_sql = ", ".join(f"'{path.as_posix()}'" for path in sources)
            temporary = target.with_suffix(".tmp.parquet")
            if temporary.is_file():
                temporary.unlink()
            con.sql(
                dedent(
                    f"""
                    COPY (
                        SELECT *
                        FROM read_parquet([{source_sql}])
                        ORDER BY query_entry, target_entry,
                                 query_chain_mapped, target_chain_mapped, source
                    ) TO '{temporary.as_posix()}'
                    (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100_000);
                    """
                )
            )
            temporary.replace(target)


def scatter_collate_partitions() -> list[list[str]]:
    partitions = [[i] for i in digits + ascii_lowercase] + [["apo"], ["pred"]]
    return partitions


def collate_partitions(*, data_dir: Path, partition: list[str]) -> None:
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
    con = duckdb.connect()
    temp_dir = data_dir / "scratch" / "duckdb" / part
    temp_dir.mkdir(exist_ok=True, parents=True)
    con.sql(f"set temp_directory='{temp_dir.as_posix()}';")

    search_db = "holo"
    src = f"*{part}.parquet"
    tgt = f"{part}.parquet"
    if part in ["apo", "pred"]:
        search_db = part
        src = "*.parquet"
    score_dir = data_dir / "scores" / f"search_db={search_db}"
    source_dir = data_dir / "dbs" / "subdbs" / f"search_db={search_db}"
    source = f"{source_dir}/{src}"
    target = f"{score_dir}/{tgt}"

    if not list(source_dir.glob(src)):
        LOG.info(f"collate_partitions: no source files matching {source}")
        return
    score_dir.mkdir(exist_ok=True, parents=True)
    if Path(target).is_file():
        Path(target).unlink()
    con.sql(
        dedent(
            f"""
                COPY
                    (select * from '{source}')
                TO
                    '{target}'
                (FORMAT PARQUET, ROW_GROUP_SIZE 500_000);
            """
        )
    )


def scatter_make_components_and_communities(
    *,
    data_dir: Path,
    metrics: list[str],
    thresholds: list[int],
    stop_on_cluster: int,
    skip_existing_clusters: bool,
) -> list[list[tuple[str, int]]]:
    values = [[(metric, threshold)] for metric in metrics for threshold in thresholds]
    if stop_on_cluster:
        values = values[:stop_on_cluster]
    if not skip_existing_clusters:
        return values

    rerun = []
    for tup in values:
        metric, threshold = tup[0]
        roots = ["clusters"]
        if is_ligand_level_metric(metric):
            roots.append("ligand_clusters")
        expected = [
            data_dir
            / root
            / f"cluster={cluster}"
            / f"directed={directed}"
            / f"metric={metric}"
            / f"threshold={threshold}.parquet"
            for root in roots
            for cluster, directed in [
                ("components", True),
                ("components", False),
                ("communities", False),
            ]
        ]
        if not all(path.is_file() for path in expected):
            rerun.append(tup)
    # Metaflow needs at least one foreach branch in order to reach the join.
    # An empty chunk is consumed as a no-op when every cluster is cached.
    return rerun or [[]]


def make_components_and_communities(
    *,
    data_dir: Path,
    metric_threshold: list[tuple[str, int]],
    skip_existing_clusters: bool,
) -> None:
    if not metric_threshold:
        LOG.info("make_components_and_communities: all clusters are cached")
        return
    [(metric, threshold)] = metric_threshold
    clusters.make_components_and_communities(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        # directed=True,
        skip_existing_clusters=skip_existing_clusters,
    )


def finalize_index(*, data_dir: Path) -> None:
    """Merge locally generated cluster IDs into the annotation indexes."""
    utils.finalize_index(data_dir=data_dir)
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
