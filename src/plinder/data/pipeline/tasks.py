# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import os
from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from pathlib import Path
from shutil import rmtree
from string import ascii_lowercase, digits
from subprocess import check_output
from textwrap import dedent
from typing import Any
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
from omegaconf import DictConfig
from tqdm import tqdm

from plinder.core.utils import gcs
from plinder.core.utils.log import setup_logger
from plinder.data import clusters, databases, leakage, splits
from plinder.data.pipeline import io, utils
from plinder.data.utils import tanimoto

LOG = setup_logger(__name__)
STAGES = [
    "download_rcsb_files",
    "download_alternative_datasets",
    "make_dbs",
    "make_entries",
    "structure_qc",
    "make_system_archives",
    "make_ligands",
    "compute_ligand_fingerprints",
    "make_ligand_scores",
    "make_sub_dbs",
    "run_batch_searches",
    "make_batch_scores",
    "collate_partitions",
    "make_mmp_index",
    "make_components_and_communities",
    "make_splits",
    "compute_ligand_leakage",
    "compute_protein_leakage",
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
    with ThreadPoolExecutor() as executor:
        futures: list[Future[Any]] = [
            executor.submit(io.download_cofactors, **kws),
            executor.submit(io.download_seqres_data, **kws),
            executor.submit(io.download_kinase_data, **kws),
            executor.submit(io.download_panther_data, **kws),
            executor.submit(io.download_components_cif, **kws),
            executor.submit(io.download_ecod_data, **kws),
            executor.submit(io.download_affinity_data, **kws),
        ]
        wait(futures, return_when=ALL_COMPLETED)
        for future in futures:
            exc = future.exception()
            if exc is not None:
                raise exc


def make_dbs(*, data_dir: Path, sub_databases: list[str]) -> None:
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
        databases.make_db(source, output_dir, db.split("_")[-1])


def scatter_make_entries(
    *,
    data_dir: Path,
    batch_size: int,
    two_char_codes: list[str],
    pdb_ids: list[str],
    force_update: bool,
) -> list[list[str]]:
    """
    Distribute annotation generation by pdb id rather than
    two character code for more symmetric distributed
    processing.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    batch_size : int
        how many codes to put in a chunk
    force_update : bool
    two_char_codes : list[str], default=[]
        only consider particular codes
    pdb_ids : list[str], default=[]
        only consider particular pdb IDs
    force_update : bool
        if True, re-process existing entries
    """
    pdb_dirs = utils.get_local_contents(
        data_dir=data_dir / "ingest",
        two_char_codes=two_char_codes,
        pdb_ids=pdb_ids,
    )
    if not force_update:
        pdb_dirs = [
            pdb_dir
            for pdb_dir in pdb_dirs
            if not utils.entry_exists(
                entry_dir=data_dir / "raw_entries",
                pdb_id=pdb_dir[-4:],
            )
        ]
    LOG.info(f"scatter_make_entries: found {len(pdb_dirs)} PDBs")
    return [
        pdb_dirs[pos : pos + batch_size] for pos in range(0, len(pdb_dirs), batch_size)
    ]


def make_entries(
    *,
    data_dir: Path,
    pdb_dirs: list[str],
    force_update: bool,
    annotation_cfg: DictConfig,
    entry_cfg: DictConfig,
    cpu: int = 1,
) -> list[str]:
    """
    Offload individual plinder annotation tasks to
    a multiprocessing queue wrapped as a subprocess.
    This is a simple way to gracefully handle seg-faults
    from C-extensions without impacting later entries in
    the list of pdb_dirs. Additionally keep a record of
    entries which failed so that they can be re-processed
    with larger resource requests.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    pdb_dirs : list[str]
        list of pdb directories to process
    force_update : bool
        if True, force re-processing
    annotation_cfg : DictConfig
        from plinder.data.pipeline.config.AnnotationConfig
    entry_cfg : DictConfig
        from plinder.data.pipeline.config.EntryConfig
    cpu : int, default=1
        number of CPUs to use

    Returns
    -------
    failed : list[str]
        list of pdb directories to re-process
    """
    input_dir = data_dir / "ingest"
    report_dir = data_dir / "reports"
    output_dir = data_dir / "raw_entries"
    scratch_dir = data_dir / "scratch" / "raw_entries"
    LOG.info(f"making {len(pdb_dirs)} entries in {output_dir}")
    output_dir.mkdir(exist_ok=True, parents=True)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    hashed_contents = utils.hash_contents(pdb_dirs)
    scratch_tasks = scratch_dir / f"tasks_{hashed_contents}.txt"
    check_finished = []
    with scratch_tasks.open("w") as f:
        for pdb_dir in pdb_dirs:
            two_char_code = pdb_dir[-3:-1]
            pdb_id = pdb_dir[-4:]
            output = output_dir / two_char_code / (pdb_id + ".json")
            if not force_update and output.is_file():
                LOG.info(f"skipping {pdb_id} since entry exists already")
                continue
            output.parent.mkdir(exist_ok=True, parents=True)
            check_finished.append((pdb_dir, output))
            # cifs are in ingest/{two_char_code}/pdb_0000{pdb_id}/
            # but vals are in reports/{two_char_code}/{pdb_id}/
            # and not all cifs have vals so gracefully handle
            [mmcif] = list((input_dir / two_char_code / pdb_dir).glob(io.CIF_GLOB))
            try:
                [report] = list((report_dir / two_char_code / pdb_id).glob(io.VAL_GLOB))
            except Exception:
                report = (
                    report_dir / two_char_code / pdb_id / (pdb_id + io.VAL_GLOB[1:])
                )
            cmd = (
                [
                    "python",
                    "-m",
                    "plinder.data.get_system_annotations",
                    f"mmcif_file={mmcif.as_posix()}",
                    f"validation_xml={report.as_posix()}",
                ]
                + [f"annotation.{k}={v}" for k, v in annotation_cfg.items()]
                + [f"entry.{k}={v}" for k, v in entry_cfg.items() if k != "save_folder"]
                + [f"entry.save_folder={output.parent.as_posix()}"]
            )
            f.write(" ".join(cmd) + "\n")
    try:
        check_output(
            [
                "python",
                "-m",
                "plinder.data.pipeline.mpqueue",
                f"{scratch_tasks.as_posix()}",
                f"--cores={max(1, cpu - 1)}",
            ],
            text=True,
            timeout=10800,
        )
        scratch_tasks.unlink()
    except Exception:
        LOG.error("beep boop mpqueue timed out")

    rerun = []
    for pdb_dir, output in check_finished:
        if output.is_file():
            continue
        rerun.append(pdb_dir)

    fail_dir = data_dir / "failed_entries"
    fail_dir.mkdir(exist_ok=True, parents=True)
    LOG.info(f"would rerun {len(rerun)} systems")
    with (fail_dir / f"fails_{hashed_contents}.txt").open("w") as f:
        for item in rerun:
            f.write(f"{item}\n")
    return rerun


def scatter_structure_qc(
    *,
    data_dir: Path,
    two_char_codes: list[str],
    batch_size: int,
) -> list[list[str]]:
    """
    Scatter two character codes for system archive generation

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
    chunks : list[list[str]]
        batches of two character codes
    """
    entry_dir = data_dir / "raw_entries"
    if len(two_char_codes):
        codes = sorted(two_char_codes)
    else:
        codes = sorted(os.listdir(entry_dir.as_posix()))
    LOG.info(f"scatter_structure_qc: found {len(codes)} two character codes")
    return [codes[pos : pos + batch_size] for pos in range(0, len(codes), batch_size)]


def structure_qc(
    *,
    data_dir: Path,
    two_char_codes: list[str],
) -> None:
    from plinder.data.final_structure_qc import (
        prepare_system_dict,
        run_structure_checks,
    )
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    entry_dir = data_dir / "raw_entries"
    err_dir = data_dir / "qc" / "logs"
    zip_dir = data_dir / "entries"
    pqt_dir = data_dir / "qc" / "index"
    err_dir.mkdir(exist_ok=True, parents=True)
    zip_dir.mkdir(exist_ok=True, parents=True)
    pqt_dir.mkdir(exist_ok=True, parents=True)
    for code in two_char_codes:
        with ZipFile(zip_dir / f"{code}.zip", "w", compression=ZIP_DEFLATED) as archive:
            with (err_dir / f"{code}_qc_fails.csv").open("w") as fails:
                fails.write("entry_json,error\n")
                LOG.info(f"structure_qc: two_char_code={code}")
                entry_dfs = []
                structure_qc = []
                for i, entry_json in enumerate((entry_dir / code).glob("*json")):
                    if not i % 25:
                        LOG.info(f"on entry={i} {entry_json}")
                    try:
                        entry = Entry.from_json(
                            entry_json, clear_non_pocket_residues=True
                        )
                    except Exception as e:
                        LOG.warn(f"failed loading {entry_json}")
                        clean = str(e).replace(",", "_").replace("\n", " ")[:50]
                        fails.write(f"{entry_json},{clean}")
                        continue
                    if len(entry.systems):
                        entry_dfs.append(entry.to_df())
                    archive.writestr(entry_json.name, entry.model_dump_json())
                    system_structure_path = entry_dir / code
                    try:
                        for system_dict in prepare_system_dict(
                            system_structure_path,
                            entry,
                        ):
                            structure_qc.extend(run_structure_checks(system_dict))
                    except Exception as e:
                        LOG.warn(f"failed structure checks for {entry_json}")
                        fails.write(f"{entry_json},{str(e).replace(',', '_')}\n")
                        continue
                if len(entry_dfs) and len(structure_qc):
                    entry_df = pd.concat(entry_dfs).reset_index(drop=True)
                    structure_df = pd.DataFrame(structure_qc)
                    pd.merge(
                        entry_df,
                        structure_df,
                        on=["system_id", "ligand_instance", "ligand_asym_id"],
                        how="left",
                    ).to_parquet(pqt_dir / f"{code}.parquet", index=False)


def scatter_make_system_archives(
    *,
    data_dir: Path,
    two_char_codes: list[str],
    batch_size: int,
) -> list[list[str]]:
    """
    Scatter two character codes for system archive generation

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
    chunks : list[list[str]]
        batches of two character codes
    """
    entry_dir = data_dir / "raw_entries"
    if len(two_char_codes):
        codes = sorted(two_char_codes)
    else:
        codes = sorted(os.listdir(entry_dir.as_posix()))
    LOG.info(f"scatter_make_system_archives: found {len(codes)} two character codes")
    return [codes[pos : pos + batch_size] for pos in range(0, len(codes), batch_size)]


def make_system_archives(
    *,
    data_dir: Path,
    two_char_codes: list[str],
) -> None:
    """
    Create a zip file of all the system files in a given
    two character code to reduce pressure on the network
    for large-scale file IO.
    """
    for code in two_char_codes:
        entry_dir = data_dir / "raw_entries" / code
        archive = data_dir / "archives" / f"{code}.zip"
        archive.parent.mkdir(exist_ok=True, parents=True)
        with ZipFile(archive.as_posix(), "w", compression=ZIP_DEFLATED) as zip_archive:
            for subdir in os.listdir(entry_dir):
                if os.path.isdir(f"{entry_dir}/{subdir}"):
                    for root, _, files in os.walk(f"{entry_dir}/{subdir}"):
                        for file in files:
                            full_path = os.path.join(root, file)
                            parts = full_path.split(os.sep)
                            arc_name = os.path.join(*parts[parts.index(code) + 1 :])
                            zip_archive.write(full_path, arc_name)


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
    entries = utils.load_entries_from_zips(data_dir=data_dir)
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
    entries = utils.load_entries_from_zips(data_dir=data_dir, pdb_ids=pdb_ids)
    hashed_contents = utils.hash_contents(pdb_ids)
    output_dir = data_dir / "ligands"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / f"{hashed_contents}.parquet"
    LOG.info("make_ligands: running save_ligand_batch")
    utils.save_ligand_batch(
        entries=entries,
        output_path=output_path,
    )


def compute_ligand_fingerprints(
    *,
    data_dir: Path,
    split_char: str = "__",
    radius: int = 2,
    nbits: int = 1024,
) -> None:
    """ """
    LOG.info("compute_ligand_fingerprints: running")
    #  data_dir / "fingerprints" / ligands_per_system.parquet
    #  data_dir / "fingerprints"  / ligands_per_inchikey.parquet
    #  data_dir / "fingerprints" /ligands_per_inchikey_ecfp4.npy
    tanimoto.compute_ligand_fingerprints(
        data_dir=data_dir,
        split_char=split_char,
        radius=radius,
        nbits=nbits,
    )


def scatter_make_ligand_scores(
    *,
    data_dir: Path,
    batch_size: int,
    number_id_col: str = "number_id_by_inchikeys",
) -> list[list[int]]:
    """ """
    ligands = pd.read_parquet(
        data_dir / "fingerprints" / "ligands_per_inchikey.parquet",
        columns=[number_id_col],
    )[number_id_col].to_list()
    LOG.info(f"scatter_make_ligand_scores: found {len(ligands)} ligands")
    return [
        ligands[pos : pos + batch_size] for pos in range(0, len(ligands), batch_size)
    ]


def make_ligand_scores(
    *,
    data_dir: Path,
    ligand_ids: list[int],
    save_top_k_similar_ligands: int = 5000,
    multiply_by: int = 100,
    number_id_col: str = "number_id_by_inchikeys",
) -> None:
    hashid = utils.hash_contents([str(i) for i in ligand_ids])
    output_path = data_dir / "ligand_scores" / f"{hashid}.parquet"
    output_path.parent.mkdir(exist_ok=True, parents=True)
    tanimoto.ligand_scores(
        ligand_ids=ligand_ids,
        data_dir=data_dir,
        output_path=output_path,
        number_id_col=number_id_col,
        save_top_k_similar_ligands=save_top_k_similar_ligands,
        multiply_by=multiply_by,
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


def scatter_collate_partitions() -> list[list[str]]:
    partitions = [[i] for i in digits + ascii_lowercase] + [["apo"], ["pred"]]
    return partitions[:1]


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
    con.sql(f"set temp_directory='/plinder/tmp/{part}';")

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

    score_dir.mkdir(exist_ok=True, parents=True)
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

    rerun = []
    for tup in values:
        metric, threshold = tup[0]
        do = False
        for cluster in ["components", "communities"]:
            for directed in [True, False]:
                if do:
                    continue
                if not (
                    data_dir
                    / f"clusters/cluster={cluster}/directed={directed}/metric={metric}/threshold={threshold}.parquet"
                ).is_file() and not skip_existing_clusters:
                    do = True
        if do:
            rerun.append(tup)
    if stop_on_cluster:
        values = values[:stop_on_cluster]
    return values


def make_components_and_communities(
    *,
    data_dir: Path,
    metric_threshold: list[tuple[str, int]],
    skip_existing_clusters: bool,
) -> None:
    [(metric, threshold)] = metric_threshold
    clusters.make_components_and_communities(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        # directed=True,
        skip_existing_clusters=skip_existing_clusters,
    )


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


def scatter_compute_ligand_leakage(
    *,
    data_dir: Path,
    test_leakage: bool,
) -> list[list[tuple[str, str, str]]]:
    paths = [path.as_posix() for path in (data_dir / "splits").glob("*parquet")]
    chunks = [
        [(path, compare_pair, metric)]
        for path in paths
        for compare_pair in [
            "train_test",
            "train_val",
            "val_test",
            "train_posebusters",
        ]
        for metric in [
            "tanimoto_similarity_max",
        ]
    ]
    if test_leakage:
        return chunks[:1]
    return chunks


def scatter_compute_protein_leakage(
    *,
    data_dir: Path,
    test_leakage: bool,
) -> list[list[tuple[str, str, str]]]:
    paths = [path.as_posix() for path in (data_dir / "splits").glob("*parquet")]
    chunks = [
        [(path, compare_pair, metric)]
        for path in paths
        for compare_pair in [
            "train_test",
            "train_val",
            "val_test",
            "train_posebusters",
        ]
        for metric in [
            "pli_qcov",
            "pocket_qcov",
            "pocket_lddt",
            "protein_seqsim_weighted_sum",
            "protein_fident_weighted_sum",
            "protein_lddt_weighted_sum",
            "protein_lddt_qcov_weighted_sum",
        ]
    ]
    if test_leakage:
        return chunks[:1]
    return chunks


def compute_ligand_leakage(
    *,
    data_dir: Path,
    inputs: list[tuple[str, str, str]],
) -> None:
    [(split_file, compare_pair, metric)] = inputs
    leakage.compute_ligand_leakage(
        data_dir=data_dir,
        split_file=split_file,
        compare_pair=compare_pair,
        metric=metric,
    )


def compute_protein_leakage(
    *,
    data_dir: Path,
    inputs: list[tuple[str, str, str]],
) -> None:
    [(split_file, compare_pair, metric)] = inputs
    leakage.compute_protein_leakage(
        data_dir=data_dir,
        split_file=split_file,
        compare_pair=compare_pair,
        metric=metric,
    )


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
                    ["kind", "reference_system_id"]
                )
            ],
        )
