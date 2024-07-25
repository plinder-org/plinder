# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import wraps
from hashlib import md5
from json import dumps
from os import listdir
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Any, Callable, Dict, Optional, TypeVar
from zipfile import ZipFile

import pandas as pd
from omegaconf import DictConfig

from plinder.core.utils.unpack import expand_config_context
from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger
from plinder.data.utils import tanimoto

if TYPE_CHECKING:
    from plinder.data.utils.annotations.aggregate_annotations import Entry
    from plinder.data.utils.annotations.get_similarity_scores import Scorer


LOG = setup_logger(__name__)
T = TypeVar("T")


def timeit(func: Callable[..., T]) -> Callable[..., T]:
    """
    Simple function timer decorator
    """

    @wraps(func)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        name = func.__name__
        mod = func.__module__
        log = setup_logger(".".join([mod, name]))
        ts = time()
        result = None
        try:
            result = func(*args, **kwargs)
            log.info(f"runtime succeeded: {time() - ts:>9.2f}s")
        except Exception as e:
            log.error(f"runtime failed: {time() - ts:>9.2f}s")
            log.error(f"{name} failed with: {repr(e)}")
            raise
        return result

    return wrapped


def entry_exists(*, entry_dir: Path, pdb_id: str) -> bool:
    """
    Check if the entry JSON file exists.

    Parameters
    ----------
    entry_dir : Path
        the directory containing entries
    pdb_id : str
        the PDB ID
    """
    two_char_code = pdb_id[-3:-1]
    output = entry_dir / two_char_code / (pdb_id + ".json")
    output.parent.mkdir(exist_ok=True, parents=True)
    return output.is_file()


@timeit
def load_entries(
    *,
    data_dir: Path,
    pdb_ids: list[str],
    clear_non_pocket_residues: bool = True,
    load_for_scoring: bool = True,
    max_protein_chains: int = 5,
    max_ligand_chains: int = 5,
) -> Dict[str, "Entry"]:
    """
    Load entries from the entries dir into a dict
    """
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    reduced = {}
    entry_dir = data_dir / "raw_entries"
    LOG.info(f"attempting to load {len(pdb_ids)} entries")
    for pdb_id in pdb_ids:
        try:
            reduced[pdb_id] = Entry.from_json(
                entry_dir / pdb_id[-3:-1] / (pdb_id + ".json"),
                clear_non_pocket_residues=clear_non_pocket_residues,
                load_for_scoring=load_for_scoring,
                max_protein_chains=max_protein_chains,
                max_ligand_chains=max_ligand_chains,
            )
        except Exception as e:
            LOG.error(f"pdb_id={pdb_id} failed with {repr(e)}")
    LOG.info(f"loaded {len(reduced)} entries from jsons")
    return reduced


@timeit
def load_entries_from_zips(
    *,
    data_dir: Path,
    two_char_codes: Optional[list[str]] = None,
    pdb_ids: Optional[list[str]] = None,
    load_for_scoring: bool = False,
) -> Dict[str, "Entry"]:
    """
    Load entries from the qc zips into a dict
    """
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    per_zip: dict[str, list[str]] | None = None
    entry_msg = "all"
    if pdb_ids is not None:
        zip_paths_set = set()
        per_zip = {}
        for pdb_id in pdb_ids:
            code = pdb_id[-3:-1]
            zip_paths_set.add(data_dir / "entries" / f"{code}.zip")
            per_zip.setdefault(code, [])
            per_zip[code].append(f"{pdb_id}.json")
        entry_msg = str(sum((len(pz) for pz in per_zip.values())))
        zip_paths = list(zip_paths_set)
    elif two_char_codes is not None:
        zip_paths = [data_dir / "entries" / f"{code}.zip" for code in two_char_codes]
    else:
        zip_paths = list((data_dir / "entries").glob("*"))
    reduced = {}
    LOG.info(f"attempting to load {entry_msg} entries from {len(zip_paths)} zips")
    for zip_path in zip_paths:
        if not zip_path.is_file():
            LOG.error(f"no archive {zip_path}, did you run structure_qc?")
            continue
        with ZipFile(zip_path) as archive:
            names = archive.namelist()
            if per_zip is not None:
                names = per_zip[zip_path.stem]
            for name in names:
                try:
                    with archive.open(name) as obj:
                        pdb_id = name.replace(".json", "")
                        reduced[pdb_id] = Entry.model_validate_json(obj.read()).prune(
                            load_for_scoring=load_for_scoring,
                        )
                except Exception as e:
                    LOG.error(f"failed to read name={name} failed with {repr(e)}")
    LOG.info(f"loaded {len(reduced)} entries from zips")
    return reduced


def get_db_sources(
    *, data_dir: Path, sub_databases: list[str] | None = None
) -> dict[str, Path]:
    if sub_databases is None:
        sub_databases = []
    dbs = {}
    if "holo" in sub_databases or not len(sub_databases):
        dbs["holo_foldseek"] = data_dir / "dbs" / "foldseek" / "foldseek"
        dbs["holo_mmseqs"] = data_dir / "dbs" / "mmseqs" / "mmseqs"
    if "apo" in sub_databases or not len(sub_databases):
        dbs["apo_foldseek"] = data_dir / "dbs" / "foldseek" / "foldseek"
        dbs["apo_mmseqs"] = data_dir / "dbs" / "mmseqs" / "mmseqs"
    if "pred" in sub_databases or not len(sub_databases):
        dbs["pred_foldseek"] = data_dir / "dbs" / "pred_foldseek" / "foldseek"
        dbs["pred_mmseqs"] = data_dir / "dbs" / "pred_mmseqs" / "mmseqs"
    return dbs


def get_scorer(
    *,
    data_dir: Path,
    pdb_ids: list[str],
    scorer_cfg: DictConfig,
    load_entries: bool,
) -> tuple["Scorer", list[str], Path]:
    from plinder.data.utils.annotations.get_similarity_scores import Scorer

    # need holo to db to compare against independently of what to score
    sub_dbs = list(set(scorer_cfg.sub_databases).union(["holo"]))
    db_sources = get_db_sources(
        data_dir=data_dir,
        sub_databases=sub_dbs,
    )
    hashed_contents = hash_contents(pdb_ids)
    if load_entries:
        entries = load_entries_from_zips(
            pdb_ids=None,  # explicitly set to None to load all entries
            data_dir=data_dir,
            load_for_scoring=True,
        )
        # TODO: bug where scatter_make_scorers uses raw ingest files
        #       instead of available entries from zips but Scorer
        #       assumes all entries are present
        entry_ids = list(set(entries.keys()).intersection(pdb_ids))
    else:
        entries = {}
        entry_ids = pdb_ids
    scores_dir = data_dir / "scores"
    sub_db_dir = data_dir / "dbs" / "subdbs"
    batch_db_dir = data_dir / "dbs" / "subdbs" / "batch_dbs" / hashed_contents
    batch_db_dir.mkdir(exist_ok=True, parents=True)
    return Scorer(
        entries=entries,
        source_to_full_db_file=db_sources,
        db_dir=sub_db_dir,
        scores_dir=scores_dir,
        minimum_threshold=scorer_cfg.minimum_threshold,
    ), entry_ids, batch_db_dir


def save_ligand_batch(
    *,
    entries: dict[str, "Entry"],
    output_path: Path,
) -> None:
    dfs = []
    for entry in entries.values():
        df = tanimoto.load_ligands_from_entry(entry=entry)
        if df is not None:
            dfs.append(df)
    LOG.info(f"save_ligand_batch: {len(dfs)} entries have usable ligands")
    df = pd.concat(dfs).drop_duplicates().reset_index(drop=True)
    for col in df.columns:
        nunique = df[col].nunique()
        LOG.info(f"save_ligand_batch: unique {col}={nunique}")
    LOG.info(f"save_ligands_batch: writing {output_path}")
    df.to_parquet(output_path, index=False)


def hash_contents(contents: list[str]) -> str:
    """
    Return a repeatable unique string identifier of a list of strings

    Parameters
    ----------
    contents : list[str]
        list of strings to identify uniquely

    Returns
    -------
    hash : str
        unique string corresponding to contents
    """
    return md5(dumps(sorted(contents)).encode("utf8")).hexdigest()


def get_local_contents(
    *,
    data_dir: Path,
    two_char_codes: Optional[list[str]] = None,
    pdb_ids: Optional[list[str]] = None,
    as_four_char_ids: bool = False,
) -> list[str]:
    """
    Starting from a root directory, assume subdirectories
    of two character codes each containing subdirectories
    of individual files. The as_ids kludge is intended to
    support both fully qualified PDB (pdb_0000{pdb_id})
    and short PDB ({pdb_id})

    Parameters
    ----------
    data_dir : Path
        directory containing two character code directories
    two_char_codes : list[str], default=None
        subset of two character codes
    pdb_ids : list[str], default=None
        subset of pdb IDs (overrides two_char_codes)
    as_four_char_ids : bool, default=False
        if True, return 4 character codes instead of nested
        subdirectories

    Returns
    -------
    contents : list[str]
        list of directory-derived metadata contents
    """
    kind, values = expand_config_context(
        pdb_ids=pdb_ids,
        two_char_codes=two_char_codes,
    )
    if kind == "pdb_id":
        return values if as_four_char_ids else [f"pdb_0000{pdb_id}" for pdb_id in values]
    codes = values if kind == "two_char_code" and len(values) else listdir(data_dir.as_posix())
    contents = []
    for code in codes:
        contents.extend(listdir((data_dir / code).as_posix()))
    if as_four_char_ids:
        return sorted([c[-4:] for c in contents])
    return sorted(contents)


def partition_batch_scores(*, partition_dir: Path, scores_dir: Path) -> None:
    """
    Consolidate individual pdb ID similarity scores parquet
    files into a pre-partitioned dataset by similarity metric
    and metric value. This partitioned dataset needs to be further
    consolidated in a join step elsewhere.

    Parameters
    ----------
    partition_dir : Path
        destination directory for consolidated scores
    scores_dir : Path
        source directory for fragmented scores
    """
    # collect the fragmented parquets
    dfs = []
    for pqt in scores_dir.glob("*.parquet"):
        df = pd.read_parquet(pqt)
        if not df.empty:
            dfs.append(df)
    if len(dfs):
        df = pd.concat(dfs).reset_index(drop=True)
        df.to_parquet(
            partition_dir,
            partition_cols=["metric", "similarity"],
            index=False,
            max_partitions=3939,
            schema=schemas.PROTEIN_SIMILARITY_SCHEMA,
        )


def should_run_stage(stage: str, run: list[str], skip: list[str]) -> bool:
    """
    Compare function name to list of whitelisted / blacklisted
    stages to determine short-circuiting behavior for pipeline
    decorator.

    Parameters
    ----------
    stage : str
        the stage in question
    run : list[str]
        list of stages to run
    skip : list[str]
        list of stages to skip

    Returns
    -------
    run : bool
        whether or not to run the function
    """
    if len(run):
        if stage in run and stage not in skip:
            return True
        return False
    elif len(skip):
        if stage in skip:
            return False
        return True
    return True


def ingest_flow_control(func: Callable[..., T]) -> Callable[..., T]:
    """
    Function decorator to apply for every stage
    in the IngestPipeline.
    """

    @wraps(func)
    def inner(pipe: Any, *args: Optional[list[str]], **kwargs: Any) -> Any:
        is_scatter = False
        is_join = False
        name = func.__name__
        if func.__name__.startswith("scatter_"):
            is_scatter = True
            name = func.__name__.replace("scatter_", "", 1)
        elif func.__name__.startswith("join_"):
            is_join = True
            name = func.__name__.replace("join_", "", 1)
        if should_run_stage(
            name,
            pipe.cfg.flow.run_specific_stages,
            pipe.cfg.flow.skip_specific_stages,
        ):
            chunks = None
            if len(args) and args[0] is not None:
                chunks = len(args[0])
            verb = "computing"
            if is_join:
                verb = "joining"
            elif is_scatter:
                verb = "producing"
            msg = f"{func.__name__} {verb}"
            if chunks is not None:
                msg += f" {chunks} parts"
            if not is_scatter:
                LOG.info(msg)
            ret = func(pipe, *args, **kwargs)
            if is_scatter and ret is not None:
                LOG.info(f"{msg} {len(ret)} chunks")  # type: ignore
            return ret
        else:
            LOG.info(f"skipping {func.__name__}")
        return [[]]

    return inner


def create_nonredundant_dataset(*, data_dir: Path) -> None:
    """
    This is called in make_mmp_index to ensure the existence of the index
    and simultaneously generates a non-redundant index for various use
    cases. Ultimately this should run as the join step of structure_qc.
    """
    if not (data_dir / "index" / "annotation_table.parquet").exists():
        dfs = []
        for path in (data_dir / "qc" / "index").glob("*"):
            dfs.append(pd.read_parquet(path))
        df = pd.concat(dfs).reset_index(drop=True)
        (data_dir / "index").mkdir(exist_ok=True, parents=True)
        df.to_parquet(data_dir / "index" / "annotation_table.parquet", index=False)
    else:
        df = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")
    df_trainable = df[
        ~(df.ligand_is_invalid | df.ligand_is_artifact | df.ligand_is_ion)
        & (df.ligand_passes_valence_checks)
        & (df.ligand_sanitization)
        & (df.ligand_passes_kekulization)
        & (df.ligand_mol_pred_loaded)
        & (df.system_num_ligand_chains == 1)
        & (df.ligand_rdkit_canonical_smiles.notna())
    ].reset_index(drop=True)

    df_trainable["uniqueness"] = (
        df_trainable["system_id_no_biounit"]
        + "_"
        + df_trainable["pli_qcov__100__strong__component"]
    )
    df_trainable_nonredundant = df_trainable.sort_values(
        "system_biounit_id"
    ).drop_duplicates("uniqueness")
    df_trainable_nonredundant.to_parquet(
        data_dir / "index" / "annotation_table_nonredundant.parquet", index=False
    )
