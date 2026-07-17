# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import multiprocessing
import shutil
from functools import wraps
from hashlib import md5
from itertools import repeat
from json import dumps, load
from os import listdir
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Any, Callable, Literal, Optional, TypeVar
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
import pyarrow.parquet as pq
from omegaconf import DictConfig

from plinder.core.structure import smallmols_similarity
from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import expand_config_context

if TYPE_CHECKING:
    from plinder.data.utils.annotations.get_similarity_scores import Scorer


LOG = setup_logger(__name__)
T = TypeVar("T")
RETIRED_ENRICHMENT_MARKERS = ("ecod", "panther", "kinase")


def _drop_retired_enrichment_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Remove annotations sourced from retired third-party enrichments."""
    retired = [
        column
        for column in df.columns
        if any(marker in column.casefold() for marker in RETIRED_ENRICHMENT_MARKERS)
    ]
    if retired:
        LOG.info(f"dropping retired enrichment columns: {sorted(retired)}")
        return df.drop(columns=retired)
    return df


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
    Check if the per-entry annotation parquet exists.

    Parameters
    ----------
    entry_dir : Path
        the directory containing entries
    pdb_id : str
        the PDB ID
    """
    two_char_code = pdb_id[-3:-1]
    output = entry_dir / two_char_code / (pdb_id + ".parquet")
    output.parent.mkdir(exist_ok=True, parents=True)
    entry_chains = entry_dir / two_char_code / pdb_id / "entry_chains.parquet"
    entry_source = entry_dir / two_char_code / pdb_id / "entry_source.parquet"
    if not (output.is_file() and entry_chains.is_file() and entry_source.is_file()):
        return False
    try:
        if "system_receptor_type" not in pq.read_schema(output).names:
            return False
        if "chain_receptor_type" not in pq.read_schema(entry_chains).names:
            return False
    except Exception:
        LOG.info(f"invalidating stale entry annotation cache for {pdb_id}")
        return False
    return True


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
        from plinder.core.scores.entries import entry_views_from_df

        annotation = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")
        entry_chains = pd.read_parquet(data_dir / "index" / "entry_chains.parquet")
        entries = entry_views_from_df(annotation, entry_chains=entry_chains)
        entry_ids = list(set(entries.keys()).intersection(pdb_ids))
    else:
        entries = {}
        entry_ids = pdb_ids
    scores_dir = data_dir / "scores"
    sub_db_dir = data_dir / "dbs" / "subdbs"
    batch_db_dir = data_dir / "dbs" / "subdbs" / "batch_dbs" / hashed_contents
    batch_db_dir.mkdir(exist_ok=True, parents=True)
    return (
        Scorer(
            entries=entries,
            source_to_full_db_file=db_sources,
            db_dir=sub_db_dir,
            scores_dir=scores_dir,
            minimum_threshold=scorer_cfg.minimum_threshold,
        ),
        entry_ids,
        batch_db_dir,
    )


def save_ligand_batch(
    *,
    annotation: pd.DataFrame,
    output_path: Path,
) -> None:
    df = smallmols_similarity.load_ligands_from_index(annotation=annotation)
    LOG.info(f"save_ligand_batch: {df['pdb_id'].nunique()} entries have usable ligands")
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
    if kind == "pdb_ids":
        return (
            values if as_four_char_ids else [f"pdb_0000{pdb_id}" for pdb_id in values]
        )
    codes = (
        values
        if kind == "two_char_codes" and len(values)
        else listdir(data_dir.as_posix())
    )
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


def get_pdb_ids_in_scoring_dataset(*, data_dir: Path) -> dict[str, list[str]]:
    """
    Get all the pdb IDs that are present in the raw scoring dataset

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    """
    found = {}
    dbs = data_dir / "dbs" / "subdbs"
    for search_db in ["holo", "apo", "pred"]:
        found[search_db] = [
            path.stem for path in (dbs / f"search_db={search_db}").glob("*parquet")
        ]
    return found


def get_alns(
    *, data_dir: Path, mapped: bool = False
) -> dict[str, dict[str, list[str]]]:
    """
    Get all the pdb IDs that are present in the raw alignment dataset

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    """
    sub = "mapped_aln" if mapped else "aln"
    found: dict[str, dict[str, list[str]]] = {}
    dbs = data_dir / "dbs" / "subdbs"
    for search_db in ["holo", "apo", "pred"]:
        found.setdefault(search_db, {})
        for aln_type in ["foldseek", "mmseqs"]:
            found[search_db].setdefault(aln_type, [])
            found[search_db][aln_type] = [
                path.stem
                for path in (dbs / f"{search_db}_{aln_type}/{sub}/").glob("*parquet")
            ]
    return found


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


def _cluster_column_name(
    *, metric: str, cluster: str, directed: bool, threshold: int, ligand: bool
) -> str:
    if cluster == "components":
        kind = "strong__component" if directed else "weak__component"
    else:
        kind = "community"
    ligand_marker = "__ligand" if ligand else ""
    return f"{metric}__{threshold}{ligand_marker}__{kind}"


def _read_local_cluster_rows(*, root: Path, node_column: str) -> pd.DataFrame:
    columns = [
        node_column,
        "label",
        "threshold",
        "metric",
        "cluster",
        "directed",
    ]
    frames = [
        pd.read_parquet(path, columns=columns)
        for path in sorted(root.rglob("*.parquet"))
    ]
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True)


def _pivot_cluster_rows(
    clusters: pd.DataFrame, *, node_column: str, ligand: bool
) -> pd.DataFrame:
    if clusters.empty:
        return pd.DataFrame(columns=[node_column])
    clusters = clusters.copy()
    clusters["directed"] = clusters["directed"].map(
        lambda value: value if isinstance(value, bool) else str(value).lower() == "true"
    )
    wide = clusters.pivot_table(
        values="label",
        index=node_column,
        columns=["metric", "cluster", "directed", "threshold"],
        aggfunc="first",
    )
    wide.columns = [
        _cluster_column_name(
            metric=str(metric),
            cluster=str(cluster),
            directed=bool(directed),
            threshold=int(threshold),
            ligand=ligand,
        )
        for metric, cluster, directed, threshold in wide.columns
    ]
    return wide.reset_index()


def add_cluster_columns(*, index: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    """Merge locally generated system and ligand clusters into the index."""
    result = index
    for root_name, node_column, ligand in [
        ("clusters", "system_id", False),
        ("ligand_clusters", "ligand_id", True),
    ]:
        clusters = _read_local_cluster_rows(
            root=data_dir / root_name,
            node_column=node_column,
        )
        if clusters.empty:
            continue
        clusters = clusters[clusters[node_column].isin(set(index[node_column]))]
        wide = _pivot_cluster_rows(
            clusters,
            node_column=node_column,
            ligand=ligand,
        )
        replacement_columns = set(wide.columns).difference({node_column})
        result = result.drop(
            columns=list(replacement_columns.intersection(result.columns))
        ).merge(wide, on=node_column, how="left", validate="many_to_one")
    return result


def add_aggregated_columns(*, index: pd.DataFrame) -> pd.DataFrame:
    """
    Add aggregated columns to the annotation table
    """
    index["biounit_num_ligands"] = index.groupby(["entry_pdb_id", "system_biounit_id"])[
        "system_id"
    ].transform("count")
    index["biounit_num_unique_ccd_codes"] = index.groupby(
        [
            "entry_pdb_id",
            "system_biounit_id",
        ]
    )["ligand_unique_ccd_code"].transform("nunique")
    index["biounit_num_proper_ligands"] = index.groupby(
        [
            "entry_pdb_id",
            "system_biounit_id",
        ]
    )["ligand_is_proper"].transform("sum")
    for n in [
        "lipinski",
        "cofactor",
        "fragment",
        "oligo",
        "artifact",
        "other",
        "covalent",
        "invalid",
        "ion",
    ]:
        index[f"system_ligand_has_{n}"] = index.groupby("system_id")[
            f"ligand_is_{n}"
        ].transform("any")
    index["system_protein_chains_total_length"] = index[
        "system_protein_chains_length"
    ].apply(sum)
    ccd_dict = (
        index.groupby("system_id")["ligand_unique_ccd_code"]
        .agg(lambda x: "-".join(sorted(set(x))))
        .to_dict()
    )
    index["system_unique_ccd_codes"] = index["system_id"].map(ccd_dict)
    ccd_proper_dict = (
        index[index["ligand_is_proper"]]
        .groupby("system_id")["ligand_unique_ccd_code"]
        .agg(lambda x: "-".join(sorted(set(x))))
        .to_dict()
    )
    index["system_proper_unique_ccd_codes"] = index["system_id"].map(ccd_proper_dict)
    return index


def finalize_index(*, data_dir: Path) -> pd.DataFrame:
    """Merge local cluster labels into the V3 annotation parquet."""
    index_path = data_dir / "index" / "annotation_table.parquet"
    if not index_path.is_file():
        raise FileNotFoundError(index_path)
    index = pd.read_parquet(index_path)
    index = add_cluster_columns(index=index, data_dir=data_dir)
    uniqueness_cluster = "pli_qcov__100__strong__component"
    if uniqueness_cluster in index.columns:
        labels = (
            index[uniqueness_cluster]
            .astype("string")
            .fillna(index["system_id"].astype("string"))
        )
        index["uniqueness"] = (
            index["system_id_no_biounit"].astype("string") + "_" + labels
        )
    index.to_parquet(index_path, index=False)
    return index


def create_entry_chain_index(
    *, data_dir: Path, force_update: bool = False
) -> pd.DataFrame:
    """Collate normalized per-entry chain metadata into one parquet."""
    output = data_dir / "index" / "entry_chains.parquet"
    output.parent.mkdir(exist_ok=True, parents=True)
    if output.exists() and not force_update:
        return pd.read_parquet(output)

    parts = sorted((data_dir / "raw_entries").glob("*/*/entry_chains.parquet"))
    frames = [pd.read_parquet(path) for path in parts]
    columns = [
        "entry_pdb_id",
        "chain_asym_id",
        "chain_auth_id",
        "chain_entity_id",
        "chain_type",
        "chain_receptor_type",
        "chain_length",
        "chain_num_unresolved_residues",
        "chain_is_holo",
        "chain_uniprot_ids",
    ]
    chains = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    chains = chains.reindex(columns=columns)
    chains.to_parquet(output, index=False)
    return chains


def create_entry_source_index(
    *, data_dir: Path, force_update: bool = False
) -> pd.DataFrame:
    """Collate one source-mmCIF revision row per PDB entry."""
    output = data_dir / "index" / "entry_sources.parquet"
    output.parent.mkdir(exist_ok=True, parents=True)
    if output.exists() and not force_update:
        return pd.read_parquet(output)

    parts = sorted((data_dir / "raw_entries").glob("*/*/entry_source.parquet"))
    columns = [
        "entry_pdb_id",
        "source_mmcif_major_revision",
        "source_mmcif_minor_revision",
    ]
    sources = (
        pd.concat([pd.read_parquet(path) for path in parts], ignore_index=True)
        if parts
        else pd.DataFrame(columns=columns)
    )
    if not sources.empty and sources["entry_pdb_id"].duplicated().any():
        duplicates = sorted(
            sources.loc[sources["entry_pdb_id"].duplicated(), "entry_pdb_id"].unique()
        )
        raise ValueError(f"duplicate entry source metadata: {duplicates}")
    sources.to_parquet(output, index=False)
    return sources


def create_index(*, data_dir: Path, force_update: bool = False) -> pd.DataFrame:
    """
    Create the index
    """
    index = data_dir / "index" / "annotation_table.parquet"
    index.parent.mkdir(exist_ok=True, parents=True)
    create_entry_chain_index(data_dir=data_dir, force_update=force_update)
    create_entry_source_index(data_dir=data_dir, force_update=force_update)

    if not index.exists() or force_update:
        dfs = []
        annotation_parts = data_dir / "raw_entries"
        for i, path in enumerate(annotation_parts.glob("*/*.parquet")):
            df = _drop_retired_enrichment_columns(pd.read_parquet(path))
            LOG.info(f"{i} {path.name} shape={df.shape}")
            if not df.empty:
                dfs.append(df)
        if not dfs:
            LOG.warning(
                f"create_index: no parquet files in {annotation_parts}, "
                "writing empty index"
            )
            pd.DataFrame().to_parquet(index, index=False)
            return pd.read_parquet(index)
        df = pd.concat(dfs).reset_index(drop=True)
        df.to_parquet(index, index=False)
    else:
        df = pd.read_parquet(index)
    old_columns = set(df.columns)
    df = _drop_retired_enrichment_columns(df)
    df = add_aggregated_columns(index=df)
    update = old_columns != set(df.columns)
    if update or force_update:
        df.to_parquet(index, index=False)
    return df


def create_nonredundant_dataset(*, data_dir: Path) -> None:
    """
    This is called in make_mmp_index to ensure the existence of the index
    and simultaneously generates a non-redundant index for various use
    cases. The initial index is collated in ``join_make_entries``.
    """
    if not (data_dir / "index" / "annotation_table.parquet").exists():
        df = create_index(data_dir=data_dir)
    else:
        df = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")
    if "uniqueness" not in df.columns:
        raise RuntimeError(
            "annotation index has no uniqueness column; run clustering and "
            "finalize_index() before creating the nonredundant dataset"
        )
    df_nonredundant = df.sort_values("system_biounit_id").drop_duplicates("uniqueness")
    df_nonredundant.to_parquet(
        data_dir / "index" / "annotation_table_nonredundant.parquet", index=False
    )


def apo_file_from_link_id(
    data_dir: Path,
    output_dir: Path,
    link_id: str,
    force_update: bool = False,
) -> dict[str, str] | None:
    import biotite.structure.io.pdbx as pdbx

    from plinder.data.utils.annotations.cif_utils import read_mmcif_file
    from plinder.data.utils.annotations.save_utils import save_cif_file

    if (output_dir / f"{link_id}.cif").exists() and not force_update:
        LOG.info(f"skipping {link_id}.cif as it already exists")
        return None

    pdb_id, chain = link_id.split("_")
    target_cif = (
        data_dir
        / "ingest"
        / pdb_id[1:3]
        / f"pdb_0000{pdb_id}"
        / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
    )
    if not target_cif.exists():
        LOG.info(f"skipping {link_id} as {target_cif} does not exist")
        return None

    cif_file_obj = read_mmcif_file(target_cif)
    atoms = pdbx.get_structure(
        cif_file_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[atoms.chain_id == chain]
    out_cif = output_dir / f"{pdb_id}_{chain}.cif"
    LOG.info(f"saving {link_id} to {out_cif}")
    save_cif_file(atoms, out_cif.stem, out_cif)
    return None


def pred_file_from_link_id(
    data_dir: Path,
    output_dir: Path,
    link_id: str,
    force_update: bool = False,
) -> None:
    import biotite.structure.io.pdbx as pdbx

    from plinder.data.utils.annotations.cif_utils import read_mmcif_file
    from plinder.data.utils.annotations.save_utils import save_cif_file

    if (output_dir / f"{link_id}.cif").exists() and not force_update:
        LOG.info(f"skipping {link_id}.cif as it already exists")
        return None

    uniprot_id, chain = link_id.split("_")
    target_cif = data_dir / "dbs" / "alphafold" / f"AF-{uniprot_id}-F1-model_v4.cif"
    if not target_cif.exists():
        LOG.info(f"skipping {link_id} as {target_cif} does not exist")
        return None

    cif_file_obj = read_mmcif_file(target_cif)
    atoms = pdbx.get_structure(
        cif_file_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[atoms.chain_id == chain]
    out_cif = output_dir / f"{uniprot_id}_{chain}.cif"
    LOG.info(f"saving {link_id} to {out_cif}")
    save_cif_file(atoms, out_cif.stem, out_cif)
    return None
    # chain_to_seqres = {c.name: c.string for c in seqres}
    # return chain_to_seqres[chain]


def pack_linked_structures(data_dir: Path, code: str, structures: bool = True) -> None:
    """
    Pack generated linked structures into a zip file for a particular
    two character code.

    Parameters
    ----------
    data_dir : Path
        plinder root dir
    code : str
        two character code
    structures : bool, default=True
        if True, make structure archives
    """
    (data_dir / "links").mkdir(exist_ok=True, parents=True)
    mode: Literal["r", "w"] = "w" if structures else "r"
    with ZipFile(
        data_dir / "links" / f"{code}.zip", mode, compression=ZIP_DEFLATED
    ) as archive:
        for search_db in ["apo", "pred"]:
            jsons = []
            root = data_dir / "linked_staging" / search_db
            system_ids = [
                system_id for system_id in listdir(root) if system_id[1:3] == code
            ]
            for system_id in system_ids:
                link_ids = listdir(f"{root}/{system_id}")
                for link_id in link_ids:
                    link = f"{root}/{system_id}/{link_id}"
                    try:
                        with open(f"{link}/scores.json") as f:
                            jsons.append(load(f))
                    except Exception:
                        pass
                    if structures:
                        try:
                            archive.write(
                                f"{link}/superposed.cif",
                                f"{search_db}/{system_id}/{link_id}/superposed.cif",
                            )
                        except Exception:
                            pass
            df = pd.DataFrame(jsons).rename(
                columns={"reference": "reference_system_id", "model": "id"}
            )
            df.to_parquet(
                data_dir / "links" / f"{search_db}_{code}.parquet", index=False
            )


def mp_pack_linked_structures(*, data_dir: Path, structures: bool = True) -> None:
    """
    Use a process pool to pack linked structures into two character code archives.

    Parameters
    ----------
    data_dir : Path
        plinder root dir
    """

    with multiprocessing.get_context("spawn").Pool() as pool:
        pool.starmap(
            pack_linked_structures,
            zip(repeat(data_dir), listdir(data_dir / "ingest"), repeat(structures)),
        )


def pack_source_structures(data_dir: Path, search_db: str) -> None:
    (data_dir / "linked_structures").mkdir(exist_ok=True, parents=True)
    with ZipFile(
        data_dir / "linked_structures" / f"{search_db}.zip",
        "w",
        compression=ZIP_DEFLATED,
    ) as archive:
        source_structures = data_dir / "linked_staging" / "source" / search_db
        for path in source_structures.rglob("*.cif"):
            archive.write(path, path.name)


def consolidate_linked_scores(*, data_dir: Path) -> None:
    """
    Consolidate linked scores into a single parquet file. Assumes
    that pack_linked_structures has been run.

    Parameters
    ----------
    data_dir : Path
        plinder root dir
    """
    for search_db in ["apo", "pred"]:
        paths = list((data_dir / "links").glob(f"{search_db}_*.parquet"))
        dfs = []
        for path in paths:
            df = pd.read_parquet(path)
            if not df.empty:
                dfs.append(df)
        ndf = pd.concat(dfs)
        odf = pd.read_parquet(
            data_dir / "linked_staging" / f"{search_db}_links.parquet"
        )
        drop = list(
            set(odf.columns.intersection(ndf.columns))
            - set(["reference_system_id", "id"])
        )
        df = pd.merge(odf.drop(columns=drop), ndf, on=["reference_system_id", "id"])
        (data_dir / "links" / f"kind={search_db}").mkdir(exist_ok=True, parents=True)
        df.to_parquet(
            data_dir / "links" / f"kind={search_db}" / "links.parquet", index=False
        )


def rename_clusters(*, data_dir: Path) -> None:
    """
    Rename cluster files to match the hive layout convention.

    Parameters
    ----------
    data_dir : Path
        plinder root dir
    """
    cluster_dir = data_dir / "clusters"
    cluster_paths = [path for path in cluster_dir.rglob("*") if path.is_file()]
    for path in cluster_paths:
        if path.name == "data.parquet":
            continue
        base = path.parent
        name = path.stem
        apath = base / name / "data.parquet"
        apath.parent.mkdir(exist_ok=True, parents=True)
        shutil.move(path, apath)
