# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import shutil
from functools import wraps
from hashlib import md5
from json import dumps, load
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Any, Callable, Optional, TypeVar

import pandas as pd
import pyarrow.parquet as pq
from omegaconf import DictConfig, OmegaConf

from plinder.core.scores.metrics import is_chemical_cluster_metric
from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger

# TODO(get_local_contents): only used by the quarantined get_local_contents below;
# restore this import if that function is confirmed live.
# from plinder.core.utils.unpack import expand_config_context

if TYPE_CHECKING:
    from plinder.data.annotations.get_similarity_scores import Scorer


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
    entry_biounit_chains = (
        entry_dir / two_char_code / pdb_id / "entry_biounit_chains.parquet"
    )
    entry_source = entry_dir / two_char_code / pdb_id / "entry_source.parquet"
    if not (
        output.is_file()
        and entry_chains.is_file()
        and entry_biounit_chains.is_file()
        and entry_source.is_file()
    ):
        return False
    try:
        if "system_receptor_type" not in pq.read_schema(output).names:
            return False
        if "chain_receptor_type" not in pq.read_schema(entry_chains).names:
            return False
        required_biounit_columns = {
            "entry_pdb_id",
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
            "chain_num_contacting_ions",
            "chain_num_contacting_artifacts",
            "chain_num_contacting_other_ligands",
        }
        if not required_biounit_columns.issubset(
            pq.read_schema(entry_biounit_chains).names
        ):
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
    foldseek_cfg: DictConfig | None = None,
    mmseqs_cfg: DictConfig | None = None,
    scratch_dir: Path | None = None,
) -> tuple["Scorer", list[str], Path]:
    from plinder.data.annotations.get_similarity_scores import Scorer
    from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig

    foldseek_config = (
        OmegaConf.to_object(foldseek_cfg)
        if foldseek_cfg is not None
        else FoldseekConfig()
    )
    mmseqs_config = (
        OmegaConf.to_object(mmseqs_cfg) if mmseqs_cfg is not None else MMSeqsConfig()
    )
    if not isinstance(foldseek_config, FoldseekConfig):
        raise TypeError("foldseek configuration did not resolve to FoldseekConfig")
    if not isinstance(mmseqs_config, MMSeqsConfig):
        raise TypeError("mmseqs configuration did not resolve to MMSeqsConfig")

    # need holo to db to compare against independently of what to score
    sub_dbs = list(set(scorer_cfg.sub_databases).union(["holo"]))
    db_sources = get_db_sources(
        data_dir=data_dir,
        sub_databases=sub_dbs,
    )
    hashed_contents = hash_contents(pdb_ids)
    if load_entries:
        from plinder.core.scores.entries import load_entry_views

        entries = load_entry_views(pdb_ids=pdb_ids, data_dir=data_dir)
        entry_ids = sorted(entries)
    else:
        entries = {}
        entry_ids = pdb_ids
    scores_dir = data_dir / "scores"
    sub_db_dir = data_dir / "dbs" / "subdbs"
    batch_db_root = scratch_dir or data_dir / "dbs" / "subdbs" / "batch_dbs"
    batch_db_dir = batch_db_root / hashed_contents
    batch_db_dir.mkdir(exist_ok=True, parents=True)
    return (
        Scorer(
            entries=entries,
            source_to_full_db_file=db_sources,
            db_dir=sub_db_dir,
            scores_dir=scores_dir,
            minimum_threshold=scorer_cfg.minimum_threshold,
            minimum_thresholds=dict(scorer_cfg.minimum_thresholds),
            max_query_protein_chains=scorer_cfg.max_query_protein_chains,
            max_query_proper_ligand_chains=(scorer_cfg.max_query_proper_ligand_chains),
            foldseek_config=foldseek_config,
            mmseqs_config=mmseqs_config,
        ),
        entry_ids,
        batch_db_dir,
    )


def save_ligand_batch(
    *,
    data_dir: Path,
    annotation: pd.DataFrame,
    output_path: Path,
) -> None:
    from plinder.data.annotations.get_similarity_scores import (
        annotate_ligand_3d_score_ability,
        load_ligands_from_index,
    )

    df = load_ligands_from_index(annotation=annotation)
    df = annotate_ligand_3d_score_ability(df, data_dir=data_dir)
    LOG.info(
        f"save_ligand_batch: collected ligands from {df['pdb_id'].nunique()} entries"
    )
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


# TODO(get_local_contents): appears unused in-repo — zero callers on HEAD, and the
# last consumer (make_sub_dbs) documents deliberately NOT routing through it. Its
# only exercise is its own tests, which are order-flaky because the None-context
# path reads the process-global cached config (leaks across the suite). Commented
# out (with its tests) to unblock the pipeline suite. Confirm whether an external
# metaflow flow references it: if stale -> delete; if used -> restore + add an
# autouse config-cache reset so the None-context path is deterministic.
# def get_local_contents(
#     *,
#     data_dir: Path,
#     two_char_codes: Optional[list[str]] = None,
#     pdb_ids: Optional[list[str]] = None,
#     as_four_char_ids: bool = False,
# ) -> list[str]:
#     """
#     Starting from a root directory, assume subdirectories
#     of two character codes each containing subdirectories
#     of individual files. The as_ids kludge is intended to
#     support both fully qualified PDB (pdb_0000{pdb_id})
#     and short PDB ({pdb_id})
#
#     Parameters
#     ----------
#     data_dir : Path
#         directory containing two character code directories
#     two_char_codes : list[str], default=None
#         subset of two character codes
#     pdb_ids : list[str], default=None
#         subset of pdb IDs (overrides two_char_codes)
#     as_four_char_ids : bool, default=False
#         if True, return 4 character codes instead of nested
#         subdirectories
#
#     Returns
#     -------
#     contents : list[str]
#         list of directory-derived metadata contents
#     """
#     kind, values = expand_config_context(
#         pdb_ids=pdb_ids,
#         two_char_codes=two_char_codes,
#     )
#     if kind == "pdb_ids":
#         return (
#             values if as_four_char_ids else [f"pdb_0000{pdb_id}" for pdb_id in values]
#         )
#     codes = (
#         values
#         if kind == "two_char_codes" and len(values)
#         else listdir(data_dir.as_posix())
#     )
#     contents = []
#     for code in codes:
#         contents.extend(listdir((data_dir / code).as_posix()))
#     if as_four_char_ids:
#         return sorted([c[-4:] for c in contents])
#     return sorted(contents)


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
        if pqt.name.endswith(".tmp.parquet"):
            continue
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
            path.stem
            for path in (dbs / f"search_db={search_db}").glob("*parquet")
            if not path.name.endswith(".tmp.parquet")
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
                if not path.name.endswith(".tmp.parquet")
                and (
                    not mapped
                    or _mapped_alignment_file_is_current(path, alignment_type=aln_type)
                )
            ]
    return found


def _mapped_alignment_file_is_current(path: Path, *, alignment_type: str) -> bool:
    """Treat corrupt and pre-compact mapped files as incomplete cache entries."""
    try:
        columns = set(pq.read_schema(path).names)
    except (OSError, ValueError) as exc:
        LOG.warning(f"ignoring unreadable mapped alignment {path}: {exc}")
        return False
    return schemas.mapped_alignment_schema_is_current(
        columns, alignment_type=alignment_type
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
            # Metaflow foreach joins need one no-op branch in order to remain
            # reachable when a stage is excluded by run_specific_stages.
            if is_scatter:
                return [[]]
        return [[]]

    return inner


def _cluster_column_name(
    *, metric: str, cluster: str, directed: bool, threshold: int, ligand: bool
) -> str:
    if directed or cluster != "set_cover":
        raise ValueError("published undirected clusters must be set covers")
    kind = "set_cover"
    ligand_marker = "__ligand" if ligand else ""
    return f"{metric}__{threshold}{ligand_marker}__{kind}"


def build_ligand_cluster_table(*, index: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    """Build one queryable cluster-assignment row per ligand."""
    node_column = "ligand_id"
    directed_cover_root = data_dir / "ligand_sampling" / "directed_set_cover"
    marker_path = data_dir / "index" / "collation.json"
    repair_started_ns: int | None = None
    if marker_path.is_file():
        with marker_path.open() as handle:
            marker = load(handle)
        if marker.get("status") == "requires_downstream_repair":
            repair_started_ns = marker_path.stat().st_mtime_ns
    set_cover_root = data_dir / "ligand_sampling" / "set_cover"
    reciprocal_paths = sorted(set_cover_root.glob("metric=*/threshold=*.parquet"))
    directed_cover_paths = sorted(
        directed_cover_root.glob("metric=*/threshold=*.parquet")
    )
    if (
        not reciprocal_paths
        and not directed_cover_paths
        and repair_started_ns is not None
    ):
        raise FileNotFoundError(
            "targeted collation repair has no rebuilt ligand clusters"
        )
    set_cover_keys = {
        (
            next(
                part.split("=", maxsplit=1)[1]
                for part in path.relative_to(set_cover_root).parts
                if part.startswith("metric=")
            ),
            int(path.stem.split("=", maxsplit=1)[1]),
        )
        for path in reciprocal_paths
    }
    directed_cover_keys = {
        (
            path.parent.name.split("=", maxsplit=1)[1],
            int(path.stem.split("=", maxsplit=1)[1]),
        )
        for path in directed_cover_paths
    }
    invalid_set_cover_metrics = sorted(
        metric
        for metric, _ in set_cover_keys
        if metric != "tanimoto_similarity_ecfp4_1024"
    )
    invalid_directed_metrics = sorted(
        metric
        for metric, _ in directed_cover_keys
        if metric == "tanimoto_similarity_ecfp4_1024"
    )
    if invalid_set_cover_metrics or invalid_directed_metrics:
        raise ValueError(
            "invalid ligand set-cover modes: "
            f"undirected={invalid_set_cover_metrics}, "
            f"directed={invalid_directed_metrics}"
        )
    node_ids = pd.Index(
        index[node_column].dropna().astype(str).unique(),
        name=node_column,
    )
    if not reciprocal_paths and not directed_cover_paths:
        return pd.DataFrame({node_column: node_ids.to_numpy()})
    proper = index["ligand_is_proper"].fillna(False).astype(bool)
    holo = index["system_type"].eq("holo")
    expected_score_nodes = set(
        index.loc[proper & holo, node_column].dropna().astype(str)
    )
    expected_fingerprint_nodes = set(
        index.loc[
            proper & index["ligand_smiles_id"].notna(),
            node_column,
        ]
        .dropna()
        .astype(str)
    )
    artifacts: list[tuple[Path, str, str, bool]] = []
    for path in reciprocal_paths:
        relative_parts = path.relative_to(set_cover_root).parts
        partitions = {
            key: value
            for key, value in (
                part.split("=", maxsplit=1) for part in relative_parts[:-1]
            )
        }
        threshold = int(path.stem.split("=", maxsplit=1)[1])
        metric = partitions["metric"]
        artifacts.append(
            (
                path,
                metric,
                _cluster_column_name(
                    metric=metric,
                    cluster="set_cover",
                    directed=False,
                    threshold=threshold,
                    ligand=True,
                ),
                False,
            )
        )
    for path in directed_cover_paths:
        metric = path.parent.name.split("=", maxsplit=1)[1]
        threshold = int(path.stem.split("=", maxsplit=1)[1])
        artifacts.append(
            (
                path,
                metric,
                f"{metric}__{threshold}__ligand__directed_set_cover",
                True,
            )
        )
    LOG.info(
        "loading %d published ligand-cluster artifacts for %d ligand IDs",
        len(artifacts),
        len(node_ids),
    )
    cluster_columns: dict[str, Any] = {}
    started = time()
    for path_index, (path, metric, column, is_directed_cover) in enumerate(
        artifacts, start=1
    ):
        artifact_threshold = int(path.stem.split("=", maxsplit=1)[1])
        if (
            repair_started_ns is not None
            and path.stat().st_mtime_ns <= repair_started_ns
        ):
            raise ValueError(
                "ligand cluster artifact predates the targeted collation "
                f"repair: {path}"
            )
        columns = [node_column, "label", "centroid_ligand_id"]
        has_coverage_centrality = False
        if is_directed_cover:
            coverage_columns = {"coverage_count", "coverage_fraction"}
            available_columns = set(pq.read_schema(path).names)
            available_coverage_columns = coverage_columns.intersection(
                available_columns
            )
            if available_coverage_columns and (
                available_coverage_columns != coverage_columns
            ):
                missing = sorted(coverage_columns.difference(available_columns))
                raise ValueError(
                    "directed ligand cover has a partial coverage-centrality "
                    f"schema: {path}; missing={missing}"
                )
            has_coverage_centrality = available_coverage_columns == coverage_columns
            if has_coverage_centrality:
                columns.extend(sorted(coverage_columns))
        labels = pd.read_parquet(path, columns=columns)
        if labels[node_column].duplicated().any():
            raise ValueError(f"duplicate ligand IDs in cluster artifact: {path}")
        labels[node_column] = labels[node_column].astype(str)
        observed_nodes = set(labels[node_column])
        expected_nodes = (
            expected_fingerprint_nodes
            if is_chemical_cluster_metric(metric)
            else expected_score_nodes
        )
        if observed_nodes != expected_nodes:
            missing = sorted(expected_nodes.difference(observed_nodes))
            extra = sorted(observed_nodes.difference(expected_nodes))
            raise ValueError(
                "ligand cluster artifact does not cover the current eligible "
                f"ligand universe: {path}; missing={missing[:10]}, "
                f"extra={extra[:10]}"
            )
        aligned = labels.set_index(node_column)["label"].reindex(node_ids)
        cluster_columns[column] = aligned.astype("string[pyarrow]").array
        if metric == "tanimoto_similarity_ecfp4_1024" and artifact_threshold == 90:
            if "entry_pdb_id" not in index.columns:
                raise ValueError(
                    "the 90-percent Tanimoto set cover requires entry_pdb_id"
                )
            cluster_column = "ligand_tanimoto_ecfp4_1024_90_cluster"
            count_column = f"{cluster_column}_num_pdb_ids"
            label_by_node = labels.set_index(node_column)["label"]
            occurrences = pd.DataFrame(
                {
                    "label": index.loc[proper & holo, node_column]
                    .astype(str)
                    .map(label_by_node),
                    "entry_pdb_id": index.loc[proper & holo, "entry_pdb_id"].astype(
                        str
                    ),
                }
            ).dropna(subset=["label"])
            pdb_counts = occurrences.groupby("label", observed=True)[
                "entry_pdb_id"
            ].nunique()
            cluster_columns[cluster_column] = aligned.astype("string[pyarrow]").array
            cluster_columns[count_column] = (
                aligned.map(pdb_counts).astype("Int32").array
            )
        labels["centroid_ligand_id"] = labels["centroid_ligand_id"].astype(str)
        labels["is_centroid"] = labels[node_column].eq(labels["centroid_ligand_id"])
        centroid_counts = labels.groupby("label", observed=True)["is_centroid"].sum()
        invalid_labels = centroid_counts[centroid_counts.ne(1)].index.tolist()
        if invalid_labels:
            raise ValueError(
                "ligand set cover must have exactly one representative row "
                f"per label: {path}; invalid={invalid_labels[:10]}"
            )
        centroid_column = f"{column}__is_centroid"
        aligned_centroids = labels.set_index(node_column)["is_centroid"].reindex(
            node_ids
        )
        cluster_columns[centroid_column] = aligned_centroids.astype("boolean").array
        if is_directed_cover:
            if has_coverage_centrality:
                if labels[["coverage_count", "coverage_fraction"]].isna().any().any():
                    raise ValueError(
                        "directed ligand cover has missing coverage centrality: "
                        f"{path}"
                    )
                if (
                    labels["coverage_count"].lt(1).any()
                    or (
                        labels["coverage_fraction"].le(0)
                        | labels["coverage_fraction"].gt(1)
                    ).any()
                ):
                    raise ValueError(
                        "directed ligand cover has invalid coverage centrality: "
                        f"{path}"
                    )
                coverage_count_column = f"{column}__coverage_count"
                coverage_fraction_column = f"{column}__coverage_fraction"
                cluster_columns[coverage_count_column] = (
                    labels.set_index(node_column)["coverage_count"]
                    .reindex(node_ids)
                    .astype("Int32")
                    .array
                )
                cluster_columns[coverage_fraction_column] = (
                    labels.set_index(node_column)["coverage_fraction"]
                    .reindex(node_ids)
                    .astype("Float32")
                    .array
                )
        if path_index % 10 == 0 or path_index == len(artifacts):
            elapsed = time() - started
            rate = path_index / elapsed
            LOG.info(
                "cluster index progress: loaded=%d/%d rate=%.2f/s " "eta_seconds=%.1f",
                path_index,
                len(artifacts),
                rate,
                (len(artifacts) - path_index) / rate,
            )
    wide = pd.DataFrame(
        {node_column: node_ids.to_numpy(), **cluster_columns},
        copy=False,
    )
    LOG.info(
        "ligand cluster table complete: rows=%d columns=%d elapsed_seconds=%.1f",
        *wide.shape,
        time() - started,
    )
    return wide


def build_interface_cluster_table(
    *, index: pd.DataFrame, data_dir: Path
) -> pd.DataFrame:
    """Build one queryable cluster-assignment row per protein interface."""
    node_column = "system_id"
    directed_cover_root = data_dir / "interface_sampling" / "directed_set_cover"
    marker_path = data_dir / "index" / "collation.json"
    repair_started_ns: int | None = None
    if marker_path.is_file():
        with marker_path.open() as handle:
            marker = load(handle)
        if marker.get("status") == "requires_downstream_repair":
            repair_started_ns = marker_path.stat().st_mtime_ns
    directed_cover_paths = sorted(
        directed_cover_root.glob("metric=*/threshold=*.parquet")
    )
    if not directed_cover_paths:
        if repair_started_ns is not None:
            raise FileNotFoundError(
                "targeted collation repair has no rebuilt interface clusters"
            )
        if index[node_column].notna().any():
            raise FileNotFoundError(
                "non-empty interface annotation has no published interface clusters"
            )
        return pd.DataFrame({node_column: pd.Series(dtype="string")})

    node_ids = pd.Index(
        index[node_column].dropna().astype(str).unique(),
        name=node_column,
    )
    membership_path = data_dir / "index/interface_membership.parquet"
    if not membership_path.is_file():
        if len(node_ids):
            raise FileNotFoundError(
                "interface cluster expansion requires representative membership: "
                f"{membership_path}"
            )
        membership = pd.DataFrame(
            columns=[
                node_column,
                "representative_system_id",
                "side_1_half_interface_id",
                "side_2_half_interface_id",
            ]
        )
    else:
        membership = pd.read_parquet(
            membership_path,
            columns=[
                node_column,
                "representative_system_id",
                "side_1_half_interface_id",
                "side_2_half_interface_id",
            ],
        )
    membership[node_column] = membership[node_column].astype(str)
    if membership[node_column].duplicated().any():
        duplicate = membership.loc[
            membership[node_column].duplicated(), node_column
        ].iloc[0]
        raise ValueError(f"duplicate interface membership for {duplicate}")
    membership_ids = set(membership[node_column])
    expected_ids = set(node_ids)
    if membership_ids != expected_ids:
        missing = sorted(expected_ids.difference(membership_ids))
        extra = sorted(membership_ids.difference(expected_ids))
        raise ValueError(
            "interface representative membership does not cover the current "
            f"interface universe: missing={missing[:10]}, extra={extra[:10]}"
        )
    membership = membership.set_index(node_column).reindex(node_ids)
    expected_representatives = set(
        membership["representative_system_id"].dropna().astype(str)
    )
    expected_half_representatives = set(
        membership[["side_1_half_interface_id", "side_2_half_interface_id"]].stack()
    )
    artifacts: list[tuple[Path, str, int, str]] = []
    for path in directed_cover_paths:
        metric = path.parent.name.split("=", maxsplit=1)[1]
        threshold = int(path.stem.split("=", maxsplit=1)[1])
        artifacts.append((path, metric, threshold, "directed_set_cover"))

    cluster_columns: dict[str, Any] = {}
    started = time()
    for path_index, (path, metric, threshold, kind) in enumerate(artifacts, start=1):
        if (
            repair_started_ns is not None
            and path.stat().st_mtime_ns <= repair_started_ns
        ):
            raise ValueError(
                "interface cluster artifact predates the targeted collation "
                f"repair: {path}"
            )
        labels = pd.read_parquet(path, columns=[node_column, "label"])
        if labels[node_column].duplicated().any():
            raise ValueError(f"duplicate interface IDs in cluster artifact: {path}")
        labels[node_column] = labels[node_column].astype(str)
        observed_nodes = set(labels[node_column])
        expected_artifact_nodes = (
            expected_half_representatives
            if metric == "interface_side_qcov"
            else expected_representatives
        )
        if observed_nodes != expected_artifact_nodes:
            missing = sorted(expected_artifact_nodes.difference(observed_nodes))
            extra = sorted(observed_nodes.difference(expected_artifact_nodes))
            raise ValueError(
                "interface cluster artifact does not cover the current interface "
                f"universe: {path}; missing={missing[:10]}, extra={extra[:10]}"
            )
        label_lookup = labels.set_index(node_column)["label"]
        if metric == "interface_side_qcov":
            for side in (1, 2):
                representative_nodes = membership[
                    f"side_{side}_half_interface_id"
                ].astype(str)
                aligned = representative_nodes.map(label_lookup)
                column = f"{metric}__{threshold}__chain_{side}_{kind}"
                cluster_columns[column] = aligned.astype("string[pyarrow]").array
        else:
            column = f"{metric}__{threshold}__{kind}"
            aligned = (
                membership["representative_system_id"].astype(str).map(label_lookup)
            )
            cluster_columns[column] = aligned.astype("string[pyarrow]").array
        if path_index % 10 == 0 or path_index == len(artifacts):
            elapsed = time() - started
            rate = path_index / elapsed
            LOG.info(
                "interface cluster index progress: loaded=%d/%d rate=%.2f/s "
                "eta_seconds=%.1f",
                path_index,
                len(artifacts),
                rate,
                (len(artifacts) - path_index) / rate,
            )

    wide = pd.DataFrame(
        {node_column: node_ids.to_numpy(), **cluster_columns},
        copy=False,
    )
    LOG.info(
        "interface cluster table complete: rows=%d columns=%d elapsed_seconds=%.1f",
        *wide.shape,
        time() - started,
    )
    return wide


def add_ligand_similarity_columns(
    *, index: pd.DataFrame, data_dir: Path
) -> pd.DataFrame:
    """Merge unique-SMILES and cofactor annotations into each ligand."""
    annotation_path = (
        data_dir / "fingerprints" / "ligand_similarity_annotations.parquet"
    )
    if not annotation_path.is_file():
        raise FileNotFoundError(
            f"missing ligand similarity annotations: {annotation_path}"
        )
    marker_path = data_dir / "index" / "collation.json"
    if marker_path.is_file():
        with marker_path.open() as handle:
            marker = load(handle)
        if (
            marker.get("status") == "requires_downstream_repair"
            and annotation_path.stat().st_mtime_ns <= marker_path.stat().st_mtime_ns
        ):
            raise ValueError(
                "ligand similarity annotations predate the targeted " "collation repair"
            )
    annotations = pd.read_parquet(annotation_path)
    artifact_smiles_column = "ligand_rdkit_canonical_smiles"
    index_smiles_column = "ligand_smiles"
    if annotations[artifact_smiles_column].duplicated().any():
        raise ValueError("ligand similarity annotations contain duplicate SMILES")
    proper_holo = index["ligand_is_proper"].fillna(False).astype(bool) & index[
        "system_type"
    ].eq("holo")
    expected_smiles = set(
        index.loc[proper_holo, index_smiles_column]
        .dropna()
        .astype(str)
        .loc[lambda values: values.ne("")]
    )
    observed_smiles = set(annotations[artifact_smiles_column].dropna().astype(str))
    if observed_smiles != expected_smiles:
        missing = sorted(expected_smiles.difference(observed_smiles))
        extra = sorted(observed_smiles.difference(expected_smiles))
        raise ValueError(
            "ligand similarity annotations do not cover the current proper "
            f"holo SMILES universe: missing={missing[:10]}, extra={extra[:10]}"
        )
    annotations = annotations.rename(
        columns={artifact_smiles_column: index_smiles_column}
    )
    replacement_columns = set(annotations.columns).difference({index_smiles_column})
    obsolete_columns = {
        "ligand_tanimoto_ecfp4_1024_90_cluster",
        "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids",
    }
    result = index.drop(
        columns=list(
            replacement_columns.intersection(index.columns)
            | obsolete_columns.intersection(index.columns)
        )
    ).merge(
        annotations,
        on=index_smiles_column,
        how="left",
        validate="many_to_one",
    )
    has_smiles = result[index_smiles_column].notna() & result[index_smiles_column].ne(
        ""
    )
    eligible = has_smiles
    if "ligand_is_proper" in result:
        eligible &= result["ligand_is_proper"].fillna(False)
    if "system_type" in result:
        eligible &= result["system_type"].eq("holo")
    if result.loc[eligible, "ligand_smiles_id"].isna().any():
        raise ValueError(
            "some proper canonical ligand SMILES lack similarity annotations"
        )
    for column in replacement_columns:
        if pd.api.types.is_bool_dtype(result[column].dtype):
            result[column] = result[column].astype("boolean")
        result.loc[~eligible, column] = pd.NA
    result["ligand_smiles_id"] = result["ligand_smiles_id"].astype("Int32")
    if "ligand_is_cofactor_like" in result:
        result["ligand_is_cofactor_like"] = result["ligand_is_cofactor_like"].astype(
            "boolean"
        )
    return result


def add_ligand_3d_score_ability_column(
    *, index: pd.DataFrame, data_dir: Path
) -> pd.DataFrame:
    """Merge occurrence-level canonical-SDF scoreability into the index."""
    expected = index["ligand_id"].notna()
    if "system_type" in index:
        expected &= index["system_type"].eq("holo")
    if "ligand_is_proper" in index:
        expected &= index["ligand_is_proper"].fillna(False)
    scoreability_column = "ligand_is_3d_score_able"
    if (
        scoreability_column in index
        and not index.loc[expected, scoreability_column].isna().any()
    ):
        abilities = index.loc[
            index["ligand_id"].notna(), ["ligand_id", scoreability_column]
        ].drop_duplicates()
        conflicts = abilities.groupby("ligand_id")[scoreability_column].nunique()
        if (conflicts > 1).any():
            raise ValueError("conflicting 3D-scoreability annotations for a ligand")
        result = index.copy()
        result.loc[~expected, scoreability_column] = False
        result[scoreability_column] = result[scoreability_column].astype("boolean")
        return result

    ligand_dataset = data_dir / "ligands"
    if not ligand_dataset.is_dir():
        raise FileNotFoundError(f"missing ligand annotation dataset: {ligand_dataset}")
    abilities = pd.read_parquet(
        ligand_dataset,
        columns=["ligand_id", "ligand_is_3d_score_able"],
    ).drop_duplicates()
    conflicts = abilities.groupby("ligand_id")["ligand_is_3d_score_able"].nunique()
    if (conflicts > 1).any():
        raise ValueError("conflicting 3D-scoreability annotations for a ligand")
    abilities = abilities.drop_duplicates(subset=["ligand_id"])
    result = index.drop(columns=["ligand_is_3d_score_able"], errors="ignore").merge(
        abilities,
        on="ligand_id",
        how="left",
        validate="many_to_one",
    )
    expected = result["ligand_id"].notna()
    if "system_type" in result:
        expected &= result["system_type"].eq("holo")
    if "ligand_is_proper" in result:
        expected &= result["ligand_is_proper"].fillna(False)
    if result.loc[expected, "ligand_is_3d_score_able"].isna().any():
        raise ValueError("some proper holo ligands lack a 3D-scoreability annotation")
    result.loc[~expected, "ligand_is_3d_score_able"] = False
    result["ligand_is_3d_score_able"] = result["ligand_is_3d_score_able"].astype(
        "boolean"
    )
    return result


def update_index_ligand_3d_score_ability(*, data_dir: Path) -> None:
    """Publish distributed 3D-scoreability annotations before protein scoring."""
    index_path = data_dir / "index" / "annotation_table.parquet"
    index = pd.read_parquet(index_path)
    add_ligand_3d_score_ability_column(index=index, data_dir=data_dir).to_parquet(
        index_path, index=False
    )


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
        # the single "oligo" flag was split into per-family mono/oligo columns;
        # aggregate each to match the ligand_is_* annotation and the documented
        # system_ligand_has_* schema (see column_descriptions/extra.tsv)
        "monosaccharide",
        "oligosaccharide",
        "mononucleotide",
        "oligonucleotide",
        "monopeptide",
        "oligopeptide",
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


def _is_ligand_cluster_column(column: str) -> bool:
    """Return whether a column belongs in the ligand-cluster sidecar."""
    return column in {
        "ligand_tanimoto_ecfp4_1024_90_cluster",
        "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids",
    } or (
        "__ligand__" in column
        and column.endswith(
            (
                "__component",
                "__community",
                "__set_cover",
                "__set_cover__is_centroid",
                "__directed_set_cover",
                "__directed_set_cover__is_centroid",
                "__directed_set_cover__coverage_count",
                "__directed_set_cover__coverage_fraction",
            )
        )
    )


def _is_interface_cluster_column(column: str) -> bool:
    """Return whether a column belongs in the interface-cluster sidecar."""
    return column.startswith(("interface_qcov__", "interface_side_qcov__")) and (
        "component" in column or "community" in column or "set_cover" in column
    )


def finalize_index(*, data_dir: Path) -> pd.DataFrame:
    """Publish enriched annotations and separate cluster-assignment tables."""
    started = time()
    index_path = data_dir / "index" / "annotation_table.parquet"
    if not index_path.is_file():
        raise FileNotFoundError(index_path)
    LOG.info("loading annotation index for final enrichment: %s", index_path)
    index = pd.read_parquet(index_path)
    index.drop(columns=["uniqueness"], errors="ignore", inplace=True)
    LOG.info("loaded annotation index: rows=%d columns=%d", *index.shape)
    index = add_ligand_3d_score_ability_column(index=index, data_dir=data_dir)
    LOG.info("merged ligand 3D-scoreability annotations")
    index = add_ligand_similarity_columns(index=index, data_dir=data_dir)
    LOG.info("merged ligand similarity annotations")
    ligand_clusters = build_ligand_cluster_table(index=index, data_dir=data_dir)
    index.drop(
        columns=[column for column in index if _is_ligand_cluster_column(column)],
        inplace=True,
    )
    interface_path = data_dir / "index" / "interface_annotation_table.parquet"
    interface_index: pd.DataFrame | None = None
    interface_clusters: pd.DataFrame | None = None
    if interface_path.is_file():
        interface_index = pd.read_parquet(interface_path)
        interface_clusters = build_interface_cluster_table(
            index=interface_index,
            data_dir=data_dir,
        )
        interface_index.drop(
            columns=[
                column
                for column in interface_index
                if _is_interface_cluster_column(column)
            ],
            inplace=True,
        )

    temporary = index_path.with_suffix(".tmp.parquet")
    temporary_interface = interface_path.with_suffix(".tmp.parquet")
    ligand_clusters_path = data_dir / "index" / "ligand_clusters.parquet"
    interface_clusters_path = data_dir / "index" / "interface_clusters.parquet"
    temporary_ligand_clusters = ligand_clusters_path.with_suffix(".tmp.parquet")
    temporary_interface_clusters = interface_clusters_path.with_suffix(".tmp.parquet")
    index_marker_removed = False
    try:
        LOG.info("staging annotation and cluster tables")
        index.to_parquet(temporary, index=False)
        ligand_clusters.to_parquet(temporary_ligand_clusters, index=False)
        if interface_index is not None:
            interface_index.to_parquet(temporary_interface, index=False)
            assert interface_clusters is not None
            interface_clusters.to_parquet(
                temporary_interface_clusters,
                index=False,
            )
        # The primary annotation table is the readability marker for this
        # generation. Install it last so an interrupted update fails closed
        # instead of exposing mismatched annotation and cluster tables.
        index_path.unlink()
        index_marker_removed = True
        if interface_index is not None:
            temporary_interface.replace(interface_path)
            temporary_interface_clusters.replace(interface_clusters_path)
        temporary_ligand_clusters.replace(ligand_clusters_path)
        temporary.replace(index_path)
    except BaseException:
        if index_marker_removed:
            index_path.unlink(missing_ok=True)
        raise
    finally:
        temporary.unlink(missing_ok=True)
        temporary_interface.unlink(missing_ok=True)
        temporary_ligand_clusters.unlink(missing_ok=True)
        temporary_interface_clusters.unlink(missing_ok=True)
    if interface_index is not None:
        LOG.info(
            "wrote interface annotations: rows=%d columns=%d",
            *interface_index.shape,
        )
    LOG.info(
        "final table publication complete: rows=%d columns=%d elapsed_seconds=%.1f",
        *index.shape,
        time() - started,
    )
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
        "chain_sequence",
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


def create_entry_biounit_chain_index(
    *, data_dir: Path, force_update: bool = False
) -> pd.DataFrame:
    """Collate biological-assembly chain membership into one parquet."""
    output = data_dir / "index" / "entry_biounit_chains.parquet"
    output.parent.mkdir(exist_ok=True, parents=True)
    if output.exists() and not force_update:
        return pd.read_parquet(output)

    parts = sorted((data_dir / "raw_entries").glob("*/*/entry_biounit_chains.parquet"))
    columns = [
        "entry_pdb_id",
        "biounit_id",
        "chain_instance",
        "chain_asym_id",
        "chain_role",
    ]
    chains = (
        pd.concat([pd.read_parquet(path) for path in parts], ignore_index=True)
        if parts
        else pd.DataFrame(columns=columns)
    ).reindex(columns=columns)
    key = ["entry_pdb_id", "biounit_id", "chain_instance"]
    if not chains.empty and chains.duplicated(key).any():
        duplicates = chains.loc[chains.duplicated(key, keep=False), key]
        raise ValueError(
            "duplicate biological-assembly chain membership: "
            f"{duplicates.drop_duplicates().to_dict(orient='records')}"
        )
    chains = chains.sort_values(key, ignore_index=True)
    chains.to_parquet(output, index=False, row_group_size=100_000)
    return chains


def create_index(*, data_dir: Path, force_update: bool = False) -> pd.DataFrame:
    """
    Create the index
    """
    index = data_dir / "index" / "annotation_table.parquet"
    index.parent.mkdir(exist_ok=True, parents=True)
    create_entry_chain_index(data_dir=data_dir, force_update=force_update)
    create_entry_biounit_chain_index(data_dir=data_dir, force_update=force_update)
    create_entry_source_index(data_dir=data_dir, force_update=force_update)

    if not index.exists() or force_update:
        dfs = []
        annotation_parts = data_dir / "raw_entries"
        # sort for deterministic collation order; glob yields filesystem order
        for i, path in enumerate(sorted(annotation_parts.glob("*/*.parquet"))):
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


def rename_clusters(*, data_dir: Path) -> None:
    """
    Rename cluster files to match the hive layout convention.

    Parameters
    ----------
    data_dir : Path
        plinder root dir
    """
    cluster_dir = data_dir / "ligand_clusters"
    cluster_paths = [path for path in cluster_dir.rglob("*") if path.is_file()]
    for path in cluster_paths:
        if path.name == "data.parquet":
            continue
        base = path.parent
        name = path.stem
        apath = base / name / "data.parquet"
        apath.parent.mkdir(exist_ok=True, parents=True)
        shutil.move(path, apath)
