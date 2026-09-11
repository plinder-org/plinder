# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import wraps
from hashlib import md5
from json import dumps, load
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Any, Callable, Optional, TypeVar

import pandas as pd
import pyarrow.parquet as pq
from omegaconf import DictConfig, OmegaConf

from plinder.core.scores.metrics import (
    CHEMICAL_CLUSTER_SUMMARY_COLUMNS,
    CHEMICAL_CLUSTER_SUMMARY_THRESHOLD,
    is_chemical_cluster_metric,
)
from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger

if TYPE_CHECKING:
    from plinder.data.annotations.get_similarity_scores import Scorer


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
        try:
            result = func(*args, **kwargs)
            log.info(f"runtime succeeded: {time() - ts:>9.2f}s")
        except Exception as e:
            log.error(f"runtime failed: {time() - ts:>9.2f}s")
            log.error(f"{name} failed with: {repr(e)}")
            raise
        return result

    return wrapped


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
    return stage not in skip and (not run or stage in run)


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
        return [[]]

    return inner


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
        metric for metric, _ in set_cover_keys if not is_chemical_cluster_metric(metric)
    )
    invalid_directed_metrics = sorted(
        metric
        for metric, _ in directed_cover_keys
        if is_chemical_cluster_metric(metric)
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
                f"{metric}__{threshold}__ligand__set_cover",
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
        summary_column = CHEMICAL_CLUSTER_SUMMARY_COLUMNS.get(metric)
        if (
            summary_column is not None
            and artifact_threshold == CHEMICAL_CLUSTER_SUMMARY_THRESHOLD
        ):
            if "entry_pdb_id" not in index.columns:
                raise ValueError(
                    f"the {artifact_threshold}-percent {metric} set cover "
                    "requires entry_pdb_id"
                )
            cluster_column = summary_column
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
    obsolete_columns = _chemical_cluster_summary_columns()
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


def _chemical_cluster_summary_columns() -> set[str]:
    """Return the 90-percent cluster sidecar columns of every chemical metric."""
    return {
        name
        for column in CHEMICAL_CLUSTER_SUMMARY_COLUMNS.values()
        for name in (column, f"{column}_num_pdb_ids")
    }


def _is_ligand_cluster_column(column: str) -> bool:
    """Return whether a column belongs in the ligand-cluster sidecar."""
    return column in _chemical_cluster_summary_columns() or (
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
