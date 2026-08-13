# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import gc
import hashlib
import heapq
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile, copytree, rmtree
from textwrap import dedent
from time import time
from typing import Any, Callable, Iterable, Literal, Mapping, Sequence, TypeVar, cast

if sys.platform == "darwin":
    # For macOS only: allow multiple OpenMP runtimes to coexist
    # (needed on macOS with conda)
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import networkit as nk
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from numpy.typing import NDArray

from plinder.core.scores.metrics import GATED_LIGAND_DIAGNOSTIC_METRICS
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

T = TypeVar("T")
ClusterEntity = Literal["ligand", "interface"]

COMPONENT_EDGE_COLUMNS = ["query_node", "target_node"]
COMPONENT_REDUCTION_VERSION = 3
COMPONENT_REDUCTION_DIRECTIONS = (False,)
SYMMETRIC_EDGE_BUCKET_COUNT = 64
SYMMETRIC_EDGE_COLUMNS = ["query_node", "target_node", "similarity"]
INTERFACE_CLUSTER_METRICS = frozenset({"interface_qcov", "interface_side_qcov"})
INTERFACE_REPRESENTATIVES = Path("index/interface_representatives.parquet")
INTERFACE_HALF_REPRESENTATIVES = Path(
    "index/interface_half_representatives.parquet"
)
INTERFACE_MEMBERSHIP = Path("index/interface_membership.parquet")


def _cluster_root(data_dir: Path, entity_type: ClusterEntity) -> Path:
    if entity_type == "ligand":
        return data_dir / "ligand_clusters"
    if entity_type == "interface":
        return data_dir / "interface_clusters"
    raise ValueError(f"unsupported cluster entity type: {entity_type}")


def _sampling_root(data_dir: Path, entity_type: ClusterEntity) -> Path:
    if entity_type == "ligand":
        return data_dir / "ligand_sampling"
    if entity_type == "interface":
        return data_dir / "interface_sampling"
    raise ValueError(f"unsupported cluster entity type: {entity_type}")


def _cluster_node_column(entity_type: ClusterEntity) -> str:
    return "ligand_id" if entity_type == "ligand" else "system_id"


def _symmetric_edge_plan_path(data_dir: Path, entity_type: ClusterEntity) -> Path:
    return _cluster_root(data_dir, entity_type) / "symmetric_edges" / "plan.json"


def _component_node_universe_paths(
    data_dir: Path, entity_type: ClusterEntity
) -> tuple[Path, Path]:
    root = _cluster_root(data_dir, entity_type) / "reductions"
    return root / "node_universe.parquet", root / "node_universe.json"


def _component_node_universe_sources(
    data_dir: Path, entity_type: ClusterEntity
) -> dict[str, dict[str, str | int]]:
    """Return every source that defines one cached clustering universe."""
    if entity_type == "ligand":
        paths = {"annotation": data_dir / "index/annotation_table.parquet"}
    else:
        paths = {
            "annotation": data_dir / "index/interface_annotation_table.parquet",
            "representatives": data_dir / INTERFACE_REPRESENTATIVES,
            "half_representatives": data_dir / INTERFACE_HALF_REPRESENTATIVES,
            "membership": data_dir / INTERFACE_MEMBERSHIP,
        }
    return {name: _component_source_signature(path) for name, path in paths.items()}


def _empty_component_edges() -> pd.DataFrame:
    """Return an empty, consistently typed component-edge table."""
    return pd.DataFrame(
        {
            "query_node": pd.Series(dtype="int64"),
            "target_node": pd.Series(dtype="int64"),
        }
    )


def _reduce_undirected_component_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Replace an undirected edge table with a connectivity-equivalent forest.

    The returned forest contains at most ``number_of_nodes - number_of_components``
    edges.  Replacing any subset of an undirected graph with such a forest preserves
    the connected components of the union with every other edge subset, which makes
    this reduction safe to apply repeatedly to independently generated shards.
    """
    if edges.empty:
        return _empty_component_edges()
    missing = sorted(set(COMPONENT_EDGE_COLUMNS).difference(edges.columns))
    if missing:
        raise ValueError(f"component edges are missing columns {missing}")

    edge_array = edges[COMPONENT_EDGE_COLUMNS].to_numpy(dtype=np.int64, copy=False)
    edge_array = edge_array[edge_array[:, 0] != edge_array[:, 1]]
    if not len(edge_array):
        return _empty_component_edges()
    edge_array = np.column_stack(
        (
            np.minimum(edge_array[:, 0], edge_array[:, 1]),
            np.maximum(edge_array[:, 0], edge_array[:, 1]),
        )
    )
    edge_array = np.unique(edge_array, axis=0)
    nodes = np.unique(edge_array)
    query_codes = np.searchsorted(nodes, edge_array[:, 0]).astype(np.uint, copy=False)
    target_codes = np.searchsorted(nodes, edge_array[:, 1]).astype(np.uint, copy=False)
    graph = nk.GraphFromCoo(
        (query_codes, target_codes),
        n=len(nodes),
        directed=False,
    )
    components = nk.components.ConnectedComponents(graph)
    components.run()

    query_parts: list[NDArray[np.int64]] = []
    target_parts: list[NDArray[np.int64]] = []
    for component in components.getComponents():
        members = np.sort(nodes[np.asarray(component, dtype=np.int64)])
        if len(members) < 2:
            continue
        query_parts.append(np.full(len(members) - 1, members[0], dtype=np.int64))
        target_parts.append(members[1:].astype(np.int64, copy=False))
    if not query_parts:
        return _empty_component_edges()
    return pd.DataFrame(
        {
            "query_node": np.concatenate(query_parts),
            "target_node": np.concatenate(target_parts),
        }
    )


def _reduce_directed_component_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Replace directed edges with an exact reachability-preserving subgraph.

    Each local strongly connected component is represented by a bidirectional
    star, while unique arcs in its condensation DAG are retained.  Reachability
    is therefore unchanged, even after independently reduced shards are united.
    Unlike an undirected forest, the condensation DAG can remain large because
    an arbitrary directed graph has no generally linear-size reachability
    representation.
    """
    if edges.empty:
        return _empty_component_edges()
    missing = sorted(set(COMPONENT_EDGE_COLUMNS).difference(edges.columns))
    if missing:
        raise ValueError(f"component edges are missing columns {missing}")
    edge_array = edges[COMPONENT_EDGE_COLUMNS].to_numpy(dtype=np.int64, copy=False)
    edge_array = edge_array[edge_array[:, 0] != edge_array[:, 1]]
    if not len(edge_array):
        return _empty_component_edges()
    edge_array = np.unique(edge_array, axis=0)
    nodes = np.unique(edge_array)
    query_codes = np.searchsorted(nodes, edge_array[:, 0]).astype(np.uint, copy=False)
    target_codes = np.searchsorted(nodes, edge_array[:, 1]).astype(np.uint, copy=False)
    graph = nk.GraphFromCoo(
        (query_codes, target_codes),
        n=len(nodes),
        directed=True,
    )
    components = nk.components.StronglyConnectedComponents(graph)
    components.run()

    component_by_node = np.empty(len(nodes), dtype=np.int64)
    representatives: list[int] = []
    query_parts: list[NDArray[np.int64]] = []
    target_parts: list[NDArray[np.int64]] = []
    for component_index, component in enumerate(components.getComponents()):
        component_codes = np.asarray(component, dtype=np.int64)
        component_by_node[component_codes] = component_index
        members = np.sort(nodes[component_codes])
        representative = int(members[0])
        representatives.append(representative)
        if len(members) < 2:
            continue
        others = members[1:].astype(np.int64, copy=False)
        roots = np.full(len(others), representative, dtype=np.int64)
        query_parts.extend([roots, others])
        target_parts.extend([others, roots])

    query_components = component_by_node[query_codes]
    target_components = component_by_node[target_codes]
    between_components = query_components != target_components
    if between_components.any():
        representatives_array = np.asarray(representatives, dtype=np.int64)
        condensation_edges = np.column_stack(
            (
                representatives_array[query_components[between_components]],
                representatives_array[target_components[between_components]],
            )
        )
        condensation_edges = np.unique(condensation_edges, axis=0)
        query_parts.append(condensation_edges[:, 0])
        target_parts.append(condensation_edges[:, 1])
    if not query_parts:
        return _empty_component_edges()
    return pd.DataFrame(
        {
            "query_node": np.concatenate(query_parts),
            "target_node": np.concatenate(target_parts),
        }
    )


class _ComponentEdgeAccumulator:
    """Bound memory while repeatedly reducing independent edge batches."""

    def __init__(
        self,
        fan_in: int,
        reducer: Callable[[pd.DataFrame], pd.DataFrame],
    ) -> None:
        if fan_in < 2:
            raise ValueError("component forest fan-in must be at least two")
        self.fan_in = fan_in
        self.reducer = reducer
        self.levels: list[list[pd.DataFrame]] = []

    def add(self, edges: pd.DataFrame) -> None:
        forest = self.reducer(edges)
        if not forest.empty:
            self._push(forest, level=0)

    def _push(self, forest: pd.DataFrame, *, level: int) -> None:
        while len(self.levels) <= level:
            self.levels.append([])
        self.levels[level].append(forest)
        if len(self.levels[level]) < self.fan_in:
            return
        merged = self.reducer(pd.concat(self.levels[level], ignore_index=True))
        self.levels[level] = []
        if not merged.empty:
            self._push(merged, level=level + 1)

    def finish(self) -> pd.DataFrame:
        forests = [forest for level in self.levels for forest in level]
        if not forests:
            return _empty_component_edges()
        while len(forests) > 1:
            reduced: list[pd.DataFrame] = []
            for start in range(0, len(forests), self.fan_in):
                group = forests[start : start + self.fan_in]
                if len(group) == 1:
                    reduced.append(group[0])
                else:
                    reduced.append(self.reducer(pd.concat(group, ignore_index=True)))
            forests = reduced
        return forests[0]


def _component_labels_from_edges(
    *, edges: pd.DataFrame, nodes: Sequence[str], directed: bool
) -> pd.DataFrame:
    """Assign deterministic size-ordered labels to a reduced edge table."""
    if not nodes:
        return pd.DataFrame(
            {
                "ligand_id": pd.Series(dtype="object"),
                "label": pd.Series(dtype="object"),
            }
        )
    if edges.empty:
        graph = nk.Graph(len(nodes), directed=directed)
    else:
        edge_array = edges[COMPONENT_EDGE_COLUMNS].to_numpy(dtype=np.int64, copy=False)
        graph = nk.GraphFromCoo(
            (
                edge_array[:, 0].astype(np.uint, copy=False),
                edge_array[:, 1].astype(np.uint, copy=False),
            ),
            n=len(nodes),
            directed=directed,
        )
    if directed:
        components = nk.components.StronglyConnectedComponents(graph)
    else:
        components = nk.components.ConnectedComponents(graph)
    components.run()
    ordered = sorted(
        components.getComponents(),
        key=lambda component: (
            -len(component),
            min(nodes[node] for node in component),
        ),
    )
    labels = [""] * len(nodes)
    for component_index, component in enumerate(ordered):
        label = f"c{component_index}"
        for node in component:
            labels[node] = label
    return pd.DataFrame({"ligand_id": list(nodes), "label": labels})


def make_exact_threshold_components(
    *,
    edge_batches: Iterable[pd.DataFrame],
    all_nodes: Iterable[str],
    thresholds: Sequence[int],
    query_column: str = "query_node",
    target_column: str = "target_node",
    similarity_column: str = "similarity",
    forest_fan_in: int = 8,
    directed: bool = False,
) -> dict[int, pd.DataFrame]:
    """Build exact weak or strong components at many thresholds in one edge scan.

    Each input edge is assigned to exactly one descending score band.  Every
    undirected batch is reduced to a spanning forest.  Directed batches retain a
    reachability-preserving SCC/condensation representation.  These reduced
    tables are recursively merged before the next threshold is labelled.
    """
    nodes, reductions = reduce_threshold_component_edges(
        edge_batches=edge_batches,
        all_nodes=all_nodes,
        thresholds=thresholds,
        query_column=query_column,
        target_column=target_column,
        similarity_column=similarity_column,
        forest_fan_in=forest_fan_in,
        directed_modes=[directed],
    )
    return _labels_from_threshold_reductions(
        nodes=nodes,
        reductions=reductions[directed],
        directed=directed,
    )


def reduce_threshold_component_edges(
    *,
    edge_batches: Iterable[pd.DataFrame],
    all_nodes: Iterable[str],
    thresholds: Sequence[int],
    query_column: str = "query_node",
    target_column: str = "target_node",
    similarity_column: str = "similarity",
    forest_fan_in: int = 8,
    directed_modes: Sequence[bool] = (False, True),
) -> tuple[list[str], dict[bool, dict[int, pd.DataFrame]]]:
    """Reduce threshold bands for weak and strong components in one edge scan."""
    ordered_thresholds = sorted(set(thresholds), reverse=True)
    if not ordered_thresholds:
        raise ValueError("at least one component threshold is required")
    if ordered_thresholds[0] > 100 or ordered_thresholds[-1] < 0:
        raise ValueError("component thresholds must be in [0, 100]")
    modes = tuple(dict.fromkeys(directed_modes))
    if not modes:
        raise ValueError("at least one component direction mode is required")
    nodes = sorted(set(map(str, all_nodes)))
    node_indexes = {node: index for index, node in enumerate(nodes)}
    accumulators = {
        directed: {
            threshold: _ComponentEdgeAccumulator(
                forest_fan_in,
                (
                    _reduce_directed_component_edges
                    if directed
                    else _reduce_undirected_component_edges
                ),
            )
            for threshold in ordered_thresholds
        }
        for directed in modes
    }
    required = {query_column, target_column, similarity_column}

    for batch in edge_batches:
        missing = sorted(required.difference(batch.columns))
        if missing:
            raise ValueError(f"component edge batch is missing columns {missing}")
        if batch.empty:
            continue
        values = batch[[query_column, target_column, similarity_column]].dropna()
        if values.empty:
            continue
        query_indexes = values[query_column].astype(str).map(node_indexes)
        target_indexes = values[target_column].astype(str).map(node_indexes)
        similarities = pd.to_numeric(values[similarity_column], errors="coerce")
        valid = query_indexes.notna() & target_indexes.notna() & similarities.notna()
        if not valid.any():
            continue
        query_array = query_indexes.loc[valid].to_numpy(dtype=np.int64, copy=False)
        target_array = target_indexes.loc[valid].to_numpy(dtype=np.int64, copy=False)
        similarity_array = similarities.loc[valid].to_numpy(dtype=float, copy=False)
        upper_bound = np.inf
        for threshold in ordered_thresholds:
            in_band = (similarity_array >= threshold) & (similarity_array < upper_bound)
            if in_band.any():
                band_edges = pd.DataFrame(
                    {
                        "query_node": query_array[in_band],
                        "target_node": target_array[in_band],
                    }
                )
                for directed in modes:
                    accumulators[directed][threshold].add(band_edges)
            upper_bound = threshold

    return nodes, {
        directed: {
            threshold: accumulators[directed][threshold].finish()
            for threshold in ordered_thresholds
        }
        for directed in modes
    }


def _labels_from_threshold_reductions(
    *,
    nodes: Sequence[str],
    reductions: Mapping[int, pd.DataFrame],
    directed: bool,
) -> dict[int, pd.DataFrame]:
    """Compose descending score-band reductions into threshold labels."""
    reducer = (
        _reduce_directed_component_edges
        if directed
        else _reduce_undirected_component_edges
    )
    result: dict[int, pd.DataFrame] = {}
    cumulative_edges = _empty_component_edges()
    for threshold in sorted(reductions, reverse=True):
        band_edges = reductions[threshold]
        cumulative_edges = reducer(
            pd.concat([cumulative_edges, band_edges], ignore_index=True)
        )
        result[threshold] = _component_labels_from_edges(
            edges=cumulative_edges,
            nodes=nodes,
            directed=directed,
        )
    return result


def make_exact_weak_and_strong_threshold_components(
    *,
    edge_batches: Iterable[pd.DataFrame],
    all_nodes: Iterable[str],
    thresholds: Sequence[int],
    query_column: str = "query_node",
    target_column: str = "target_node",
    similarity_column: str = "similarity",
    forest_fan_in: int = 8,
) -> dict[bool, dict[int, pd.DataFrame]]:
    """Build weak and strong component labels from a single edge scan."""
    nodes, reductions = reduce_threshold_component_edges(
        edge_batches=edge_batches,
        all_nodes=all_nodes,
        thresholds=thresholds,
        query_column=query_column,
        target_column=target_column,
        similarity_column=similarity_column,
        forest_fan_in=forest_fan_in,
        directed_modes=[False, True],
    )
    return {
        directed: _labels_from_threshold_reductions(
            nodes=nodes,
            reductions=reductions[directed],
            directed=directed,
        )
        for directed in [False, True]
    }


def _component_node_hash(nodes: Sequence[str]) -> str:
    digest = hashlib.sha256()
    for node in nodes:
        encoded = node.encode()
        digest.update(len(encoded).to_bytes(8, byteorder="little"))
        digest.update(encoded)
    return digest.hexdigest()


def _component_selection_hash(
    *, nodes: Sequence[str], eligible_systems: set[str] | None
) -> str:
    selected = list(nodes)
    if eligible_systems is not None:
        selected.extend(f"system:{system_id}" for system_id in sorted(eligible_systems))
    return _component_node_hash(selected)


def _component_source_signature(path: Path) -> dict[str, str | int]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _raw_component_score_sources(
    *,
    data_dir: Path,
    chemical: bool,
    entity_type: ClusterEntity = "ligand",
) -> list[Path]:
    """Return raw score sources before reciprocal edge aggregation."""
    if entity_type == "interface":
        if chemical:
            raise ValueError("interface clustering does not use chemical scores")
        return sorted((data_dir / "interface_scores").glob("shard=*.parquet"))
    if chemical:
        return sorted((data_dir / "ligand_scores").glob("*.parquet"))
    return sorted((data_dir / "scores" / "search_db=holo").rglob("*.parquet"))


def _symmetric_edge_plan_hash(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def prepare_symmetric_edge_plan(
    *,
    data_dir: Path,
    metrics: Sequence[str],
    source_batch_size: int = 20,
    bucket_count: int = SYMMETRIC_EDGE_BUCKET_COUNT,
    entity_type: ClusterEntity = "ligand",
) -> dict[str, Any]:
    """Plan sharded reciprocal-minimum edges without scanning score contents."""
    if source_batch_size < 1:
        raise ValueError("symmetric-edge source batch size must be positive")
    if bucket_count < 1:
        raise ValueError("symmetric-edge bucket count must be positive")
    selected_metrics = list(dict.fromkeys(metrics))
    chemical_metric = "tanimoto_similarity_ecfp4_1024"
    if entity_type == "interface" and (
        not selected_metrics
        or not set(selected_metrics).issubset(INTERFACE_CLUSTER_METRICS)
    ):
        raise ValueError(
            "interface clustering supports interface_qcov and " "interface_side_qcov"
        )
    batches: list[dict[str, Any]] = []
    source_groups = (
        [("interface", selected_metrics)]
        if entity_type == "interface"
        else [
            (
                "score",
                [metric for metric in selected_metrics if metric != chemical_metric],
            ),
            (
                "chemical",
                [metric for metric in selected_metrics if metric == chemical_metric],
            ),
        ]
    )
    for kind, source_metrics in source_groups:
        if not source_metrics:
            continue
        sources = _raw_component_score_sources(
            data_dir=data_dir,
            chemical=kind == "chemical",
            entity_type=entity_type,
        )
        if not sources:
            if entity_type != "interface":
                raise FileNotFoundError(f"no raw {kind} score sources found")
            nodes, _ = component_node_universe(
                data_dir=data_dir,
                metric=source_metrics[0],
                entity_type=entity_type,
            )
            if nodes:
                raise FileNotFoundError(
                    "no raw interface score sources found for a non-empty "
                    "interface universe"
                )
        for start in range(0, len(sources), source_batch_size):
            selected_sources = sources[start : start + source_batch_size]
            relative_sources = [
                source.resolve().relative_to(data_dir.resolve()).as_posix()
                for source in selected_sources
            ]
            key = hashlib.sha256(
                json.dumps(
                    [kind, relative_sources],
                    separators=(",", ":"),
                ).encode()
            ).hexdigest()[:16]
            batches.append(
                {
                    "key": key,
                    "kind": kind,
                    "metrics": source_metrics,
                    "sources": [
                        _component_source_signature(source)
                        for source in selected_sources
                    ],
                }
            )
    payload: dict[str, Any] = {
        "version": 1,
        "entity_type": entity_type,
        "metrics": selected_metrics,
        "source_batch_size": source_batch_size,
        "bucket_count": bucket_count,
        "batches": batches,
    }
    payload["plan_hash"] = _symmetric_edge_plan_hash(payload)
    plan_path = _symmetric_edge_plan_path(data_dir, entity_type)
    if plan_path.is_file():
        try:
            if json.loads(plan_path.read_text()) == payload:
                return payload
        except (OSError, ValueError, json.JSONDecodeError):
            pass
    plan_path.parent.mkdir(exist_ok=True, parents=True)
    temporary = plan_path.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(plan_path)
    return payload


def load_symmetric_edge_plan(
    data_dir: Path,
    *,
    entity_type: ClusterEntity = "ligand",
    metrics: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Load and validate the exact raw-score inputs for symmetrization."""
    plan_path = _symmetric_edge_plan_path(data_dir, entity_type)
    if not plan_path.is_file():
        raise FileNotFoundError(
            f"symmetric-edge planning must run before clustering: {plan_path}"
        )
    plan = cast(dict[str, Any], json.loads(plan_path.read_text()))
    payload = {key: value for key, value in plan.items() if key != "plan_hash"}
    if (
        plan.get("version") != 1
        or plan.get("entity_type", "ligand") != entity_type
        or plan.get("plan_hash") != (_symmetric_edge_plan_hash(payload))
    ):
        raise ValueError(f"invalid symmetric-edge plan: {plan_path}")
    selected_metrics = set(metrics or plan["metrics"])
    unknown_metrics = selected_metrics.difference(plan["metrics"])
    if unknown_metrics:
        raise ValueError(
            f"metrics are not in the symmetric-edge plan: {sorted(unknown_metrics)}"
        )
    for batch in plan["batches"]:
        if selected_metrics.isdisjoint(batch["metrics"]):
            continue
        for signature in batch["sources"]:
            source = Path(str(signature["path"]))
            if not source.is_file() or _component_source_signature(source) != signature:
                raise ValueError(
                    f"raw score source changed after symmetric-edge planning: {source}"
                )
    return plan


def _symmetric_fragment_dir(
    *,
    data_dir: Path,
    batch_key: str,
    entity_type: ClusterEntity = "ligand",
) -> Path:
    return (
        _cluster_root(data_dir, entity_type)
        / "symmetric_edges"
        / "fragments"
        / f"batch={batch_key}"
    )


def _completed_symmetric_fragment_batch(
    *,
    data_dir: Path,
    plan: Mapping[str, Any],
    batch: Mapping[str, Any],
    output_prefix: str | None = None,
    entity_type: ClusterEntity = "ligand",
) -> dict[str, Any] | None:
    output_dir = _symmetric_fragment_dir(
        data_dir=data_dir,
        batch_key=str(batch["key"]),
        entity_type=entity_type,
    )
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("version") != 1:
            return None
        if manifest.get("entity_type", "ligand") != entity_type:
            return None
        if manifest.get("plan_hash") != plan["plan_hash"]:
            return None
        if manifest.get("batch") != batch:
            return None
        outputs = manifest["outputs"]
        if output_prefix is not None:
            outputs = [
                output
                for output in outputs
                if str(output["path"]).startswith(output_prefix)
            ]
        for output in outputs:
            path = output_dir / str(output["path"])
            if not path.is_file() or path.stat().st_size != int(output["size"]):
                return None
        return cast(dict[str, Any], manifest)
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def symmetric_edge_fragment_batch_is_complete(
    *,
    data_dir: Path,
    batch: Mapping[str, Any],
    entity_type: ClusterEntity = "ligand",
) -> bool:
    """Return whether a planned fragment batch is complete and current."""
    plan = load_symmetric_edge_plan(data_dir, entity_type=entity_type)
    planned_batch = next(
        (value for value in plan["batches"] if str(value["key"]) == str(batch["key"])),
        None,
    )
    return (
        planned_batch == batch
        and _completed_symmetric_fragment_batch(
            data_dir=data_dir,
            plan=plan,
            batch=batch,
            entity_type=entity_type,
        )
        is not None
    )


def write_symmetric_edge_fragment_batch(
    *,
    data_dir: Path,
    batch: Mapping[str, Any],
    scratch_dir: Path,
    threads: int = 1,
    force_update: bool = False,
    read_paths: Sequence[Path] | None = None,
    entity_type: ClusterEntity = "ligand",
) -> dict[str, Any]:
    """Aggregate one raw-source batch into canonical-pair hash fragments."""
    if threads < 1:
        raise ValueError("symmetric-edge threads must be positive")
    plan = load_symmetric_edge_plan(data_dir, entity_type=entity_type)
    planned_batch = next(
        (value for value in plan["batches"] if str(value["key"]) == str(batch["key"])),
        None,
    )
    if planned_batch != batch:
        raise ValueError(
            f"batch is not part of the current symmetric-edge plan: {batch}"
        )
    if not force_update:
        completed = _completed_symmetric_fragment_batch(
            data_dir=data_dir,
            plan=plan,
            batch=batch,
            entity_type=entity_type,
        )
        if completed is not None:
            return completed

    source_paths = [Path(str(value["path"])) for value in batch["sources"]]
    scan_paths = list(read_paths or source_paths)
    if len(scan_paths) != len(source_paths):
        raise ValueError("symmetric-edge read paths do not match planned sources")
    scratch_output = scratch_dir / f"symmetric-fragments-{batch['key']}"
    rmtree(scratch_output, ignore_errors=True)
    scratch_output.mkdir(exist_ok=True, parents=True)

    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    temporary_root = scratch_dir / "duckdb"
    temporary_root.mkdir(exist_ok=True, parents=True)
    connection.sql(f"SET temp_directory='{temporary_root.as_posix()}'")
    paths_sql = ", ".join(
        f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'" for path in scan_paths
    )
    metrics = [str(metric) for metric in batch["metrics"]]
    if batch["kind"] == "chemical":
        metric = metrics[0]
        escaped_metric = metric.replace('"', '""')
        selected_sql = dedent(
            f"""
            SELECT
                '{metric.replace("'", "''")}'::VARCHAR AS metric,
                cast(query_ligand_id AS VARCHAR) AS query_node,
                cast(target_ligand_id AS VARCHAR) AS target_node,
                cast("{escaped_metric}" AS DOUBLE) AS similarity
            FROM read_parquet([{paths_sql}], union_by_name=true)
            """
        )
    elif batch["kind"] == "interface":
        metrics_sql = ", ".join(
            f"'{metric.replace(chr(39), chr(39) * 2)}'" for metric in metrics
        )
        selected_sql = dedent(
            f"""
            SELECT
                cast(metric AS VARCHAR) AS metric,
                cast(query_system AS VARCHAR) AS query_node,
                cast(target_system AS VARCHAR) AS target_node,
                cast(similarity AS DOUBLE) AS similarity
            FROM read_parquet([{paths_sql}], union_by_name=true)
            WHERE cast(metric AS VARCHAR) IN ({metrics_sql})
            """
        )
    else:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metrics[0],
            entity_type=entity_type,
        )
        connection.register(
            "eligible_systems",
            pd.DataFrame({"system_id": sorted(eligible_systems or set())}),
        )
        connection.register(
            "eligible_ligands",
            pd.DataFrame({"ligand_id": sorted(nodes)}),
        )
        metrics_sql = ", ".join(
            f"'{metric.replace(chr(39), chr(39) * 2)}'" for metric in metrics
        )
        selected_sql = dedent(
            f"""
            SELECT
                cast(scores.metric AS VARCHAR) AS metric,
                cast(scores.query_ligand_id AS VARCHAR) AS query_node,
                cast(scores.target_ligand_id AS VARCHAR) AS target_node,
                cast(scores.similarity AS DOUBLE) AS similarity
            FROM read_parquet([{paths_sql}], union_by_name=true) AS scores
            INNER JOIN eligible_systems AS query_system
              ON cast(scores.query_system AS VARCHAR) = query_system.system_id
            INNER JOIN eligible_systems AS target_system
              ON cast(scores.target_system AS VARCHAR) = target_system.system_id
            INNER JOIN eligible_ligands AS query_ligand
              ON cast(scores.query_ligand_id AS VARCHAR) = query_ligand.ligand_id
            INNER JOIN eligible_ligands AS target_ligand
              ON cast(scores.target_ligand_id AS VARCHAR) = target_ligand.ligand_id
            WHERE cast(scores.metric AS VARCHAR) IN ({metrics_sql})
            """
        )
    if batch["kind"] == "chemical":
        aggregate_sql = """
            max(similarity)::DOUBLE AS forward_similarity,
            max(similarity)::DOUBLE AS reverse_similarity
        """
    else:
        aggregate_sql = """
            max(similarity) FILTER (WHERE query_node = low_node)::DOUBLE
                AS forward_similarity,
            max(similarity) FILTER (WHERE query_node = high_node)::DOUBLE
                AS reverse_similarity
        """
    query = dedent(
        f"""
        COPY (
            WITH selected AS (
                {selected_sql}
            ), canonical AS (
                SELECT
                    metric,
                    query_node,
                    target_node,
                    least(query_node, target_node) AS low_node,
                    greatest(query_node, target_node) AS high_node,
                    similarity
                FROM selected
                WHERE query_node != target_node
                  AND similarity IS NOT NULL
            ), bucketed AS (
                SELECT
                    *,
                    (
                        hash(low_node, high_node)
                        % {int(plan["bucket_count"])}
                    )::INTEGER AS bucket
                FROM canonical
            )
            SELECT
                metric,
                bucket,
                low_node,
                high_node,
                {aggregate_sql}
            FROM bucketed
            GROUP BY metric, bucket, low_node, high_node
        ) TO '{scratch_output.as_posix()}'
        (
            FORMAT PARQUET,
            COMPRESSION ZSTD,
            ROW_GROUP_SIZE 500000,
            PARTITION_BY (metric, bucket)
        )
        """
    )
    started = time()
    LOG.info(
        "symmetric fragment start: batch=%s kind=%s sources=%d metrics=%s",
        batch["key"],
        batch["kind"],
        len(source_paths),
        metrics,
    )
    connection.sql(query)
    connection.close()
    output_files = sorted(scratch_output.rglob("*.parquet"))
    manifest: dict[str, Any] = {
        "version": 1,
        "entity_type": entity_type,
        "status": "complete",
        "plan_hash": plan["plan_hash"],
        "batch": dict(batch),
        "outputs": [
            {
                "path": path.relative_to(scratch_output).as_posix(),
                "size": path.stat().st_size,
            }
            for path in output_files
        ],
        "elapsed_seconds": time() - started,
    }
    output_dir = _symmetric_fragment_dir(
        data_dir=data_dir,
        batch_key=str(batch["key"]),
        entity_type=entity_type,
    )
    temporary_output = output_dir.with_name(f"{output_dir.name}.tmp")
    rmtree(temporary_output, ignore_errors=True)
    copytree(scratch_output, temporary_output)
    rmtree(output_dir, ignore_errors=True)
    temporary_output.replace(output_dir)
    temporary_manifest = output_dir / "manifest.tmp.json"
    temporary_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary_manifest.replace(output_dir / "manifest.json")
    LOG.info(
        "symmetric fragment complete: batch=%s files=%d elapsed_seconds=%.1f",
        batch["key"],
        len(output_files),
        time() - started,
    )
    return manifest


def _symmetric_edge_path(
    *,
    data_dir: Path,
    metric: str,
    bucket: int,
    entity_type: ClusterEntity = "ligand",
) -> Path:
    return (
        _cluster_root(data_dir, entity_type)
        / "symmetric_edges"
        / f"metric={metric}"
        / f"bucket={bucket:03d}.parquet"
    )


def write_symmetric_edge_shard(
    *,
    data_dir: Path,
    metric: str,
    bucket: int,
    scratch_dir: Path,
    threads: int = 1,
    force_update: bool = False,
    entity_type: ClusterEntity = "ligand",
) -> dict[str, Any]:
    """Merge directional maxima into one reciprocal-minimum edge shard."""
    if threads < 1:
        raise ValueError("symmetric-edge threads must be positive")
    plan = load_symmetric_edge_plan(data_dir, entity_type=entity_type)
    if metric not in plan["metrics"]:
        raise ValueError(f"metric is not in the symmetric-edge plan: {metric}")
    if bucket < 0 or bucket >= int(plan["bucket_count"]):
        raise ValueError(f"symmetric-edge bucket is out of range: {bucket}")
    fragments: list[Path] = []
    for batch in plan["batches"]:
        if metric not in batch["metrics"]:
            continue
        completed = _completed_symmetric_fragment_batch(
            data_dir=data_dir,
            plan=plan,
            batch=batch,
            output_prefix=f"metric={metric}/bucket={bucket}/",
            entity_type=entity_type,
        )
        if completed is None:
            raise ValueError(
                f"symmetric-edge fragment batch is incomplete: {batch['key']}"
            )
        output_dir = _symmetric_fragment_dir(
            data_dir=data_dir,
            batch_key=str(batch["key"]),
            entity_type=entity_type,
        )
        fragments.extend(output_dir.glob(f"metric={metric}/bucket={bucket}/*.parquet"))
    fragments = sorted(fragments)
    output = _symmetric_edge_path(
        data_dir=data_dir,
        metric=metric,
        bucket=bucket,
        entity_type=entity_type,
    )
    manifest_path = output.with_suffix(".json")
    source_signatures = [_component_source_signature(path) for path in fragments]
    if not force_update and output.is_file() and manifest_path.is_file():
        try:
            cached_manifest = json.loads(manifest_path.read_text())
            if (
                cached_manifest.get("version") == 1
                and cached_manifest.get("plan_hash") == plan["plan_hash"]
                and cached_manifest.get("fragments") == source_signatures
                and output.stat().st_size == int(cached_manifest["size"])
            ):
                return cast(dict[str, Any], cached_manifest)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass

    import duckdb

    scratch_dir.mkdir(exist_ok=True, parents=True)
    local_output = scratch_dir / f"{metric}-{bucket:03d}.parquet"
    local_output.unlink(missing_ok=True)
    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
    if fragments:
        local_fragment_dir = scratch_dir / "fragments"
        rmtree(local_fragment_dir, ignore_errors=True)
        local_fragment_dir.mkdir(exist_ok=True, parents=True)
        scan_fragments: list[Path] = []
        for index, fragment in enumerate(fragments):
            local_fragment = local_fragment_dir / f"{index:05d}.parquet"
            copyfile(fragment, local_fragment)
            scan_fragments.append(local_fragment)
        paths_sql = ", ".join(
            f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'"
            for path in scan_fragments
        )
        selected_sql = dedent(
            f"""
            WITH combined AS (
                SELECT
                    cast(low_node AS VARCHAR) AS query_node,
                    cast(high_node AS VARCHAR) AS target_node,
                    max(forward_similarity)::DOUBLE AS forward_similarity,
                    max(reverse_similarity)::DOUBLE AS reverse_similarity
                FROM read_parquet([{paths_sql}], union_by_name=true)
                GROUP BY query_node, target_node
            )
            SELECT
                query_node,
                target_node,
                forward_similarity,
                reverse_similarity,
                CASE
                    WHEN forward_similarity IS NOT NULL
                     AND reverse_similarity IS NOT NULL
                    THEN least(forward_similarity, reverse_similarity)
                END::DOUBLE AS similarity,
                CASE
                    WHEN forward_similarity IS NULL THEN reverse_similarity
                    WHEN reverse_similarity IS NULL THEN forward_similarity
                    ELSE greatest(forward_similarity, reverse_similarity)
                END::DOUBLE AS maximum_similarity
            FROM combined
            WHERE forward_similarity IS NOT NULL
               OR reverse_similarity IS NOT NULL
            """
        )
    else:
        selected_sql = """
            SELECT
                NULL::VARCHAR AS query_node,
                NULL::VARCHAR AS target_node,
                NULL::DOUBLE AS forward_similarity,
                NULL::DOUBLE AS reverse_similarity,
                NULL::DOUBLE AS similarity,
                NULL::DOUBLE AS maximum_similarity
            WHERE false
        """
    started = time()
    connection.sql(
        dedent(
            f"""
            COPY ({selected_sql})
            TO '{local_output.as_posix()}'
            (
                FORMAT PARQUET,
                COMPRESSION ZSTD,
                ROW_GROUP_SIZE 500000
            )
            """
        )
    )
    row_count = int(
        connection.sql(
            f"SELECT count(*) FROM read_parquet('{local_output.as_posix()}')"
        ).fetchone()[0]
    )
    connection.close()
    output.parent.mkdir(exist_ok=True, parents=True)
    temporary_output = output.with_suffix(".tmp.parquet")
    copyfile(local_output, temporary_output)
    temporary_output.replace(output)
    manifest: dict[str, Any] = {
        "version": 1,
        "status": "complete",
        "plan_hash": plan["plan_hash"],
        "metric": metric,
        "bucket": bucket,
        "fragments": source_signatures,
        "rows": row_count,
        "size": output.stat().st_size,
        "elapsed_seconds": time() - started,
    }
    temporary_manifest = manifest_path.with_suffix(".tmp.json")
    temporary_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary_manifest.replace(manifest_path)
    LOG.info(
        "symmetric edge shard complete: metric=%s bucket=%d rows=%d "
        "fragments=%d elapsed_seconds=%.1f",
        metric,
        bucket,
        row_count,
        len(fragments),
        time() - started,
    )
    return manifest


def _completed_component_reduction(
    *,
    output_dir: Path,
    source_signature: Mapping[str, str | int],
    node_hash: str,
    selection_hash: str,
    thresholds: Sequence[int],
) -> dict[str, Any] | None:
    manifest_path = output_dir / "reduction.json"
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("version") != COMPONENT_REDUCTION_VERSION:
            return None
        if manifest.get("directed_modes") != list(COMPONENT_REDUCTION_DIRECTIONS):
            return None
        if manifest.get("source") != source_signature:
            return None
        if manifest.get("node_hash") != node_hash:
            return None
        if manifest.get("selection_hash") != selection_hash:
            return None
        if manifest.get("thresholds") != list(thresholds):
            return None
        outputs = manifest["outputs"]
        for output in outputs:
            path = output_dir / str(output["path"])
            if not path.is_file() or path.stat().st_size != int(output["size"]):
                return None
            schema = pq.read_schema(path)
            if not set(COMPONENT_EDGE_COLUMNS).issubset(schema.names):
                return None
        return cast(dict[str, Any], manifest)
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def write_component_reduction_shard(
    *,
    edge_batches: Iterable[pd.DataFrame],
    all_nodes: Iterable[str],
    thresholds: Sequence[int],
    source_path: Path,
    output_dir: Path,
    query_column: str = "query_node",
    target_column: str = "target_node",
    similarity_column: str = "similarity",
    forest_fan_in: int = 8,
    force_update: bool = False,
    selection_hash: str | None = None,
) -> dict[str, Any]:
    """Persist resumable directed-component reductions for one score shard."""
    ordered_thresholds = sorted(set(thresholds), reverse=True)
    nodes = sorted(set(map(str, all_nodes)))
    source_signature = _component_source_signature(source_path)
    node_hash = _component_node_hash(nodes)
    effective_selection_hash = selection_hash or node_hash
    if not force_update:
        completed = _completed_component_reduction(
            output_dir=output_dir,
            source_signature=source_signature,
            node_hash=node_hash,
            selection_hash=effective_selection_hash,
            thresholds=ordered_thresholds,
        )
        if completed is not None:
            return completed

    _, reductions = reduce_threshold_component_edges(
        edge_batches=edge_batches,
        all_nodes=nodes,
        thresholds=ordered_thresholds,
        query_column=query_column,
        target_column=target_column,
        similarity_column=similarity_column,
        forest_fan_in=forest_fan_in,
        directed_modes=COMPONENT_REDUCTION_DIRECTIONS,
    )
    outputs: list[dict[str, str | int | bool]] = []
    for directed in COMPONENT_REDUCTION_DIRECTIONS:
        for threshold in ordered_thresholds:
            relative = (
                Path(f"directed={str(directed).lower()}")
                / f"threshold={threshold}.parquet"
            )
            output = output_dir / relative
            output.parent.mkdir(exist_ok=True, parents=True)
            temporary = output.with_suffix(".tmp.parquet")
            reductions[directed][threshold].to_parquet(
                temporary,
                index=False,
                compression="zstd",
            )
            temporary.replace(output)
            outputs.append(
                {
                    "directed": directed,
                    "threshold": threshold,
                    "path": relative.as_posix(),
                    "rows": len(reductions[directed][threshold]),
                    "size": output.stat().st_size,
                }
            )
    manifest: dict[str, Any] = {
        "version": COMPONENT_REDUCTION_VERSION,
        "status": "complete",
        "directed_modes": list(COMPONENT_REDUCTION_DIRECTIONS),
        "source": source_signature,
        "node_hash": node_hash,
        "selection_hash": effective_selection_hash,
        "node_count": len(nodes),
        "thresholds": ordered_thresholds,
        "outputs": outputs,
    }
    output_dir.mkdir(exist_ok=True, parents=True)
    manifest_path = output_dir / "reduction.json"
    temporary_manifest = manifest_path.with_suffix(".tmp.json")
    temporary_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary_manifest.replace(manifest_path)
    return manifest


def merge_component_reduction_shards(
    *,
    reduction_dirs: Sequence[Path],
    all_nodes: Iterable[str],
    thresholds: Sequence[int],
    forest_fan_in: int = 8,
    parquet_batch_size: int = 250_000,
    selection_hash: str | None = None,
) -> dict[bool, dict[int, pd.DataFrame]]:
    """Merge completed shard reductions into exact directed-component labels."""
    if parquet_batch_size < 1:
        raise ValueError("component Parquet batch size must be positive")
    ordered_thresholds = sorted(set(thresholds), reverse=True)
    nodes = sorted(set(map(str, all_nodes)))
    node_hash = _component_node_hash(nodes)
    effective_selection_hash = selection_hash or node_hash
    reduction_count = len(reduction_dirs)
    started = time()
    LOG.info(
        "component merge: loading %d source reductions for %d nodes at thresholds=%s",
        reduction_count,
        len(nodes),
        ordered_thresholds,
    )
    accumulators = {
        directed: {
            threshold: _ComponentEdgeAccumulator(
                forest_fan_in,
                (
                    _reduce_directed_component_edges
                    if directed
                    else _reduce_undirected_component_edges
                ),
            )
            for threshold in ordered_thresholds
        }
        for directed in COMPONENT_REDUCTION_DIRECTIONS
    }
    for index, reduction_dir in enumerate(reduction_dirs, start=1):
        manifest_path = reduction_dir / "reduction.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        source = cast(dict[str, str | int], manifest["source"])
        source_path = Path(str(source["path"]))
        completed = _completed_component_reduction(
            output_dir=reduction_dir,
            source_signature=_component_source_signature(source_path),
            node_hash=node_hash,
            selection_hash=effective_selection_hash,
            thresholds=ordered_thresholds,
        )
        if completed is None:
            raise ValueError(
                f"stale or incomplete component reduction: {reduction_dir}"
            )
        output_by_key = {
            (bool(output["directed"]), int(output["threshold"])): output
            for output in completed["outputs"]
        }
        for directed in COMPONENT_REDUCTION_DIRECTIONS:
            for threshold in ordered_thresholds:
                output = output_by_key[(directed, threshold)]
                path = reduction_dir / str(output["path"])
                parquet = pq.ParquetFile(path)
                for batch in parquet.iter_batches(
                    batch_size=parquet_batch_size,
                    columns=COMPONENT_EDGE_COLUMNS,
                ):
                    accumulators[directed][threshold].add(batch.to_pandas())
        if index % 25 == 0 or index == reduction_count:
            elapsed = time() - started
            rate = index / elapsed
            LOG.info(
                "component merge source progress: loaded=%d/%d rate=%.1f/s "
                "eta_seconds=%.1f",
                index,
                reduction_count,
                rate,
                (reduction_count - index) / rate,
            )

    result: dict[bool, dict[int, pd.DataFrame]] = {}
    for directed in COMPONENT_REDUCTION_DIRECTIONS:
        reductions: dict[int, pd.DataFrame] = {}
        for threshold in ordered_thresholds:
            phase_started = time()
            LOG.info(
                "component merge: reducing directed=%s threshold=%d",
                directed,
                threshold,
            )
            reductions[threshold] = accumulators[directed][threshold].finish()
            LOG.info(
                "component merge: reduced directed=%s threshold=%d edges=%d "
                "elapsed_seconds=%.1f",
                directed,
                threshold,
                len(reductions[threshold]),
                time() - phase_started,
            )
        phase_started = time()
        result[directed] = _labels_from_threshold_reductions(
            nodes=nodes,
            reductions=reductions,
            directed=directed,
        )
        LOG.info(
            "component merge: labeled directed=%s thresholds=%d "
            "elapsed_seconds=%.1f",
            directed,
            len(ordered_thresholds),
            time() - phase_started,
        )
    LOG.info("component merge complete: elapsed_seconds=%.1f", time() - started)
    return result


def component_score_sources(
    *,
    data_dir: Path,
    metric: str,
    entity_type: ClusterEntity = "ligand",
) -> list[Path]:
    """Return compact reciprocal-minimum shards for one clustering metric."""
    plan = load_symmetric_edge_plan(
        data_dir,
        entity_type=entity_type,
        metrics=[metric],
    )
    if metric not in plan["metrics"]:
        raise ValueError(f"metric is not in the symmetric-edge plan: {metric}")
    sources = [
        _symmetric_edge_path(
            data_dir=data_dir,
            metric=metric,
            bucket=bucket,
            entity_type=entity_type,
        )
        for bucket in range(int(plan["bucket_count"]))
    ]
    incomplete: list[Path] = []
    required_columns = {
        "query_node",
        "target_node",
        "forward_similarity",
        "reverse_similarity",
        "similarity",
        "maximum_similarity",
    }
    for bucket, path in enumerate(sources):
        manifest_path = path.with_suffix(".json")
        try:
            manifest = json.loads(manifest_path.read_text())
            valid = (
                path.is_file()
                and manifest.get("version") == 1
                and manifest.get("status") == "complete"
                and manifest.get("plan_hash") == plan["plan_hash"]
                and manifest.get("metric") == metric
                and int(manifest.get("bucket", -1)) == bucket
                and int(manifest.get("size", -1)) == path.stat().st_size
                and required_columns.issubset(pq.read_schema(path).names)
            )
        except (OSError, TypeError, ValueError, KeyError, json.JSONDecodeError):
            valid = False
        if not valid:
            incomplete.append(path)
    if incomplete:
        raise FileNotFoundError(
            f"symmetric-edge shards are incomplete for {metric}: "
            f"{[path.name for path in incomplete[:10]]}"
        )
    return sources


def iter_component_score_batches(
    *,
    source_path: Path,
    metric: str,
    batch_size: int = 250_000,
    eligible_systems: set[str] | None = None,
) -> Iterable[pd.DataFrame]:
    """Read one compact reciprocal-minimum ligand-edge shard in bounded batches."""
    if batch_size < 1:
        raise ValueError("component score batch size must be positive")
    parquet = pq.ParquetFile(source_path)
    columns = SYMMETRIC_EDGE_COLUMNS
    missing = sorted(set(columns).difference(parquet.schema.names))
    if missing:
        raise ValueError(f"component score shard {source_path} is missing {missing}")
    for batch in parquet.iter_batches(
        batch_size=batch_size,
        columns=columns,
    ):
        frame = batch.to_pandas()
        if not frame.empty:
            yield frame[SYMMETRIC_EDGE_COLUMNS]


def _cached_component_node_universe(
    data_dir: Path,
    *,
    entity_type: ClusterEntity = "ligand",
) -> tuple[list[str], set[str] | None] | None:
    cache_path, manifest_path = _component_node_universe_paths(data_dir, entity_type)
    if not cache_path.is_file() or not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("version") != 3 or not manifest.get("universe_hash"):
            return None
        if manifest.get("entity_type", "ligand") != entity_type:
            return None
        if manifest.get("sources") != _component_node_universe_sources(
            data_dir, entity_type
        ):
            return None
        if cache_path.stat().st_size != int(manifest["cache_size"]):
            return None
        if set(pq.read_schema(cache_path).names) != {"kind", "id"}:
            return None
        frame = pd.read_parquet(cache_path, columns=["kind", "id"])
        if entity_type == "ligand":
            nodes = frame.loc[frame["kind"].eq("ligand"), "id"].astype(str).tolist()
            ligand_systems = set(
                frame.loc[frame["kind"].eq("system"), "id"].astype(str).tolist()
            )
            if len(nodes) != int(manifest["ligand_count"]):
                return None
            if len(ligand_systems) != int(manifest["system_count"]):
                return None
            systems: set[str] | None = ligand_systems
        else:
            nodes = frame.loc[frame["kind"].eq("interface"), "id"].astype(str).tolist()
            side_nodes = (
                frame.loc[frame["kind"].eq("interface_side"), "id"].astype(str).tolist()
            )
            systems = None
            if len(nodes) != int(manifest["interface_count"]):
                return None
            if len(side_nodes) != int(manifest["interface_side_count"]):
                return None
            nodes.extend(side_nodes)
        return nodes, systems
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def prepare_component_node_universe(
    data_dir: Path, *, entity_type: ClusterEntity = "ligand"
) -> dict[str, Any]:
    """Materialize the shared clustering-node universe before scattering."""
    cached = _cached_component_node_universe(data_dir, entity_type=entity_type)
    if cached is not None:
        nodes, systems = cached
        if entity_type == "interface":
            interface_side_count = sum("::side=" in node for node in nodes)
            return {
                "status": "cached",
                "interface_count": len(nodes) - interface_side_count,
                "interface_side_count": interface_side_count,
            }
        return {
            "status": "cached",
            "ligand_count": len(nodes),
            "system_count": len(systems or set()),
        }

    if entity_type == "interface":
        representative_path = data_dir / INTERFACE_REPRESENTATIVES
        half_representative_path = data_dir / INTERFACE_HALF_REPRESENTATIVES
        nodes = sorted(
            set(
                pd.read_parquet(
                    representative_path,
                    columns=["representative_system_id"],
                )["representative_system_id"]
                .dropna()
                .astype(str)
            )
        )
        side_nodes = sorted(
            set(
                pd.read_parquet(
                    half_representative_path,
                    columns=["half_interface_id"],
                )["half_interface_id"]
                .dropna()
                .astype(str)
            )
        )
        system_ids: list[str] = []
        frame = pd.concat(
            [
                pd.DataFrame({"kind": "interface", "id": nodes}),
                pd.DataFrame({"kind": "interface_side", "id": side_nodes}),
            ],
            ignore_index=True,
        )
        counts = {
            "interface_count": len(nodes),
            "interface_side_count": len(side_nodes),
        }
    else:
        annotation = _eligible_annotation(data_dir)
        nodes = sorted(
            set(
                annotation.loc[annotation["ligand_is_proper"], "ligand_id"]
                .dropna()
                .astype(str)
            )
        )
        system_ids = sorted(set(annotation["system_id"].dropna().astype(str)))
        frame = pd.concat(
            [
                pd.DataFrame({"kind": "ligand", "id": nodes}),
                pd.DataFrame({"kind": "system", "id": system_ids}),
            ],
            ignore_index=True,
        )
        counts = {
            "ligand_count": len(nodes),
            "system_count": len(system_ids),
        }
    cache_path, manifest_path = _component_node_universe_paths(data_dir, entity_type)
    cache_path.parent.mkdir(exist_ok=True, parents=True)
    temporary_cache = cache_path.with_suffix(".tmp.parquet")
    frame.to_parquet(temporary_cache, index=False, compression="zstd")
    temporary_cache.replace(cache_path)
    universe_digest = hashlib.sha256()
    for kind, node_id in frame[["kind", "id"]].itertuples(index=False, name=None):
        universe_digest.update(str(kind).encode())
        universe_digest.update(b"\0")
        universe_digest.update(str(node_id).encode())
        universe_digest.update(b"\n")
    manifest = {
        "version": 3,
        "entity_type": entity_type,
        "sources": _component_node_universe_sources(data_dir, entity_type),
        "cache_size": cache_path.stat().st_size,
        "universe_hash": universe_digest.hexdigest(),
        **counts,
    }
    temporary_manifest = manifest_path.with_suffix(".tmp.json")
    temporary_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary_manifest.replace(manifest_path)
    return {"status": "complete", **manifest}


def component_node_universe(
    *,
    data_dir: Path,
    metric: str,
    entity_type: ClusterEntity = "ligand",
) -> tuple[list[str], set[str] | None]:
    if entity_type == "interface":
        if metric not in INTERFACE_CLUSTER_METRICS:
            raise ValueError(f"unsupported interface clustering metric: {metric}")
        cached = _cached_component_node_universe(data_dir, entity_type=entity_type)
        if cached is not None:
            cached_nodes, systems = cached
            is_side = metric == "interface_side_qcov"
            side_nodes = [node for node in cached_nodes if "::side=" in node]
            return (
                side_nodes
                if is_side
                else [node for node in cached_nodes if "::side=" not in node],
                systems,
            )
        if metric == "interface_side_qcov":
            source = data_dir / INTERFACE_HALF_REPRESENTATIVES
            column = "half_interface_id"
        else:
            source = data_dir / INTERFACE_REPRESENTATIVES
            column = "representative_system_id"
        nodes = sorted(
            set(pd.read_parquet(source, columns=[column])[column].dropna().astype(str))
        )
        return nodes, None
    if metric == "tanimoto_similarity_ecfp4_1024":
        nodes = (
            pd.read_parquet(
                data_dir / "fingerprints" / "ligands_per_smiles.parquet",
                columns=["ligand_smiles_id"],
            )["ligand_smiles_id"]
            .astype(str)
            .tolist()
        )
        return sorted(set(nodes)), None
    cached = _cached_component_node_universe(data_dir, entity_type=entity_type)
    if cached is not None:
        return cached
    annotation = _eligible_annotation(data_dir)
    systems = set(annotation["system_id"].dropna().astype(str))
    node_set = set(
        annotation.loc[annotation["ligand_is_proper"], "ligand_id"].dropna().astype(str)
    )
    return sorted(node_set), systems


def component_reduction_dir(
    *,
    data_dir: Path,
    metric: str,
    source_path: Path,
    entity_type: ClusterEntity = "ligand",
) -> Path:
    """Return a stable local reduction path without embedding an absolute root."""
    try:
        source_name = source_path.resolve().relative_to(data_dir.resolve()).as_posix()
    except ValueError:
        source_name = source_path.name
    source_key = hashlib.sha256(source_name.encode()).hexdigest()[:16]
    return (
        _cluster_root(data_dir, entity_type)
        / "reductions"
        / f"metric={metric}"
        / f"source={source_key}"
    )


def directed_cover_reduction_dir(
    *,
    data_dir: Path,
    metric: str,
    source_path: Path,
    entity_type: ClusterEntity = "ligand",
) -> Path:
    """Return the reduction path for any-direction cover connectivity."""
    try:
        source_name = source_path.resolve().relative_to(data_dir.resolve()).as_posix()
    except ValueError:
        source_name = source_path.name
    source_key = hashlib.sha256(source_name.encode()).hexdigest()[:16]
    return (
        _sampling_root(data_dir, entity_type)
        / "directed_set_cover"
        / "reductions"
        / f"metric={metric}"
        / f"source={source_key}"
    )


def _directed_cover_connectivity_thresholds(
    thresholds: Sequence[int],
) -> list[int]:
    """Include the threshold-50 graph needed by stricter cover passes."""
    result = set(thresholds)
    if any(threshold > 50 for threshold in result):
        result.add(50)
    return sorted(result, reverse=True)


def _iter_directed_cover_edge_batches(
    *, source_path: Path, batch_size: int = 250_000
) -> Iterable[pd.DataFrame]:
    """Read maximum-direction edges for directed-cover weak connectivity."""
    parquet = pq.ParquetFile(source_path)
    columns = ["query_node", "target_node", "maximum_similarity"]
    missing = sorted(set(columns).difference(parquet.schema.names))
    if missing:
        raise ValueError(
            f"directed-cover score shard {source_path} is missing {missing}"
        )
    for batch in parquet.iter_batches(batch_size=batch_size, columns=columns):
        frame = batch.to_pandas().rename(columns={"maximum_similarity": "similarity"})
        if not frame.empty:
            yield frame[SYMMETRIC_EDGE_COLUMNS]


def make_directed_cover_component_reduction(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    source_path: Path,
    batch_size: int = 250_000,
    forest_fan_in: int = 8,
    force_update: bool = False,
    read_path: Path | None = None,
    all_nodes: Sequence[str] | None = None,
    eligible_systems: set[str] | None = None,
    entity_type: ClusterEntity = "ligand",
) -> dict[str, Any]:
    """Reduce one shard by edges present in either score direction."""
    if all_nodes is None:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
            entity_type=entity_type,
        )
    else:
        nodes = sorted(set(map(str, all_nodes)))
    selection_hash = _component_selection_hash(
        nodes=nodes,
        eligible_systems=eligible_systems,
    )
    return write_component_reduction_shard(
        edge_batches=_iter_directed_cover_edge_batches(
            source_path=read_path or source_path,
            batch_size=batch_size,
        ),
        all_nodes=nodes,
        thresholds=_directed_cover_connectivity_thresholds(thresholds),
        source_path=source_path,
        output_dir=directed_cover_reduction_dir(
            data_dir=data_dir,
            metric=metric,
            source_path=source_path,
            entity_type=entity_type,
        ),
        forest_fan_in=forest_fan_in,
        force_update=force_update,
        selection_hash=selection_hash,
    )


def merge_directed_cover_component_reductions(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    forest_fan_in: int = 8,
    parquet_batch_size: int = 250_000,
    entity_type: ClusterEntity = "ligand",
) -> dict[int, pd.DataFrame]:
    """Merge exact weak components used to bound directed set cover."""
    thresholds = _directed_cover_connectivity_thresholds(thresholds)
    sources = component_score_sources(
        data_dir=data_dir,
        metric=metric,
        entity_type=entity_type,
    )
    nodes, eligible_systems = component_node_universe(
        data_dir=data_dir,
        metric=metric,
        entity_type=entity_type,
    )
    selection_hash = _component_selection_hash(
        nodes=nodes,
        eligible_systems=eligible_systems,
    )
    labels = merge_component_reduction_shards(
        reduction_dirs=[
            directed_cover_reduction_dir(
                data_dir=data_dir,
                metric=metric,
                source_path=source,
                entity_type=entity_type,
            )
            for source in sources
        ],
        all_nodes=nodes,
        thresholds=thresholds,
        forest_fan_in=forest_fan_in,
        parquet_batch_size=parquet_batch_size,
        selection_hash=selection_hash,
    )[False]
    for threshold, label_frame in labels.items():
        output = (
            _sampling_root(data_dir, entity_type)
            / "directed_set_cover"
            / "reductions"
            / f"metric={metric}"
            / "labels"
            / f"threshold={threshold}.parquet"
        )
        output.parent.mkdir(exist_ok=True, parents=True)
        temporary = output.with_suffix(".tmp.parquet")
        label_frame.rename(
            columns={"ligand_id": _cluster_node_column(entity_type)}
        ).to_parquet(temporary, index=False, compression="zstd")
        temporary.replace(output)
    return labels


def score_component_reduction_is_complete(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    source_path: Path,
    all_nodes: Sequence[str] | None = None,
    eligible_systems: set[str] | None = None,
    entity_type: ClusterEntity = "ligand",
) -> bool:
    """Return whether one source/metric reduction matches all current inputs."""
    if all_nodes is None:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
            entity_type=entity_type,
        )
    else:
        nodes = sorted(set(map(str, all_nodes)))
    node_hash = _component_node_hash(nodes)
    selection_hash = _component_selection_hash(
        nodes=nodes,
        eligible_systems=eligible_systems,
    )
    return (
        _completed_component_reduction(
            output_dir=component_reduction_dir(
                data_dir=data_dir,
                metric=metric,
                source_path=source_path,
                entity_type=entity_type,
            ),
            source_signature=_component_source_signature(source_path),
            node_hash=node_hash,
            selection_hash=selection_hash,
            thresholds=sorted(set(thresholds), reverse=True),
        )
        is not None
    )


def directed_cover_component_reduction_is_complete(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    source_path: Path,
    all_nodes: Sequence[str] | None = None,
    eligible_systems: set[str] | None = None,
    entity_type: ClusterEntity = "ligand",
) -> bool:
    """Return whether one any-direction connectivity reduction is current."""
    if all_nodes is None:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
            entity_type=entity_type,
        )
    else:
        nodes = sorted(set(map(str, all_nodes)))
    node_hash = _component_node_hash(nodes)
    selection_hash = _component_selection_hash(
        nodes=nodes,
        eligible_systems=eligible_systems,
    )
    return (
        _completed_component_reduction(
            output_dir=directed_cover_reduction_dir(
                data_dir=data_dir,
                metric=metric,
                source_path=source_path,
                entity_type=entity_type,
            ),
            source_signature=_component_source_signature(source_path),
            node_hash=node_hash,
            selection_hash=selection_hash,
            thresholds=_directed_cover_connectivity_thresholds(thresholds),
        )
        is not None
    )


def make_score_component_reduction(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    source_path: Path,
    batch_size: int = 250_000,
    forest_fan_in: int = 8,
    force_update: bool = False,
    read_path: Path | None = None,
    all_nodes: Sequence[str] | None = None,
    eligible_systems: set[str] | None = None,
    entity_type: ClusterEntity = "ligand",
) -> dict[str, Any]:
    """Map one final score shard into resumable component reductions."""
    if all_nodes is None:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
            entity_type=entity_type,
        )
    else:
        nodes = sorted(set(map(str, all_nodes)))
    selection_hash = _component_selection_hash(
        nodes=nodes,
        eligible_systems=eligible_systems,
    )
    return write_component_reduction_shard(
        edge_batches=iter_component_score_batches(
            source_path=read_path or source_path,
            metric=metric,
            batch_size=batch_size,
            eligible_systems=eligible_systems,
        ),
        all_nodes=nodes,
        thresholds=thresholds,
        source_path=source_path,
        output_dir=component_reduction_dir(
            data_dir=data_dir,
            metric=metric,
            source_path=source_path,
            entity_type=entity_type,
        ),
        forest_fan_in=forest_fan_in,
        force_update=force_update,
        selection_hash=selection_hash,
    )


def merge_score_component_reductions(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    forest_fan_in: int = 8,
    parquet_batch_size: int = 250_000,
    entity_type: ClusterEntity = "ligand",
) -> dict[bool, dict[int, pd.DataFrame]]:
    """Reduce reciprocal-minimum shards into internal connectivity labels."""
    started = time()
    sources = component_score_sources(
        data_dir=data_dir,
        metric=metric,
        entity_type=entity_type,
    )
    if not sources:
        raise FileNotFoundError(f"no component score sources found for {metric}")
    nodes, eligible_systems = component_node_universe(
        data_dir=data_dir,
        metric=metric,
        entity_type=entity_type,
    )
    selection_hash = _component_selection_hash(
        nodes=nodes,
        eligible_systems=eligible_systems,
    )
    LOG.info(
        "component merge metric=%s: sources=%d nodes=%d thresholds=%s",
        metric,
        len(sources),
        len(nodes),
        sorted(set(thresholds), reverse=True),
    )
    labels = merge_component_reduction_shards(
        reduction_dirs=[
            component_reduction_dir(
                data_dir=data_dir,
                metric=metric,
                source_path=source,
                entity_type=entity_type,
            )
            for source in sources
        ],
        all_nodes=nodes,
        thresholds=thresholds,
        forest_fan_in=forest_fan_in,
        parquet_batch_size=parquet_batch_size,
        selection_hash=selection_hash,
    )
    for directed in COMPONENT_REDUCTION_DIRECTIONS:
        for threshold, label_frame in labels[directed].items():
            node_column = _cluster_node_column(entity_type)
            internal_frame = label_frame.rename(columns={"ligand_id": node_column})
            internal = (
                _cluster_root(data_dir, entity_type)
                / "reductions"
                / f"metric={metric}"
                / "labels"
                / f"directed={str(directed).lower()}"
                / f"threshold={threshold}.parquet"
            )
            internal.parent.mkdir(exist_ok=True, parents=True)
            temporary_internal = internal.with_suffix(".tmp.parquet")
            internal_frame.to_parquet(
                temporary_internal,
                index=False,
                compression="zstd",
            )
            temporary_internal.replace(internal)
    LOG.info(
        "component merge metric=%s complete: elapsed_seconds=%.1f",
        metric,
        time() - started,
    )
    return labels


def count_crossing_component_edges(
    *,
    edge_batches: Iterable[pd.DataFrame],
    labels_by_threshold: Mapping[int, pd.DataFrame],
    query_column: str = "query_node",
    target_column: str = "target_node",
    similarity_column: str = "similarity",
) -> dict[int, int]:
    """Count qualifying edges that cross generated component labels."""
    if not labels_by_threshold:
        raise ValueError("component labels are required for validation")
    label_maps: dict[int, dict[str, str]] = {}
    for threshold, labels in labels_by_threshold.items():
        missing = sorted({"ligand_id", "label"}.difference(labels.columns))
        if missing:
            raise ValueError(f"component labels are missing columns {missing}")
        if labels["ligand_id"].duplicated().any():
            raise ValueError(f"component labels at {threshold} contain duplicate nodes")
        label_maps[threshold] = dict(
            zip(labels["ligand_id"].astype(str), labels["label"].astype(str))
        )
    crossings = dict.fromkeys(label_maps, 0)
    required = {query_column, target_column, similarity_column}
    for batch in edge_batches:
        missing = sorted(required.difference(batch.columns))
        if missing:
            raise ValueError(f"component edge batch is missing columns {missing}")
        if batch.empty:
            continue
        values = batch[[query_column, target_column, similarity_column]].dropna()
        similarities = pd.to_numeric(values[similarity_column], errors="coerce")
        for threshold, labels in label_maps.items():
            qualifies = similarities >= threshold
            if not qualifies.any():
                continue
            query_labels = values.loc[qualifies, query_column].astype(str).map(labels)
            target_labels = values.loc[qualifies, target_column].astype(str).map(labels)
            eligible = query_labels.notna() & target_labels.notna()
            crossings[threshold] += int(
                (query_labels.loc[eligible] != target_labels.loc[eligible]).sum()
            )
    return crossings


def make_nk_graph(
    df: pd.DataFrame,
    num_systems: int,
    system_ids_cat: pd.CategoricalDtype,
    directed: bool,
    weighted: bool,
    query_col: str = "query_system",
    target_col: str = "target_system",
    similarity_col: str = "similarity",
) -> tuple[nk.graph.Graph, pd.CategoricalDtype]:
    """
    Make networkit graphs

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe of edges
    num_systems : int
        Number of systems
    system_ids_cat : pd.CategoricalDtype
        system ids, must match exactly the entries in df
    directed : bool
        Whether to make directed graph

    Returns
    -------
    graph
    """
    if df.empty:
        return nk.Graph(
            num_systems, weighted=weighted, directed=directed
        ), system_ids_cat

    df[[query_col, target_col]] = df[[query_col, target_col]].astype(system_ids_cat)
    if weighted:
        graph = nk.GraphFromCoo(
            (
                df[similarity_col].values / 100.0,
                (
                    df[query_col].cat.codes.to_numpy(dtype=np.uint, copy=False),
                    df[target_col].cat.codes.to_numpy(dtype=np.uint, copy=False),
                ),
            ),
            weighted=True,
            n=num_systems,
            directed=directed,
        )
    else:
        graph = nk.GraphFromCoo(
            (
                df[query_col].cat.codes.to_numpy(dtype=np.uint, copy=False),
                df[target_col].cat.codes.to_numpy(dtype=np.uint, copy=False),
            ),
            n=num_systems,
            directed=directed,
        )
    return graph, system_ids_cat


def _eligible_annotation(data_dir: Path) -> pd.DataFrame:
    columns = [
        "system_id",
        "ligand_id",
        "system_type",
        "ligand_is_proper",
    ]
    annotation = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
        columns=columns,
    )
    return annotation[
        (annotation["system_type"] == "holo")
        & annotation["ligand_is_proper"].fillna(False).astype(bool)
    ]


def _greedy_set_cover(
    graph: nk.graph.Graph,
    nodes: Sequence[str],
) -> list[tuple[str, list[str]]]:
    """Partition a threshold graph by greedy cover and best-representative assignment.

    The first pass repeatedly selects the candidate whose closed neighborhood
    covers the most uncovered nodes, with stable node ID as the tie-break. The
    selected representatives stay fixed; every other node is then assigned to
    the connected representative with its highest edge weight. Consequently
    every member has a direct threshold-qualified edge to its representative.
    """
    if graph.isDirected():
        raise ValueError("set covering requires an undirected graph")
    if not graph.isWeighted():
        raise ValueError("set covering requires edge weights")
    if graph.numberOfNodes() != len(nodes):
        raise ValueError("graph node count does not match ligand IDs")
    size = len(nodes)
    uncovered = set(range(size))
    centroids: list[int] = []
    selected: set[int] = set()
    heap = [
        (-(graph.degree(node) + 1), str(nodes[node]), node)
        for node in range(size)
    ]
    heapq.heapify(heap)
    while uncovered:
        while heap:
            negative_count, node_id, centroid = heapq.heappop(heap)
            if centroid in selected:
                continue
            covered = ({centroid} if centroid in uncovered else set())
            covered.update(
                int(neighbor)
                for neighbor in graph.iterNeighbors(centroid)
                if int(neighbor) in uncovered
            )
            actual_count = len(covered)
            if actual_count == -negative_count:
                break
            if actual_count:
                heapq.heappush(
                    heap,
                    (-actual_count, node_id, centroid),
                )
        else:
            raise RuntimeError("greedy set cover left nodes uncovered")
        selected.add(centroid)
        centroids.append(centroid)
        uncovered.difference_update(covered)
        if len(centroids) % 100_000 == 0:
            LOG.info(
                "set-cover progress: nodes=%d representatives=%d covered=%d "
                "uncovered=%d",
                len(nodes),
                len(centroids),
                len(nodes) - len(uncovered),
                len(uncovered),
            )

    centroid_set = set(centroids)
    groups: dict[int, list[str]] = {centroid: [] for centroid in centroids}
    for node, ligand_id in enumerate(nodes):
        if node in centroid_set:
            representative = node
        else:
            candidates = [
                int(neighbor)
                for neighbor in graph.iterNeighbors(node)
                if int(neighbor) in centroid_set
            ]
            if not candidates:
                raise RuntimeError(
                    f"greedy cover left {ligand_id} without a connected centroid"
                )
            representative = min(
                candidates,
                key=lambda centroid: (
                    -float(graph.weight(node, centroid)),
                    str(nodes[centroid]),
                ),
            )
        groups[representative].append(str(ligand_id))
    return [
        (str(nodes[centroid]), groups[centroid]) for centroid in centroids
    ]


@dataclass(frozen=True)
class _RepresentativeSelection:
    representative: int
    order: int
    threshold: int | None
    marginal_gain: int


@dataclass(frozen=True)
class _RepresentativeAssignment:
    query: int
    representative: int
    similarity: float
    threshold: int | None


def _selection_mask(
    values: NDArray[np.bool_] | None,
    *,
    size: int,
    name: str,
) -> NDArray[np.bool_]:
    if values is None:
        return np.ones(size, dtype=bool)
    result = np.asarray(values, dtype=bool)
    if result.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},), got {result.shape}")
    return result.copy()


def _incoming_cover(
    graph: nk.graph.Graph,
    representative: int,
    *,
    minimum_weight: float,
    uncovered: NDArray[np.bool_],
) -> NDArray[np.int64]:
    covered = [representative] if uncovered[representative] else []
    covered.extend(
        int(query)
        for query in graph.iterInNeighbors(representative)
        if int(query) != representative
        and uncovered[int(query)]
        and float(graph.weight(int(query), representative)) >= minimum_weight
    )
    return np.asarray(covered, dtype=np.int64)


def _greedy_directed_set_cover(
    graph: nk.graph.Graph,
    nodes: Sequence[str],
    *,
    primary_threshold: int,
    fallback_threshold: int = 50,
    fallback_graph: nk.graph.Graph | None = None,
    candidate_mask: NDArray[np.bool_] | None = None,
    target_mask: NDArray[np.bool_] | None = None,
) -> tuple[list[_RepresentativeSelection], list[_RepresentativeAssignment]]:
    """Select representatives by residual gain at two directed thresholds.

    An edge ``query -> representative`` means that selecting the target covers
    the query. At each threshold, the candidate covering the most uncovered
    queries is selected, with stable node ID as the only tie-break. Existing
    representatives are expanded when the threshold drops. Any target still
    uncovered after the fallback pass represents itself.
    """
    if not graph.isDirected() or not graph.isWeighted():
        raise ValueError("directed representative covering requires a weighted graph")
    if graph.numberOfNodes() != len(nodes):
        raise ValueError("graph node count does not match ligand IDs")
    if fallback_graph is None:
        fallback_graph = graph
    if not fallback_graph.isDirected() or not fallback_graph.isWeighted():
        raise ValueError("fallback covering requires a weighted directed graph")
    if fallback_graph.numberOfNodes() != len(nodes):
        raise ValueError("fallback graph node count does not match ligand IDs")
    if not 0 <= fallback_threshold <= primary_threshold <= 100:
        raise ValueError(
            "representative thresholds must satisfy "
            "0 <= fallback_threshold <= primary_threshold <= 100"
        )
    size = len(nodes)
    candidates = _selection_mask(
        candidate_mask,
        size=size,
        name="candidate_mask",
    )
    targets = _selection_mask(target_mask, size=size, name="target_mask")

    primary_weight = primary_threshold / 100.0
    fallback_weight = fallback_threshold / 100.0
    uncovered = targets.copy()
    selected = np.zeros(size, dtype=bool)

    def make_heap(
        active_graph: nk.graph.Graph,
        minimum_weight: float,
        heap_candidates: NDArray[np.bool_],
    ) -> list[tuple[int, str, int]]:
        heap = [
            (
                -len(
                    _incoming_cover(
                        active_graph,
                        node,
                        minimum_weight=minimum_weight,
                        uncovered=uncovered,
                    )
                ),
                str(nodes[node]),
                node,
            )
            for node in np.flatnonzero(heap_candidates)
        ]
        heapq.heapify(heap)
        return heap

    def best_candidate(
        heap: list[tuple[int, str, int]],
        active_graph: nk.graph.Graph,
        minimum_weight: float,
    ) -> tuple[int | None, NDArray[np.int64]]:
        while heap:
            negative_gain, node_id, candidate = heapq.heappop(heap)
            if selected[candidate]:
                continue
            covered = _incoming_cover(
                active_graph,
                candidate,
                minimum_weight=minimum_weight,
                uncovered=uncovered,
            )
            gain = len(covered)
            if gain == -negative_gain:
                # A self-only choice has no edge-derived gain. Try the lower
                # threshold before assigning any remaining node to itself.
                if len(covered) and np.any(covered != candidate):
                    return candidate, covered
                continue
            heapq.heappush(
                heap,
                (-gain, node_id, candidate),
            )
        return None, np.asarray([], dtype=np.int64)

    selections: list[_RepresentativeSelection] = []

    def select_from_heap(
        heap: list[tuple[int, str, int]],
        *,
        active_graph: nk.graph.Graph,
        minimum_weight: float,
        selection_threshold: int,
    ) -> None:
        while uncovered.any():
            representative, newly_covered = best_candidate(
                heap,
                active_graph,
                minimum_weight,
            )
            if representative is None:
                return
            selected[representative] = True
            uncovered[newly_covered] = False
            selections.append(
                _RepresentativeSelection(
                    representative=representative,
                    order=len(selections),
                    threshold=selection_threshold,
                    marginal_gain=len(newly_covered),
                )
            )
            if len(selections) % 100_000 == 0:
                LOG.info(
                    "directed representative cover progress: nodes=%d "
                    "selected=%d covered=%d uncovered=%d",
                    size,
                    len(selections),
                    int(targets.sum() - uncovered.sum()),
                    int(uncovered.sum()),
                )

    select_from_heap(
        make_heap(graph, primary_weight, candidates & ~selected),
        active_graph=graph,
        minimum_weight=primary_weight,
        selection_threshold=primary_threshold,
    )

    if fallback_threshold < primary_threshold and uncovered.any():
        # Representatives selected from stricter edges remain valid choices.
        # Before adding any new representative, let them cover all remaining
        # queries that reach them at the fallback threshold.
        for representative in np.flatnonzero(selected):
            newly_covered = _incoming_cover(
                fallback_graph,
                int(representative),
                minimum_weight=fallback_weight,
                uncovered=uncovered,
            )
            uncovered[newly_covered] = False
        select_from_heap(
            make_heap(
                fallback_graph,
                fallback_weight,
                candidates & ~selected,
            ),
            active_graph=fallback_graph,
            minimum_weight=fallback_weight,
            selection_threshold=fallback_threshold,
        )

    while uncovered.any():
        remaining = np.flatnonzero(uncovered)
        representative = min(
            map(int, remaining),
            key=lambda node: str(nodes[node]),
        )
        selected[representative] = True
        uncovered[representative] = False
        selections.append(
            _RepresentativeSelection(
                representative=representative,
                order=len(selections),
                threshold=None,
                marginal_gain=1,
            )
        )

    selection_by_node = {
        selection.representative: selection for selection in selections
    }
    assignments: list[_RepresentativeAssignment] = []
    for query in np.flatnonzero(targets):
        query = int(query)
        if query in selection_by_node:
            representative = query
            similarity = 100.0
            assignment_threshold = primary_threshold
        else:
            strict: list[int] = []
            for target in graph.iterNeighbors(query):
                representative = int(target)
                if (
                    representative in selection_by_node
                    and float(graph.weight(query, representative)) >= primary_weight
                ):
                    strict.append(representative)
            options = strict
            assignment_threshold = primary_threshold
            if not options and fallback_threshold < primary_threshold:
                options = []
                for target in fallback_graph.iterNeighbors(query):
                    representative = int(target)
                    if (
                        representative in selection_by_node
                        and float(fallback_graph.weight(query, representative))
                        >= fallback_weight
                    ):
                        options.append(representative)
                assignment_threshold = fallback_threshold
            if not options:
                raise RuntimeError(
                    f"directed cover left {nodes[query]} without a representative"
                )
            representative = min(
                options,
                key=lambda candidate: (
                    -float(
                        (
                            graph
                            if assignment_threshold == primary_threshold
                            else fallback_graph
                        ).weight(query, candidate)
                    ),
                    str(nodes[candidate]),
                ),
            )
            assignment_graph = (
                graph
                if assignment_threshold == primary_threshold
                else fallback_graph
            )
            similarity = 100.0 * float(
                assignment_graph.weight(query, representative)
            )
        assignments.append(
            _RepresentativeAssignment(
                query=query,
                representative=representative,
                similarity=similarity,
                threshold=assignment_threshold,
            )
        )
    return selections, assignments


def _greedy_directed_centroid_cover(
    graph: nk.graph.Graph,
    nodes: Sequence[str],
) -> list[tuple[str, str, float, int, float]]:
    """Compatibility wrapper for a full node-level directed cover."""
    selections, assignments = _greedy_directed_set_cover(
        graph,
        nodes,
        primary_threshold=0,
        fallback_threshold=0,
    )
    del selections
    target_mask = np.ones(len(nodes), dtype=bool)
    coverage_counts = [
        len(
            _incoming_cover(
                graph,
                node,
                minimum_weight=0.0,
                uncovered=target_mask,
            )
        )
        for node in range(len(nodes))
    ]
    return [
        (
            str(nodes[assignment.query]),
            str(nodes[assignment.representative]),
            assignment.similarity,
            coverage_counts[assignment.query],
            coverage_counts[assignment.query] / len(nodes),
        )
        for assignment in assignments
    ]


def _expand_fingerprint_cover(
    *, data_dir: Path, assignments: pd.DataFrame
) -> pd.DataFrame:
    """Expand SMILES-node assignments to ligand IDs and concrete centroids."""
    ligands_per_smiles = pd.read_parquet(
        data_dir / "fingerprints/ligands_per_smiles.parquet",
        columns=["ligand_rdkit_canonical_smiles", "ligand_smiles_id"],
    )
    annotation = pd.read_parquet(
        data_dir / "index/annotation_table.parquet",
        columns=[
            "ligand_id",
            "ligand_rdkit_canonical_smiles",
            "ligand_is_proper",
        ],
        filters=[("ligand_is_proper", "==", True)],
    )
    smiles_to_id = dict(
        zip(
            ligands_per_smiles["ligand_rdkit_canonical_smiles"],
            ligands_per_smiles["ligand_smiles_id"].astype(str),
        )
    )
    annotation["node"] = annotation["ligand_rdkit_canonical_smiles"].map(smiles_to_id)
    annotation.dropna(subset=["node"], inplace=True)
    ligands_by_node = {
        str(node): tuple(sorted(set(group["ligand_id"].astype(str))))
        for node, group in annotation.groupby("node", sort=False)
    }
    representative_by_node = {
        str(node): min(group["ligand_id"].astype(str))
        for node, group in annotation.groupby("node", sort=False)
    }
    expanded = assignments.copy()
    expanded["ligand_id"] = expanded["ligand_id"].map(ligands_by_node)
    expanded["centroid_ligand_id"] = expanded["centroid_node"].map(
        representative_by_node
    )
    return (
        expanded.dropna(subset=["ligand_id", "centroid_ligand_id"])
        .explode("ligand_id")
        .reset_index(drop=True)
    )


def directed_set_cover_is_complete(
    path: Path,
    *,
    entity_type: ClusterEntity = "ligand",
) -> bool:
    """Return whether a cached directed-cover artifact has the current schema."""
    centroid_column = (
        "centroid_ligand_id"
        if entity_type == "ligand"
        else "centroid_system_id"
    )
    required_columns = {
        _cluster_node_column(entity_type),
        centroid_column,
        "similarity_to_centroid",
        "coverage_count",
        "coverage_fraction",
        "representative_selection_order",
        "representative_selection_threshold",
        "representative_marginal_gain",
        "assignment_threshold",
        "label",
        "metric",
        "threshold",
        "directed",
    }
    try:
        return path.is_file() and required_columns.issubset(
            pq.read_schema(path).names
        )
    except (OSError, ValueError):
        return False


def set_cover_is_complete(path: Path) -> bool:
    """Return whether a cached ligand-Tanimoto cover has the current schema."""
    required_columns = {
        "ligand_id",
        "centroid_ligand_id",
        "label",
        "metric",
        "cluster",
        "threshold",
        "directed",
    }
    try:
        return path.is_file() and required_columns.issubset(
            pq.read_schema(path).names
        )
    except (OSError, ValueError):
        return False


def make_directed_set_cover(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
    skip_existing: bool = False,
    scratch_dir: Path | None = None,
    threads: int = 1,
    entity_type: ClusterEntity = "ligand",
) -> Path:
    """Build a directed greedy set cover from score shards."""
    started = time()
    if entity_type == "ligand" and metric in GATED_LIGAND_DIAGNOSTIC_METRICS:
        raise ValueError(
            f"{metric} cannot be clustered directly; use " "sucos_shape_pocket_qcov"
        )
    if threads < 1:
        raise ValueError("directed-cover threads must be positive")
    if not 0 <= threshold <= 100:
        raise ValueError("directed-cover threshold must be in [0, 100]")
    node_column = _cluster_node_column(entity_type)
    output = (
        _sampling_root(data_dir, entity_type)
        / "directed_set_cover"
        / f"metric={metric}"
        / f"threshold={threshold}.parquet"
    )
    if skip_existing and directed_set_cover_is_complete(
        output,
        entity_type=entity_type,
    ):
        return output
    fallback_threshold = min(threshold, 50)
    component_path = (
        _sampling_root(data_dir, entity_type)
        / "directed_set_cover"
        / "reductions"
        / f"metric={metric}"
        / "labels"
        / f"threshold={fallback_threshold}.parquet"
    )
    if not component_path.is_file():
        raise FileNotFoundError(
            "directed-cover connectivity must be merged before covering: "
            f"{component_path}"
        )
    component_labels = pd.read_parquet(component_path)
    component_labels[node_column] = component_labels[node_column].astype(str)
    component_labels["label"] = component_labels["label"].astype(str)
    component_labels["component"] = pd.factorize(component_labels["label"], sort=False)[
        0
    ].astype(np.uint32)
    component_labels = component_labels.sort_values(
        ["component", node_column], kind="stable"
    ).reset_index(drop=True)
    component_labels["component_node"] = (
        component_labels.groupby("component", sort=False).cumcount().astype(np.uint32)
    )
    nodes_by_component = {
        int(component): group[node_column].tolist()
        for component, group in component_labels.groupby("component", sort=False)
    }
    component_lookup = component_labels[[node_column, "component", "component_node"]]
    sources = component_score_sources(
        data_dir=data_dir,
        metric=metric,
        entity_type=entity_type,
    )

    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    temporary_root = scratch_dir or data_dir / "scratch/duckdb/directed_set_cover"
    temporary_root.mkdir(exist_ok=True, parents=True)
    connection.sql(f"SET temp_directory='{temporary_root.as_posix()}'")
    connection.register("component_labels", component_lookup)
    paths_sql = ", ".join(
        f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'" for path in sources
    )
    query = dedent(
        f"""
        WITH directed_edges AS (
            SELECT
                cast(query_node AS VARCHAR) AS query_node,
                cast(target_node AS VARCHAR) AS target_node,
                forward_similarity::DOUBLE AS similarity
            FROM read_parquet([{paths_sql}])
            WHERE forward_similarity >= {fallback_threshold}
            UNION ALL
            SELECT
                cast(target_node AS VARCHAR) AS query_node,
                cast(query_node AS VARCHAR) AS target_node,
                reverse_similarity::DOUBLE AS similarity
            FROM read_parquet([{paths_sql}])
            WHERE reverse_similarity >= {fallback_threshold}
        ), labeled AS (
            SELECT
                query_labels.component AS query_component,
                target_labels.component AS target_component,
                query_labels.component_node::UINTEGER AS query_node,
                target_labels.component_node::UINTEGER AS target_node,
                CAST(round(directed_edges.similarity) AS UTINYINT) AS similarity,
                query_labels.component != target_labels.component AS crossing_edge
            FROM directed_edges
            INNER JOIN component_labels AS query_labels
              ON directed_edges.query_node = query_labels.{node_column}
            INNER JOIN component_labels AS target_labels
              ON directed_edges.target_node = target_labels.{node_column}
        )
        SELECT
            query_component AS component,
            query_node,
            target_node,
            similarity,
            crossing_edge
        FROM labeled
        ORDER BY crossing_edge DESC, component, query_node, target_node
        """
    )
    LOG.info(
        "directed cover edge query start: metric=%s threshold=%d "
        "fallback_threshold=%d sources=%d nodes=%d components=%d",
        metric,
        threshold,
        fallback_threshold,
        len(sources),
        len(component_labels),
        len(nodes_by_component),
    )
    staged_edges = temporary_root / (
        f"{metric}-{threshold}-{fallback_threshold}-directed-edges.parquet"
    )
    staged_edges.unlink(missing_ok=True)
    staged_schema = pa.schema(
        [
            ("component", pa.uint32()),
            ("query_node", pa.uint32()),
            ("target_node", pa.uint32()),
            ("similarity", pa.uint8()),
        ]
    )
    staged_rows = 0
    stage_started = time()
    writer = pq.ParquetWriter(
        staged_edges,
        staged_schema,
        compression="zstd",
        use_dictionary=False,
        write_statistics=False,
    )
    try:
        reader = connection.execute(query).fetch_record_batch(rows_per_batch=250_000)
        for batch_index, record_batch in enumerate(reader, start=1):
            crossing_index = record_batch.schema.get_field_index("crossing_edge")
            crossing_edges = int(
                np.asarray(
                    record_batch.column(crossing_index).to_numpy(zero_copy_only=False),
                    dtype=bool,
                ).sum()
            )
            if crossing_edges:
                raise ValueError(
                    f"directed-cover component validation failed for {metric} at "
                    f"{threshold}: {crossing_edges} crossing edges"
                )
            compact = pa.Table.from_batches([record_batch]).select(
                staged_schema.names
            ).cast(staged_schema, safe=False)
            writer.write_table(compact, row_group_size=1_000_000)
            staged_rows += len(record_batch)
            if batch_index % 20 == 0:
                elapsed = time() - stage_started
                LOG.info(
                    "directed cover edge staging progress: metric=%s threshold=%d "
                    "batches=%d rows=%d rate=%.1f_rows/s elapsed_seconds=%.1f",
                    metric,
                    threshold,
                    batch_index,
                    staged_rows,
                    staged_rows / elapsed,
                    elapsed,
                )
    except BaseException:
        writer.close()
        connection.close()
        staged_edges.unlink(missing_ok=True)
        raise
    else:
        writer.close()
        connection.close()
    LOG.info(
        "directed cover edge staging complete: metric=%s threshold=%d rows=%d "
        "size_gb=%.2f elapsed_seconds=%.1f",
        metric,
        threshold,
        staged_rows,
        staged_edges.stat().st_size / (1024**3),
        time() - stage_started,
    )
    del reader
    gc.collect()

    assignments: list[dict[str, Any]] = []
    processed_nodes: set[str] = set()
    representative_count = 0
    current_component: int | None = None
    current_nodes: list[str] = []
    current_graph: nk.graph.Graph | None = None
    current_primary_graph: nk.graph.Graph | None = None

    def finish_component() -> None:
        nonlocal current_component, current_graph, current_primary_graph
        nonlocal representative_count
        if (
            current_component is None
            or current_graph is None
            or current_primary_graph is None
        ):
            return
        if len(current_nodes) > 1 and current_graph.numberOfEdges() == 0:
            raise ValueError(
                f"non-singleton directed-cover component {current_component} "
                "has no edges"
            )
        selections, component_assignments = _greedy_directed_set_cover(
            current_primary_graph,
            current_nodes,
            primary_threshold=threshold,
            fallback_threshold=fallback_threshold,
            fallback_graph=current_graph,
        )
        selection_by_node = {
            selection.representative: selection for selection in selections
        }
        all_targets = np.ones(len(current_nodes), dtype=bool)
        coverage_counts = [
            len(
                _incoming_cover(
                    current_primary_graph,
                    node,
                    minimum_weight=threshold / 100.0,
                    uncovered=all_targets,
                )
            )
            for node in range(len(current_nodes))
        ]
        for assignment in component_assignments:
            selection = selection_by_node[assignment.representative]
            assignments.append(
                {
                    node_column: str(current_nodes[assignment.query]),
                    "centroid_node": str(
                        current_nodes[assignment.representative]
                    ),
                    "similarity_to_centroid": assignment.similarity,
                    "coverage_count": coverage_counts[assignment.query],
                    "coverage_fraction": coverage_counts[assignment.query]
                    / len(current_nodes),
                    "representative_selection_order": representative_count
                    + selection.order,
                    "representative_selection_threshold": selection.threshold,
                    "representative_marginal_gain": selection.marginal_gain,
                    "assignment_threshold": assignment.threshold,
                }
            )
        representative_count += len(selections)
        processed_nodes.update(current_nodes)

    streamed_rows = 0
    stream_started = time()
    try:
        staged_parquet = pq.ParquetFile(staged_edges)
        for batch_index, record_batch in enumerate(
            staged_parquet.iter_batches(batch_size=250_000), start=1
        ):
            frame = record_batch.to_pandas()
            streamed_rows += len(frame)
            if batch_index % 20 == 0:
                elapsed = time() - stream_started
                LOG.info(
                    "directed cover graph progress: metric=%s threshold=%d "
                    "batches=%d rows=%d rate=%.1f_rows/s elapsed_seconds=%.1f",
                    metric,
                    threshold,
                    batch_index,
                    streamed_rows,
                    streamed_rows / elapsed,
                    elapsed,
                )
            for component, group in frame.groupby("component", sort=False):
                component = int(component)
                if component != current_component:
                    finish_component()
                    current_component = component
                    current_nodes = nodes_by_component[component]
                    current_graph = nk.Graph(
                        len(current_nodes),
                        weighted=True,
                        directed=True,
                    )
                    current_primary_graph = (
                        current_graph
                        if threshold == fallback_threshold
                        else nk.Graph(
                            len(current_nodes),
                            weighted=True,
                            directed=True,
                        )
                    )
                assert current_graph is not None
                assert current_primary_graph is not None
                similarities = group["similarity"].to_numpy(
                    dtype=float, copy=False
                )
                query_nodes = group["query_node"].to_numpy(
                    dtype=np.uint, copy=False
                )
                target_nodes = group["target_node"].to_numpy(
                    dtype=np.uint, copy=False
                )
                weighted_edges = (
                    similarities / 100.0,
                    (query_nodes, target_nodes),
                )
                current_graph.addEdges(weighted_edges)
                if current_primary_graph is not current_graph:
                    primary_edges = similarities >= threshold
                    if primary_edges.any():
                        current_primary_graph.addEdges(
                            (
                                similarities[primary_edges] / 100.0,
                                (
                                    query_nodes[primary_edges],
                                    target_nodes[primary_edges],
                                ),
                            )
                        )
        finish_component()
    finally:
        staged_edges.unlink(missing_ok=True)
    for component, nodes in nodes_by_component.items():
        missing = [node for node in nodes if node not in processed_nodes]
        if len(missing) > 1:
            raise ValueError(
                f"non-singleton directed-cover component {component} was "
                "absent from edges"
            )
        for node in missing:
            assignments.append(
                {
                    node_column: node,
                    "centroid_node": node,
                    "similarity_to_centroid": 100.0,
                    "coverage_count": 1,
                    "coverage_fraction": 1.0,
                    "representative_selection_order": representative_count,
                    "representative_selection_threshold": None,
                    "representative_marginal_gain": 1,
                    "assignment_threshold": threshold,
                }
            )
            representative_count += 1

    published = pd.DataFrame.from_records(
        assignments,
        columns=[
            node_column,
            "centroid_node",
            "similarity_to_centroid",
            "coverage_count",
            "coverage_fraction",
            "representative_selection_order",
            "representative_selection_threshold",
            "representative_marginal_gain",
            "assignment_threshold",
        ],
    )
    published["coverage_count"] = published["coverage_count"].astype("Int32")
    published["coverage_fraction"] = published["coverage_fraction"].astype(
        "Float32"
    )
    published["representative_selection_order"] = published[
        "representative_selection_order"
    ].astype("Int32")
    published["representative_selection_threshold"] = published[
        "representative_selection_threshold"
    ].astype("Int16")
    published["representative_marginal_gain"] = published[
        "representative_marginal_gain"
    ].astype("Int32")
    published["assignment_threshold"] = published["assignment_threshold"].astype(
        "Int16"
    )
    if metric == "tanimoto_similarity_ecfp4_1024":
        published = _expand_fingerprint_cover(
            data_dir=data_dir,
            assignments=published,
        )
    else:
        published[
            "centroid_ligand_id" if entity_type == "ligand" else "centroid_system_id"
        ] = published["centroid_node"]
    group_sizes = published.groupby("centroid_node", sort=False).size()
    ordered_centroids = sorted(
        group_sizes.index.astype(str),
        key=lambda centroid: (-int(group_sizes.loc[centroid]), centroid),
    )
    published["label"] = published["centroid_node"].map(
        {centroid: f"c{index}" for index, centroid in enumerate(ordered_centroids)}
    )
    published["metric"] = metric
    published["threshold"] = threshold
    published["directed"] = True
    published.drop(columns=["centroid_node"], inplace=True)
    published.sort_values(node_column, inplace=True, kind="stable")
    output.parent.mkdir(exist_ok=True, parents=True)
    temporary = output.with_suffix(".tmp.parquet")
    published.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(output)
    LOG.info(
        "directed cover complete: metric=%s threshold=%d assignments=%d "
        "centroids=%d elapsed_seconds=%.1f",
        metric,
        threshold,
        len(published),
        published[
            "centroid_ligand_id" if entity_type == "ligand" else "centroid_system_id"
        ].nunique(),
        time() - started,
    )
    return output


def make_set_cover(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
    skip_existing_clusters: bool = False,
    scratch_dir: Path | None = None,
    threads: int = 1,
    entity_type: ClusterEntity = "ligand",
) -> Path:
    """Build an undirected greedy set cover on reciprocal-minimum edges."""
    started = time()
    if entity_type != "ligand" or metric != "tanimoto_similarity_ecfp4_1024":
        raise ValueError("undirected set cover is supported only for ligand Tanimoto")
    if entity_type == "ligand" and metric in GATED_LIGAND_DIAGNOSTIC_METRICS:
        raise ValueError(
            f"{metric} is evaluated only for ligand pairs with positive pocket "
            "coverage and cannot be clustered directly; use "
            "sucos_shape_pocket_qcov"
        )
    if threads < 1:
        raise ValueError("set-cover threads must be positive")
    node_column = _cluster_node_column(entity_type)
    output = (
        _sampling_root(data_dir, entity_type)
        / "set_cover"
        / f"metric={metric}"
        / f"threshold={threshold}.parquet"
    )
    if skip_existing_clusters and set_cover_is_complete(output):
        return output
    component_path = (
        _cluster_root(data_dir, entity_type)
        / "reductions"
        / f"metric={metric}"
        / "labels"
        / "directed=false"
        / f"threshold={threshold}.parquet"
    )
    if not component_path.is_file():
        raise FileNotFoundError(
            "internal connectivity labels must be merged before set cover: "
            f"{component_path}"
        )
    component_labels = pd.read_parquet(component_path)
    component_labels[node_column] = component_labels[node_column].astype(str)
    if component_labels[node_column].duplicated().any():
        raise ValueError(f"duplicate node IDs in component labels: {component_path}")
    component_labels["label"] = component_labels["label"].astype(str)
    component_labels["component"] = pd.factorize(component_labels["label"], sort=False)[
        0
    ].astype(np.uint32)
    component_labels = component_labels.sort_values(
        ["component", node_column], kind="stable"
    ).reset_index(drop=True)
    component_labels["component_node"] = (
        component_labels.groupby("component", sort=False).cumcount().astype(np.uint32)
    )
    nodes_by_component = {
        int(component): group[node_column].tolist()
        for component, group in component_labels.groupby("component", sort=False)
    }
    component_lookup = component_labels[[node_column, "component", "component_node"]]
    sources = component_score_sources(
        data_dir=data_dir,
        metric=metric,
        entity_type=entity_type,
    )
    if not sources:
        raise FileNotFoundError(f"no set-cover score sources found for {metric}")

    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    temporary_root = scratch_dir or data_dir / "scratch" / "duckdb" / "set_cover"
    temporary_root.mkdir(exist_ok=True, parents=True)
    connection.sql(f"SET temp_directory='{temporary_root.as_posix()}'")
    connection.register("component_labels", component_lookup)
    paths_sql = ", ".join(
        f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'" for path in sources
    )
    query = dedent(
        f"""
        WITH selected AS (
            SELECT
                cast(query_node AS VARCHAR) AS query_node,
                cast(target_node AS VARCHAR) AS target_node,
                cast(similarity AS DOUBLE) AS similarity
            FROM read_parquet([{paths_sql}])
            WHERE similarity >= {threshold}
        ), labeled AS (
            SELECT
                query_labels.component AS query_component,
                target_labels.component AS target_component,
                query_labels.component_node::UINTEGER AS query_node,
                target_labels.component_node::UINTEGER AS target_node,
                selected.similarity,
                query_labels.component != target_labels.component AS crossing_edge
            FROM selected
            INNER JOIN component_labels AS query_labels
                ON selected.query_node = query_labels.{node_column}
            INNER JOIN component_labels AS target_labels
                ON selected.target_node = target_labels.{node_column}
        )
        SELECT
            query_component AS component,
            query_node,
            target_node,
            similarity,
            crossing_edge
        FROM labeled
        ORDER BY crossing_edge DESC, component, query_node, target_node
        """
    )
    LOG.info(
        "set-cover edge query start: metric=%s threshold=%d sources=%d "
        "nodes=%d components=%d elapsed_seconds=%.1f",
        metric,
        threshold,
        len(sources),
        len(component_labels),
        len(nodes_by_component),
        time() - started,
    )
    query_started = time()
    reader = connection.execute(query).fetch_record_batch(rows_per_batch=250_000)
    LOG.info(
        "set-cover edge query ready: metric=%s threshold=%d elapsed_seconds=%.1f",
        metric,
        threshold,
        time() - query_started,
    )
    cover_groups: list[tuple[str, list[str]]] = []
    processed_nodes: set[str] = set()
    current_component: int | None = None
    current_nodes: list[str] = []
    current_graph: nk.graph.Graph | None = None

    def finish_component() -> None:
        nonlocal current_component, current_graph
        if current_component is None or current_graph is None:
            return
        if len(current_nodes) > 1 and current_graph.numberOfEdges() == 0:
            raise ValueError(
                f"non-singleton component {current_component} has no edges"
            )
        if current_graph.numberOfEdges() == 0:
            groups = [(node, [node]) for node in current_nodes]
        else:
            groups = _greedy_set_cover(
                current_graph,
                current_nodes,
            )
        cover_groups.extend(groups)
        processed_nodes.update(current_nodes)

    streamed_rows = 0
    stream_started = time()
    for batch_index, record_batch in enumerate(reader, start=1):
        frame = record_batch.to_pandas()
        if frame.empty:
            continue
        crossing_edges = int(frame["crossing_edge"].fillna(False).sum())
        if crossing_edges:
            raise ValueError(
                f"connectivity partition validation failed for {metric} at "
                f"{threshold}: {crossing_edges} crossing edges"
            )
        frame = frame.dropna(subset=["component"])
        streamed_rows += len(frame)
        if batch_index % 20 == 0:
            elapsed = time() - stream_started
            LOG.info(
                "set-cover edge stream progress: metric=%s threshold=%d "
                "batches=%d rows=%d rate=%.1f_rows/s elapsed_seconds=%.1f",
                metric,
                threshold,
                batch_index,
                streamed_rows,
                streamed_rows / elapsed,
                elapsed,
            )
        for component, group in frame.groupby("component", sort=False):
            component = int(component)
            if component != current_component:
                finish_component()
                current_component = component
                current_nodes = nodes_by_component[component]
                current_graph = nk.Graph(
                    len(current_nodes),
                    weighted=True,
                    directed=False,
                )
            assert current_graph is not None
            query_nodes = group["query_node"].to_numpy(dtype=np.uint, copy=False)
            target_nodes = group["target_node"].to_numpy(dtype=np.uint, copy=False)
            similarities = group["similarity"].to_numpy(dtype=float, copy=False)
            current_graph.addEdges((similarities / 100.0, (query_nodes, target_nodes)))
    finish_component()
    connection.close()
    LOG.info(
        "set-cover edge stream complete: metric=%s threshold=%d rows=%d "
        "elapsed_seconds=%.1f",
        metric,
        threshold,
        streamed_rows,
        time() - stream_started,
    )

    unprocessed = set(component_labels[node_column]) - processed_nodes
    for component, nodes in nodes_by_component.items():
        missing = [node for node in nodes if node in unprocessed]
        if len(missing) > 1:
            raise ValueError(
                f"non-singleton component {component} was absent from edges"
            )
        cover_groups.extend((node, [node]) for node in missing)
    ordered_groups = sorted(
        cover_groups,
        key=lambda item: (-len(item[1]), item[0]),
    )
    published = pd.DataFrame.from_records(
        [
            {
                node_column: node,
                "centroid_node": representative,
                "label": f"c{cover_index}",
            }
            for cover_index, (representative, group) in enumerate(ordered_groups)
            for node in sorted(group)
        ],
        columns=[node_column, "centroid_node", "label"],
    )
    published = _expand_fingerprint_cover(
        data_dir=data_dir,
        assignments=published,
    )
    published.drop(columns=["centroid_node"], inplace=True)
    published.sort_values(node_column, inplace=True, kind="stable")
    if published[node_column].duplicated().any():
        raise ValueError(
            "ligand Tanimoto set cover produced duplicate ligand assignments"
        )
    published["metric"] = metric
    published["directed"] = False
    published["threshold"] = threshold
    published["cluster"] = "set_cover"
    output.parent.mkdir(exist_ok=True, parents=True)
    temporary = output.with_suffix(".tmp.parquet")
    published.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(output)
    return output
