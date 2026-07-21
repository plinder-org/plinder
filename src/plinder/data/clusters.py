# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import hashlib
import json
import os
import sys
from pathlib import Path
from textwrap import dedent
from time import time
from typing import Any, Callable, Iterable, Mapping, Sequence, TypeVar, cast

if sys.platform == "darwin":
    # For macOS only: allow multiple OpenMP runtimes to coexist
    # (needed on macOS with conda)
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import networkit as nk
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from numpy.typing import NDArray

from plinder.core.scores.metrics import GATED_LIGAND_DIAGNOSTIC_METRICS
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import LIGAND_CLUSTER_SCHEMA

LOG = setup_logger(__name__)

T = TypeVar("T")

COMPONENT_EDGE_COLUMNS = ["query_node", "target_node"]
COMPONENT_NODE_UNIVERSE_RELATIVE = Path(
    "ligand_clusters/reductions/node_universe.parquet"
)
COMPONENT_NODE_UNIVERSE_MANIFEST_RELATIVE = Path(
    "ligand_clusters/reductions/node_universe.json"
)


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
    """Persist resumable weak/strong reductions for one immutable score shard."""
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
        directed_modes=[False, True],
    )
    outputs: list[dict[str, str | int | bool]] = []
    for directed in [False, True]:
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
        "status": "complete",
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
    """Merge completed shard reductions into exact weak/strong labels."""
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
        for directed in [False, True]
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
        for directed in [False, True]:
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
    for directed in [False, True]:
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


def component_score_sources(*, data_dir: Path, metric: str) -> list[Path]:
    """Return the immutable Parquet shards that contain one clustering metric."""
    if metric == "tanimoto_similarity_ecfp4_1024":
        return sorted((data_dir / "ligand_scores").glob("*.parquet"))
    return sorted((data_dir / "scores" / "search_db=holo").rglob("*.parquet"))


def iter_component_score_batches(
    *,
    source_path: Path,
    metric: str,
    batch_size: int = 250_000,
    eligible_systems: set[str] | None = None,
) -> Iterable[pd.DataFrame]:
    """Read one score shard as bounded, normalized ligand-edge batches."""
    if batch_size < 1:
        raise ValueError("component score batch size must be positive")
    parquet = pq.ParquetFile(source_path)
    chemical = metric == "tanimoto_similarity_ecfp4_1024"
    if chemical:
        columns = ["query_ligand_id", "target_ligand_id", metric]
    else:
        columns = [
            "query_system",
            "query_ligand_id",
            "target_system",
            "target_ligand_id",
            "metric",
            "similarity",
        ]
    missing = sorted(set(columns).difference(parquet.schema.names))
    if missing:
        raise ValueError(f"component score shard {source_path} is missing {missing}")
    for batch in parquet.iter_batches(
        batch_size=batch_size,
        columns=columns,
    ):
        frame = batch.to_pandas()
        if chemical:
            frame = frame.rename(
                columns={
                    "query_ligand_id": "query_node",
                    "target_ligand_id": "target_node",
                    metric: "similarity",
                }
            )
        else:
            frame = frame[frame["metric"].astype(str).eq(metric)]
            if eligible_systems is not None:
                frame = frame[
                    frame["query_system"].astype(str).isin(eligible_systems)
                    & frame["target_system"].astype(str).isin(eligible_systems)
                ]
            frame = frame.rename(
                columns={
                    "query_ligand_id": "query_node",
                    "target_ligand_id": "target_node",
                }
            )
        if not frame.empty:
            yield frame[["query_node", "target_node", "similarity"]]


def _cached_component_node_universe(
    data_dir: Path,
) -> tuple[list[str], set[str]] | None:
    annotation_path = data_dir / "index" / "annotation_table.parquet"
    cache_path = data_dir / COMPONENT_NODE_UNIVERSE_RELATIVE
    manifest_path = data_dir / COMPONENT_NODE_UNIVERSE_MANIFEST_RELATIVE
    if (
        not annotation_path.is_file()
        or not cache_path.is_file()
        or not manifest_path.is_file()
    ):
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("version") != 1:
            return None
        if manifest.get("annotation") != _component_source_signature(annotation_path):
            return None
        if cache_path.stat().st_size != int(manifest["cache_size"]):
            return None
        if set(pq.read_schema(cache_path).names) != {"kind", "id"}:
            return None
        frame = pd.read_parquet(cache_path, columns=["kind", "id"])
        nodes = frame.loc[frame["kind"].eq("ligand"), "id"].astype(str).tolist()
        systems = set(frame.loc[frame["kind"].eq("system"), "id"].astype(str).tolist())
        if len(nodes) != int(manifest["ligand_count"]):
            return None
        if len(systems) != int(manifest["system_count"]):
            return None
        return nodes, systems
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def prepare_component_node_universe(data_dir: Path) -> dict[str, Any]:
    """Materialize the shared ligand/system universe once before scattering."""
    cached = _cached_component_node_universe(data_dir)
    if cached is not None:
        nodes, systems = cached
        return {
            "status": "cached",
            "ligand_count": len(nodes),
            "system_count": len(systems),
        }

    annotation_path = data_dir / "index" / "annotation_table.parquet"
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
    cache_path = data_dir / COMPONENT_NODE_UNIVERSE_RELATIVE
    manifest_path = data_dir / COMPONENT_NODE_UNIVERSE_MANIFEST_RELATIVE
    cache_path.parent.mkdir(exist_ok=True, parents=True)
    temporary_cache = cache_path.with_suffix(".tmp.parquet")
    frame.to_parquet(temporary_cache, index=False, compression="zstd")
    temporary_cache.replace(cache_path)
    manifest = {
        "version": 1,
        "annotation": _component_source_signature(annotation_path),
        "cache_size": cache_path.stat().st_size,
        "ligand_count": len(nodes),
        "system_count": len(system_ids),
    }
    temporary_manifest = manifest_path.with_suffix(".tmp.json")
    temporary_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary_manifest.replace(manifest_path)
    return {"status": "complete", **manifest}


def component_node_universe(
    *, data_dir: Path, metric: str
) -> tuple[list[str], set[str] | None]:
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
    cached = _cached_component_node_universe(data_dir)
    if cached is not None:
        return cached
    annotation = _eligible_annotation(data_dir)
    systems = set(annotation["system_id"].dropna().astype(str))
    nodes = set(
        annotation.loc[annotation["ligand_is_proper"], "ligand_id"].dropna().astype(str)
    )
    return sorted(nodes), systems


def component_reduction_dir(*, data_dir: Path, metric: str, source_path: Path) -> Path:
    """Return a stable local reduction path without embedding an absolute root."""
    try:
        source_name = source_path.resolve().relative_to(data_dir.resolve()).as_posix()
    except ValueError:
        source_name = source_path.name
    source_key = hashlib.sha256(source_name.encode()).hexdigest()[:16]
    return (
        data_dir
        / "ligand_clusters"
        / "reductions"
        / f"metric={metric}"
        / f"source={source_key}"
    )


def score_component_reduction_is_complete(
    *,
    data_dir: Path,
    metric: str,
    thresholds: Sequence[int],
    source_path: Path,
    all_nodes: Sequence[str] | None = None,
    eligible_systems: set[str] | None = None,
) -> bool:
    """Return whether one source/metric reduction matches all current inputs."""
    if all_nodes is None:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
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
            ),
            source_signature=_component_source_signature(source_path),
            node_hash=node_hash,
            selection_hash=selection_hash,
            thresholds=sorted(set(thresholds), reverse=True),
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
) -> dict[str, Any]:
    """Map one final score shard into resumable component reductions."""
    if all_nodes is None:
        nodes, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
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
) -> dict[bool, dict[int, pd.DataFrame]]:
    """Reduce all expected score-shard maps and publish component labels."""
    started = time()
    sources = component_score_sources(data_dir=data_dir, metric=metric)
    if not sources:
        raise FileNotFoundError(f"no component score sources found for {metric}")
    nodes, eligible_systems = component_node_universe(
        data_dir=data_dir,
        metric=metric,
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
            )
            for source in sources
        ],
        all_nodes=nodes,
        thresholds=thresholds,
        forest_fan_in=forest_fan_in,
        parquet_batch_size=parquet_batch_size,
        selection_hash=selection_hash,
    )
    for directed in [False, True]:
        for threshold, label_frame in labels[directed].items():
            phase_started = time()
            internal = (
                data_dir
                / "ligand_clusters"
                / "reductions"
                / f"metric={metric}"
                / "labels"
                / f"directed={str(directed).lower()}"
                / f"threshold={threshold}.parquet"
            )
            internal.parent.mkdir(exist_ok=True, parents=True)
            temporary_internal = internal.with_suffix(".tmp.parquet")
            label_frame.to_parquet(
                temporary_internal,
                index=False,
                compression="zstd",
            )
            temporary_internal.replace(internal)
            output = (
                data_dir
                / "ligand_clusters"
                / "cluster=components"
                / f"directed={directed}"
                / f"metric={metric}"
                / f"threshold={threshold}.parquet"
            )
            published = label_frame.copy()
            if metric == "tanimoto_similarity_ecfp4_1024":
                published = expand_fingerprint_clusters_to_ligands(
                    data_dir=data_dir,
                    labeldf=published,
                )
            published["metric"] = metric
            published["directed"] = directed
            published["threshold"] = threshold
            published["cluster"] = "components"
            output.parent.mkdir(exist_ok=True, parents=True)
            temporary = output.with_suffix(".tmp.parquet")
            published.to_parquet(
                temporary,
                schema=LIGAND_CLUSTER_SCHEMA,
                index=False,
            )
            temporary.replace(output)
            LOG.info(
                "component merge published: metric=%s directed=%s threshold=%d "
                "ligands=%d elapsed_seconds=%.1f",
                metric,
                directed,
                threshold,
                len(published),
                time() - phase_started,
            )
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


def expand_fingerprint_clusters_to_ligands(
    *,
    data_dir: Path,
    labeldf: pd.DataFrame,
) -> pd.DataFrame:
    ligands_per_smiles = pd.read_parquet(
        data_dir / "fingerprints/ligands_per_smiles.parquet",
        columns=["ligand_rdkit_canonical_smiles", "ligand_smiles_id"],
    )
    annotation_df = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
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
            ligands_per_smiles["ligand_smiles_id"],
        )
    )
    annotation_df["ligand_smiles_id"] = annotation_df[
        "ligand_rdkit_canonical_smiles"
    ].map(smiles_to_id)
    annotation_df.dropna(subset=["ligand_smiles_id"], inplace=True)
    annotation_df["ligand_smiles_id"] = annotation_df["ligand_smiles_id"].astype(int)
    smiles_to_ligands = {
        int(smiles_id): tuple(sorted(set(group["ligand_id"].astype(str))))
        for smiles_id, group in annotation_df.groupby("ligand_smiles_id")
    }
    labeldf["ligand_id"] = labeldf["ligand_id"].astype(int).map(smiles_to_ligands)
    return (
        labeldf.dropna(subset=["ligand_id"])
        .explode("ligand_id")
        .drop_duplicates(subset=["ligand_id"], keep="first")
    )


def _eligible_annotation(data_dir: Path) -> pd.DataFrame:
    columns = [
        "system_id",
        "ligand_id",
        "system_type",
        "system_num_protein_chains",
        "system_num_ligand_chains",
        "ligand_is_proper",
    ]
    annotation = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
        columns=columns,
    )
    return annotation[
        (annotation["system_type"] == "holo")
        & (annotation["system_num_protein_chains"] <= 5)
        & (annotation["system_num_ligand_chains"] <= 5)
    ]


def make_communities(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
    skip_existing_clusters: bool = False,
    scratch_dir: Path | None = None,
    threads: int = 1,
) -> None:
    """Stream one metric into independent weak-component PLM graphs."""
    started = time()
    if metric in GATED_LIGAND_DIAGNOSTIC_METRICS:
        raise ValueError(
            f"{metric} is evaluated only for ligand pairs with positive pocket "
            "coverage and cannot be clustered directly; use "
            "sucos_shape_pocket_qcov"
        )
    if threads < 1:
        raise ValueError("community threads must be positive")
    output = (
        data_dir
        / "ligand_clusters"
        / "cluster=communities"
        / "directed=False"
        / f"metric={metric}"
        / f"threshold={threshold}.parquet"
    )
    if skip_existing_clusters and output.is_file():
        return
    component_path = (
        data_dir
        / "ligand_clusters"
        / "reductions"
        / f"metric={metric}"
        / "labels"
        / "directed=false"
        / f"threshold={threshold}.parquet"
    )
    if not component_path.is_file():
        raise FileNotFoundError(
            f"weak component labels must be merged before communities: "
            f"{component_path}"
        )
    component_labels = pd.read_parquet(component_path)
    component_labels["ligand_id"] = component_labels["ligand_id"].astype(str)
    if component_labels["ligand_id"].duplicated().any():
        raise ValueError(f"duplicate ligand IDs in component labels: {component_path}")
    component_labels["label"] = component_labels["label"].astype(str)
    component_labels["component"] = pd.factorize(component_labels["label"], sort=False)[
        0
    ].astype(np.uint32)
    component_labels = component_labels.sort_values(
        ["component", "ligand_id"], kind="stable"
    ).reset_index(drop=True)
    component_labels["component_node"] = (
        component_labels.groupby("component", sort=False).cumcount().astype(np.uint32)
    )
    nodes_by_component = {
        int(component): group["ligand_id"].tolist()
        for component, group in component_labels.groupby("component", sort=False)
    }
    component_lookup = component_labels[["ligand_id", "component", "component_node"]]
    sources = component_score_sources(data_dir=data_dir, metric=metric)
    if not sources:
        raise FileNotFoundError(f"no community score sources found for {metric}")

    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    temporary_root = scratch_dir or data_dir / "scratch" / "duckdb" / "communities"
    temporary_root.mkdir(exist_ok=True, parents=True)
    connection.sql(f"SET temp_directory='{temporary_root.as_posix()}'")
    connection.register("component_labels", component_lookup)
    paths_sql = ", ".join(
        f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'" for path in sources
    )
    chemical = metric == "tanimoto_similarity_ecfp4_1024"
    if chemical:
        selected_sql = f"""
            SELECT
                cast(query_ligand_id AS VARCHAR) AS query_node,
                cast(target_ligand_id AS VARCHAR) AS target_node,
                cast({metric} AS DOUBLE) AS similarity
            FROM read_parquet([{paths_sql}])
            WHERE {metric} >= {threshold}
        """
    else:
        _, eligible_systems = component_node_universe(
            data_dir=data_dir,
            metric=metric,
        )
        systems = pd.DataFrame({"system_id": sorted(eligible_systems or set())})
        connection.register("eligible_systems", systems)
        escaped_metric = metric.replace("'", "''")
        selected_sql = f"""
            SELECT
                cast(scores.query_ligand_id AS VARCHAR) AS query_node,
                cast(scores.target_ligand_id AS VARCHAR) AS target_node,
                cast(scores.similarity AS DOUBLE) AS similarity
            FROM read_parquet([{paths_sql}]) AS scores
            INNER JOIN eligible_systems AS query_system
                ON cast(scores.query_system AS VARCHAR) = query_system.system_id
            INNER JOIN eligible_systems AS target_system
                ON cast(scores.target_system AS VARCHAR) = target_system.system_id
            WHERE scores.metric = '{escaped_metric}'
              AND scores.similarity >= {threshold}
        """
    query = dedent(
        f"""
        WITH selected AS (
            {selected_sql}
        ), labeled AS (
            SELECT
                query_labels.component AS query_component,
                target_labels.component AS target_component,
                least(
                    query_labels.component_node,
                    target_labels.component_node
                )::UINTEGER AS query_node,
                greatest(
                    query_labels.component_node,
                    target_labels.component_node
                )::UINTEGER AS target_node,
                selected.similarity,
                query_labels.component != target_labels.component AS crossing_edge
            FROM selected
            INNER JOIN component_labels AS query_labels
                ON selected.query_node = query_labels.ligand_id
            INNER JOIN component_labels AS target_labels
                ON selected.target_node = target_labels.ligand_id
        ), symmetrized_edges AS (
            SELECT
                query_component,
                target_component,
                query_node,
                target_node,
                crossing_edge,
                max(similarity)::DOUBLE AS similarity
            FROM labeled
            WHERE crossing_edge OR query_node != target_node
            GROUP BY
                query_component,
                target_component,
                query_node,
                target_node,
                crossing_edge
        ), weighted AS (
            SELECT
                symmetrized_edges.*,
                sum(similarity) FILTER (
                    WHERE NOT crossing_edge
                ) OVER ()::DOUBLE AS total_similarity
            FROM symmetrized_edges
        )
        SELECT
            query_component AS component,
            query_node,
            target_node,
            similarity,
            total_similarity,
            crossing_edge
        FROM weighted
        ORDER BY crossing_edge, component, query_node, target_node
        """
    )
    LOG.info(
        "community edge query start: metric=%s threshold=%d sources=%d "
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
        "community edge query ready: metric=%s threshold=%d elapsed_seconds=%.1f",
        metric,
        threshold,
        time() - query_started,
    )
    nk.setNumberOfThreads(threads)
    community_groups: list[list[str]] = []
    processed_nodes: set[str] = set()
    current_component: int | None = None
    current_nodes: list[str] = []
    current_graph: nk.graph.Graph | None = None
    current_similarity = 0.0
    total_similarity = 0.0

    def finish_component() -> None:
        nonlocal current_component, current_graph, current_similarity
        if current_component is None or current_graph is None:
            return
        if len(current_nodes) > 1 and current_graph.numberOfEdges() == 0:
            raise ValueError(
                f"non-singleton weak component {current_component} has no edges"
            )
        if current_graph.numberOfEdges() == 0:
            groups = [[node] for node in current_nodes]
        else:
            gamma = current_similarity / total_similarity
            communities = nk.community.detectCommunities(
                current_graph,
                nk.community.PLM(current_graph, gamma=gamma),
            )
            groups = [
                [current_nodes[node] for node in communities.getMembers(index)]
                for index in range(communities.numberOfSubsets())
            ]
        community_groups.extend(groups)
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
                f"weak component validation failed for {metric} at {threshold}: "
                f"{crossing_edges} crossing edges"
            )
        observed_total = frame["total_similarity"].dropna()
        if not observed_total.empty:
            total_similarity = float(observed_total.iloc[0])
        frame = frame.dropna(subset=["component"])
        streamed_rows += len(frame)
        if batch_index % 20 == 0:
            elapsed = time() - stream_started
            LOG.info(
                "community edge stream progress: metric=%s threshold=%d "
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
                current_similarity = 0.0
            assert current_graph is not None
            query_nodes = group["query_node"].to_numpy(dtype=np.uint, copy=False)
            target_nodes = group["target_node"].to_numpy(dtype=np.uint, copy=False)
            similarities = group["similarity"].to_numpy(dtype=float, copy=False)
            current_graph.addEdges((similarities / 100.0, (query_nodes, target_nodes)))
            current_similarity += float(similarities.sum())
    finish_component()
    connection.close()
    LOG.info(
        "community edge stream complete: metric=%s threshold=%d rows=%d "
        "elapsed_seconds=%.1f",
        metric,
        threshold,
        streamed_rows,
        time() - stream_started,
    )

    unprocessed = set(component_labels["ligand_id"]) - processed_nodes
    for component, nodes in nodes_by_component.items():
        missing = [node for node in nodes if node in unprocessed]
        if len(missing) > 1:
            raise ValueError(
                f"non-singleton weak component {component} was absent from edges"
            )
        community_groups.extend([[node] for node in missing])
    ordered_groups = sorted(
        community_groups,
        key=lambda group: (-len(group), min(group)),
    )
    published = pd.DataFrame(
        [
            {"ligand_id": node, "label": f"c{community_index}"}
            for community_index, group in enumerate(ordered_groups)
            for node in sorted(group)
        ]
    )
    if chemical:
        published = expand_fingerprint_clusters_to_ligands(
            data_dir=data_dir,
            labeldf=published,
        )
    published["metric"] = metric
    published["directed"] = False
    published["threshold"] = threshold
    published["cluster"] = "communities"
    output.parent.mkdir(exist_ok=True, parents=True)
    temporary = output.with_suffix(".tmp.parquet")
    published.to_parquet(
        temporary,
        schema=LIGAND_CLUSTER_SCHEMA,
        index=False,
    )
    temporary.replace(output)
