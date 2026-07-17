# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os
import sys
from pathlib import Path
from time import time
from typing import Any, Callable, TypeVar, cast

if sys.platform == "darwin":
    # For macOS only: allow multiple OpenMP runtimes to coexist
    # (needed on macOS with conda)
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import networkit as nk
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from plinder.core.scores.metrics import (
    GATED_LIGAND_DIAGNOSTIC_METRICS,
    is_ligand_level_metric,
)
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import CLUSTER_SCHEMA, LIGAND_CLUSTER_SCHEMA

LOG = setup_logger(__name__)

T = TypeVar("T")


def make_nk_components(
    graph: nk.graph.Graph,
    is_directed: bool,
) -> tuple[list[tuple[int, str]], int]:
    """
    Get connected component clusters

    Parameters
    ----------
    graph : nk.graph.Graph
        Input graph
    is_directed : bool
        is directed graph

    Returns
    -------
    tuple[list[tuple[int, str]], int]
    """
    if is_directed:
        cc = nk.components.StronglyConnectedComponents(graph)
    else:
        cc = nk.components.ConnectedComponents(graph)
    cc.run()
    components = cc.getComponents()
    return (
        [
            (node, f"c{idx}")
            for idx, component in enumerate(sorted(components, key=len, reverse=True))
            for node in component
        ],
        len(components),
    )


def make_nk_communities(
    graph: nk.graph.Graph, directed: bool
) -> tuple[list[tuple[int, str]], int]:
    """
    Get community clusters

    Parameters
    ----------
    graph : nk.graph.Graph

    Returns
    -------
    tuple[list[tuple[int, str]], int]
    """
    assert not directed
    if sys.platform == "darwin":
        # For macOS only: limit to 1 thread to avoid segfault in PLM with multiple OMP runtimes
        nk.setNumberOfThreads(1)
    communities = nk.community.detectCommunities(graph, nk.community.PLM(graph))
    community_list = [
        communities.getMembers(i) for i in range(communities.numberOfSubsets())
    ]
    return (
        [
            (node, f"c{idx}")
            for idx, component in enumerate(
                sorted(community_list, key=len, reverse=True)
            )
            for node in component
        ],
        len(community_list),
    )


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


def get_labels(
    graph: nk.graph.Graph,
    system_ids_and_singletons: set[str],
    system_ids_cat: pd.CategoricalDtype,
    clustering_fn: Callable[[nk.graph.Graph, bool], tuple[list[tuple[int, str]], int]],
    directed: bool,
    node_column: str = "system_id",
) -> list[dict[str, str]]:
    """
    Get community clusters

    Parameters
    ----------
    graph : nk.graph.Graph
        Input graph
    system_ids_and_singletons : set[str]
        All system ids, both the ones in graph \
            and singletons not in graph
    system_ids_cat :  pd.CategoricalDtype
        system ids categories
    clustering_fn : Callable
        clustering function, could be a
        component or community detector
    directed : bool
        Is directed?

    Returns
    -------
    list[dict[str, str]]
    """
    if graph.numberOfNodes() == 0:
        node_labels: list[tuple[int, str]] = []
        max_cluster = 0
    else:
        node_labels, max_cluster = clustering_fn(graph, directed)
    labels: list[dict[str, str]] = [
        {
            node_column: system_ids_cat.categories[sys_int_id],
            "label": cluster_id,
        }
        for sys_int_id, cluster_id in node_labels
    ]
    # Add singletons
    singletons = set(system_ids_and_singletons) - {
        system_ids_cat.categories[node] for node in graph.iterNodes()
    }
    for i, system_id in enumerate(sorted(singletons)):
        labels.append(
            {
                node_column: system_id,
                "label": f"c{max_cluster + i}",
            }
        )
    return labels


def explode_ligand_clusters(
    *,
    data_dir: Path,
    labeldf: pd.DataFrame,
) -> pd.DataFrame:
    ligands_per_system = pd.read_parquet(
        data_dir / "fingerprints/ligands_per_system.parquet"
    )
    annotation_df = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
        columns=[
            "system_id",
            "ligand_molecular_weight",
            "ligand_rdkit_canonical_smiles",
        ],
        filters=[
            ("system_id", "in", set(ligands_per_system["system_id"])),
            ("ligand_is_ion", "==", False),
            ("ligand_is_artifact", "==", False),
        ],
    )
    mapr = dict(
        zip(
            ligands_per_system["ligand_rdkit_canonical_smiles"],
            ligands_per_system["number_id_by_inchikeys"],
        )
    )
    annotation_df["number_id_by_inchikeys"] = annotation_df[
        "ligand_rdkit_canonical_smiles"
    ].map(mapr)
    annotation_df.dropna(subset=["number_id_by_inchikeys"], inplace=True)
    annotation_df["number_id_by_inchikeys"] = annotation_df[
        "number_id_by_inchikeys"
    ].astype(int)
    annotation_df = annotation_df.sort_values(
        by=["system_id", "ligand_molecular_weight"], ascending=[True, False]
    ).drop_duplicates(subset=["system_id"], keep="first")
    ligand_to_system: dict[int, set[str]] = {}
    for ligand_id, group in annotation_df.groupby("number_id_by_inchikeys"):
        ligand_to_system[int(ligand_id)] = set(group["system_id"])
    labeldf["system_id"] = labeldf["system_id"].astype(int).map(ligand_to_system)
    labeldf = labeldf.dropna(subset=["system_id"]).explode("system_id")
    return labeldf


def make_cluster_file(
    *,
    graph: nk.graph.Graph,
    system_ids_and_singletons: set[str],
    system_ids_cat: pd.CategoricalDtype,
    data_dir: Path,
    metric: str,
    threshold: int,
    directed: bool,
    cluster: str,
    skip_existing_clusters: bool,
    node_column: str = "system_id",
    output_root: str = "clusters",
    explode_fingerprint_ligands: bool = False,
) -> None:
    cluster_file = (
        data_dir
        / output_root
        / f"cluster={cluster}"
        / f"directed={directed}"
        / f"metric={metric}"
        / f"threshold={threshold}.parquet"
    )
    LOG.info(f"make_cluster_file: {cluster_file}")
    if cluster_file.is_file() and skip_existing_clusters:
        LOG.info(f"skipping {cluster_file} because it exists")
        return
    func = make_nk_components if cluster == "components" else make_nk_communities
    labels = get_labels(
        graph,
        system_ids_and_singletons,
        system_ids_cat,
        func,
        directed,
        node_column,
    )
    if len(labels):
        labeldf = pd.DataFrame(labels)
        labeldf["metric"] = metric
        labeldf["directed"] = directed
        labeldf["threshold"] = threshold
        labeldf["cluster"] = cluster
        if explode_fingerprint_ligands:
            labeldf = explode_ligand_clusters(data_dir=data_dir, labeldf=labeldf)
        LOG.info(f"saving {cluster_file}")
        t0 = time()
        cluster_file.parent.mkdir(exist_ok=True, parents=True)
        schema = CLUSTER_SCHEMA if node_column == "system_id" else LIGAND_CLUSTER_SCHEMA
        labeldf.to_parquet(cluster_file, schema=schema, index=False)
        t1 = time()
        LOG.info(f"make_cluster_file: saving took {t1 - t0:.2f}s")


def _read_local_score_rows(
    *, data_dir: Path, metric: str, threshold: int
) -> pd.DataFrame:
    """Read one metric from score files generated in this ingest directory."""
    score_paths = sorted((data_dir / "scores" / "search_db=holo").rglob("*.parquet"))
    frames: list[pd.DataFrame] = []
    wanted = [
        "query_system",
        "query_ligand_id",
        "target_system",
        "target_ligand_id",
        "metric",
        "similarity",
    ]
    for path in score_paths:
        read_schema = cast(Any, pq.read_schema)
        available = set(read_schema(path).names)
        if not {"query_system", "target_system", "metric", "similarity"}.issubset(
            available
        ):
            continue
        frame = pd.read_parquet(
            path,
            columns=[column for column in wanted if column in available],
            filters=[("metric", "==", metric), ("similarity", ">=", threshold)],
        )
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=wanted)
    return pd.concat(frames, ignore_index=True)


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


def _collapse_edges(
    rows: pd.DataFrame, *, query_column: str, target_column: str
) -> pd.DataFrame:
    if rows.empty or not {query_column, target_column}.issubset(rows.columns):
        return pd.DataFrame(columns=["query_node", "target_node", "similarity"])
    edges = rows.dropna(subset=[query_column, target_column])[
        [query_column, target_column, "similarity"]
    ].copy()
    if edges.empty:
        return pd.DataFrame(columns=["query_node", "target_node", "similarity"])
    edges[[query_column, target_column]] = edges[[query_column, target_column]].astype(
        str
    )
    return (
        edges.groupby([query_column, target_column], as_index=False)["similarity"]
        .max()
        .rename(columns={query_column: "query_node", target_column: "target_node"})
    )


def prepare_score_edges(
    *, data_dir: Path, metric: str, threshold: int
) -> tuple[pd.DataFrame, set[str], pd.DataFrame | None, set[str]]:
    """Prepare system edges and, for ligand metrics, ligand edges locally."""
    annotation = _eligible_annotation(data_dir)
    system_nodes = set(annotation["system_id"].dropna().astype(str))
    ligand_nodes = set(
        annotation.loc[annotation["ligand_is_proper"], "ligand_id"].dropna().astype(str)
    )
    rows = _read_local_score_rows(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
    )
    rows = rows[
        rows["query_system"].astype(str).isin(system_nodes)
        & rows["target_system"].astype(str).isin(system_nodes)
    ]
    ligand_metric = is_ligand_level_metric(metric)
    has_ligand_ids = {"query_ligand_id", "target_ligand_id"}.issubset(rows.columns)
    if ligand_metric and has_ligand_ids:
        rows = rows[
            rows["query_ligand_id"].astype(str).isin(ligand_nodes)
            & rows["target_ligand_id"].astype(str).isin(ligand_nodes)
        ]
    system_edges = _collapse_edges(
        rows,
        query_column="query_system",
        target_column="target_system",
    )
    ligand_edges: pd.DataFrame | None = None
    if ligand_metric:
        if has_ligand_ids:
            ligand_edges = _collapse_edges(
                rows,
                query_column="query_ligand_id",
                target_column="target_ligand_id",
            )
        elif not rows.empty:
            LOG.warning(
                f"metric={metric} has no ligand identifiers; "
                "creating only legacy system-level clusters"
            )
    LOG.info(
        f"metric={metric} threshold={threshold}: {len(system_edges)} system edges, "
        f"{0 if ligand_edges is None else len(ligand_edges)} ligand edges"
    )
    return system_edges, system_nodes, ligand_edges, ligand_nodes


def prepare_df_ligand(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
) -> tuple[pd.DataFrame, set[str]]:
    LOG.info(f"threshold={threshold} metric={metric} getting ligand_ids")
    t0 = time()
    system_ids_and_singletons = set(
        pd.read_parquet(data_dir / "fingerprints/ligands_per_system.parquet")[
            "number_id_by_inchikeys"
        ].astype(str)
    )
    t1 = time()
    LOG.info(f"getting {len(system_ids_and_singletons)} ligand_ids took {t1 - t0:.2f}s")
    if not len(system_ids_and_singletons):
        LOG.info("no ligand_ids found, returning")
        return (
            pd.DataFrame(columns=["query_node", "target_node", "similarity"]),
            set(),
        )
    frames = []
    for path in sorted((data_dir / "ligand_scores").glob("*.parquet")):
        frame = pd.read_parquet(
            path,
            columns=[
                "query_ligand_id",
                "target_ligand_id",
                "tanimoto_similarity_max",
            ],
            filters=[("tanimoto_similarity_max", ">=", threshold)],
        )
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return (
            pd.DataFrame(columns=["query_node", "target_node", "similarity"]),
            system_ids_and_singletons,
        )
    df = pd.concat(frames, ignore_index=True)
    df[["query_ligand_id", "target_ligand_id"]] = df[
        ["query_ligand_id", "target_ligand_id"]
    ].astype(str)
    LOG.info(f"found {len(df.index)} similarity scores")
    df.rename(
        columns={
            "query_ligand_id": "query_node",
            "target_ligand_id": "target_node",
            "tanimoto_similarity_max": "similarity",
        },
        inplace=True,
    )
    return df, system_ids_and_singletons


def _make_clusters_for_edges(
    *,
    data_dir: Path,
    edges: pd.DataFrame,
    all_nodes: set[str],
    metric: str,
    threshold: int,
    skip_existing_clusters: bool,
    node_column: str,
    output_root: str,
    explode_fingerprint_ligands: bool = False,
) -> None:
    edge_nodes = sorted(set(edges["query_node"]).union(edges["target_node"]))
    node_categories = pd.CategoricalDtype(categories=edge_nodes)
    directed_graph, directed_categories = make_nk_graph(
        edges.copy(),
        len(edge_nodes),
        node_categories,
        directed=True,
        weighted=False,
        query_col="query_node",
        target_col="target_node",
    )
    make_cluster_file(
        graph=directed_graph,
        system_ids_cat=directed_categories,
        directed=True,
        cluster="components",
        system_ids_and_singletons=all_nodes,
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        skip_existing_clusters=skip_existing_clusters,
        node_column=node_column,
        output_root=output_root,
        explode_fingerprint_ligands=explode_fingerprint_ligands,
    )
    undirected_graph, undirected_categories = make_nk_graph(
        edges.copy(),
        len(edge_nodes),
        node_categories,
        directed=False,
        weighted=True,
        query_col="query_node",
        target_col="target_node",
    )
    for cluster in ["components", "communities"]:
        make_cluster_file(
            graph=undirected_graph,
            system_ids_cat=undirected_categories,
            directed=False,
            cluster=cluster,
            system_ids_and_singletons=all_nodes,
            data_dir=data_dir,
            metric=metric,
            threshold=threshold,
            skip_existing_clusters=skip_existing_clusters,
            node_column=node_column,
            output_root=output_root,
            explode_fingerprint_ligands=explode_fingerprint_ligands,
        )


def make_components_and_communities(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
    skip_existing_clusters: bool = False,
) -> None:
    if metric in GATED_LIGAND_DIAGNOSTIC_METRICS:
        raise ValueError(
            f"{metric} is evaluated only for ligand pairs with positive pocket "
            "coverage and cannot be clustered directly; use "
            "sucos_shape_pocket_qcov"
        )
    if metric == "tanimoto_similarity_max":
        edges, nodes = prepare_df_ligand(
            data_dir=data_dir,
            metric=metric,
            threshold=threshold,
        )
        _make_clusters_for_edges(
            data_dir=data_dir,
            edges=edges,
            all_nodes=nodes,
            metric=metric,
            threshold=threshold,
            skip_existing_clusters=skip_existing_clusters,
            node_column="system_id",
            output_root="clusters",
            explode_fingerprint_ligands=True,
        )
        return

    system_edges, system_nodes, ligand_edges, ligand_nodes = prepare_score_edges(
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
    )
    _make_clusters_for_edges(
        data_dir=data_dir,
        edges=system_edges,
        all_nodes=system_nodes,
        metric=metric,
        threshold=threshold,
        skip_existing_clusters=skip_existing_clusters,
        node_column="system_id",
        output_root="clusters",
    )
    if ligand_edges is not None:
        _make_clusters_for_edges(
            data_dir=data_dir,
            edges=ligand_edges,
            all_nodes=ligand_nodes,
            metric=metric,
            threshold=threshold,
            skip_existing_clusters=skip_existing_clusters,
            node_column="ligand_id",
            output_root="ligand_clusters",
        )
