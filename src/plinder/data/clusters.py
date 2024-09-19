# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from time import time
from typing import Callable, TypeVar

import networkit as nk
import numpy as np
import pandas as pd

from plinder.core import scores
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import CLUSTER_SCHEMA

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
    node_labels, max_cluster = clustering_fn(graph, directed)
    labels: list[dict[str, str]] = [
        {
            "system_id": system_ids_cat.categories[sys_int_id],
            "label": cluster_id,
        }
        for sys_int_id, cluster_id in node_labels
    ]
    # Add singletons
    singletons = set(system_ids_and_singletons) - {
        system_ids_cat.categories[node] for node in graph.iterNodes()
    }
    for i, system_id in enumerate(singletons):
        labels.append(
            {
                "system_id": system_id,
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
) -> None:
    cluster_file = (
        data_dir
        / "clusters"
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
    )
    if len(labels):
        labeldf = pd.DataFrame(labels)
        labeldf["metric"] = metric
        labeldf["directed"] = False
        labeldf["threshold"] = threshold
        labeldf["cluster"] = cluster
        if metric == "tanimoto_similarity_max":
            labeldf = explode_ligand_clusters(data_dir=data_dir, labeldf=labeldf)
        LOG.info(f"saving {cluster_file}")
        t0 = time()
        cluster_file.parent.mkdir(exist_ok=True, parents=True)
        labeldf.to_parquet(cluster_file, schema=CLUSTER_SCHEMA)
        t1 = time()
        LOG.info(f"make_cluster_file: saving took {t1-t0:.2f}s")


def prepare_df_protein(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
) -> pd.DataFrame:
    LOG.info(f"threshold={threshold} metric={metric} getting system_ids")
    t0 = time()
    system_ids_and_singletons = set(
        pd.read_parquet(
            data_dir / "index" / "annotation_table.parquet",
            columns=["system_id"],
            filters=[
                ("system_type", "==", "holo"),
                ("system_num_protein_chains", "<=", 5),
                ("system_num_ligand_chains", "<=", 5),
            ],
        )["system_id"]
    )
    t1 = time()
    LOG.info(f"getting {len(system_ids_and_singletons)} system_ids took {t1-t0:.2f}s")
    if not len(system_ids_and_singletons):
        LOG.info("no system_ids found, returning")
        return
    df = scores.query_protein_similarity(
        search_db="holo",
        columns=[
            "query_system",
            "target_system",
            "similarity",
        ],
        filters=[
            ("similarity", ">=", threshold),
            ("metric", "==", metric),
        ],
    )
    if df is None or df.empty:
        LOG.info("no protein similarity scores found, returning")
        return
    LOG.info(f"found {len(df.index)} similarity scores")
    return df, system_ids_and_singletons


def prepare_df_ligand(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
) -> pd.DataFrame:
    LOG.info(f"threshold={threshold} metric={metric} getting ligand_ids")
    t0 = time()
    system_ids_and_singletons = set(
        pd.read_parquet(data_dir / "fingerprints/ligands_per_system.parquet")[
            "number_id_by_inchikeys"
        ].astype(str)
    )
    t1 = time()
    LOG.info(f"getting {len(system_ids_and_singletons)} ligand_ids took {t1-t0:.2f}s")
    if not len(system_ids_and_singletons):
        LOG.info("no ligand_ids found, returning")
        return
    df = scores.query_ligand_similarity(
        columns=[
            "query_ligand_id",
            "target_ligand_id",
            "tanimoto_similarity_max",
        ],
        filters=[
            ("tanimoto_similarity_max", ">=", threshold),
        ],
    )
    if df is None or df.empty:
        LOG.info("no ligand similarity scores found, returning")
        return
    df[["query_ligand_id", "target_ligand_id"]] = df[
        ["query_ligand_id", "target_ligand_id"]
    ].astype(str)
    LOG.info(f"found {len(df.index)} similarity scores")
    df.rename(
        columns={
            "query_ligand_id": "query_system",
            "target_ligand_id": "target_system",
            "tanimoto_similarity_max": "similarity",
        },
        inplace=True,
    )
    return df, system_ids_and_singletons


def make_components_and_communities(
    *,
    data_dir: Path,
    metric: str,
    threshold: int,
    skip_existing_clusters: bool = False,
) -> None:
    if metric == "tanimoto_similarity_max":
        df, system_ids_and_singletons = prepare_df_ligand(
            data_dir=data_dir,
            metric=metric,
            threshold=threshold,
        )
    else:
        df, system_ids_and_singletons = prepare_df_protein(
            data_dir=data_dir,
            metric=metric,
            threshold=threshold,
        )
    system_ids_df = set(df["query_system"]).union(set(df["target_system"]))
    system_ids_cat = pd.CategoricalDtype(categories=list(system_ids_df))
    df[["query_system", "target_system"]] = df[
        ["query_system", "target_system"]
    ].astype(system_ids_cat)
    directed_graph, system_ids_directed_cat = make_nk_graph(
        df, len(system_ids_df), system_ids_cat, directed=True, weighted=False
    )
    make_cluster_file(
        graph=directed_graph,
        system_ids_cat=system_ids_directed_cat,
        directed=True,
        cluster="components",
        system_ids_and_singletons=system_ids_and_singletons,
        data_dir=data_dir,
        metric=metric,
        threshold=threshold,
        skip_existing_clusters=skip_existing_clusters,
    )
    del directed_graph
    del system_ids_directed_cat
    undirected_graph, system_ids_undirected_cat = make_nk_graph(
        df, len(system_ids_df), system_ids_cat, directed=False, weighted=True
    )
    for cluster in ["components", "communities"]:
        make_cluster_file(
            graph=undirected_graph,
            system_ids_cat=system_ids_undirected_cat,
            directed=False,
            cluster=cluster,
            system_ids_and_singletons=system_ids_and_singletons,
            data_dir=data_dir,
            metric=metric,
            threshold=threshold,
            skip_existing_clusters=skip_existing_clusters,
        )
