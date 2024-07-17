# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import networkx as nx
import numpy as np
import pandas as pd


def add_singletons(
    key_to_cluster: dict[str, str], max_cluster: int, all_system_ids: set[str]
) -> dict[str, str]:
    singletons = all_system_ids - set(key_to_cluster)
    for i, system_id in enumerate(singletons):
        key_to_cluster[system_id] = f"c{max_cluster + i}"
    return key_to_cluster


def get_key_to_cluster_nx(clusters: list[set[str]]) -> dict[str, str]:
    key_to_cluster = {}
    for i, c in enumerate(clusters):
        for k in c:
            key_to_cluster[k] = f"c{i}"
    return key_to_cluster


def make_communities_nx(
    graph: nx.Graph,
    system_ids: set[str],
) -> dict[str, str]:
    communities = sorted(
        nx.community.asyn_lpa_communities(graph, weight="similarity"),
        key=len,
        reverse=True,
    )
    key_to_community = get_key_to_cluster_nx(communities)
    key_to_community = add_singletons(key_to_community, len(communities), system_ids)
    return key_to_community


def make_components_and_communities_nx(
    df: pd.DataFrame, system_ids: set[str], strong: bool
) -> list[dict[str, str]]:
    graph = nx.from_pandas_edgelist(
        df,
        source="query_system",
        target="target_system",
        edge_key="similarity",
        edge_attr="similarity",
        create_using=nx.DiGraph() if strong else nx.Graph(),
    )
    if strong:
        components = sorted(
            nx.strongly_connected_components(graph), key=len, reverse=True
        )
    else:
        components = sorted(nx.connected_components(graph), key=len, reverse=True)
    key_to_component = get_key_to_cluster_nx(components)
    key_to_component = add_singletons(key_to_component, len(components), system_ids)
    key_to_community = make_communities_nx(graph, system_ids)
    labels = []
    for system_id in system_ids:
        labels.append(
            {
                "system_id": system_id,
                "component": key_to_component[system_id],
                "community": key_to_community[system_id],
            }
        )
    return labels


def make_components_and_communities_gt(
    df: pd.DataFrame, system_ids: set[str], strong: bool
) -> list[dict[str, str]]:
    import graph_tool.all as gt

    graph = gt.Graph(
        list(df.itertuples(index=False, name=None)),
        directed=strong,
        hashed=True,
        hash_type="string",
        eprops=[("similarity", "int")],
    )
    key_to_component, component_sizes = gt.label_components(
        graph, directed=strong, attractors=False
    )
    # sort component name largest first
    component_indices = np.argsort(component_sizes)[::-1]
    mapping = dict({v: i for i, v in enumerate(component_indices)})
    key_to_component = dict(
        zip(graph.vp["ids"], [f"c{mapping[x]}" for x in key_to_component])
    )
    key_to_component = add_singletons(key_to_component, len(mapping), system_ids)
    labels = []
    for system_id in system_ids:
        labels.append(
            {
                "system_id": system_id,
                "component": key_to_component[system_id],
            }
        )
    return labels
