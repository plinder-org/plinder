# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from collections import Counter, defaultdict
from copy import copy
from dataclasses import dataclass, field
from hashlib import md5
from io import StringIO
from json import dumps
from pathlib import Path
from typing import Any

import networkit as nk
import pandas as pd
from omegaconf import DictConfig, ListConfig, OmegaConf
from tqdm import tqdm

from plinder.core import scores
from plinder.core.utils.log import setup_logger
from plinder.data import clusters
from plinder.data.databases import get_ids_in_db

LOG = setup_logger(__name__)


@dataclass
class GraphConfig:
    metric: str
    threshold: int  # edges above this threshold are kept
    depth: int  # neighbors at this depth are counted as leakage


@dataclass
class SplitConfig:
    graph_configs: list[GraphConfig] = field(
        default_factory=lambda: [
            GraphConfig("pli_qcov", 30, 1),
            GraphConfig("pocket_qcov", 50, 1),
        ]
    )
    # how many unique congeneric IDs passing quality to consider as MMS
    mms_unique_quality_count: int = 3
    # what kind of cluster to use for sampling test
    test_cluster_cluster: str = "communities"
    # metric to use for sampling representatives from each test cluster
    test_cluster_metric: str = "pli_qcov"
    # threshold to use for sampling representatives from each test cluster
    test_cluster_threshold: int = 50
    # directed to use for sampling representatives from each test cluster
    test_cluster_directed: bool = False
    # ?
    cluster_column: str = "cluster"
    # max number of representatives from each test cluster
    num_test_representatives: int = 3
    # test should not be singletons
    min_test_cluster_size: int = 5
    # test should not be in too big communities or cause too many train cases to be removed
    max_test_leakage_count: int = 300
    # test should not have too few or too many interactions
    min_max_test_pli: tuple[int, int] = (3, 50)
    # test should not have too few or too many pocket residues
    min_max_test_pocket: tuple[int, int] = (5, 100)
    # fraction of systems to choose for test
    test_fraction: float = 0.01
    # what kind of cluster to use for sampling val
    val_cluster_cluster: str = "components"
    # metric to use for splitting train and val
    val_cluster_metric: str = "pocket_qcov"
    # threshold to use for splitting train and val
    val_cluster_threshold: int = 50
    # directed to use for splitting train and val
    val_cluster_directed: bool = False
    # max number of representatives from each val cluster
    num_val_representatives: int = 3
    # val should not be singletons
    min_val_cluster_size: int = 5
    # val should not have too few or too many interactions
    min_max_val_pli: tuple[int, int] = (2, 50)
    # val should not have too few or too many pocket residues
    min_max_val_pocket: tuple[int, int] = (5, 100)
    # fraction of systems to choose for val
    val_fraction: float = 0.01


_map = {
    "split": SplitConfig,
}


def _clean_to_python(cfg: Any) -> Any:
    """
    Recursively convert an omegaconf object to a plain
    python dict. Account for ListConfig and DictConfig
    in the nesting.

    Parameters
    ----------
    cfg : Any
        omegaconf container or primitive
    """
    if isinstance(cfg, DictConfig):
        return {k: _clean_to_python(v) for k, v in cfg.items()}
    elif isinstance(cfg, ListConfig):
        return [_clean_to_python(v) for v in cfg]
    else:
        return cfg


def get_config_hash(config_obj: Any) -> str:
    """
    Get a unique string representation of the config
    passed to a splitting run.

    Parameters
    ----------
    config_obj : dict | dataclass | DictConfig
        the configuration to hash

    Returns
    -------
    config_hash : str
        the hash of the configuration
    """
    if isinstance(config_obj, dict):
        contents = config_obj
    elif isinstance(config_obj, DictConfig):
        contents = _clean_to_python(config_obj)
    else:
        contents = config_obj.__dict__
    config_hash = md5(dumps(contents, sort_keys=True).encode("utf-8")).hexdigest()
    return config_hash


def get_default_config() -> DictConfig:
    default = DictConfig(
        OmegaConf.merge(
            {
                "split": OmegaConf.structured(SplitConfig()),
            }
        )
    )
    return default


def get_config(config_contents: str) -> DictConfig:
    default = get_default_config()
    from_file = DictConfig(OmegaConf.load(StringIO(config_contents)))
    cfg = OmegaConf.merge(default, from_file)
    return DictConfig({str(k): _map[str(k)](**v) for k, v in cfg.items()})


def find_neighbors_upto_specific_depth(
    graph: nk.graph.Graph, depth: int, target_nodes: set[int]
) -> set[int]:
    """
    Find neighbors of a set of nodes

    Parameters
    ----------
    graph : nk.graph.Graph
        graph of systems
    depth : int
        search distance
    target_nodes: set[int]
        nodes of interest

    Returns
    -------
    Neighbors of all nodes
    """
    nodes_to_deleak = copy(target_nodes)
    neighbors = set()
    for depth in range(depth):
        next_nodes_to_deleak = {
            j for i in nodes_to_deleak for j in graph.iterNeighbors(i)
        } - nodes_to_deleak
        neighbors |= next_nodes_to_deleak
        nodes_to_deleak = next_nodes_to_deleak
    return neighbors


def deleak_entry_nk_single(
    system_id: str,
    cluster_members: set[str],
    graphs: list[nk.graph.Graph],
    graph_configs: list[GraphConfig],
    vertex_mappings: list[dict[str, int]],
    vertex_ids: list[list[str]],
) -> set[str]:
    leaked_systems = set(cluster_members) - {system_id}
    for graph, graph_config, system_id_to_vertex, vertex_ids_x in zip(
        graphs, graph_configs, vertex_mappings, vertex_ids
    ):
        if system_id not in system_id_to_vertex:
            continue
        leaked_neighbors = find_neighbors_upto_specific_depth(
            graph, graph_config.depth, {system_id_to_vertex[system_id]}
        )
        # those neighbors are leaky
        # leaked_systems.update(leaked_neighbors)
        leaked_systems.update({vertex_ids_x[x] for x in leaked_neighbors})
    # the cluster + neighbors needs to be removed if that system is chosen
    # return leaked_systems.union(cluster_members)
    return leaked_systems


def prep_data_for_desired_properties(
    data_dir: Path, cfg: DictConfig
) -> tuple[pd.DataFrame, pd.DataFrame, set[str], dict[str, set[str]]]:
    """
    Deleak a specific fraction of systems id

    Parameters
    ----------
    cfg : DictConfig
        splitting config

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame,  set]
        mms dataframe, entries annotations,  high quality systems
    """
    LOG.info(OmegaConf.to_yaml(cfg))

    entries = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")

    # sanity check to only split entries that were successfully scored
    foldseek_ids = get_ids_in_db(
        data_dir=data_dir, search_db="holo", aln_type="foldseek"
    )
    mmseqs_ids = get_ids_in_db(data_dir=data_dir, search_db="holo", aln_type="mmseqs")
    foldseek_entries_present = foldseek_ids[1].str.replace("pdb_0000", "").str[:4]
    mmseqs_entries_present = mmseqs_ids[1].str.replace("pdb_0000", "").str[:4]
    all_entries_present = set(foldseek_entries_present).intersection(
        mmseqs_entries_present
    )

    entries = entries[
        (entries["system_type"] == "holo")
        & (entries["system_num_interacting_protein_chains"] <= 5)
        & (entries["system_num_ligand_chains"] <= 5)
        & (entries["entry_pdb_id"].isin(all_entries_present))
    ].reset_index(drop=True)
    entries["uniqueness"] = (
        entries["system_id_no_biounit"]
        + "_"
        + entries["pli_qcov__100__strong__component"]
    )
    nonredundant_to_all = defaultdict(set)
    for system_id, unique_id in zip(entries["system_id"], entries["uniqueness"]):
        nonredundant_to_all[unique_id].add(system_id)
    entries = entries.sort_values("system_biounit_id").drop_duplicates("uniqueness")
    all_system_ids = set(entries["system_id"])

    LOG.info(f"loaded {len(all_system_ids)} from annotation table")

    entries["proto_test"] = entries["system_pass_validation_criteria"].fillna(False)
    quality = set(entries[entries["proto_test"]]["system_id"])
    LOG.info(f"{len(quality)} systems are of high quality")

    apo_links = scores.query_protein_similarity(
        search_db="apo",
        columns=["query_system", "target_system"],
        filters=[
            ("similarity", ">=", 95),
            ("metric", "==", "pocket_fident"),
        ],
    )
    assert apo_links is not None, "apo_links is None"
    num_apo_links = apo_links.groupby("query_system").size()
    LOG.info(f"num_apo_links has {len(num_apo_links)} elements")

    pred_links = scores.query_protein_similarity(
        search_db="pred",
        columns=["query_system", "target_system"],
        filters=[
            ("similarity", ">=", 100),
            ("search_db", "==", "pred"),
            ("metric", "==", "pocket_fident"),
        ],
    )
    assert pred_links is not None, "pred_links is None"
    num_pred_links = pred_links.groupby("query_system").size()
    LOG.info(f"num_pred_links has {len(num_pred_links)} elements")

    entries["has_apo"] = entries["system_id"].map(lambda x: num_apo_links.get(x, 0) > 0)
    entries["has_pred"] = entries["system_id"].map(
        lambda x: num_pred_links.get(x, 0) > 0
    )

    mms_df = pd.read_parquet(data_dir / "mmp/plinder_mmp_series.parquet")
    LOG.info(f"read mmp series has {len(mms_df.index)} records")
    mms_count = (
        mms_df[mms_df["system_id"].isin(quality)].groupby("congeneric_id").size()
    )
    LOG.info(f"mmp series after quality control {len(mms_count.index)} records")
    mms_df["mms_unique_quality_count"] = mms_df["congeneric_id"].map(
        lambda x: mms_count.get(x, 0)
    )
    good_mms_systems = set(
        mms_df[
            mms_df["mms_unique_quality_count"] >= cfg.split.mms_unique_quality_count
        ]["system_id"]
    )
    LOG.info(f"good mms systems: {len(good_mms_systems)}")
    entries["system_has_mms"] = entries["system_id"].isin(good_mms_systems)
    LOG.info(
        f"entries with good mms: {len(entries[entries['system_has_mms']].drop_duplicates('system_id').index)}"
    )
    return mms_df, entries, quality, nonredundant_to_all


def load_graphs(
    data_dir: Path,
    systems: set[str],
    cfg: DictConfig,
) -> tuple[list[nk.graph.Graph], list[dict[str, int]], list[list[str]]]:
    """
    Load graph into networkit graph

    Parameters
    ----------
    data_dir : Path,
    graph_config : GraphConfig

    Returns
    -------
    tuple[list[nk.graph.Graph], list[dict[str, int]]]
        similarity graph, vertex mapping
    """
    graphs = []
    vertex_mappings = []
    vertex_ids = []
    graph_dir = data_dir / "graph_queries"
    graph_dir.mkdir(exist_ok=True, parents=True)

    for graph_config in tqdm(cfg.split.graph_configs):
        graph_file = (
            graph_dir / f"{graph_config.metric}_{graph_config.threshold}.parquet"
        )
        if graph_file.is_file():
            df = pd.read_parquet(graph_file)
        else:
            df = scores.query_protein_similarity(
                search_db="holo",
                columns=["query_system", "target_system"],
                filters=[
                    ("similarity", ">=", graph_config.threshold),
                    ("metric", "==", graph_config.metric),
                ],
            )
            df.to_parquet(graph_file, index=False)
        df = df[
            (df["query_system"].isin(systems)) & (df["target_system"].isin(systems))
        ].reset_index(drop=True)
        LOG.info(f"Loaded {graph_config.metric} similarity dataframe")
        system_ids = set(df.query_system.unique()).union(set(df.target_system))
        system_ids_cat = pd.CategoricalDtype(categories=list(system_ids))
        graph, system_ids_cat = clusters.make_nk_graph(
            df, len(system_ids), system_ids_cat, directed=False, weighted=False
        )

        LOG.info(f"Built graph for {graph_config.metric}")
        vertex_mapping = {v: idx for idx, v in enumerate(system_ids_cat.categories)}
        del df
        graphs.append(graph)
        vertex_mappings.append(vertex_mapping)
        vertex_ids.append(list(system_ids_cat.categories))
        LOG.info(f"Built graph for {graph_config.metric}")
        LOG.info(f"Graph has {graph.numberOfNodes()} nodes")
    return graphs, vertex_mappings, vertex_ids


def get_sampling_clusters(
    cfg: DictConfig, data_dir: Path, entries: pd.DataFrame
) -> tuple[dict[str, set[str]], dict[str, str], pd.DataFrame, set[str], dict[str, str]]:
    """
    Get sampling clusters

    Parameters
    ----------
    data_dir : Path
    entries : pd.DataFrame

    Returns
    -------
    tuple[
        dict[str, set], dict[str, str],
        pd.DataFrame, set[str], dict[str, str]]
        (test-sampling-cluster-to-systems dictionary,
        test-systems-to-sampling-cluster dictionary,
        entries annotations,
        test-system-ids,
        validation-system-to-sampling-cluster dictionary
    """
    sampling_cluster_file = (
        data_dir
        / "clusters"
        / f"cluster={cfg.split.test_cluster_cluster}"
        / f"directed={cfg.split.test_cluster_directed}"
        / f"metric={cfg.split.test_cluster_metric}"
        / f"threshold={cfg.split.test_cluster_threshold}.parquet"
    )
    cluster_pq = pd.read_parquet(sampling_cluster_file)
    labels = dict(zip(cluster_pq["system_id"], cluster_pq["label"]))
    LOG.info(f"loaded {len(set(labels.values()))} communities for sampling")
    entries["cluster"] = entries["system_id"].map(labels)

    cluster_to_size = Counter(entries.drop_duplicates("system_id")["cluster"])
    system_to_cluster = dict(zip(entries["system_id"], entries["cluster"]))
    cluster_to_systems = defaultdict(set)
    for s in system_to_cluster:
        cluster_to_systems[system_to_cluster[s]].add(s)

    val_cluster_file = (
        data_dir
        / "clusters"
        / f"cluster={cfg.split.val_cluster_cluster}"
        / f"directed={cfg.split.val_cluster_directed}"
        / f"metric={cfg.split.val_cluster_metric}"
        / f"threshold={cfg.split.val_cluster_threshold}.parquet"
    )
    cluster_pq = pd.read_parquet(val_cluster_file)
    labels = dict(zip(cluster_pq["system_id"], cluster_pq["label"]))
    LOG.info(
        f"loaded {len(set(labels.values()))} communities for splitting train and val"
    )
    entries["cluster_for_val_split"] = entries["system_id"].map(labels)

    system_to_cluster_val = dict(
        zip(entries["system_id"], entries["cluster_for_val_split"])
    )
    cluster_to_systems_val = defaultdict(set)
    for s in system_to_cluster_val:
        cluster_to_systems_val[system_to_cluster_val[s]].add(s)

    entries["split"] = "train"
    entries.loc[entries["proto_test"], "split"] = "proto-test"

    test_systems = set(
        row.system_id
        for _, row in tqdm(
            entries[entries["proto_test"]][
                ["system_id", cfg.split.cluster_column]
            ].iterrows()
        )
        if cluster_to_size[row[cfg.split.cluster_column]]
        >= cfg.split.min_test_cluster_size
    )
    LOG.info(f"found {len(test_systems)} test systems")

    return (
        cluster_to_systems,
        system_to_cluster,
        entries,
        test_systems,
        system_to_cluster_val,
    )


def remove_high_degree_systems(
    graphs: nk.graph.Graph,
    max_test_leakage_count: int,
    test_systems: set[str],
    vertex_mappings: list[dict[str, int]],
) -> set[str]:
    """
    Prune system id list to remove high degree nodes

    Parameters
    ----------
    graphs : nk.graph.Graph
        Similarity graph
    max_test_leakage_count : int
        Maximum leakge allowed in proto-test systems
    test_systems : set[str]
        test systems selected for deleaking
    vertex_mappings : dict[str, str]
        systems id to vertex dictionary

    Returns
    -------
    set[str]
        test systems selected for deleaking
    """
    ignore_systems = set()
    for system_id in tqdm(test_systems):
        if (
            max(
                [
                    graphs[idx].degreeOut(vertex_mapping[system_id])
                    if system_id in vertex_mapping
                    else 0
                    for idx, vertex_mapping in enumerate(vertex_mappings)
                ]
            )
            > max_test_leakage_count
        ):
            ignore_systems.add(system_id)
    test_systems = test_systems.difference(ignore_systems)
    LOG.info(
        f"found {len(test_systems)} test systems after removing by first degree neighbor leakage count"
    )
    return test_systems


def prioritize_test_sample(
    cfg: DictConfig,
    entries: pd.DataFrame,
    system_id_to_leakage: dict[str, set[str]],
    test_systems: set[str],
    mms_df: pd.DataFrame,
    cluster_to_systems: dict[str, set[str]],
    system_to_cluster: dict[str, str],
    quality: set[str],
) -> tuple[dict[str, set[str]], set[str]]:
    """
    Prioritize test systems to ensure selected clusters met a certain criteria.
    - Minimum cluster size and doesn't have more than a \
    specified number of potential leaky systems.
    - We are prioritizing representatives that have have:
       - high apo links
       - has binding affinity
       - high mms links
       - low leaky count

    Parameters
    ----------
    cfg : DictConfig
        Config dictionary
    entries : pd.DataFrame
        entries annotation
    system_id_to_leakage : dict[str, set[str]]
        test system id mapped to leakage
    test_systems : set[str]
        test system ids
    mms_df : pd.DataFrame
        matched-molecular series dataframe
    cluster_to_systems : dict[str, set[str]]
        cluster ids mapped to sets of  system ids
    system_to_cluster : dict[str, str]
         system id mapped to cluster ids
    quality : set[set]
         set of high quality test systems

    Returns
    -------
    tuple[dict[str, set[str]], set[str]]
        new_system_id_cluster_dict, test_system_ids
        system ids mapped to cluster of nodes marked for removal, test system ids
    """
    entries["leakage_count"] = entries["system_id"].map(
        lambda x: len(
            system_id_to_leakage.get(x, cluster_to_systems[system_to_cluster[x]])
        )
    )
    test_data = entries[entries["system_id"].isin(test_systems)].drop_duplicates(
        "system_id"
    )

    test_data = test_data[
        (test_data["leakage_count"] >= cfg.split.min_test_cluster_size)
        & (test_data["leakage_count"] <= cfg.split.max_test_leakage_count)
        & (test_data["system_num_pocket_residues"] >= cfg.split.min_max_test_pocket[0])
        & (test_data["system_num_pocket_residues"] <= cfg.split.min_max_test_pocket[1])
        & (test_data["system_num_interactions"] >= cfg.split.min_max_test_pli[0])
        & (test_data["system_num_interactions"] <= cfg.split.min_max_test_pli[1])
    ]
    LOG.info(
        f"found {len(set(test_data['system_id']))} test systems after removing by second pass neighbor leakage count"
    )
    test_data = (
        test_data.sort_values(
            [
                "system_has_mms",
                "has_apo",
                "system_has_binding_affinity",
                "leakage_count",
            ],
            ascending=[False, False, False, True],
        )
        .drop_duplicates("system_id")
        .groupby(cfg.split.cluster_column)
        .head(cfg.split.num_test_representatives)
    )
    LOG.info(
        f"test_data after groupby.head of num_test_representatives: {len(set(test_data['system_id']))}"
    )
    max_test = int(cfg.split.test_fraction * entries["system_id"].nunique())
    all_test_system_ids = list(
        test_data.sort_values("cluster")["system_id"].drop_duplicates()
    )
    test_system_ids = set(all_test_system_ids[:max_test])
    remaining = set(all_test_system_ids).difference(test_system_ids)
    test_clusters = set(system_to_cluster[x] for x in test_system_ids)
    for x in remaining:
        if system_to_cluster[x] in test_clusters:
            test_system_ids.add(x)
    LOG.info(
        f"test_data after max_test of {cfg.split.test_fraction} ({max_test}): {len(test_system_ids)}"
    )

    test_congeneric_ids = set(
        mms_df[mms_df["system_id"].isin(test_system_ids)]["congeneric_id"]
    )
    test_system_ids |= set(
        mms_df[mms_df["congeneric_id"].isin(test_congeneric_ids)]["system_id"]
    ).intersection(quality)
    LOG.info(f"test_data after adding congeneric series: {len(test_system_ids)}")
    return entries, test_system_ids


def assign_split_membership(
    cfg: DictConfig,
    entries: pd.DataFrame,
    nonredundant_to_all: dict[str, set[str]],
    test_system_ids: set[str],
    system_id_to_leakage: dict[str, set[str]],
    system_to_cluster_val: dict[str, str],
    split_path: Path,
) -> pd.DataFrame:
    """
    Assign system ids to train, val, test and \
    marked others for removal.


    Parameters
    ----------
    data_dir : Path
        Data directory
    entries : pd.DataFrame
        entries annotations
    test_system_ids : set[str]
        set of test systems ids
    system_id_to_leakage : dict[str, set[str]]
        Dictionary of test systems mapped to leaked system
    system_to_cluster_val : dict[str, str]
        System ids mapped to validation cluster
    relpath : str
        Relative path to save splits

    Returns
    -------
    pd.DataFrame
        Dataframe of splits
    """
    LOG.info(f"have {len(test_system_ids)} test system IDs")
    entries["split"] = "train"
    entries.loc[entries["system_id"].isin(test_system_ids), "split"] = "test"
    remove_from_train = set()

    for system_id in entries[entries["split"] == "test"]["system_id"]:
        remove_from_train |= system_id_to_leakage.get(system_id, set())

    LOG.info(f"found {len(remove_from_train)} entries to remove from train")
    entries.loc[
        (entries["split"] == "train") & entries["system_id"].isin(remove_from_train),
        "split",
    ] = "removed"

    train_val_systems = list(
        set(
            entries[
                (entries["split"] == "train")
                & (
                    entries["system_num_pocket_residues"]
                    >= cfg.split.min_max_val_pocket[0]
                )
                & (
                    entries["system_num_pocket_residues"]
                    <= cfg.split.min_max_val_pocket[1]
                )
                & (entries["system_num_interactions"] >= cfg.split.min_max_val_pli[0])
                & (entries["system_num_interactions"] <= cfg.split.min_max_val_pli[1])
            ]["system_id"]
        )
    )
    train_val_clusters = defaultdict(list)

    for x in train_val_systems:
        train_val_clusters[system_to_cluster_val[x]].append(x)

    cluster_order = sorted(train_val_clusters, key=lambda x: len(train_val_clusters[x]))
    for start_cluster_index, x in enumerate(cluster_order):
        if len(train_val_clusters[x]) > cfg.split.min_val_cluster_size:
            break
    val_system_ids: list[str] = []
    remove_ids: list[str] = []
    max_val = int(cfg.split.val_fraction * entries["system_id"].nunique())
    system_to_cluster = dict(zip(entries["system_id"], entries["cluster"]))
    while len(val_system_ids) < max_val:
        clusters = defaultdict(list)
        for x in train_val_clusters[cluster_order[start_cluster_index]]:
            clusters[system_to_cluster[x]].append(x)
        to_add = []
        to_remove = []
        for c in clusters:
            to_add.extend(clusters[c][: cfg.split.num_val_representatives])
            if len(clusters[c]) > cfg.split.num_val_representatives:
                to_remove.extend(clusters[c][cfg.split.num_val_representatives :])
        val_system_ids.extend(to_add)
        remove_ids.extend(to_remove)
        start_cluster_index += 1
    # train_clusters = list(
    #     set(
    #         entries[entries["split"] == "train"]["cluster_for_val_split"].sample(
    #             frac=1.0
    #         )
    #     )
    # )
    # val_clusters = set(train_clusters[: len(train_clusters) // 10])
    # train_system_ids = set(entries[entries["split"] == "train"]["system_id"])
    # val_system_ids = set(
    #     [x for x in train_system_ids if system_to_cluster_val[x] in val_clusters]
    # )
    entries.loc[entries["system_id"].isin(set(val_system_ids)), "split"] = "val"
    LOG.info(f"Additional {len(set(remove_ids))} to be removed due to val")
    entries.loc[entries["system_id"].isin(set(remove_ids)), "split"] = "removed"
    final = entries[
        [
            "system_id",
            "uniqueness",
            "split",
            "cluster",
            "cluster_for_val_split",
        ]
    ].drop_duplicates("system_id")
    final["system_id"] = final["uniqueness"].map(nonredundant_to_all)
    final = final.explode("system_id")
    for k, v in final.value_counts("split").items():
        LOG.info(f"final system splits {k} {v}")
    LOG.info(f"writing out final system splits {len(final.index)} records")
    final.to_parquet(split_path, index=False)
    return final


def split(*, data_dir: Path, cfg: DictConfig, relpath: str) -> pd.DataFrame:
    hash_str = get_config_hash(cfg)
    LOG.info(f"the config path produced this hash: {hash_str}")
    mms_df, entries, quality, nonredundant_to_all = prep_data_for_desired_properties(
        data_dir, cfg
    )

    (
        cluster_to_systems,
        system_to_cluster,
        entries,
        test_systems,
        system_to_cluster_val,
    ) = get_sampling_clusters(cfg, data_dir, entries)

    graphs, vertex_mappings, vertex_ids = load_graphs(
        data_dir, set(entries["system_id"]), cfg
    )
    low_degree_test_systems = remove_high_degree_systems(
        graphs, cfg.split.max_test_leakage_count, test_systems, vertex_mappings
    )
    system_id_to_leakage = {}
    for system_id in tqdm(low_degree_test_systems):
        system_id_to_leakage[system_id] = deleak_entry_nk_single(
            system_id,
            cluster_to_systems[system_to_cluster[system_id]],
            graphs,
            cfg.split.graph_configs,
            vertex_mappings,
            vertex_ids,
        )

    entries, test_system_ids = prioritize_test_sample(
        cfg,
        entries,
        system_id_to_leakage,
        low_degree_test_systems,
        mms_df,
        cluster_to_systems,
        system_to_cluster,
        quality,
    )

    for system_id in test_system_ids:
        if system_id not in system_id_to_leakage:
            system_id_to_leakage[system_id] = deleak_entry_nk_single(
                system_id,
                cluster_to_systems[system_to_cluster[system_id]],
                graphs,
                cfg.split.graph_configs,
                vertex_mappings,
                vertex_ids,
            )
    (data_dir / "splits").mkdir(exist_ok=True)
    breadcrumb = "_".join(Path(relpath).parts[2:]).replace(".", "_")
    split_path = data_dir / "splits" / f"splits_{breadcrumb}_{hash_str}.parquet"
    return assign_split_membership(
        cfg,
        entries,
        nonredundant_to_all,
        test_system_ids,
        system_id_to_leakage,
        system_to_cluster_val,
        split_path,
    )
