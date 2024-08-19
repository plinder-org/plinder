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
from plinder.data.compare_split_props import (
    label_has_mms_within_subset,
    reset_lipinski_and_other,
)

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
            GraphConfig("pli_unique_qcov", 50, 1),
            # GraphConfig("pocket_qcov", 50, 1),
            GraphConfig("protein_seqsim_weighted_sum", 30, 1),
        ]
    )
    # how many unique congeneric IDs passing quality to consider as MMS
    mms_unique_quality_count: int = 5
    # what kind of cluster to use for sampling test
    test_cluster_cluster: str = "communities"
    # metric to use for sampling representatives from each test cluster
    test_cluster_metric: str = "pli_unique_qcov"
    # threshold to use for sampling representatives from each test cluster
    test_cluster_threshold: int = 50
    # directed to use for sampling representatives from each test cluster
    test_cluster_directed: bool = False
    # max number of representatives from each test cluster
    num_test_representatives: int = 1
    # test should not be singletons
    min_test_cluster_size: int = 5
    # test should not be too unique
    min_test_leakage_count: int = 30
    # test should not be in too big communities or cause too many train cases to be removed
    max_test_leakage_count: int = 1000
    # maximum fraction of systems that can be removed due to test set selection
    max_removed_fraction: float = 0.15
    # fraction of systems to choose for test
    test_fraction: float = 0.0025
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
    min_val_cluster_size: int = 30
    # fraction of systems to choose for val
    val_fraction: float = 0.0025
    # test/val should not have too few or too many interactions
    min_max_pli: tuple[int, int] = (3, 50)
    # test/val should not have too few or too many pocket residues
    min_max_pocket: tuple[int, int] = (5, 100)
    # test/val should not have too small or too large ligands
    min_max_ligand: tuple[int, int] = (200, 800)
    # priority columns to use for scoring systems with a weight attached to each column
    priority_columns: dict[str, float] = field(
        default_factory=lambda: {
            "system_has_mms": 0.5,
            "system_has_binding_affinity": 1.0,
            "leakage_count": 0,
            "system_has_cofactor": -0.5,
        }
    )


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
        leaked_systems.update({vertex_ids_x[x] for x in leaked_neighbors})
    # the cluster + neighbors needs to be removed if that system is chosen
    return leaked_systems


def prep_data_for_desired_properties(
    data_dir: Path,
    cfg: DictConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    """
    Load data and add annotations relevant for splitting

    Parameters
    ----------
    cfg : DictConfig
        splitting config

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]
        mms dataframe, entries annotations, mapping of nonredundant systems to all systems
    """
    LOG.info(OmegaConf.to_yaml(cfg))

    entries = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")

    # sanity check to only split entries that were successfully scored
    # foldseek_ids = get_ids_in_db(
    #     data_dir=data_dir, search_db="holo", aln_type="foldseek"
    # )
    # mmseqs_ids = get_ids_in_db(data_dir=data_dir, search_db="holo", aln_type="mmseqs")
    # foldseek_entries_present = foldseek_ids[1].str.replace("pdb_0000", "").str[:4]
    # mmseqs_entries_present = mmseqs_ids[1].str.replace("pdb_0000", "").str[:4]
    # all_entries_present = set(foldseek_entries_present).intersection(
    #     mmseqs_entries_present
    # )
    all_entries_present = set(
        pd.read_csv(data_dir / "dbs" / "subdbs" / "holo_ids.csv")["pdb_id"]
    )

    entries = entries[
        (entries["system_type"] == "holo")
        & (entries["system_num_interacting_protein_chains"] <= 5)
        & (entries["system_num_ligand_chains"] <= 5)
        & (entries["entry_pdb_id"].isin(all_entries_present))
    ].reset_index(drop=True)
    LOG.info(f"loaded {entries['system_id'].nunique()} systems from annotation table")
    entries = reset_lipinski_and_other(entries)
    entries["system_has_lipinski"] = entries.groupby("system_id")[
        "ligand_is_lipinski"
    ].transform("any")
    entries["system_has_cofactor"] = entries.groupby("system_id")[
        "ligand_is_cofactor"
    ].transform("any")
    entries["ligand_is_proper"] = (
        ~entries["ligand_is_ion"] & ~entries["ligand_is_artifact"]
    )
    entries["system_proper_num_ligand_chains"] = entries.groupby("system_id")[
        "ligand_is_proper"
    ].transform("sum")
    entries["system_proper_num_interactions"] = (
        entries["ligand_num_interactions"]
        .where(entries["ligand_is_proper"], other=0)
        .groupby(entries["system_id"])
        .transform("sum")
    )
    entries["system_proper_pocket_num_residues"] = (
        entries["ligand_num_neighboring_residues"]
        .where(entries["ligand_is_proper"], other=0)
        .groupby(entries["system_id"])
        .transform("sum")
    )
    entries["system_proper_ligand_max_molecular_weight"] = (
        entries["ligand_molecular_weight"]
        .where(entries["ligand_is_proper"], other=float("-inf"))
        .groupby(entries["system_id"])
        .transform("max")
    )

    entries["uniqueness"] = (
        entries["system_id_no_biounit"]
        + "_"
        + entries["pli_qcov__100__strong__component"]
    )
    nonredundant_to_all = defaultdict(set)
    for system_id, unique_id in zip(entries["system_id"], entries["uniqueness"]):
        nonredundant_to_all[unique_id].add(system_id)
    entries = entries.sort_values("system_biounit_id").drop_duplicates(
        "uniqueness"
    )  # after this point there's only one row per system
    all_system_ids = set(entries["system_id"])

    LOG.info(f"loaded {len(all_system_ids)} nonredundant systems from annotation table")

    entries["system_pass_statistics_criteria"] = (
        (entries["system_proper_pocket_num_residues"] >= cfg.split.min_max_pocket[0])
        & (entries["system_proper_pocket_num_residues"] <= cfg.split.min_max_pocket[1])
        & (entries["system_proper_num_interactions"] >= cfg.split.min_max_pli[0])
        & (entries["system_proper_num_interactions"] <= cfg.split.min_max_pli[1])
        & (
            entries["system_proper_ligand_max_molecular_weight"]
            >= cfg.split.min_max_ligand[0]
        )
        & (
            entries["system_proper_ligand_max_molecular_weight"]
            <= cfg.split.min_max_ligand[1]
        )
    )
    entries["system_pass_validation_criteria"] = (
        entries["system_pass_validation_criteria"].fillna(False).astype(bool)
    )
    LOG.info(
        f"{entries['system_pass_validation_criteria'].sum()} systems pass validation criteria"
    )
    LOG.info(
        f"{entries['system_pass_statistics_criteria'].sum()} systems pass statistics criteria"
    )
    entries["proto_test"] = (
        entries["system_pass_validation_criteria"]
        & entries["system_pass_statistics_criteria"]
    )
    LOG.info(
        f"{len(entries[entries['proto_test']].index)} systems pass quality and statistics filters"
    )

    apo_links = pd.read_parquet(
        data_dir / "links" / "apo_links.parquet", columns=["reference_system_id"]
    )
    apo_links["system_id"] = apo_links["reference_system_id"]
    assert apo_links is not None, "apo_links is None"
    num_apo_links = apo_links.groupby("system_id").size()
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

    entries["system_has_apo"] = entries["system_id"].map(
        lambda x: num_apo_links.get(x, 0) > 0
    )
    entries["system_has_pred"] = entries["system_id"].map(
        lambda x: num_pred_links.get(x, 0) > 0
    )
    entries["system_has_apo_or_pred"] = (
        entries["system_has_apo"] | entries["system_has_pred"]
    )

    mms_df = pd.read_parquet(data_dir / "mmp/plinder_mmp_series.parquet")
    LOG.info(f"read mmp series has {len(mms_df.index)} records")
    quality = set(entries[entries["system_pass_validation_criteria"]]["system_id"])
    good_mms_systems = label_has_mms_within_subset(
        mms_df, quality, cfg.split.mms_unique_quality_count
    )
    # mms_df = mms_df[mms_df["system_id"].isin(quality)]
    # LOG.info(
    #     f"mmp series after quality control {len(mms_df.index)} records for {mms_df['system_id'].nunique()} systems"
    # )
    # mms_count = mms_df.groupby("congeneric_id").size()
    # mms_df["mms_unique_quality_count"] = mms_df["congeneric_id"].map(mms_count)
    # good_mms_systems = set(
    #     mms_df[
    #         mms_df["mms_unique_quality_count"] >= cfg.split.mms_unique_quality_count
    #     ]["system_id"]
    # )
    entries["system_has_mms"] = entries["system_id"].isin(good_mms_systems)
    LOG.info(f"entries with good mms: {len(entries[entries['system_has_mms']].index)}")
    return mms_df, entries, nonredundant_to_all


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
    LOG.info(
        f"loaded {len(set(labels.values()))} {cfg.split.test_cluster_cluster} clusters for sampling"
    )
    entries["cluster"] = entries["system_id"].map(labels)

    cluster_to_size = Counter(entries["cluster"])
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
        f"loaded {len(set(labels.values()))} {cfg.split.val_cluster_cluster} clusters for splitting train and val"
    )
    entries["cluster_for_val_split"] = entries["system_id"].map(labels)

    entries["split"] = "train"
    entries.loc[entries["proto_test"], "split"] = "proto-test"

    test_systems = set(
        row.system_id
        for _, row in tqdm(
            entries[entries["proto_test"]][["system_id", "cluster"]].iterrows()
        )
        if cluster_to_size[row["cluster"]] >= cfg.split.min_test_cluster_size
    )
    LOG.info(f"found {len(test_systems)} test systems")

    return (
        cluster_to_systems,
        system_to_cluster,
        entries,
        test_systems,
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
    test_data = entries[entries["system_id"].isin(test_systems)]

    test_data = test_data[
        (test_data["leakage_count"] >= cfg.split.min_test_leakage_count)
        & (test_data["leakage_count"] <= cfg.split.max_test_leakage_count)
    ]
    LOG.info(
        f"found {len(set(test_data['system_id']))} test systems after removing by second pass neighbor leakage count"
    )
    priority_column_order = sorted(
        cfg.split.priority_columns,
        key=lambda x: cfg.split.priority_columns[x],
        reverse=True,
    )
    ascending = [cfg.split.priority_columns[x] <= 0 for x in priority_column_order]
    test_data = (
        test_data.sort_values(
            priority_column_order,
            ascending=ascending,
        )
        .groupby("cluster")
        .head(cfg.split.num_test_representatives)
        .reset_index(drop=True)
    )
    LOG.info(
        f"test_data after groupby.head of num_test_representatives: {len(set(test_data['system_id']))}"
    )
    max_test = int(cfg.split.test_fraction * entries["system_id"].nunique())

    test_data["system_score"] = (
        test_data[priority_column_order].mul(cfg.split.priority_columns).sum(axis=1)
    )
    test_data["cluster_score"] = test_data.groupby("cluster")["system_score"].transform(
        "max"
    )
    all_test_system_ids = (
        test_data.sort_values(
            ["cluster_score", "system_score", "leakage_count"],
            ascending=[False, False, True],
        )
        .groupby("cluster", sort=False)
        .apply(
            lambda x: x.sort_values(
                ["system_score", "leakage_count"], ascending=[False, True]
            )
        )
        .reset_index(drop=True)["system_id"]
        .tolist()
    )
    max_removed = int(cfg.split.max_removed_fraction * entries["system_id"].nunique())
    test_system_ids = set()
    to_remove = set()
    for x in all_test_system_ids:
        if len(test_system_ids) < max_test:
            test_system_ids.add(x)
            to_remove |= system_id_to_leakage.get(x, set())
        if len(to_remove) > max_removed:
            break
    # test_system_ids = set(all_test_system_ids[:max_test])
    remaining = set(all_test_system_ids).difference(test_system_ids)
    test_clusters = set(system_to_cluster[x] for x in test_system_ids)
    for x in remaining:
        if system_to_cluster[x] in test_clusters:
            test_system_ids.add(x)
    LOG.info(
        f"test_data after max_test of {cfg.split.test_fraction} ({max_test}): {len(test_system_ids)}"
    )

    test_congeneric_ids = set(
        mms_df[
            mms_df["system_id"].isin(test_system_ids)
            & (mms_df["mms_unique_quality_count"] >= cfg.split.mms_unique_quality_count)
        ]["congeneric_id"]
    )
    test_system_ids |= set(
        mms_df[mms_df["congeneric_id"].isin(test_congeneric_ids)]["system_id"]
    )
    LOG.info(f"test_data after adding congeneric series: {len(test_system_ids)}")
    return entries, test_system_ids


def assign_split_membership(
    cfg: DictConfig,
    entries: pd.DataFrame,
    nonredundant_to_all: dict[str, set[str]],
    test_system_ids: set[str],
    system_id_to_leakage: dict[str, set[str]],
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
    relpath : str
        Relative path to save splits

    Returns
    -------
    pd.DataFrame
        Dataframe of splits
    """
    LOG.info(f"have {len(test_system_ids)} test system IDs")
    pdb_id_to_systems = (
        entries.groupby("entry_pdb_id")["system_id"].apply(set).to_dict()
    )

    entries["split"] = "train"
    entries.loc[entries["system_id"].isin(test_system_ids), "split"] = "test"
    remove_from_train = (
        set.union(
            *(
                system_id_to_leakage.get(system_id, set())
                for system_id in test_system_ids
            )
        )
        - test_system_ids
    )
    LOG.info(f"found {len(remove_from_train)} entries to remove from train")
    entries.loc[
        (entries["split"] == "train") & entries["system_id"].isin(remove_from_train),
        "split",
    ] = "removed"

    # Group all train systems by their validation clusters
    entries["system_score_val"] = (
        entries["system_pass_validation_criteria"]
        + entries["system_pass_statistics_criteria"]
    )
    train_val_clusters = (
        entries[entries["split"] == "train"]
        .sort_values("system_score_val", ascending=False)
        .groupby("cluster_for_val_split")["system_id"]
        .apply(list)
        .to_dict()
    )
    passing_statistics_criteria = set(
        entries[entries["system_pass_statistics_criteria"]]["system_id"]
    )
    passing_validation_criteria = set(
        entries[entries["system_pass_validation_criteria"]]["system_id"]
    )
    # Sort clusters by size
    cluster_order = sorted(train_val_clusters, key=lambda x: len(train_val_clusters[x]))

    # Find the starting index for clusters that meet the minimum size
    start_cluster_index = next(
        (
            i
            for i, x in enumerate(cluster_order)
            if len(train_val_clusters[x]) >= cfg.split.min_val_cluster_size
        ),
        len(cluster_order),
    )

    val_system_ids = []
    remove_ids = set()
    max_val = int(cfg.split.val_fraction * entries["system_id"].nunique())

    for cluster in cluster_order[start_cluster_index:]:
        systems = train_val_clusters[cluster]

        if all(x not in passing_validation_criteria for x in systems) or all(
            x not in passing_statistics_criteria for x in systems
        ):
            # If we don't have enough valid systems to meet the minimum size, skip this cluster
            continue

        val_system_ids.extend(systems[: cfg.split.num_val_representatives])

        # Remove all systems with this cluster ID from train
        remove_ids.update(systems)

        if len(val_system_ids) >= max_val:
            break
    remove_ids -= set(val_system_ids)
    LOG.info(
        f"val_system_ids after max_val of {cfg.split.val_fraction} ({max_val}): {len(val_system_ids)}"
    )

    # Remove val PDB IDs
    for x in val_system_ids:
        remove_ids |= pdb_id_to_systems[x.split("__")[0]]
    remove_ids -= set(val_system_ids)

    entries.loc[entries["system_id"].isin(set(val_system_ids)), "split"] = "val"
    LOG.info(f"Additional {len(remove_ids)} to be removed due to val")
    entries.loc[entries["system_id"].isin(remove_ids), "split"] = "removed"
    final = entries[
        [
            "system_id",
            "uniqueness",
            "split",
            "cluster",
            "cluster_for_val_split",
            "system_pass_validation_criteria",
            "system_pass_statistics_criteria",
            "system_proper_num_ligand_chains",
            "system_proper_pocket_num_residues",
            "system_proper_num_interactions",
            "system_proper_ligand_max_molecular_weight",
            "system_has_binding_affinity",
            "system_has_lipinski",
            "system_has_apo_or_pred",
            "system_has_mms",
        ]
    ].reset_index(drop=True)
    final["system_id"] = final["uniqueness"].map(nonredundant_to_all)
    final = final.explode("system_id")
    for k, v in final.value_counts("split").items():
        LOG.info(f"final system splits {k} {v}")
    LOG.info(f"writing out final system splits {len(final.index)} records")
    final.to_parquet(split_path, index=False)
    return final


def split(*, data_dir: Path, cfg: DictConfig, relpath: str) -> pd.DataFrame:
    OmegaConf.save(config=cfg, f=data_dir / "splits" / f"split_{relpath}.yaml")

    mms_df, entries, nonredundant_to_all = prep_data_for_desired_properties(
        data_dir=data_dir,
        cfg=cfg,
    )

    (
        cluster_to_systems,
        system_to_cluster,
        entries,
        test_systems,
    ) = get_sampling_clusters(cfg=cfg, data_dir=data_dir, entries=entries)

    pdb_id_to_systems = (
        entries.groupby("entry_pdb_id")["system_id"].apply(set).to_dict()
    )
    LOG.info(f"loaded {len(pdb_id_to_systems)} pdb_id to system_id mappings")

    graphs, vertex_mappings, vertex_ids = load_graphs(
        data_dir=data_dir, systems=set(entries["system_id"]), cfg=cfg
    )
    low_degree_test_systems = remove_high_degree_systems(
        graphs=graphs,
        max_test_leakage_count=cfg.split.max_test_leakage_count,
        test_systems=test_systems,
        vertex_mappings=vertex_mappings,
    )
    system_id_to_leakage = {}
    for system_id in tqdm(low_degree_test_systems):
        system_id_to_leakage[system_id] = deleak_entry_nk_single(
            system_id=system_id,
            cluster_members=cluster_to_systems[system_to_cluster[system_id]],
            graphs=graphs,
            graph_configs=cfg.split.graph_configs,
            vertex_mappings=vertex_mappings,
            vertex_ids=vertex_ids,
        ).union(pdb_id_to_systems[system_id.split("__")[0]]) - set([system_id])

    entries, test_system_ids = prioritize_test_sample(
        cfg=cfg,
        entries=entries,
        system_id_to_leakage=system_id_to_leakage,
        test_systems=low_degree_test_systems,
        mms_df=mms_df,
        cluster_to_systems=cluster_to_systems,
        system_to_cluster=system_to_cluster,
    )

    for system_id in test_system_ids:
        if system_id not in system_id_to_leakage:
            system_id_to_leakage[system_id] = deleak_entry_nk_single(
                system_id=system_id,
                cluster_members=cluster_to_systems[system_to_cluster[system_id]],
                graphs=graphs,
                graph_configs=cfg.split.graph_configs,
                vertex_mappings=vertex_mappings,
                vertex_ids=vertex_ids,
            ).union(pdb_id_to_systems[system_id.split("__")[0]]) - set([system_id])
    (data_dir / "splits").mkdir(exist_ok=True)
    split_path = data_dir / "splits" / f"split_{relpath}.parquet"
    return assign_split_membership(
        cfg=cfg,
        entries=entries,
        nonredundant_to_all=nonredundant_to_all,
        test_system_ids=test_system_ids,
        system_id_to_leakage=system_id_to_leakage,
        split_path=split_path,
    )
