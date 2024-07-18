# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path

import pandas as pd

from plinder.core.utils import gcs
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def get_workshop_splits() -> list[str]:
    remote_dir = "gs://plinder-collab-bucket/2024-04"
    return [
        f"{remote_dir}/v1/splits/v0/pdbbind/pdbbind.csv",
        f"{remote_dir}/v1/splits/v0/pdbbind_lp/pdbbind_lp.csv",
        f"{remote_dir}/v1/splits/v0/equibind/equibind.csv",
        f"{remote_dir}/v1/splits/v0/dockgen/dockgen.csv",
        f"{remote_dir}/v1/workshop_split/ecod_split.csv",
        f"{remote_dir}/v1/workshop_split/plinder_v0_no_posbuster.csv",
        f"{remote_dir}/v1/workshop_split/time_split.csv",
    ]


def get_missing_pdb_ids(data_dir: Path) -> set[str]:
    entries = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")
    seqsim_cluster = pd.read_parquet(
        data_dir
        / "clusters"
        / "cluster=components"
        / "directed=False"
        / "metric=protein_seqsim_max"
        / "threshold=50.parquet"
    )

    LOG.info(f"entries {entries.shape}")
    LOG.info(f"seqsim cluster {seqsim_cluster.shape}")
    s = seqsim_cluster.groupby("label").size()
    missing = set(
        seqsim_cluster[seqsim_cluster["label"].isin(s[s == 1].index.values)][
            "system_id"
        ]
    )
    LOG.info(f"singular system_ids {len(set(missing))}")
    pdb_ids = entries[entries["system_id"].isin(missing)][["system_id", "entry_pdb_id"]]
    LOG.info(f"missing entry_pdb_ids {pdb_ids['entry_pdb_id'].nunique()}")
    pdb_ids.to_csv("missing_pdb_ids.csv", index=False)
    LOG.info("missing the following pdb ids in the scores dataset:")
    LOG.info(",".join(set(pdb_ids["entry_pdb_id"])))
    return set(pdb_ids["entry_pdb_id"])


def get_index_and_split_sets(
    data_dir: Path,
    split_file: str,
    compare_pair: str,
) -> tuple[pd.DataFrame, str, pd.DataFrame, set[str], set[str]]:
    index = pd.read_parquet(data_dir / "index/annotation_table.parquet")
    missing = get_missing_pdb_ids(data_dir)
    if split_file.endswith(".csv"):
        split = pd.read_csv(split_file)
    else:
        split = pd.read_parquet(split_file)
    if compare_pair == "train_posebusters":
        left = set(split[split["split"] == "train"]["system_id"]).difference(missing)
        text = gcs.download_as_str(
            gcs_path="gs://plinder-collab-bucket/2024-04/v1/posebuster_system_ids.csv"
        )
        right = set(line.strip() for line in text.splitlines()).difference(missing)
    else:
        left_target, right_target = compare_pair.split("_")
        left = set(split[split["split"] == left_target]["system_id"]).difference(
            missing
        )
        right = set(split[split["split"] == right_target]["system_id"]).difference(
            missing
        )
    replace_keys = {
        "ecod_split": "plinder-ECOD",
        "plinder_v0_no_posbuster": "plinder-v0",
        "time_split": "plinder-Time",
    }
    split_name = Path(split_file).stem
    split_name = replace_keys.get(split_name, split_name)
    LOG.info(
        f"split_file={split_file}, split_name={split_name}, compare_pair={compare_pair}"
    )
    LOG.info(f"num_left={len(left)}, num_right={len(right)}")
    return index, split_name, split, left, right


def compute_ligand_leakage(
    *,
    data_dir: Path,
    split_file: str,
    compare_pair: str,
    metric: str,
) -> None:
    df, split_name, split, left, right = get_index_and_split_sets(
        data_dir, split_file, compare_pair
    )

    ligands_per_system = pd.read_parquet(
        data_dir / "fingerprints/ligands_per_system.parquet"
    )
    mapr = dict(
        zip(
            ligands_per_system["ligand_rdkit_canonical_smiles"],
            ligands_per_system["number_id_by_inchikeys"],
        )
    )
    df["number_id_by_inchikeys"] = df["ligand_rdkit_canonical_smiles"].map(mapr)
    LOG.info(f"number_id_by_inchikeys={df['number_id_by_inchikeys'].nunique()}")

    left_ligand_ids = set(df[df["system_id"].isin(left)]["number_id_by_inchikeys"])
    right_ligand_ids = set(df[df["system_id"].isin(right)]["number_id_by_inchikeys"])
    LOG.info(f"left_ligand_ids={len(left_ligand_ids)}")
    LOG.info(f"right_ligand_ids={len(right_ligand_ids)}")

    ligands = pd.read_parquet(
        data_dir / "ligand_scores",
        columns=["query_ligand_id", "target_ligand_id", metric],
        filters=[
            [
                ("query_ligand_id", "in", left_ligand_ids),
                ("target_ligand_id", "in", right_ligand_ids),
            ],
            [
                ("query_ligand_id", "in", right_ligand_ids),
                ("target_ligand_id", "in", left_ligand_ids),
            ],
        ],
    )

    LOG.info(f"ligands={len(ligands.index)}")
    LOG.info(f"unique_query_ligand_id={ligands['query_ligand_id'].nunique()}")
    LOG.info(f"unique_target_ligand_id={ligands['target_ligand_id'].nunique()}")

    updated_query_ligand_ids = []
    for q, t in zip(ligands["query_ligand_id"], ligands["target_ligand_id"]):
        if t in right_ligand_ids:
            updated_query_ligand_ids.append(t)
        else:
            updated_query_ligand_ids.append(q)
    ligands["updated_query_ligand_id"] = updated_query_ligand_ids
    LOG.info(
        f"unique_updated_query_ligand_id={ligands['updated_query_ligand_id'].nunique()}"
    )

    idx = ligands.groupby("updated_query_ligand_id")[metric].idxmax()
    ligands = ligands.loc[idx]
    ligand_to_system = {}
    for ligand_id, group in ligands_per_system.groupby("number_id_by_inchikeys"):
        ligand_to_system[ligand_id] = set(group["system_id"])
    ligands["query_system"] = ligands["query_ligand_id"].map(ligand_to_system)
    leakage_dir = data_dir / "leakage"
    leakage_dir.mkdir(exist_ok=True, parents=True)
    exploded = ligands.explode("query_system")
    LOG.info(f"exploded={len(exploded.index)}")
    filename = leakage_dir / f"{split_name}__{metric}__{compare_pair}.parquet"
    LOG.info(f"writing {filename}")
    exploded[["query_system", metric]].rename(
        columns={
            metric: "similarity",
        }
    ).to_parquet(filename, index=False)


def compute_protein_leakage(
    *,
    data_dir: Path,
    split_file: str,
    compare_pair: str,
    metric: str,
) -> None:
    df, split_name, split, left, right = get_index_and_split_sets(
        data_dir, split_file, compare_pair
    )
    # quality_config = TestCriteria()
    # df["system_num_covalent_ligands"] = df.groupby("system_id")["ligand_is_covalent"].transform("sum")
    # df["passes_quality"] = df.apply(lambda row: get_high_quality_systems(row, criteria=quality_config), axis=1)
    # passes_quality = set(df[df["passes_quality"]]["system_id"])

    leakage_dir = data_dir / "leakage"
    leakage_dir.mkdir(exist_ok=True, parents=True)
    # for metric in [
    #     "pli_qcov",
    #     "pocket_qcov",
    #     "pocket_lddt",
    #     "protein_seqsim_weighted_sum",
    #     "protein_fident_weighted_sum",
    #     "protein_lddt_weighted_sum",
    #     "protein_lddt_qcov_weighted_sum",
    # ]:
    scores = pd.read_parquet(
        data_dir / "scores" / "search_db=holo",
        columns=["query_system", "target_system", "similarity"],
        filters=[
            [
                ("query_system", "in", left),
                ("target_system", "in", right),
                ("metric", "==", metric),
            ],
            [
                ("query_system", "in", right),
                ("target_system", "in", left),
                ("metric", "==", metric),
            ],
        ],
    )
    updated_query_systems = []
    for q, t in zip(scores["query_system"], scores["target_system"]):
        if t in right:
            updated_query_systems.append(t)
        else:
            updated_query_systems.append(q)

    scores["updated_query_system"] = updated_query_systems

    LOG.info(f"scores={len(scores.index)}")
    idx = scores.groupby("updated_query_system")["similarity"].idxmax()
    scores = scores.loc[idx][["updated_query_system", "similarity"]].rename(
        columns={"updated_query_system": "query_system"}
    )
    LOG.info(f"scores_after_groupby={len(scores.index)}")
    filename = leakage_dir / f"{split_name}__{metric}__{compare_pair}.parquet"
    LOG.info(f"writing {filename}")
    scores.to_parquet(filename, index=False)
    # if compare_pair in ["train_test", "val_test"]:
    #     scores = pd.read_parquet(
    #         data_dir / "scores" / "search_db=holo",
    #         columns=["query_system", "similarity"],
    #         filters=[
    #             [
    #                 ("query_system", "in", left.intersection(passes_quality)),
    #                 ("target_system", "in", right.intersection(passes_quality)),
    #                 ("metric", "==", metric),
    #             ],
    #             [
    #                 ("query_system", "in", right.intersection(passes_quality)),
    #                 ("target_system", "in", left.intersection(passes_quality)),
    #                 ("metric", "==", metric),
    #             ],
    #         ],
    #     )
    #     LOG.info(f"scores_hq={len(scores.index)}")
    #     idx = scores.groupby("query_system")["similarity"].idxmax()
    #     scores = scores.loc[idx]
    #     LOG.info(f"scores_after_groupby_hq={len(scores.index)}")
    #     filename = leakage_dir / f"{split_name}__{metric}__{compare_pair}_hq.parquet"
    #     LOG.info(f"writing {filename}")
    #     scores.to_parquet(filename, index=False)
