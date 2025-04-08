# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from cloudpathlib import AnyPath
from tqdm import tqdm

import plinder.core.utils.config
from plinder.core.scores.index import query_index
from plinder.core.scores.protein import cross_similarity as protein_cross_similarity
from plinder.core.structure import smallmols_similarity
from plinder.core.utils.log import setup_logger

cfg = plinder.core.get_config()

LOG = setup_logger(__name__)

SIMILARITY_METRICS = (
    # pli
    "pli_unique_qcov",
    # protein
    "protein_seqsim_weighted_sum",
    "protein_fident_weighted_sum",
    "protein_lddt_weighted_sum",
    # pocket
    "pocket_fident_qcov",
    "pocket_lddt_qcov",
    "pocket_lddt",
    "pocket_qcov",
    # ligand
    "tanimoto_similarity_max",
    "mmp_similarity_max",
)

STRATIFICATION_LABELS = [
    "novel_pocket_pli",
    "novel_ligand",
    "novel_protein",
    "novel_all",
]


def compute_protein_max_similarities(
    left: set[str], right: set[str], metric: str, output_file: Path
) -> None:
    LOG.info(
        f"compute_protein_max_similarities: Computing max similarities for {metric}"
    )
    protein_cross_similarity(
        query_systems=left,
        target_systems=right,
        metric=metric,
    ).rename(
        columns={"query_system": "system_id", "target_system": "train_system_id"}
    ).to_parquet(output_file, index=False)
    LOG.info(
        f"compute_protein_max_similarities: Done computing max similarities for {metric}"
    )


def compute_ligand_ecfp_max_similarities(
    df: pd.DataFrame,
    split_label: str,
    train_label: str,
    test_label: str,
    output_file: Path,
) -> None:
    if "fp" not in df.columns:
        smiles_fp_dict = {
            smi: smallmols_similarity.mol2morgan_fp(smi)
            for smi in df["ligand_rdkit_canonical_smiles"].drop_duplicates().to_list()
        }
        df["fp"] = df["ligand_rdkit_canonical_smiles"].map(smiles_fp_dict)

    df_test = df.loc[df[split_label] == test_label][
        ["system_id", "ligand_rdkit_canonical_smiles", "fp"]
    ].copy()

    (
        df_test["tanimoto_similarity_max"],
        argmax_array,
    ) = smallmols_similarity.tanimoto_maxsim_and_argmax(
        df.loc[df[split_label] == train_label, "fp"].to_list(),
        df_test["fp"].to_list(),
    )
    # get most similar smiles in train
    train_smiles = df.loc[
        df[split_label] == train_label, "ligand_rdkit_canonical_smiles"
    ].to_list()
    df_test[f"tanimoto_most_similar_{train_label}_smiles"] = [
        train_smiles[idx] for idx in argmax_array
    ]

    df_test.drop("fp", axis=1).groupby(
        ["system_id", "ligand_rdkit_canonical_smiles"]
    ).agg("max").reset_index().to_parquet(output_file, index=False)


def compute_ligand_mmp_max_similarities(
    df: pd.DataFrame,
    split_label: str,
    train_label: str,
    test_label: str,
    output_file: Path,
) -> None:
    mmp_path = AnyPath(f"{cfg.data.plinder_remote}/mmp/plinder_mms.csv.gz")
    mmp_sim_dict: dict[
        str, dict[str, float]
    ] = smallmols_similarity.get_mmp_similarity_dict(mmp_path=mmp_path)

    if "inchikey" not in df.columns:
        smi_inchikey_map = {
            smi: smallmols_similarity.smiles2inchikey(smi, remove_stereo=True)
            for smi in df.ligand_rdkit_canonical_smiles.unique()
        }
        df["inchikey"] = df.ligand_rdkit_canonical_smiles.map(smi_inchikey_map)
    df_test = df.loc[df[split_label] == test_label][
        ["system_id", "ligand_rdkit_canonical_smiles", "inchikey"]
    ].copy()
    train_inchikeys = (
        df[df.split == train_label]["inchikey"].drop_duplicates().to_list()
    )
    test_inchikeys = df_test["inchikey"].drop_duplicates().to_list()
    test_train_sims = {}
    for test_ink in test_inchikeys:
        if test_ink in train_inchikeys:
            # by definition 1 if in train set
            max_similarity = 1.0
            max_sim_inchikey = test_ink
        else:
            # get similarity dictionary
            sim_dict_i = mmp_sim_dict.get(test_ink, {})
            if sim_dict_i:
                similarities = [sim_dict_i.get(ink, 0) for ink in train_inchikeys]
                max_similarity = np.max(similarities)
                max_sim_inchikey = train_inchikeys[np.argmax(similarities)]
            else:
                max_similarity = 0.0

        # get first SMILES in train set that match most similar inchikey
        if max_similarity > 0:
            most_similar_smiles_in_train = df.loc[
                (df[split_label] == train_label) & (df["inchikey"] == max_sim_inchikey),
                "ligand_rdkit_canonical_smiles",
            ].iloc[0]
        else:
            most_similar_smiles_in_train = None
        # store
        test_train_sims[test_ink] = [max_similarity, most_similar_smiles_in_train]

    df_test["mmp_similarity_max"] = df_test["inchikey"].apply(
        lambda x: test_train_sims[x][0]
    )
    df_test[f"mmp_most_similar_{train_label}_smiles"] = df_test["inchikey"].apply(
        lambda x: test_train_sims[x][1]
    )
    # save
    df_test.groupby(["system_id", "ligand_rdkit_canonical_smiles"]).agg(
        "max"
    ).reset_index().to_parquet(output_file, index=False)


@dataclass
class StratifiedTestSet:
    split_df: pd.DataFrame
    output_dir: Path
    split_label: str = "split"
    train_label: str = "train"
    test_label: str = "test"
    similarity_thresholds: dict[str, int] = field(
        default_factory=lambda: dict(
            pli_unique_qcov=50,
            pocket_lddt_qcov=50,
            pocket_lddt=50,
            pocket_qcov=50,
            protein_seqsim_weighted_sum=30,
            protein_lddt_weighted_sum=50,
            tanimoto_similarity_max=30,
            mmp_similarity_max=50,
        )
    )
    similarity_combinations: dict[str, list[str]] = field(
        default_factory=lambda: {
            "novel_pocket_pli": ["pli_unique_qcov", "pocket_qcov", "pocket_lddt_qcov"],
            "novel_protein": [
                "protein_seqsim_weighted_sum",
                "protein_lddt_weighted_sum",
            ],
            "novel_ligand": [
                "tanimoto_similarity_max",
                "mmp_similarity_max",
            ],
            "novel_all": [
                "pli_unique_qcov",
                "pocket_qcov",
                "pocket_lddt",
                "protein_seqsim_weighted_sum",
                "protein_lddt_weighted_sum",
                "tanimoto_similarity_max",
            ],
        }
    )
    max_similarities: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(  # type: ignore
            columns=["system_id"] + list(SIMILARITY_METRICS)
        )
    )

    @classmethod
    def from_split(
        cls,
        split_file: Path,
        output_dir: Path,
        split_label: str = "split",
        train_label: str = "train",
        test_label: str = "test",
        overwrite: bool = False,
    ) -> "StratifiedTestSet":
        if split_file.name.endswith(".csv"):
            split_df = pd.read_csv(split_file)
        else:
            split_df = pd.read_parquet(split_file)
        assert all(x in split_df.columns for x in [split_label, "system_id"])
        split_df = split_df[split_df[split_label].isin([train_label, test_label])][
            [split_label, "system_id"]
        ].reset_index(drop=True)
        data = cls(
            split_df=split_df,
            output_dir=output_dir,
            split_label=split_label,
            train_label=train_label,
            test_label=test_label,
        )
        data.output_dir.mkdir(exist_ok=True)
        data.compute_train_test_max_similarity(overwrite=overwrite)
        data.assign_test_set_quality()
        data.stratify_test_set()
        data.max_similarities.to_parquet(data.output_dir / f"{test_label}_set.parquet")
        return data

    def stratify_test_set(self) -> None:
        for label, metric_list in self.similarity_combinations.items():
            self.max_similarities[label] = np.logical_and.reduce(
                [
                    self.max_similarities[metric] < self.similarity_thresholds[metric]
                    for metric in metric_list
                ]
            )
            LOG.info(
                f'stratify_test_set: Found {self.max_similarities[self.max_similarities[label]]["system_id"].nunique()} systems labelled {label} ({self.max_similarities[self.max_similarities[label] & self.max_similarities["passes_quality"]]["system_id"].nunique()} passing quality)'
            )
        self.max_similarities["not_novel"] = np.logical_and.reduce(
            [~self.max_similarities[label] for label in self.similarity_combinations]
        )
        LOG.info(
            f'stratify_test_set: Found {self.max_similarities[self.max_similarities["not_novel"]]["system_id"].nunique()} systems labelled not_novel ({self.max_similarities[self.max_similarities["not_novel"] & self.max_similarities["passes_quality"]]["system_id"].nunique()} passing quality)'
        )

    def get_filename(self, metric: str) -> Path:
        return (
            self.output_dir
            / f"max_similarities__{self.test_label}_vs_{self.train_label}__{metric}.parquet"
        )

    def compute_train_test_max_similarity(
        self, overwrite: bool = False
    ) -> pd.DataFrame:
        left, right = (
            set(
                self.split_df[self.split_df[self.split_label] == self.train_label][
                    "system_id"
                ]
            ),
            set(
                self.split_df[self.split_df[self.split_label] == self.test_label][
                    "system_id"
                ]
            ),
        )
        LOG.info(
            f"compute_train_test_max_similarity: Found {len(left)} train and {len(right)} test systems"
        )
        for metric in tqdm(SIMILARITY_METRICS):
            if overwrite or not (self.get_filename(metric)).exists():
                if metric in ["tanimoto_similarity_max", "mmp_similarity_max"]:
                    # avoid repetitive local assignment
                    if "df" not in locals():
                        df = query_index(
                            columns=[
                                "system_id",
                                "ligand_rdkit_canonical_smiles",
                            ],
                            filters=[  # type: ignore
                                ("system_id", "in", left.union(right)),
                                ("ligand_is_ion", "==", False),
                                ("ligand_is_artifact", "==", False),
                            ],
                            splits=["*"],
                        ).drop(columns=["split"])
                        df = df.merge(self.split_df, on="system_id", how="left")
                    if metric == "tanimoto_similarity_max":
                        compute_ligand_ecfp_max_similarities(
                            df,
                            self.split_label,
                            self.train_label,
                            self.test_label,
                            self.get_filename(metric),
                        )
                    elif metric == "mmp_similarity_max":
                        compute_ligand_mmp_max_similarities(
                            df,
                            self.split_label,
                            self.train_label,
                            self.test_label,
                            self.get_filename(metric),
                        )
                else:
                    compute_protein_max_similarities(
                        left, right, metric, self.get_filename(metric)
                    )

        per_metric_similarities = []
        for metric in SIMILARITY_METRICS:
            df = pd.read_parquet(self.get_filename(metric))
            df = df.loc[df.groupby("system_id")[metric].idxmax()]
            df = df[df["system_id"].isin(right)].reset_index(drop=True)
            if "train_system_id" in df.columns:
                df = df.drop(columns="train_system_id")
            if "ligand_rdkit_canonical_smiles" in df.columns:
                df = df.drop(columns="ligand_rdkit_canonical_smiles")
            per_metric_similarities.append(df.set_index("system_id"))
        self.max_similarities = pd.concat(
            per_metric_similarities, join="outer", axis=1
        ).reset_index()
        LOG.info(
            f'compute_train_test_max_similarity: Got max similarities for {self.max_similarities["system_id"].nunique()} systems'
        )
        systems_with_similarities = set(self.max_similarities["system_id"])
        extra_rows = []
        for system in right.difference(systems_with_similarities):
            extra_row = {"system_id": system}
            for metric in SIMILARITY_METRICS:
                extra_row[metric] = np.nan
            extra_rows.append(extra_row)
        if len(extra_rows):
            LOG.info(
                f"compute_train_test_max_similarity: Adding nan similarities for {len(extra_rows)} systems"
            )
            self.max_similarities = pd.concat(
                [self.max_similarities, pd.DataFrame(extra_rows)]
            )
        with pd.option_context("future.no_silent_downcasting", True):
            # FutureWarning: Downcasting object dtype arrays on .fillna, .ffill, .bfill is deprecated and will change in a future version.
            self.max_similarities = self.max_similarities.fillna(0)

    def assign_test_set_quality(self) -> None:
        df = query_index(
            filters=[
                (
                    "system_id",
                    "in",
                    set(
                        self.split_df[
                            self.split_df[self.split_label] == self.test_label
                        ]["system_id"]
                    ),
                )
            ],  # type: ignore
            columns=[
                "system_id",
                "system_pass_validation_criteria",
            ],
            splits=["*"],
        ).drop(columns=["split"])
        quality = dict(
            zip(df["system_id"], df["system_pass_validation_criteria"].fillna(False))
        )
        missing_systems = set(
            self.max_similarities[~self.max_similarities["system_id"].isin(quality)][
                "system_id"
            ]
        )
        if len(missing_systems):
            LOG.info(
                f"Discarding {len(missing_systems)} as they are not in the plindex"
            )
            self.max_similarities = self.max_similarities[
                self.max_similarities["system_id"].isin(quality)
            ].reset_index(drop=True)
        self.max_similarities["passes_quality"] = self.max_similarities[
            "system_id"
        ].map(lambda x: quality.get(x, False))
        LOG.info(
            f'assign_test_set_quality: Found {self.max_similarities[self.max_similarities["passes_quality"]]["system_id"].nunique()} '
            f'out of {self.max_similarities["system_id"].nunique()} systems passing quality'
        )


def stratify_cmd(args: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Annotate and stratify a test set")
    parser.add_argument(
        "--split_file",
        type=Path,
        help="Path to split file with [system_id, split] as columns",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output folder where similarity and stratification data are saved",
    )
    parser.add_argument(
        "--split_label",
        type=str,
        default="split",
        help="split=<split_label> is used to get split systems",
    )
    parser.add_argument(
        "--train_label",
        type=str,
        default="train",
        help="split=<train_label> is used to get train systems",
    )
    parser.add_argument(
        "--test_label",
        type=str,
        default="test",
        help="split=<test_label> is used to get test systems",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite max similarity files",
    )

    ns, unknown_args = parser.parse_known_args(args=args)
    if len(unknown_args):
        LOG.warning(f"ignoring arguments {unknown_args}")

    StratifiedTestSet.from_split(
        split_file=Path(ns.split_file),
        output_dir=Path(ns.output_dir),
        split_label=ns.split_label,
        train_label=ns.train_label,
        test_label=ns.test_label,
        overwrite=ns.overwrite,
    )


if __name__ == "__main__":
    stratify_cmd()
