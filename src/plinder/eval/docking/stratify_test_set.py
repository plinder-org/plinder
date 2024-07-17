# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

SIMILARITY_METRICS = (
    # pli
    "pli_qcov",
    # protein
    "protein_seqsim_qcov_weighted_sum",
    "protein_seqsim_weighted_sum",
    "protein_fident_qcov_weighted_sum",
    "protein_fident_weighted_sum",
    "protein_lddt_qcov_weighted_sum",
    "protein_lddt_weighted_sum",
    "protein_qcov_weighted_sum",
    # pocket
    "pocket_fident_qcov",
    "pocket_fident",
    "pocket_lddt_qcov",
    "pocket_lddt",
    "pocket_qcov",
    # ligand
    "tanimoto_similarity_max",
)


@dataclass
class TestCriteria:
    max_entry_resolution: float = 3.5
    max_entry_r: float = 0.4
    max_entry_rfree: float = 0.45
    max_entry_r_minus_rfree: float = 0.05
    ligand_max_num_unresolved_heavy_atoms: int = 0
    ligand_max_alt_count: int = 1
    ligand_min_average_occupancy: float = 0.8
    ligand_min_average_rscc: float = 0.8
    ligand_max_average_rsr: float = 0.3
    ligand_max_percent_outliers_clashes: float = 0
    pocket_max_num_unresolved_heavy_atoms: int = 0
    pocket_max_alt_count: int = 1
    pocket_min_average_occupancy: float = 0.8
    pocket_min_average_rscc: float = 0.8
    pocket_max_average_rsr: float = 0.3
    pocket_max_percent_outliers_clashes: int = 100


def get_high_quality_systems(row: pd.Series, criteria: TestCriteria) -> bool:
    if row.system_type != "holo":
        return False
    if row.entry_r is not None and row.system_ligand_average_rscc is not None:
        quality = [
            # ENTRY
            row.entry_resolution <= criteria.max_entry_resolution,
            row.entry_r <= criteria.max_entry_r,
            row.entry_rfree <= criteria.max_entry_rfree,
            row.entry_r_minus_rfree <= criteria.max_entry_r_minus_rfree,
            # LIGAND
            row.system_ligand_num_unresolved_heavy_atoms
            <= row.system_num_covalent_ligands
            + criteria.ligand_max_num_unresolved_heavy_atoms,
            row.system_ligand_max_alt_count
            <= criteria.ligand_max_alt_count,  # NOTE: max_alt_count is misnomer - this counts number of total conformers!
            row.system_ligand_average_occupancy
            >= criteria.ligand_min_average_occupancy,
            row.system_ligand_average_rscc >= criteria.ligand_min_average_rscc,
            row.system_ligand_average_rsr <= criteria.ligand_max_average_rsr,
            row.system_ligand_percent_outliers_clashes
            <= criteria.ligand_max_percent_outliers_clashes,
            # POCKET
            row.system_pocket_num_unresolved_heavy_atoms
            <= criteria.pocket_max_num_unresolved_heavy_atoms,
            row.system_pocket_max_alt_count <= criteria.pocket_max_alt_count,
            row.system_pocket_average_occupancy
            >= criteria.pocket_min_average_occupancy,
            row.system_pocket_average_rscc >= criteria.pocket_min_average_rscc,
            row.system_pocket_average_rsr <= criteria.pocket_max_average_rsr,
            row.system_pocket_percent_outliers_clashes
            <= criteria.pocket_max_percent_outliers_clashes,
        ]
        if np.logical_and.reduce(quality):
            return True
    return False


def compute_ligand_max_similarities(
    data_dir: Path, left: set[str], right: set[str], metric: str, output_file: Path
) -> None:
    ligands_per_system = pd.read_parquet(
        data_dir / "fingerprints/ligands_per_system.parquet"
    )
    annotation_df = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
        filters=[
            ("system_type", "==", "holo"),
            ("system_num_interacting_protein_chains", "<=", 5),
            ("system_num_ligand_chains", "<=", 5),
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

    left_ligand_ids = set(
        annotation_df[annotation_df["system_id"].isin(left)]["number_id_by_inchikeys"]
    )
    right_ligand_ids = set(
        annotation_df[annotation_df["system_id"].isin(right)]["number_id_by_inchikeys"]
    )

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
    updated_query_ligand_ids = []
    for q, t in zip(ligands["query_ligand_id"], ligands["target_ligand_id"]):
        if t in right_ligand_ids:
            updated_query_ligand_ids.append(t)
        else:
            updated_query_ligand_ids.append(q)
    ligands["updated_query_ligand_id"] = updated_query_ligand_ids

    idx = ligands.groupby("updated_query_ligand_id")[metric].idxmax()
    ligands = ligands.loc[idx]
    ligand_to_system = {}
    for ligand_id, group in ligands_per_system.groupby("number_id_by_inchikeys"):
        ligand_to_system[ligand_id] = set(group["system_id"])
    ligands["query_system"] = ligands["query_ligand_id"].map(ligand_to_system)
    (
        ligands.explode("query_system")
        .rename(
            columns={
                "query_system": "system_id",
            }
        )
        .drop_duplicates("system_id")[["system_id", metric]]
        .reset_index(drop=True)
    ).to_parquet(output_file)


def compute_max_similarities(
    data_dir: Path, left: set[str], right: set[str], metric: str, output_file: Path
) -> None:
    if metric == "tanimoto_similarity_max":
        compute_ligand_max_similarities(data_dir, left, right, metric, output_file)
    else:
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
        updated_target_systems = []
        for q, t in zip(scores["query_system"], scores["target_system"]):
            if t in right:
                updated_query_systems.append(t)
                updated_target_systems.append(q)
            else:
                updated_query_systems.append(q)
                updated_target_systems.append(t)

        scores["updated_query_system"] = updated_query_systems
        scores["updated_target_system"] = updated_target_systems
        idx = scores.groupby("updated_query_system")["similarity"].idxmax()
        scores = scores.loc[idx][
            ["updated_query_system", "similarity", "updated_target_system"]
        ].rename(
            columns={
                "updated_query_system": "system_id",
                "updated_target_system": "train_system_id",
                "similarity": metric,
            }
        )
        scores.reset_index(drop=True).to_parquet(output_file)


@dataclass
class StratifiedTestSet:
    split_df: pd.DataFrame
    data_dir: Path
    output_dir: Path
    train_label: str = "train"
    test_label: str = "test"
    similarity_thresholds: dict[str, int] = field(
        default_factory=lambda: dict(
            pli_qcov=50,
            pocket_lddt_qcov=50,
            pocket_lddt=50,
            pocket_qcov=50,
            protein_seqsim_weighted_sum=30,
            protein_lddt_qcov_weighted_sum=70,
            tanimoto_similarity_max=30,
        )
    )
    similarity_combinations: dict[str, list[str]] = field(
        default_factory=lambda: {
            "novel_pocket_pli": ["pli_qcov", "pocket_qcov", "pocket_lddt_qcov"],
            "novel_pocket_ligand": [
                "pli_qcov",
                "pocket_qcov",
                "pocket_lddt_qcov",
                "tanimoto_similarity_max",
            ],
            "novel_protein": [
                "protein_seqsim_weighted_sum",
                "protein_lddt_qcov_weighted_sum",
            ],
            "novel_all": [
                "pli_qcov",
                "pocket_qcov",
                "pocket_lddt",
                "protein_seqsim_weighted_sum",
                "protein_lddt_qcov_weighted_sum",
                "tanimoto_similarity_max",
            ],
        }
    )
    max_similarities: pd.DataFrame = pd.DataFrame(
        columns=["system_id"] + list(SIMILARITY_METRICS)
    )

    @classmethod
    def from_split(
        cls,
        split_file: Path,
        data_dir: Path,
        output_dir: Path,
        train_label: str = "train",
        test_label: str = "test",
        num_processes: int = 8,
        overwrite: bool = False,
    ) -> "StratifiedTestSet":
        if split_file.name.endswith(".csv"):
            split_df = pd.read_csv(split_file)
        else:
            split_df = pd.read_parquet(split_file)
        assert all(x in split_df.columns for x in ["split", "system_id"])
        data = cls(
            split_df=split_df,
            data_dir=data_dir,
            output_dir=output_dir,
            train_label=train_label,
            test_label=test_label,
        )
        data.output_dir.mkdir(exist_ok=True)
        data.compute_train_test_max_similarity(
            num_processes=num_processes, overwrite=overwrite
        )
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
        self, num_processes: int, overwrite: bool = False
    ) -> pd.DataFrame:
        left, right = (
            set(self.split_df[self.split_df["split"] == self.train_label]["system_id"]),
            set(self.split_df[self.split_df["split"] == self.test_label]["system_id"]),
        )
        LOG.info(
            f"compute_train_test_max_similarity: Found {len(left)} train and {len(right)} test systems"
        )
        multiprocessing.set_start_method("spawn")
        with multiprocessing.Pool(num_processes) as p:
            p.starmap(
                compute_max_similarities,
                [
                    (
                        self.data_dir,
                        left,
                        right,
                        metric,
                        self.get_filename(metric),
                    )
                    for metric in SIMILARITY_METRICS
                    if overwrite or not (self.get_filename(metric)).exists()
                ],
            )
        per_metric_similarities = []
        for metric in SIMILARITY_METRICS:
            df = pd.read_parquet(self.get_filename(metric))
            df = df.loc[df.groupby("system_id")[metric].idxmax()]
            df = df[df["system_id"].isin(right)].reset_index(drop=True)
            if "train_system_id" in df.columns:
                df = df.drop(columns="train_system_id")
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
        self.max_similarities = self.max_similarities.fillna(0)

    def assign_test_set_quality(self) -> None:
        df = pd.read_parquet(self.data_dir / "index" / "annotation_table.parquet")
        df = df[
            df["system_id"].isin(
                self.split_df[self.split_df["split"] == self.test_label]["system_id"]
            )
        ]
        df["system_num_covalent_ligands"] = df["system_id"].map(
            df.groupby("system_id")["ligand_is_covalent"].sum()
        )
        criteria = TestCriteria()
        df["quality"] = df.apply(get_high_quality_systems, criteria=criteria, axis=1)
        quality = dict(zip(df["system_id"], df["quality"]))
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


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Annotate and stratify a test set")
    parser.add_argument(
        "--split_file",
        type=Path,
        help="Path to split file with [system_id, split] as columns",
    )
    parser.add_argument(
        "--data_dir",
        type=Path,
        help="Path to plinder data",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output folder where similarity and stratification data are saved",
    )
    parser.add_argument(
        "--num_processes",
        type=int,
        default=8,
        help="Number of processes to use",
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

    args = parser.parse_args()

    StratifiedTestSet.from_split(
        split_file=args.split_file,
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        train_label=args.train_label,
        test_label=args.test_label,
        num_processes=args.num_processes,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
