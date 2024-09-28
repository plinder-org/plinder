# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from plinder.core.index.utils import get_plindex

try:
    import matplotlib.pyplot as plt
    import mols2grid
    import seaborn as sns
    from matplotlib_venn import venn3
    from scipy.stats import ks_2samp
except ImportError:
    raise ImportError("Please run: pip install plinder[plots] to use this module")


from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def get_ks_test(
    df1: pd.DataFrame, df2: pd.DataFrame, metric: str
) -> tuple[float, float]:
    res = ks_2samp(
        df1[metric],
        df2[metric],
        method="auto",
        axis=0,
        nan_policy="omit",
        keepdims=False,
    )
    return res.statistic, res.pvalue


def label_has_mms_within_subset(
    mms_df: pd.DataFrame, subset: set[str], minimum_count: int
) -> tuple[pd.DataFrame, set[str]]:
    mms_df = mms_df[mms_df["system_id"].isin(subset)].reset_index(drop=True)
    LOG.info(
        f"mmp series after subset {len(mms_df.index)} records for {mms_df['system_id'].nunique()} systems"
    )
    # count number of unique ligands (not counting repeats)
    mms_count = mms_df.groupby("congeneric_id")["congeneric_ligand_ccd_code"].agg(
        "nunique"
    )
    mms_df["mms_unique_count"] = mms_df["congeneric_id"].map(mms_count)
    good_mms_systems = set(
        mms_df[mms_df["mms_unique_count"] >= minimum_count]["system_id"]
    )
    return mms_df, good_mms_systems


@dataclass
class SplitPropertiesPlotter:
    data_dir: Path
    split_file: Path
    output_dir: Path
    stratified_files: dict[str, Path] = field(default_factory=dict)
    plindex: pd.DataFrame = field(init=False)
    system_plindex: pd.DataFrame = field(init=False)
    mms_count: int = field(default=3)
    colors: dict[str, str] = field(
        default_factory=lambda: {
            "train": "#ff9999",  # Light red
            "test": "#66b3ff",  # Light blue
            "val": "#99ff99",  # Light green
            "removed": "#d3d3d3",  # Light grey
            "train_text": "red",
            "test_text": "blue",
            "val_text": "green",
            "False": "#D3D3D3",
            "True": "#99ff99",
        }
    )
    system_descriptors: dict[str, int | None] = field(
        default_factory=lambda: {
            # ligand
            "ligand_molecular_weight": 50,
            "ligand_tpsa": 10,
            "ligand_num_heavy_atoms": 5,
            "ligand_num_rot_bonds": 2,
            "ligand_num_rings": 1,
            "ligand_crippen_clogp": 1,
            # pocket
            "system_proper_num_pocket_residues": 5,
            # protein
            "system_protein_chains_total_length": None,
        }
    )
    priority_columns: list[str] = field(
        default_factory=lambda: [
            "system_has_binding_affinity",
            "system_has_mms",
            "system_has_apo_or_pred",
            "system_ligand_has_lipinski",
            "system_pass_validation_criteria",
            "system_pass_statistics_criteria",
        ]
    )
    domain_columns: list[str] = field(
        default_factory=lambda: [
            "system_pocket_ECOD_t_name",
            "system_pocket_CATH",
            "system_pocket_SCOP2B",
            "system_pocket_Pfam",
        ]
    )
    ligand_types: list[str] = field(
        default_factory=lambda: [
            f"system_ligand_has_{x}"
            for x in [
                "lipinski",
                "covalent",
                "cofactor",
                "oligo",
                "ion",
                "fragment",
                "artifact",
                "other",
            ]
        ]
    )
    diversity_columns: dict[str, str] = field(
        default_factory=lambda: {
            "system_id": "Systems",
            "entry_pdb_id": "PDB IDs",
            "system_proper_ligand_unique_ccd_codes": "Ligand CCD codes",
            "pli_unique_qcov__50__communities": "PLI communities",
            "pocket_qcov__50__communities": "Pocket communities",
            "tanimoto_similarity_max__50__communities": "Ligand communities",
        }
    )

    @classmethod
    def from_files(
        cls,
        data_dir: Path,
        split_file: Path,
        output_dir: Path = Path("split_plots"),
        plindex_file: Path | None = None,
        stratified_train_test_file: Path | None = None,
        stratified_train_val_file: Path | None = None,
        stratified_val_test_file: Path | None = None,
        mms_count: int = 3,
        make_plots: bool = True,
    ) -> "SplitPropertiesPlotter":
        stratified_files = {}
        if stratified_train_test_file is not None:
            stratified_files["train_vs_test"] = stratified_train_test_file
        if stratified_train_val_file is not None:
            stratified_files["train_vs_val"] = stratified_train_val_file
        if stratified_val_test_file is not None:
            stratified_files["val_vs_test"] = stratified_val_test_file
        plotter = cls(
            data_dir=data_dir,
            split_file=split_file,
            stratified_files=stratified_files,
            output_dir=output_dir,
            mms_count=mms_count,
        )
        if plindex_file is None:
            plindex = get_plindex().drop(columns=["split"])
        else:
            plindex = pd.read_parquet(plindex_file)
        plotter.plindex = plotter.merge_splits_and_plindex(plindex)
        plotter.system_plindex = plotter.plindex.drop_duplicates("system_id")
        plotter.merge_stratification()
        if make_plots:
            plotter.plot_all()
        return plotter

    def __post_init__(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def save_ligand_report_html(
        self,
        split_name: str,
        smiles_col: str = "ligand_rdkit_canonical_smiles",
        bg_color_col: str = "system_has_mms",
        sort_col: str = "system_id",
    ) -> None:
        color_col = f"tanimoto_similarity_max__train_vs_{split_name}"
        df = self.plindex[
            (self.plindex["split"] == split_name)
            & (~self.plindex["ligand_is_ion"])
            & (~self.plindex["ligand_is_artifact"])
        ]
        system_df = df.drop_duplicates("system_id")
        grouped_smiles_df = (
            df.groupby("system_id")[smiles_col].agg(list).apply(".".join).reset_index()
        )
        df = system_df.merge(grouped_smiles_df, on="system_id", suffixes=["_", ""])

        grid = mols2grid.MolGrid(
            df,
            smiles_col=smiles_col,
        )
        grid.save(
            self.output_dir / f"{split_name}.html",
            # rename columns for easier interpretation and formatting
            rename={
                "pli_unique_qcov__50__communities": "PLI community ID",
                "tanimoto_similarity_max__50__communities": "Ligand community ID",
                "tanimoto_similarity_max__30__communities": "Ligand community ID",
                "system_num_protein_chains": "Receptor chain count",
                "system_proper_num_ligand_chains": "Ligand count",
            },
            # set what's displayed on the grid
            subset=[
                "system_id",
                "img",
                "Ligand community ID",
                "PLI community ID",
                color_col,
            ],
            # set what's displayed on the hover tooltip
            tooltip=[
                "Receptor chain count",
                "Ligand count",
                "system_has_binding_affinity",
                "system_has_mms",
                "system_has_apo_or_pred",
                "system_ligand_has_lipinski",
                "system_ligand_has_cofactor",
            ],
            n_items_per_page=48,
            # style for the grid labels and tooltips
            style={
                color_col: lambda x: "color: red; font-weight: bold;" if x > 30 else "",
                "__all__": lambda x: "background-color: azure;"
                if x[bg_color_col]
                else "",
            },
            transform={color_col: lambda x: round(x, 0)},
            # sort the grid in a different order by default
            sort_by=sort_col,
        )

    def merge_splits_and_plindex(self, plindex: pd.DataFrame) -> pd.DataFrame:
        split = pd.read_parquet(self.split_file)
        plindex = plindex[plindex["system_id"].isin(split["system_id"])].reset_index(
            drop=True
        )
        plindex = pd.merge(plindex, split, on="system_id", suffixes=("", "_y"))
        mms_df = pd.read_parquet(self.data_dir / "mmp/plinder_mmp_series.parquet")
        split_subdf = plindex[
            plindex["system_pass_validation_criteria"].fillna(False)
        ].drop_duplicates("system_id")
        good_mms_systems = set()
        for split in plindex["split"].unique():
            split_systems = set(
                split_subdf[split_subdf["split"] == split]["system_id"].unique()
            )
            good_mms_systems |= label_has_mms_within_subset(
                mms_df,
                split_systems,
                self.mms_count,
            )[1]
        plindex["system_has_mms"] = plindex["system_id"].isin(good_mms_systems)

        return plindex

    def merge_stratification(self) -> None:
        self.stratified = {
            split: pd.read_parquet(self.stratified_files[split])
            .drop_duplicates("system_id")
            .rename(
                mapper=lambda x: f"{x}__{split}"
                if x != "system_id" and "novel" not in x
                else x,
                axis=1,
            )
            for split in self.stratified_files
        }
        for split in self.stratified:
            self.system_plindex = self.system_plindex.merge(
                self.stratified[split],
                on="system_id",
                how="left",
            )
            self.plindex = self.plindex.merge(
                self.stratified[split],
                on="system_id",
                how="left",
            )

    def plot_ligands(self, split_names: list[str] = ["test", "val"]) -> None:
        plindex = self.plindex[
            self.plindex["ligand_is_proper"] & self.plindex["split"].isin(split_names)
        ]
        for split in split_names:
            ligand_df = plindex[plindex["split"] == split][
                ["system_id", "ligand_rdkit_canonical_smiles"]
            ].reset_index(drop=True)
            if split in self.stratified_files:
                stratified = pd.read_parquet(self.stratified_files[split])
                stratified["novel_ligand"] = stratified["tanimoto_similarity_max"] < 30
                ligand_df = ligand_df.merge(
                    stratified[
                        ["system_id", "tanimoto_similarity_max", "novel_ligand"]
                    ],
                    on="system_id",
                    how="left",
                )

            print(f"{split}: {len(ligand_df)}")

    def plot_split_proportions(self) -> None:
        fig, ax = plt.subplots()
        counts = self.system_plindex["split"].value_counts()
        ax.pie(
            counts,
            labels=counts.index,
            autopct="%1.1f%%",
            colors=[self.colors[x] for x in counts.index],
        )
        ax.set_title("Split proportions")
        plt.savefig(self.output_dir / "split_proportions.png")

    def print_stratification_table(self) -> None:
        if len(self.stratified) == 0:
            LOG.info("No stratified data found, skipping stratification table")
            return
        stratified_dfs = {
            split: self.stratified[split][
                [c for c in self.stratified[split].columns if "novel" in c]
            ]
            .apply(lambda x: x.value_counts())
            .fillna(0)
            .astype(int)
            for split in self.stratified
        }
        df_combined = pd.concat(stratified_dfs.values(), keys=stratified_dfs.keys())
        df_combined = df_combined.reset_index().rename(columns={"level_0": "split"})
        df_combined = (
            df_combined[df_combined["level_1"]]
            .drop(columns=["level_1"])
            .reset_index(drop=True)
        )
        df_combined["total"] = df_combined["split"].apply(
            lambda x: self.system_plindex[
                self.system_plindex["split"] == x.split("_")[-1]
            ].shape[0]
        )
        df_combined["split"] = df_combined["split"].apply(
            lambda x: x.replace("_", " ").upper()
        )
        df_combined.to_csv(self.output_dir / "stratification_table.csv", index=False)

    def print_overall_diversity(self) -> None:
        cluster_df = pd.DataFrame(
            {
                col: self.system_plindex.groupby("split")[col].nunique()
                for col in self.diversity_columns
            }
        ).rename(columns=self.diversity_columns)
        cluster_df.to_csv(self.output_dir / "overall_diversity.csv", index=True)

    def plot_molecular_descriptors(
        self, val_label: str = "val", test_label: str = "test"
    ) -> None:
        fig, axes = self.get_axes(len(self.system_descriptors))
        plindex_to_use = self.plindex[self.plindex["ligand_is_proper"]].reset_index(
            drop=True
        )
        for i, prop in enumerate(self.system_descriptors):
            stat, pval = get_ks_test(
                plindex_to_use[plindex_to_use.split == val_label][[prop]],
                plindex_to_use[plindex_to_use.split == test_label][[prop]],
                metric=prop,
            )
            self.plot_hist_plot(
                plindex_to_use,
                prop,
                axes[i],
            )
            ks_info = f"KS test ({val_label}/{test_label}):\nstatistic = {stat:.3f}\np-value = {pval:.3e}"
            axes[i].text(
                0.95,
                0.05,
                ks_info,
                transform=axes[i].transAxes,
                verticalalignment="bottom",
                horizontalalignment="right",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                fontsize=8,
            )
        self.clear_unused_axes(axes, len(self.system_descriptors))
        plt.tight_layout()
        plt.savefig(self.output_dir / "molecular_descriptors.png")

    def plot_hist_plot(self, df: pd.DataFrame, plot_metric: str, ax: plt.Axes) -> None:
        if plot_metric.startswith("system"):
            df = df.drop(columns=["system_id"])
        sns.histplot(
            df,
            x=plot_metric,
            hue="split",
            hue_order=["train", "val", "test"],
            element="step",
            stat="density",
            common_norm=False,
            binwidth=self.system_descriptors[plot_metric],
            bins=None
            if self.system_descriptors[plot_metric]
            else np.linspace(0, 4, 50),
            log_scale=plot_metric.startswith("system_protein_chains_total_length"),
            ax=ax,
            palette=self.colors,
        )
        ax.set_xlabel("")
        name = (
            plot_metric.replace("system_", "")
            .replace("proper", "")
            .replace("ligand_", "")
            .replace("protein_chains_asym_id", "protein")
            .replace("protein_chains_", "")
            .replace("_", " ")
            .upper()
        )
        ax.set_title(name)

    def plot_priorities(self) -> None:
        fig, axes = self.get_axes(len(self.priority_columns))

        for idx, column in enumerate(self.priority_columns):
            column_data = self.system_plindex[["split", column]].astype(str)
            data = column_data.groupby("split")[column].value_counts()
            percentages = column_data.groupby("split")[column].value_counts(
                normalize=True
            )
            split_order = ["train", "val", "test", "removed"]
            percentages = percentages.reindex(split_order, level=0)
            unstacked = percentages.unstack().loc[split_order].fillna(0) * 100
            unstacked = unstacked[["True", "False"]]

            unstacked.plot(
                kind="bar",
                stacked=True,
                ax=axes[idx],
                color=[self.colors[x] for x in unstacked.columns],
                rot=0,
                legend=False,
            )

            axes[idx].set_title(column.replace("system_", "").replace("_", " ").title())
            axes[idx].set_ylabel("Systems (%)")
            axes[idx].set_ylim(0, 100)
            axes[idx].set_xlabel("")

            names = percentages.index.get_level_values(0).unique()
            for c in axes[idx].containers:
                labels = []
                for v in names:
                    key = (v, c.get_label())
                    if key in data.index:
                        labels.append(int(data[key]))
                    else:
                        labels.append(0)
                axes[idx].bar_label(c, labels=labels, label_type="center", fontsize=12)
        self.clear_unused_axes(axes, len(self.priority_columns))
        fig.suptitle("Priorities across splits", fontsize=16)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper left", ncol=2)
        plt.tight_layout()
        plt.savefig(self.output_dir / "priorities.png")

    def plot_venn(self, column: str, ax: plt.Axes, label: str) -> None:
        per_split = {
            split: self.system_plindex[
                self.system_plindex[column].notna()
                & (self.system_plindex["split"] == split)
            ][["system_id", column]]
            for split in ["train", "val", "test"]
        }

        per_split_groups = {split: set(per_split[split][column]) for split in per_split}

        v = venn3(
            (
                per_split_groups["train"],
                per_split_groups["val"],
                per_split_groups["test"],
            ),
            set_labels=None,
            ax=ax,
            alpha=0.6,
            set_colors=(
                self.colors["train"],
                self.colors["val"],
                self.colors["test"],
            ),
        )

        for idx in ["100", "010", "001", "110", "101", "011", "111"]:
            if v.get_patch_by_id(idx):
                if idx == "100":
                    groups = per_split_groups["train"] - (
                        per_split_groups["val"] | per_split_groups["test"]
                    )
                    systems = {
                        "train": per_split["train"][
                            per_split["train"][column].isin(groups)
                        ]
                    }
                elif idx == "010":
                    groups = per_split_groups["val"] - (
                        per_split_groups["train"] | per_split_groups["test"]
                    )
                    systems = {
                        "val": per_split["val"][per_split["val"][column].isin(groups)]
                    }
                elif idx == "001":
                    groups = per_split_groups["test"] - (
                        per_split_groups["train"] | per_split_groups["val"]
                    )
                    systems = {
                        "test": per_split["test"][
                            per_split["test"][column].isin(groups)
                        ]
                    }
                elif idx == "110":
                    groups = (
                        per_split_groups["train"]
                        & per_split_groups["val"] - per_split_groups["test"]
                    )
                    systems = {
                        "train": per_split["train"][
                            per_split["train"][column].isin(groups)
                        ],
                        "val": per_split["val"][per_split["val"][column].isin(groups)],
                    }
                elif idx == "101":
                    groups = (
                        per_split_groups["train"]
                        & per_split_groups["test"] - per_split_groups["val"]
                    )
                    systems = {
                        "train": per_split["train"][
                            per_split["train"][column].isin(groups)
                        ],
                        "test": per_split["test"][
                            per_split["test"][column].isin(groups)
                        ],
                    }
                elif idx == "011":
                    groups = (
                        per_split_groups["val"]
                        & per_split_groups["test"] - per_split_groups["train"]
                    )
                    systems = {
                        "val": per_split["val"][per_split["val"][column].isin(groups)],
                        "test": per_split["test"][
                            per_split["test"][column].isin(groups)
                        ],
                    }
                elif idx == "111":
                    groups = (
                        per_split_groups["train"]
                        & per_split_groups["val"]
                        & per_split_groups["test"]
                    )
                    systems = {
                        "train": per_split["train"][
                            per_split["train"][column].isin(groups)
                        ],
                        "val": per_split["val"][per_split["val"][column].isin(groups)],
                        "test": per_split["test"][
                            per_split["test"][column].isin(groups)
                        ],
                    }

                n = len(groups)
                if n == 0:
                    v.get_label_by_id(idx).set_text("")
                    continue
                v.get_label_by_id(idx).set_text(str(n))

                x, y = v.get_label_by_id(idx).get_position()
                y += 0.06
                v.get_label_by_id(idx).set_position((x, y))
                for split in systems:
                    n = len(systems[split])
                    if n == 0:
                        continue
                    y -= 0.04
                    ax.text(
                        x,
                        y,
                        str(n),
                        color=self.colors[f"{split}_text"],
                        ha="center",
                        va="center",
                    )

        ax.set_title(label, fontsize=15)

    def plot_domain_classifications(self) -> None:
        fig, axes = self.get_axes(len(self.domain_columns))
        for i, c in enumerate(self.domain_columns):
            self.plot_venn(
                c, axes[i], c.replace("system_pocket_", "").replace("_", " ").upper()
            )
        self.clear_unused_axes(axes, len(self.domain_columns))
        fig.suptitle("Domain Classifications")
        plt.tight_layout()
        plt.savefig(self.output_dir / "domain_classifications.png")

    def get_axes(self, num_plots: int) -> tuple[plt.Figure, np.ndarray]:
        num_cols = int(np.ceil(np.sqrt(num_plots)))
        num_rows = int(np.ceil(num_plots / num_cols))
        fig, axes = plt.subplots(
            num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows)
        )
        return fig, axes.flatten()

    def clear_unused_axes(self, axes: np.ndarray, num_plots: int) -> None:
        for j in range(num_plots, len(axes)):
            axes[j].axis("off")

    def plot_clusters(self, cluster: str = "communities", threshold: int = 50) -> None:
        cluster_names = [
            c
            for c in self.system_plindex.columns
            if c.endswith(f"__{threshold}__{cluster}")
            and "weighted_max" not in c
            and ("tanimoto" in c or "max" not in c)
        ]

        fig, axes = self.get_axes(len(cluster_names))

        for i, c in enumerate(cluster_names):
            name = (
                c.replace(f"__{threshold}__{cluster}", "")
                .replace("weighted_sum", "")
                .replace("max", "")
            )
            if name.startswith("protein"):
                name = name.replace("qcov", "global")
            else:
                name = name.replace("qcov", "shared")
            name = name.replace("_", " ").upper()
            self.plot_venn(
                c,
                axes[i],
                name,
            )

        self.clear_unused_axes(axes, len(cluster_names))

        fig.suptitle(f"Clusters ({cluster.capitalize()}, {threshold}%)", fontsize=16)
        plt.tight_layout()
        plt.savefig(self.output_dir / "plinder_clusters.png")

    def plot_ligand_types(self) -> None:
        counts_per_split = (
            self.system_plindex.groupby("split")[self.ligand_types]
            .mean()
            .mul(100)
            .to_dict()
        )
        labels = [
            c.replace("system_ligand_has_", "").capitalize() for c in self.ligand_types
        ]
        split_names = ["train", "val", "test", "removed"]
        fig, axes = self.get_axes(len(split_names))

        bar_colors = plt.cm.Pastel2.colors
        for i, split in enumerate(split_names):
            split_data = self.system_plindex[self.system_plindex["split"] == split]
            total_systems = len(split_data)
            percentages = [counts_per_split[c][split] for c in self.ligand_types]
            counts = [int(p / 100 * total_systems) for p in percentages]

            bars = axes[i].bar(
                np.arange(len(labels)),
                percentages,
                width=1,
                color=bar_colors,
                edgecolor="black",
                label=split,
                linewidth=1,
            )
            axes[i].set_xticks(np.arange(len(labels)))
            axes[i].set_xticklabels(labels, rotation=70)
            axes[i].set_ylim(0, 100)
            axes[i].set_title(f"{split.capitalize()} (Total: {total_systems})")

            # Add count labels to bars
            for bar, count in zip(bars, counts):
                axes[i].text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bar.get_height() + 2,
                    f"{count}",
                    ha="center",
                    va="bottom",
                    rotation=70,
                    fontsize=10,
                )
        fig.suptitle("Distribution of Ligand Types Contained in Systems")
        self.clear_unused_axes(axes, len(split_names))
        plt.tight_layout()
        plt.savefig(self.output_dir / "ligand_types.png")

    def plot_chain_composition(self) -> None:
        split_names = ["train", "val", "test", "removed"]
        fig, axes = self.get_axes(len(split_names))

        for i, split in enumerate(split_names):
            split_stats = self.system_plindex[self.system_plindex["split"] == split][
                [
                    "system_num_protein_chains",
                    "system_proper_num_ligand_chains",
                ]
            ]
            counts = {
                "Single Protein & Single Ligand": split_stats[
                    (split_stats["system_num_protein_chains"] == 1)
                    & (split_stats["system_proper_num_ligand_chains"] == 1)
                ].shape[0],
                "Single Protein & Multiple Ligands": split_stats[
                    (split_stats["system_num_protein_chains"] == 1)
                    & (split_stats["system_proper_num_ligand_chains"] > 1)
                ].shape[0],
                "Multiple Proteins & Single Ligand": split_stats[
                    (split_stats["system_proper_num_ligand_chains"] == 1)
                    & (split_stats["system_num_protein_chains"] > 1)
                ].shape[0],
                "Multiple Proteins & Multiple Ligands": split_stats[
                    (split_stats["system_proper_num_ligand_chains"] > 1)
                    & (split_stats["system_num_protein_chains"] > 1)
                ].shape[0],
            }

            wedges, texts, autotexts = axes[i].pie(
                list(counts.values()),
                colors=plt.cm.Pastel2.colors,
                autopct=lambda pct: f"{pct:.1f}%\n{int(pct/100.*sum(counts.values())):d}",
                textprops={"fontsize": 8},
                wedgeprops={"linewidth": 0.5, "edgecolor": "black"},
            )
            axes[i].set_title(f"{split.capitalize()} ({sum(counts.values())})")

        self.clear_unused_axes(axes, len(split_names))
        fig.suptitle("System Chain Composition")
        fig.legend(wedges, list(counts.keys()), loc="lower center", ncol=2)
        plt.tight_layout()
        plt.savefig(self.output_dir / "chain_composition.png")

    def plot_all(self) -> None:
        self.plot_split_proportions()
        self.print_stratification_table()
        try:
            if "train_vs_test" in self.stratified.keys():
                self.save_ligand_report_html("test")
            if "train_vs_val" in self.stratified.keys():
                self.save_ligand_report_html("val")
        except Exception as e:
            LOG.error(f"error saving ligand report: {e}")
        self.print_overall_diversity()
        self.plot_molecular_descriptors()
        self.plot_priorities()
        self.plot_domain_classifications()
        self.plot_clusters()
        self.plot_ligand_types()
        self.plot_chain_composition()


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Compare split properties")
    parser.add_argument(
        "--split_file",
        type=Path,
        help="Path to split file",
    )
    parser.add_argument(
        "--data_dir",
        type=Path,
        help="Path to plinder data",
    )
    parser.add_argument(
        "--stratified_train_test_file",
        type=Path,
        help="Path to stratified train vs test file",
    )
    parser.add_argument(
        "--stratified_train_val_file",
        type=Path,
        help="Path to stratified train vs val file",
    )
    parser.add_argument(
        "--stratified_val_test_file",
        type=Path,
        help="Path to stratified val vs test file",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output directory",
        default=Path("split_plots"),
    )
    parser.add_argument(
        "--mms_count",
        type=int,
        help="Number of systems to consider as an MMS",
        default=3,
    )

    args = parser.parse_args()

    SplitPropertiesPlotter.from_files(
        split_file=args.split_file,
        data_dir=args.data_dir,
        stratified_train_test_file=args.stratified_train_test_file,
        stratified_train_val_file=args.stratified_train_val_file,
        stratified_val_test_file=args.stratified_val_test_file,
        mms_count=args.mms_count,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
