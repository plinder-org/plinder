# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ks_2samp

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


SYSTEM_PROPS = {
    # ligand
    "ligand_molecular_weight": 50,
    # would be nice to add tpsa, but currently not calculated in annotations
    # 'ligand_tpsa' : 10,
    "ligand_num_heavy_atoms": 5,
    "ligand_num_rot_bonds": 2,
    "ligand_num_rings": 1,
    "ligand_crippen_clogp": 1,
    # pocket
    "system_num_pocket_residues": 20,
    # protein
    "system_interacting_protein_chains_length": None,
}


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


def plot_hist_plot(
    df: pd.DataFrame,
    plot_metric: str,
    output_png_file: Path,
    overwrite: bool = False,
) -> None:
    g = sns.histplot(
        df,
        x=plot_metric,
        hue="split",
        element="step",
        stat="density",
        common_norm=False,
        binwidth=SYSTEM_PROPS[plot_metric],
        bins=None if SYSTEM_PROPS[plot_metric] else np.linspace(0, 4, 50),
        log_scale=plot_metric.startswith("system_interacting_protein_chains_length"),
    )
    if overwrite or not output_png_file.exists():
        g.figure.savefig(output_png_file)
    g.figure.clf()


def plot_molecular_descriptor_splits(
    split_file: Path,
    data_dir: Path,
    output_dir: Path,
    train_label: str = "train",
    val_label: str = "val",
    test_label: str = "test",
    overwrite: bool = False,
) -> None:
    if split_file.name.endswith(".csv"):
        split_df = pd.read_csv(split_file)
    else:
        split_df = pd.read_parquet(split_file)
        assert all(x in split_df.columns for x in ["split", "system_id"])

    # remove other splits that are not used
    split_df = split_df[split_df.split.isin([train_label, val_label, test_label])]
    # add annotations
    annotation_df = pd.read_parquet(
        data_dir / "index" / "annotation_table.parquet",
        filters=[
            ("system_type", "==", "holo"),
            ("system_num_interacting_protein_chains", "<=", 5),
            ("system_num_ligand_chains", "<=", 5),
        ],
    )
    # sum all residues
    annotation_df["system_interacting_protein_chains_length"] = annotation_df[
        "system_interacting_protein_chains_length"
    ].apply(lambda x: sum(int(i) for i in x.split(";")))

    split_df = split_df.merge(annotation_df, on="system_id", how="left")

    for prop in SYSTEM_PROPS:
        stat, pval = get_ks_test(
            split_df[split_df.split == val_label],
            split_df[split_df.split == test_label],
            metric=prop,
        )
        LOG.info(
            f"KS test statistic for {prop} {val_label}/{test_label}: {stat} (p-value {pval})"
        )
        output_png_file = output_dir / f"hist1d_{prop}.png"
        plot_hist_plot(split_df, prop, output_png_file, overwrite=overwrite)


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
        "--train_label",
        type=str,
        default="train",
        help="split=<train_label> is used to get train systems",
    )
    parser.add_argument(
        "--val_label",
        type=str,
        default="val",
        help="split=<val_label> is used to get val systems",
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

    plot_molecular_descriptor_splits(
        split_file=args.split_file,
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        train_label=args.train_label,
        val_label=args.val_label,
        test_label=args.test_label,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
