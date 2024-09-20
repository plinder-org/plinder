# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

from plinder.core.utils.log import setup_logger
from plinder.eval.docking.stratify_test_set import STRATIFICATION_LABELS

LOG = setup_logger(__name__)

# TODO: need a better way to keep this in sync with stratify_test_set.py SIMILARITY_METRICS
METRICS_DICT = {
    "pocket_qcov": "POCKET SHARED",
    "pocket_lddt": "POCKET LDDT",
    "protein_lddt_weighted_sum": "PROTEIN LDDT",
    "protein_seqsim_weighted_sum": "PROTEIN SEQSIM",
    "pli_unique_qcov": "PLI SHARED",
    "tanimoto_similarity_max": "LIGAND SIMILARITY",
}

ANALYSES = {
    "S": "Success Rate",
    "XS": "Excess Success Rate",
    "EF": "Enrichment Factor",
    "delta_RMSD": "Change in mean RMSD",
    "delta_lDDT_PLI": "Change in mean lDDT-PLI",
    "CDF": "Cumulative Leakage Fraction",
}


@dataclass
class EvaluationResults:
    data: pd.DataFrame
    output_dir: Path

    @classmethod
    def from_scores_and_data_files(
        cls, score_file: Path, data_file: Path, output_dir: Path, top_n: int = 10
    ) -> "EvaluationResults":
        scores_df = pd.read_parquet(score_file)
        data_df = pd.read_parquet(data_file)
        merged_df = scores_df[scores_df["reference"].isin(data_df["system_id"])].merge(
            data_df, left_on="reference", right_on="system_id", how="left"
        )
        merged_df["bisy_rmsd_wave"] = merged_df["bisy_rmsd_wave"].fillna(
            np.nanmax(merged_df["bisy_rmsd_wave"])
        )
        merged_df["lddt_pli_wave"] = merged_df["lddt_pli_wave"].fillna(0)
        merged_df["success"] = merged_df["bisy_rmsd_wave"] <= 2
        merged_df["rank"] = merged_df["rank"].astype(int)
        output_dir.mkdir(exist_ok=True)
        merged_df.to_parquet(output_dir / "merged.parquet")
        data = cls(merged_df, output_dir)
        data.write_results_table(True, top_n)
        data.write_results_table(False, top_n)
        data.write_leakage_plots(True, top_n)
        data.write_leakage_plots(False, top_n)
        return data

    def write_results_table(
        self, only_passes_quality: bool = True, top_n: int = 10
    ) -> None:
        results = []
        sub_df = self.data[self.data["rank"] <= top_n].reset_index(drop=True)
        output_suffix = ""
        if only_passes_quality:
            sub_df = sub_df[sub_df["passes_quality"]].reset_index(drop=True)
            output_suffix += "_quality"
            if not len(sub_df):
                LOG.warning("No systems passed quality control")
                return
        for label in ["all"] + STRATIFICATION_LABELS:
            if label == "all":
                sub_df_label = sub_df.copy(deep=True)
            else:
                sub_df_label = sub_df[sub_df[label]].reset_index(drop=True)
                if not len(sub_df_label):
                    LOG.warning(f"No systems in {label} category at top {top_n}")
                    continue
            row = {
                "Subset": label,
                "No. systems": sub_df_label["reference"].nunique(),
                "Top n": top_n,
            }
            sub_df_top_n = sub_df_label.loc[
                sub_df_label[sub_df_label["rank"] <= top_n]
                .groupby("reference")["bisy_rmsd_wave"]
                .idxmin()
            ]
            row["Success Rate (%)"] = (
                100
                * sub_df_top_n[sub_df_top_n["success"]]["reference"].nunique()
                / sub_df_top_n["reference"].nunique()
            )
            row["Median RMSD"] = np.median(sub_df_top_n["bisy_rmsd_wave"])
            row["Stdev RMSD"] = np.std(sub_df_top_n["bisy_rmsd_wave"])
            row["Mean lDDT-PLI"] = np.median(sub_df_top_n["lddt_pli_wave"])
            row["Stdev lDDT-PLI"] = np.std(sub_df_top_n["lddt_pli_wave"])
            results.append(row)
        result_df = pd.DataFrame(results)
        result_df.to_csv(self.output_dir / f"results{output_suffix}.csv", index=False)

    def write_leakage_plots(
        self, only_passes_quality: bool = True, top_n: int = 1
    ) -> None:
        sub_df = self.data[self.data["rank"] <= top_n].reset_index(drop=True)
        output_suffix = f"_topn{top_n}"
        if only_passes_quality:
            sub_df = sub_df[sub_df["passes_quality"]].reset_index(drop=True)
            if not len(sub_df):
                LOG.warning(f"No systems passed quality control at top {top_n}")
                return
            output_suffix += "_quality"
        df_analysis = perf_vs_traindist_all(sub_df)
        for plot_type in ANALYSES:
            plot_leakage_vs_perfomance(
                df_analysis,
                analysis_Y=plot_type,
                output_html_file=self.output_dir / f"{plot_type}{output_suffix}.html",
            )


def perf_vs_traindist(
    df: pd.DataFrame,
    metric: str = "pli_qcov",
    bins: np.ndarray = np.linspace(0, 100, 11),
) -> dict[str, np.ndarray]:
    sr_total = sum(df["success"]) / len(df)
    mean_rmsd_total = sum(df["bisy_rmsd_wave"]) / len(df)
    mean_lddt_pli_total = sum(df["lddt_pli_wave"]) / len(df)
    y: dict[str, list[float]] = {x: [] for x in ANALYSES}
    for dist in bins:
        train_leaked = df[metric] >= 100 - dist
        sr = sum(df[train_leaked]["success"]) / np.max([sum(train_leaked), 1])
        mean_rmsd = sum(df["bisy_rmsd_wave"][train_leaked]) / np.max(
            [sum(train_leaked), 1]
        )
        mean_lddt_pli = sum(df["lddt_pli_wave"][train_leaked]) / np.max(
            [sum(train_leaked), 1]
        )
        fraction_leaked = sum(train_leaked) / len(train_leaked)
        y["S"].append(sr)
        y["XS"].append(sr - sr_total)
        y["EF"].append(sr / sr_total)
        y["delta_RMSD"].append(mean_rmsd - mean_rmsd_total)
        y["delta_lDDT_PLI"].append(mean_lddt_pli_total - mean_lddt_pli)
        y["CDF"].append(fraction_leaked)
    Y = [y]
    trends = {k: np.array([yi[k] for yi in Y]) for k in Y[0].keys()}
    return trends


def perf_vs_traindist_all(
    df: pd.DataFrame, xbins: np.ndarray = np.linspace(0, 100, 11)
) -> pd.DataFrame:
    data = []
    for metric in METRICS_DICT:
        trends = perf_vs_traindist(df, metric=metric, bins=xbins)
        for k, v in trends.items():
            data += [
                {
                    "metric": metric,
                    "cutoff": x,
                    "mean": np.mean(v[:, ix]),
                    "err": np.std(v[:, ix]) / np.sqrt(len(v[:, ix])),
                    "analysis": k,
                }
                for ix, x in enumerate(xbins)
            ]
    return pd.DataFrame.from_records(data)


def plot_leakage_vs_perfomance(
    df_analysis: pd.DataFrame,
    output_html_file: Path,
    analysis_Y: str,
    plot_dimension: tuple[int, int] = (400, 650),
) -> None:
    xaxis_label = "Distance to Train Set"

    symbols = ["circle", "square", "star", "diamond", "hourglass", "pentagon", "cross"]
    colors = ["lightblue", "pink", "lightgreen", "purple", "goldenrod", "red", "blue"]

    df_plot = df_analysis[df_analysis.analysis == analysis_Y].reset_index(drop=True)
    df_plot.metric = df_plot.metric.map(METRICS_DICT)
    df_plot = df_plot.rename(
        columns={
            "cutoff": xaxis_label,
            "mean": ANALYSES[analysis_Y],
            "metric": "Metric",
        }
    )
    no_of_colors = len(df_plot.Metric.unique())

    fig = px.scatter(
        df_plot,
        x=xaxis_label,
        y=ANALYSES[analysis_Y],
        color="Metric",
        error_y="err",
        color_discrete_sequence=colors[:no_of_colors],
        symbol="Metric",
        symbol_sequence=symbols[:no_of_colors],
        height=plot_dimension[0],
        width=plot_dimension[1],
    )

    fig.update_layout(
        {"plot_bgcolor": "rgba(0, 0, 0, 0)"},
        legend=dict(
            traceorder="normal",
            font=dict(size=16, color="black"),
        ),
        font=dict(size=18),
        yaxis=dict(tickfont=dict(size=18)),
        xaxis=dict(tickfont=dict(size=18)),
        legend_title=None,
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor="black", color="black")
    fig.update_yaxes(showline=True, linewidth=2, linecolor="black", color="black")
    fig.update_traces(
        marker=dict(size=12, line=dict(width=1, color="DarkSlateGrey")),
        selector=dict(mode="markers"),
    )
    fig.write_html(output_html_file)


def plot_cmd(args: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Write evaluation results for a stratified test set"
    )
    parser.add_argument(
        "--score_file",
        type=Path,
        help="Path to scores file generated by plinder_eval",
    )
    parser.add_argument(
        "--data_file",
        type=Path,
        help="Path to test data file generated by plinder_stratify",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output folder where result tables and plots are saved",
    )
    parser.add_argument(
        "--max_top_n",
        type=int,
        default=1,
        help="Maximum rank to use for Top n calculations",
    )

    ns, unknown_args = parser.parse_known_args(args=args)
    if len(unknown_args):
        LOG.warning(f"ignoring arguments {unknown_args}")
    Path(ns.output_dir).mkdir(parents=True, exist_ok=True)
    EvaluationResults.from_scores_and_data_files(
        score_file=ns.score_file,
        data_file=ns.data_file,
        output_dir=ns.output_dir,
        top_n=ns.max_top_n,
    )


if __name__ == "__main__":
    plot_cmd()
