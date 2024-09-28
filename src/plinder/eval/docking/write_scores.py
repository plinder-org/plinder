# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import multiprocessing
from pathlib import Path
from typing import Any

import numpy as np
import ost
import pandas as pd

from plinder.core import PlinderSystem
from plinder.core.utils.log import setup_logger
from plinder.eval.docking import utils

ost.PushVerbosityLevel(-1)
LOG = setup_logger(__name__)


def evaluate(
    model_system_id: str,
    reference_system_id: str,
    receptor_file: Path | str,
    ligand_files: list[Path | str],
    predictions_dir: Path | None = None,
    flexible: bool = False,
    posebusters: bool = False,
    posebusters_full: bool = False,
) -> dict[str, Any]:
    """
    Evaluate a single receptor - ligand pair

    Parameters
    ----------
    model_system_id: str
        The ID of the model system (e.g "prediction_0_rank1")
    reference_system_id: str
        The PLINDER systemID of the reference system
    receptor_file: Path
        The path to the receptor CIF/PDB file
    ligand_files: list[Path]
        The path to the ligand SDF files
    predictions_dir: Path | None
        The path to the directory containing the predictions (used if receptor/ligand files are not provided)
    flexible: bool
        Run protein scoring, for flexible docking and cofolding
    posebusters: bool
        Run posebuster scoring
    posebusters_full: bool
        Run posebuster scoring and return full report
    """
    reference_system = PlinderSystem(system_id=reference_system_id)
    receptor_file = Path(receptor_file)
    ligand_file_paths = [Path(ligand_file) for ligand_file in ligand_files]
    if not receptor_file.exists():
        if predictions_dir is not None and (predictions_dir / receptor_file).exists():
            receptor_file = predictions_dir / receptor_file
    if not receptor_file.exists():
        raise FileNotFoundError(f"Receptor file {receptor_file} could not be found")
    if not all(ligand_file.exists() for ligand_file in ligand_file_paths):
        if predictions_dir is not None:
            ligand_file_paths = [
                predictions_dir / ligand_file for ligand_file in ligand_file_paths
            ]
    assert ligand_file_paths is not None and all(
        ligand_file.exists() for ligand_file in ligand_file_paths
    ), f"Ligand files {ligand_file_paths} could not be found"
    return utils.ModelScores.from_model_files(
        model_system_id,
        receptor_file,
        ligand_file_paths,
        reference_system,
        score_protein=flexible,
        score_posebusters=posebusters,
        score_posebusters_full_report=posebusters_full,
    ).summarize_scores()


def write_scores_as_json(
    scorer_input: pd.Series,
    output_file: Path,
    overwrite: bool = False,
    predictions_dir: Path | None = None,
    flexible: bool = False,
    posebusters: bool = False,
    posebusters_full: bool = False,
) -> None:
    """
    Makes score file for a single reference - model pair

    Parameters
    ----------
    scorer_input: pd.Series
        row with [id, reference_system_id, receptor_file, ligand_file, rank, confidence] attributes
    output_file: Path
        Path to save score json file
    overwrite: bool
        Skips existing if False
    predictions_dir: Path | None
        Path to predictions directory
    flexible: bool
        Run protein scoring, for flexible docking and cofolding
    posebusters: bool
        Run posebuster scoring
    posebusters_full: bool
        Run posebuster scoring and return full report
    """
    try:
        if not overwrite and output_file.exists():
            LOG.warning(f"get_scores: {output_file} exists")
            return
        reference_system = PlinderSystem(system_id=scorer_input.reference_system_id)
        receptor_file = None
        if scorer_input.receptor_file is not None and not np.isnan(
            scorer_input.receptor_file
        ):
            receptor_file = Path(scorer_input.receptor_file)
        ligand_file = Path(scorer_input.ligand_file)

        if receptor_file is not None and not receptor_file.exists():
            if (
                predictions_dir is not None
                and (predictions_dir / receptor_file).exists()
            ):
                receptor_file = predictions_dir / receptor_file
            else:
                assert reference_system.receptor_cif is not None
                receptor_file = Path(reference_system.receptor_cif)
        else:
            receptor_file = Path(reference_system.receptor_cif)
        if ligand_file is not None and not ligand_file.exists():
            if predictions_dir is not None and (predictions_dir / ligand_file).exists():
                ligand_file = predictions_dir / ligand_file
        assert receptor_file is not None
        assert ligand_file is not None
        scores = evaluate(
            scorer_input.id,
            scorer_input.reference_system_id,
            receptor_file,
            [
                ligand_file,
            ],  # TODO: change to accept multi-ligand prediction input
            predictions_dir,
            flexible,
            posebusters,
            posebusters_full,
        )
        output_file.parent.mkdir(exist_ok=True, parents=True)
        with open(output_file, "w") as f:
            json.dump(scores, f)
    except Exception as e:
        LOG.error(f"get_scores: {e}")
        return


def score_test_set(
    prediction_file: Path,
    output_dir: Path,
    output_file: Path,
    num_processes: int = 8,
    overwrite: bool = False,
    posebusters: bool = False,
    flexible: bool = False,
    posebusters_full: bool = False,
) -> None:
    """
    Makes score files for a test set

    Parameters
    ----------
    prediction_file: Path
        A CSV file with the columns [id, reference_system_id, receptor_file, ligand_file, rank, confidence]
    system_dir: Path
        Path to *extracted* system files of the test set
    output_dir: Path
        Path to save score files, in the form <output_dir>/<id>/<rank>.json
    output_file: Path
        Path to save all scores
    num_processes: int
        Number of processes to use
    overwrite: bool
        Skips existing files if False
    posebusters: bool
        Run posebuster scoring
    flexible: bool
        Run protein scoring, for flexible docking and cofolding
    posebusters_full: bool
        Run posebuster scoring and return full report
    """
    predictions = pd.read_csv(prediction_file)
    predictions["output_file"] = predictions.apply(
        lambda row: output_dir / row["id"] / f"{row['rank']}.json", axis=1
    )
    predictions["exists"] = predictions["output_file"].map(lambda x: x.exists())
    if not overwrite:
        LOG.info(
            f"Skipping {predictions['exists'].sum()} case as score files exists and overwrite=False"
        )
        predictions = predictions[~predictions["exists"]]
    LOG.info(f"Running evaluation on {predictions.shape[0]} cases")
    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(num_processes) as p:
        p.starmap(
            write_scores_as_json,
            [
                (
                    row,
                    row.output_file,
                    overwrite,
                    prediction_file.parent,
                    flexible,
                    posebusters,
                    posebusters_full,
                )
                for _, row in predictions.iterrows()
            ],
        )
    scores = []
    for model_dir in output_dir.iterdir():
        for json_file in model_dir.iterdir():
            with open(json_file) as f:
                try:
                    s = json.load(f)
                except Exception as e:
                    LOG.error(
                        f"score_test_set: Error loading scores file {json_file}: {e}"
                    )
                s["rank"] = json_file.stem
                scores.append(s)
    pd.DataFrame(scores).to_parquet(output_file, index=False)


def evaluate_cmd(args: list[str] | None = None) -> None:
    """
    Evaluate a set of docking/cofolding predictions against a set of reference Plinder systems.
    """
    import argparse

    parser = argparse.ArgumentParser(usage=evaluate_cmd.__doc__)
    parser.add_argument(
        "--prediction_file",
        type=Path,
        help="Path to prediction file with [id, reference_system_id, receptor_file, ligand_file, rank, confidence] as columns",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output folder where score JSON files are saved",
    )
    parser.add_argument(
        "--num_processes",
        type=int,
        help="Number of processes to use",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite max similarity files",
    )
    parser.add_argument(
        "--posebusters",
        action="store_true",
        help="Run posebuster scoring",
    )
    parser.add_argument(
        "--flexible",
        action="store_true",
        help="Run protein scoring, for flexible docking and cofolding",
    )
    parser.add_argument(
        "--posebusters_full",
        action="store_true",
        help="Run posebuster scoring and return full report",
    )
    ns, unknown_args = parser.parse_known_args(args=args)
    if len(unknown_args):
        LOG.warning(f"ignoring arguments {unknown_args}")
    LOG.info(f"Evaluating {ns.prediction_file} and saving to {ns.output_dir}")
    Path(ns.output_dir).mkdir(parents=True, exist_ok=True)
    Path(ns.output_dir / "scores").mkdir(parents=True, exist_ok=True)

    score_test_set(
        prediction_file=ns.prediction_file,
        output_dir=Path(ns.output_dir) / "scores",
        output_file=Path(ns.output_dir) / "scores.parquet",
        num_processes=ns.num_processes,
        overwrite=ns.overwrite,
        posebusters=ns.posebusters,
        flexible=ns.flexible,
        posebusters_full=ns.posebusters_full,
    )


if __name__ == "__main__":
    evaluate_cmd()
