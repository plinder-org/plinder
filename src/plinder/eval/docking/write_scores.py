# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import multiprocessing
from pathlib import Path

import numpy as np
import ost
import pandas as pd

from plinder.core.system.system import PlinderSystem
from plinder.core.utils.log import setup_logger
from plinder.eval.docking import utils

ost.PushVerbosityLevel(-1)
LOG = setup_logger(__name__)


def write_scores_as_json(
    scorer_input: pd.Series,
    output_file: Path,
    overwrite: bool = False,
    predictions_dir: Path | None = None,
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
    """
    try:
        if not overwrite and output_file.exists():
            LOG.warning(f"get_scores: {output_file} exists")
            return
        reference_system = PlinderSystem(system_id=scorer_input.reference_system_id)
        import sys

        print(scorer_input, file=sys.stderr, flush=True)
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
        scores = utils.ModelScores.from_files(
            scorer_input.id,
            receptor_file,
            [
                ligand_file,
            ],  # TODO: change to accept multi-ligand prediction input
            reference_system,
            score_protein=False,  # TODO: add to ScorerConfig class
            score_posebusters=False,  # TODO: add to ScorerConfig class
        ).summarize_scores()
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
    )


if __name__ == "__main__":
    evaluate_cmd()
