# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
import multiprocessing
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import ost
import pandas as pd

from plinder.core.utils.log import setup_logger
from plinder.eval.docking import utils

ost.PushVerbosityLevel(-1)
LOG = setup_logger(__name__)


def write_scores_as_json(
    scorer_input: pd.Series,
    system_dir: Path,
    output_file: Path,
    overwrite: bool = False,
) -> None:
    """
    Makes score file for a single reference - model pair

    Parameters
    ----------
    scorer_input: pd.Series
        row with [id, reference_system_id, receptor_file, ligand_file, rank, confidence] attributes
    system_dir: Path
        Path to *extracted* system files of the test set
    output_file: Path
        Path to save score json file
    overwrite: bool
        Skips existing if False
    """
    try:
        if not overwrite and output_file.exists():
            LOG.warning(f"get_scores: {output_file} exists")
            return
        reference_system = utils.ReferenceSystem.from_reference_system(
            system_dir, scorer_input.reference_system_id
        )
        receptor_file = scorer_input.receptor_file
        if np.isnan(receptor_file):
            # Rigid docking, use reference system's receptor cif file
            receptor_file = reference_system.receptor_cif_file
        scores = utils.ModelScores.from_files(
            scorer_input.id,
            Path(receptor_file),
            [
                Path(scorer_input.ligand_file)
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


def extract_test_systems(
    system_dir: Path, systems: list[str], output_dir: Path, overwrite: bool = False
) -> None:
    per_two_char = defaultdict(set)
    for system in systems:
        per_two_char[system[1:3]].add(system)
    LOG.info(
        f"Extracting {len(systems)} systems from {len(per_two_char)} two character codes"
    )
    for two_char in per_two_char:
        with ZipFile(system_dir / f"{two_char}.zip") as archive:
            for name in archive.namelist():
                if name.split("/")[0] in per_two_char[two_char]:
                    output_file = output_dir / name
                    if not overwrite and output_file.exists():
                        continue
                    archive.extract(name, output_dir)


def score_test_set(
    prediction_file: Path,
    system_dir: Path,
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
                    system_dir,
                    row.output_file,
                    overwrite,
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


def extract_and_score_test_set(
    prediction_file: Path,
    data_dir: Path,
    output_dir: Path,
    num_processes: int,
    overwrite: bool = False,
) -> None:
    predictions = pd.read_csv(prediction_file)
    test_systems = set(predictions["reference_system_id"])
    system_dir = output_dir / "test_systems"
    system_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "scores").mkdir(parents=True, exist_ok=True)
    if not overwrite:
        test_systems = test_systems - set(x.name for x in system_dir.iterdir())
    system_dir.mkdir(exist_ok=True)
    extract_test_systems(
        data_dir / "systems", list(test_systems), system_dir, overwrite=overwrite
    )
    score_test_set(
        prediction_file,
        system_dir,
        output_dir / "scores",
        output_dir / "scores.parquet",
        num_processes=num_processes,
        overwrite=overwrite,
    )


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Evaluate a set of docking predictions"
    )
    parser.add_argument(
        "--prediction_file",
        type=Path,
        help="Path to prediction file with [id, reference_system_id, receptor_file, ligand_file, rank, confidence] as columns",
    )
    parser.add_argument(
        "--data_dir",
        type=Path,
        help="Path to plinder data",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output folder where test systems and score JSON files are saved",
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

    args = parser.parse_args()

    extract_and_score_test_set(
        args.prediction_file,
        args.data_dir,
        args.output_dir,
        num_processes=args.num_processes,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
