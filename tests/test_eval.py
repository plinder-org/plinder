# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from plinder.core.index.system import PlinderSystem
from plinder.eval.docking import utils
from plinder.eval.docking.make_plots import plot_cmd
from plinder.eval.docking.stratify_test_set import stratify_cmd
from plinder.eval.docking.write_scores import evaluate_cmd


def mock_path_eval(
    *, rel: str = "", download: bool = False, force_progress: bool = False
):
    obj = Path(
        "/".join(
            [
                str(os.getenv("PLINDER_MOUNT")),
                str(os.getenv("PLINDER_BUCKET")),
                str(os.getenv("PLINDER_RELEASE")),
            ]
        )
    )
    return obj / rel if rel else obj


def mock_tanimoto(x, y):
    return np.random.rand(max(len(x), len(y)))


@pytest.fixture
def mock_cpl_eval(read_plinder_eval_mount, monkeypatch):
    # patch cpl at core.utils not core.index.utils because of unpack
    monkeypatch.setattr(
        "plinder.core.utils.cpl.get_plinder_path",
        mock_path_eval,
    )
    monkeypatch.setattr(
        "plinder.core.utils.cpl.download_paths",
        lambda **kws: None,
    )
    monkeypatch.setattr(
        "plinder.eval.docking.stratify_test_set.smallmolecules.tanimoto_maxsim_matrix",
        mock_tanimoto,
    )


def test_single_protein_single_ligand_scoring(
    system_1a3b, predicted_pose_1a3b, mock_cpl_eval
):
    reference_system = PlinderSystem(system_id=system_1a3b)
    scores = utils.ModelScores.from_model_files(
        predicted_pose_1a3b.parent.name,
        Path(reference_system.receptor_cif),
        [predicted_pose_1a3b],
        reference_system,
        score_protein=True,
        score_posebusters=True,
    ).summarize_scores()
    true_scores = {
        "model": "1a3b__1__1.B__1.D",
        "reference": "1a3b__1__1.B__1.D",
        "num_reference_ligands": 1,
        "num_model_ligands": 1,
        "num_reference_proteins": 1,
        "num_model_proteins": 1,
        "fraction_reference_ligands_mapped": 1.0,
        "fraction_model_ligands_mapped": 1.0,
        "lddt_pli_ave": 0.8581504702194357,
        "lddt_pli_wave": 0.8581504702194357,
        "bisy_rmsd_ave": 1.6171839155715722,
        "bisy_rmsd_wave": 1.6171839155715722,
        "lddt_lp_ave": 1.0,
        "lddt_lp_wave": 1.0,
        "posebusters_mol_pred_loaded": True,
        "posebusters_mol_cond_loaded": True,
        "posebusters_sanitization": True,
        "posebusters_all_atoms_connected": True,
        "posebusters_bond_lengths": True,
        "posebusters_bond_angles": True,
        "posebusters_internal_steric_clash": True,
        "posebusters_aromatic_ring_flatness": True,
        "posebusters_double_bond_flatness": True,
        "posebusters_internal_energy": True,
        "posebusters_protein-ligand_maximum_distance": True,
        "posebusters_minimum_distance_to_protein": False,
        "posebusters_minimum_distance_to_organic_cofactors": True,
        "posebusters_minimum_distance_to_inorganic_cofactors": True,
        "posebusters_minimum_distance_to_waters": True,
        "posebusters_volume_overlap_with_protein": True,
        "posebusters_volume_overlap_with_organic_cofactors": True,
        "posebusters_volume_overlap_with_inorganic_cofactors": True,
        "posebusters_volume_overlap_with_waters": True,
        "fraction_reference_proteins_mapped": 1.0,
        "fraction_model_proteins_mapped": 1.0,
        "lddt": 1.0,
        "bb_lddt": 1.0,
        "per_chain_lddt_ave": 0.9960159362549801,
        "per_chain_bb_lddt_ave": 1.0,
    }
    for k in true_scores:
        assert k in scores
        if type(true_scores[k]) == float:
            assert np.isclose(
                scores[k], true_scores[k]
            ), f"{k}: {scores[k]} != {true_scores[k]}"
        else:
            assert scores[k] == true_scores[k], f"{k}: {scores[k]} != {true_scores[k]}"


def test_multi_protein_single_ligand_scoring(
    system_1ai5, predicted_pose_1ai5, mock_cpl_eval
):
    reference_system = PlinderSystem(system_id=system_1ai5)
    scores = utils.ModelScores.from_model_files(
        predicted_pose_1ai5.parent.name,
        Path(reference_system.receptor_cif),
        [predicted_pose_1ai5],
        reference_system,
        score_protein=True,
        score_posebusters=True,
    ).summarize_scores()
    true_scores = {
        "model": "1ai5__1__1.A_1.B__1.D",
        "reference": "1ai5__1__1.A_1.B__1.D",
        "num_reference_ligands": 1,
        "num_model_ligands": 1,
        "num_reference_proteins": 2,
        "num_model_proteins": 2,
        "fraction_reference_ligands_mapped": 1.0,
        "fraction_model_ligands_mapped": 1.0,
        "lddt_pli_ave": 0.5106951871657754,
        "lddt_pli_wave": 0.5106951871657754,
        "bisy_rmsd_ave": 3.6651428915654645,
        "bisy_rmsd_wave": 3.6651428915654645,
        "lddt_lp_ave": 1.0,
        "lddt_lp_wave": 1.0,
        "posebusters_mol_pred_loaded": True,
        "posebusters_mol_cond_loaded": True,
        "posebusters_sanitization": True,
        "posebusters_all_atoms_connected": True,
        "posebusters_bond_lengths": True,
        "posebusters_bond_angles": True,
        "posebusters_internal_steric_clash": True,
        "posebusters_aromatic_ring_flatness": True,
        "posebusters_double_bond_flatness": True,
        "posebusters_internal_energy": True,
        "posebusters_protein-ligand_maximum_distance": True,
        "posebusters_minimum_distance_to_protein": False,
        "posebusters_minimum_distance_to_organic_cofactors": True,
        "posebusters_minimum_distance_to_inorganic_cofactors": True,
        "posebusters_minimum_distance_to_waters": True,
        "posebusters_volume_overlap_with_protein": False,
        "posebusters_volume_overlap_with_organic_cofactors": True,
        "posebusters_volume_overlap_with_inorganic_cofactors": True,
        "posebusters_volume_overlap_with_waters": True,
        "fraction_reference_proteins_mapped": 1.0,
        "fraction_model_proteins_mapped": 1.0,
        "lddt": 1.0,
        "bb_lddt": 1.0,
        "per_chain_lddt_ave": 0.9991023339317774,
        "per_chain_bb_lddt_ave": 1.0,
        "qs_global": 1.0,
        "qs_best": 1.0,
        "dockq_wave": 1.0,
        "dockq_ave": 1.0,
    }
    for k in true_scores:
        assert k in scores
        if type(true_scores[k]) == float:
            assert np.isclose(
                scores[k], true_scores[k]
            ), f"{k}: {scores[k]} != {true_scores[k]}"
        else:
            assert scores[k] == true_scores[k], f"{k}: {scores[k]} != {true_scores[k]}"


@pytest.fixture
def prediction_csv(read_plinder_eval_mount, mock_cpl_eval, tmp_path):
    csv = f"""\
id,reference_system_id,receptor_file,rank,confidence,ligand_file
1ai5__1__1.A_1.B__1.D,1ai5__1__1.A_1.B__1.D,,1,1.0,{read_plinder_eval_mount}/predicted_poses/1ai5__1__1.A_1.B__1.D/rank1.sdf
1a3b__1__1.B__1.D,1a3b__1__1.B__1.D,,1,1.0,{read_plinder_eval_mount}/predicted_poses/1a3b__1__1.B__1.D/rank1.sdf
"""
    fn = tmp_path / "prediction.csv"
    with fn.open("w") as f:
        f.write(csv)
    return fn


def test_evaluate_stratify_plot_cmds(prediction_csv, mock_cpl_eval):
    from plinder.core.utils import config

    config._config._clear()
    cfg = config.get_config()

    args = [
        "--prediction_file",
        f"{prediction_csv}",
        "--output_dir",
        f"{prediction_csv.parent}",
        "--num_processes",
        "8",
    ]
    evaluate_cmd(args=args)

    score_df = pd.read_parquet(f"{prediction_csv.parent}/scores.parquet")
    assert np.allclose(
        score_df.sort_values(by="reference").bisy_rmsd_wave.to_list(),
        [1.617184, 3.665143],
    )

    args = [
        "--split_file",
        f"{Path(cfg.data.plinder_dir) / 'test_split.parquet'}",
        "--output_dir",
        f"{prediction_csv.parent}",
    ]
    stratify_cmd(args=args)

    stratify_df = pd.read_parquet(f"{prediction_csv.parent}/test_set.parquet")
    assert np.allclose(stratify_df.novel_all.to_list(), [True, True])

    args = [
        "--score_file",
        f"{prediction_csv.parent}/scores.parquet",
        "--data_file",
        f"{prediction_csv.parent}/test_set.parquet",
        "--output_dir",
        f"{prediction_csv.parent}/plots",
        "--max_top_n",
        "1",
    ]
    plot_cmd(args=args)
    result_df = pd.read_csv(Path(prediction_csv.parent) / "plots" / "results.csv")
    truth = pd.read_csv(Path(cfg.data.plinder_dir) / "results.csv")
    assert result_df.equals(truth)
    assert (Path(prediction_csv.parent) / "plots" / "merged.parquet").exists()
    assert (
        Path(prediction_csv.parent) / "plots" / "delta_lDDT_PLI_topn1.html"
    ).exists()
