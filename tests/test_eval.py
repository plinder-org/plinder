# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from plinder.eval.docking import utils
import numpy as np


def test_single_protein_single_ligand_scoring(system_1a3b, predicted_pose_1a3b):
    reference_system = utils.ReferenceSystem.from_reference_system(
        system_1a3b.parent, system_1a3b.name
    )
    scores = utils.ModelScores.from_files(
        predicted_pose_1a3b.parent.name,
        reference_system.receptor_cif_file,
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
        "lddt_pli_ave": 0.8895061728395062,
        "lddt_pli_wave": 0.8895061728395062,
        "lddt_pli_amd_ave": 0.8581504702194357,
        "lddt_pli_amd_wave": 0.8581504702194357,
        "scrmsd_ave": 1.6171839155715722,
        "scrmsd_wave": 1.6171839155715722,
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
            assert np.isclose(scores[k], true_scores[k])
        else:
            assert scores[k] == true_scores[k]


def test_multi_protein_single_ligand_scoring(system_1ai5, predicted_pose_1ai5):
    reference_system = utils.ReferenceSystem.from_reference_system(
        system_1ai5.parent, system_1ai5.name
    )
    scores = utils.ModelScores.from_files(
        predicted_pose_1ai5.parent.name,
        reference_system.receptor_cif_file,
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
        "lddt_pli_ave": 0.557840616966581,
        "lddt_pli_wave": 0.557840616966581,
        "lddt_pli_amd_ave": 0.5106951871657754,
        "lddt_pli_amd_wave": 0.5106951871657754,
        "scrmsd_ave": 3.6651428915654645,
        "scrmsd_wave": 3.6651428915654645,
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
            assert np.isclose(scores[k], true_scores[k])
        else:
            assert scores[k] == true_scores[k]
