from pathlib import Path

import numpy as np
import pytest
from plinder.core import index
from plinder.core.loader.featurizer import system_featurizer


@pytest.mark.parametrize(
    "system_id",
    [
        "19hc__1__1.A_1.B__1.D_1.L_1.Q_1.S_1.U",
        "19hc__1__1.A_1.B__1.E_1.F_1.H_1.J_1.O",
        "19hc__1__1.A_1.B__1.K_1.M_1.N",
        "19hc__1__1.A_1.B__1.V_1.X_1.Y",
        "19hc__1__1.A__1.G",
        "19hc__1__1.A__1.I",
    ],
)
def test_plinder_system(system_id, read_plinder_mount):
    index.PlinderSystem(system_id=system_id)


@pytest.mark.parametrize(
    "system_id",
    [
        "19hc__1__1.A__1.C",
        "19hc__1__1.B__1.P",
    ],
)
def test_plinder_system_fails(system_id, read_plinder_mount):
    with pytest.raises(ValueError):
        index.PlinderSystem(system_id=system_id).system


def test_plinder_system_system_files(read_plinder_mount):
    system_id = "19hc__1__1.A_1.B__1.V_1.X_1.Y"
    s = index.PlinderSystem(system_id=system_id)
    assert len(s.structures) == 10
    assert len(s.ligand_sdfs) == 3
    assert len(s.system_cif)
    assert len(s.receptor_cif)
    assert len(s.receptor_pdb)
    assert len(s.sequences)
    assert s.chain_mapping is not None and len(s.chain_mapping)
    assert s.water_mapping is not None and len(s.water_mapping)
    assert Path(s.system_cif).is_file()
    assert Path(s.receptor_cif).is_file()
    assert Path(s.receptor_pdb).is_file()
    assert Path(s.sequences_fasta).is_file()
    assert isinstance(s.chain_mapping, dict)
    assert isinstance(s.water_mapping, dict)


def test_plinder_structure(read_plinder_mount):
    system_id = "1avd__1__1.A__1.C"
    s = index.PlinderSystem(system_id=system_id)
    holo_struc = s.holo_structure
    ligand_mols = holo_struc.ligand_mols
    # test the mask order for smiles
    assert np.all(
        ligand_mols["1.C"][3][0]
        == np.array([[13, 4, 5, 7, 9, 10, 1, 0, 3, 6, 8, 12, 11, 2]])
    )
    assert holo_struc.protein_sequence is not None
    assert len(holo_struc.protein_sequence)
    assert holo_struc.protein_atom_array is not None
    assert len(holo_struc.protein_atom_array)
    assert holo_struc.ligand_sdfs is not None
    assert len(holo_struc.ligand_sdfs)
    feat = system_featurizer(s)
    # Check the selected apo structure
    assert feat[1] == "1vyo_B"
    # Check seqres to holo mask
    seqres_holo_mask = feat[0]["holo_features"][
        "holo_input_sequence_residue_mask_stacked"
    ]
    assert np.allclose(np.array(seqres_holo_mask[:, :2]), np.array([[0, 0]]))
    # Check seqres to apo mask
    seqres_apo_mask = feat[0]["apo_features"][
        "apo_input_sequence_residue_mask_stacked"
    ]
    assert np.allclose(np.array(seqres_apo_mask[:, -3:]), np.array([[0, 0, 0]]))
    # Check apo to holo mask
    apo_holo_mask = feat[0]["apo_features"]["apo_input_sequence_residue_mask_stacked"]
    assert np.allclose(np.array(apo_holo_mask[:, -3:]), np.array([[0, 0, 0]]))

    assert np.allclose(np.array(
        feat[0]["apo_features"]["apo_protein_coordinates_stacked"].shape
    ), np.array([1, 958, 3]))
    assert np.allclose(np.array(
        feat[0]["apo_features"]["apo_protein_calpha_coordinates_stacked"].shape
    ), np.array([1, 122, 3]))


    assert np.allclose(np.array(
        feat[0]["holo_features"]["holo_protein_calpha_coordinates_stacked"].shape
    ), np.array([1, 123, 3]))

    assert np.allclose(np.array(
        feat[0]["holo_features"]["holo_protein_coordinates_stacked"].shape
    ), np.array([1, 964, 3]))

    assert np.allclose(np.array(
        feat[0]["sequence_features"]["input_sequence_residue_feat_stack"].shape
    ), np.array([1, 128, 21]))

    assert np.allclose(np.array(
        feat[0]["sequence_features"]["input_sequence_full_atom_feat_stack"].shape
    ), np.array([1, 1007, 12]))
