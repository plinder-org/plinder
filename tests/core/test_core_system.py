from pathlib import Path

import numpy as np
import pytest
from plinder.core import index


@pytest.mark.parametrize(
    "system_id",
    [
        "1avd__1__1.A__1.C",
        "19hc__1__1.A__1.I",
        "19hc__1__1.B__1.T",
    ],
)
def test_plinder_system(system_id, read_plinder_mount):
    index.PlinderSystem(system_id=system_id).system


def test_plinder_system_system_files(read_plinder_mount):
    system_id = "1avd__1__1.A__1.C"
    s = index.PlinderSystem(system_id=system_id)
    assert len(s.ligand_sdfs) >= 1
    assert len(s.system_cif)
    assert len(s.receptor_cif)
    assert len(s.sequences)
    assert Path(s.system_cif).is_file()
    assert Path(s.receptor_cif).is_file()
    assert Path(s.sequences_fasta).is_file()


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
