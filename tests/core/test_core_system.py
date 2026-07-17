from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
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


def test_plinder_system_receptor_type():
    system = index.PlinderSystem.__new__(index.PlinderSystem)
    system.system_id = "1abc__1__1.A_1.N__1.L"
    system._system = pd.DataFrame(
        {
            "system_receptor_type": ["protein+dna"],
            "system_protein_chains_asym_id": [["1.A", "1.N"]],
        }
    )
    system._entry_chains = pd.DataFrame(
        {
            "chain_asym_id": ["A", "N"],
            "chain_receptor_type": ["protein", "dna"],
        }
    )

    assert system.receptor_type == "protein+dna"
    assert system.receptor_chain_types == {"1.A": "protein", "1.N": "dna"}


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
    assert s.receptor_structure.array_length() > 0
    assert set(s.ligand_structures) == set(s.ligand_sdfs)
    assert all(
        structure.array_length() > 0 for structure in s.ligand_structures.values()
    )


def test_plinder_system_unpacks_canonical_ligand_archive(
    write_plinder_mount, monkeypatch
):
    from plinder.core.utils import config, cpl

    monkeypatch.setenv("PLINDER_OFFLINE", "true")
    config._config._clear()
    monkeypatch.setattr(cpl, "_CLIENTS", {})
    archive = write_plinder_mount / "ligand_archives" / "av.zip"
    archive.parent.mkdir(parents=True)
    with ZipFile(archive, "w", compression=ZIP_DEFLATED) as zip_file:
        zip_file.writestr("1avd/ligand_files/C.sdf", "canonical ASU ligand")

    system = index.PlinderSystem(system_id="1avd__1__1.A__1.C")

    assert system.canonical_ligand_sdfs == {
        "1.C": (
            write_plinder_mount / "ligand_archives" / "1avd/ligand_files/C.sdf"
        ).as_posix()
    }


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
