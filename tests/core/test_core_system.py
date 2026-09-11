from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from plinder.core import index
from plinder.core.release import PlinderRelease


@pytest.mark.usefixtures("read_plinder_mount")
@pytest.mark.parametrize(
    "system_id",
    [
        "1avd__1__1.A__1.C",
        "19hc__1__1.A__1.I",
        "19hc__1__1.B__1.T",
    ],
)
def test_plinder_system(system_id):
    index.PlinderSystem(system_id=system_id).system


def test_plinder_system_uses_explicit_release(cif_2y4i, tmp_path):
    release_root = tmp_path / "release"
    index_dir = release_root / "index"
    index_dir.mkdir(parents=True)
    system_id = "2y4i__1__1.A__1.B"
    pd.DataFrame({
        "entry_pdb_id": ["2y4i"],
        "system_id": [system_id],
        "system_biounit_id": ["1"],
        "system_receptor_type": ["protein"],
        "system_protein_chains_asym_id": [["1.A"]],
        "ligand_instance": [1],
        "ligand_asym_id": ["B"],
    }).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame({
        "entry_pdb_id": ["2y4i"],
        "entry_release_date": ["2000-01-01"],
    }).to_parquet(index_dir / "entry_metadata.parquet", index=False)
    pd.DataFrame({
        "entry_pdb_id": ["2y4i"],
        "chain_asym_id": ["A"],
        "chain_receptor_type": ["protein"],
    }).to_parquet(index_dir / "entry_chains.parquet", index=False)
    pd.DataFrame({
        "entry_pdb_id": ["2y4i", "2y4i"],
        "biounit_id": ["1", "1"],
        "chain_instance": ["1.A", "1.B"],
        "chain_asym_id": ["A", "B"],
        "chain_role": ["receptor", "ligand"],
    }).to_parquet(index_dir / "entry_biounit_chains.parquet", index=False)

    release = PlinderRelease(data_dir=release_root)
    system = index.PlinderSystem(
        system_id=system_id,
        release=release,
        source_mmcif=cif_2y4i,
    )

    assert system.reconstruction_dir == (
        release_root / "reconstructed_systems" / system_id
    )
    assert system.system["entry_release_date"].tolist() == ["2000-01-01"]
    assert system.entry["system_id"].tolist() == [system_id]
    assert system.receptor_chain_types == {"1.A": "protein"}
    assert system.biounit_chains["chain_instance"].tolist() == ["1.A", "1.B"]


def test_plinder_system_receptor_type():
    system = index.PlinderSystem.__new__(index.PlinderSystem)
    system.system_id = "1abc__1__1.A_1.N__1.L"
    system._system = pd.DataFrame({
        "system_receptor_type": ["protein+dna"],
        "system_protein_chains_asym_id": [["1.A", "1.N"]],
    })
    system._entry_chains = pd.DataFrame({
        "chain_asym_id": ["A", "N"],
        "chain_receptor_type": ["protein", "dna"],
    })

    assert system.receptor_type == "protein+dna"
    assert system.receptor_chain_types == {"1.A": "protein", "1.N": "dna"}


def test_plinder_system_smiles_uses_published_ligand_smiles():
    system = index.PlinderSystem.__new__(index.PlinderSystem)
    system._system = pd.DataFrame({
        "ligand_instance": [1],
        "ligand_asym_id": ["L"],
        "ligand_smiles": ["CCO"],
        "ligand_resolved_smiles": ["CC"],
    })

    assert system.smiles == {"1.L": "CCO"}


def test_plinder_system_system_files(cached_plinder_system):
    system_id = "1avd__1__1.A__1.C"
    s = cached_plinder_system(system_id)
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


def test_plinder_system_extracts_canonical_ligand_archive(
    write_plinder_mount, monkeypatch
):
    from plinder.core.utils import config, cpl

    monkeypatch.setenv("PLINDER_OFFLINE", "true")
    config._config._clear()
    monkeypatch.setattr(cpl, "_CLIENTS", {})
    archive = write_plinder_mount / "ligand_archives" / "av.parquet"
    archive.parent.mkdir(parents=True)
    pd.DataFrame({
        "pdb_id": ["1avd"],
        "ligand_asym_id": ["C"],
        "sdf": [b"canonical ASU ligand"],
    }).to_parquet(archive, index=False)

    system = index.PlinderSystem(system_id="1avd__1__1.A__1.C")

    assert system.canonical_ligand_sdfs == {
        "1.C": (
            write_plinder_mount / "ligand_archives" / "1avd/ligand_files/C.sdf"
        ).as_posix()
    }


def test_plinder_structure(cached_plinder_system):
    system_id = "1avd__1__1.A__1.C"
    s = cached_plinder_system(system_id)
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
