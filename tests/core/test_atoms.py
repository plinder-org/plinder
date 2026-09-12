# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from zipfile import ZipFile

import biotite.structure as struc
import numpy as np
import pytest
from biotite.structure.atoms import AtomArray, AtomArrayStack
from rdkit import Chem

from plinder.core.structure import atoms
from plinder.core.structure.smallmols_utils import generate_input_conformer
from plinder.core.structure.structure import Structure


def test_cif_loader(cif_1qz5_unzipped):
    arr = atoms.atom_array_from_cif_file(cif_1qz5_unzipped, use_author_fields=False)

    assert isinstance(arr, AtomArray)
    assert arr.shape == (3256,)
    assert sorted(set(arr.chain_id)) == ["A", "B", "C", "D", "E"]


def test_apply_mask(cif_atom_array):
    arr = cif_atom_array.copy()
    mask = arr.atom_name == "CA"

    assert mask.sum() == 360
    assert atoms.apply_mask(arr, mask).shape == (360,)
    stacked = atoms.apply_mask(struc.stack([arr, arr]), mask)
    assert isinstance(stacked, AtomArrayStack)
    assert stacked.shape == (2, 360)
    with pytest.raises(TypeError, match="AtomArray or AtomArrayStack"):
        atoms.apply_mask(mask, mask)


def test_resn2seq():
    # MSE keeps its methionine letter; SEC and PYL fall back to C and K because
    # biotite's protein alphabet lacks U and O; non-residues become X.
    names = ["ALA", "MSE", "SEC", "PYL", "HOH", "GLY", "UNK"]

    assert atoms.resn2seq(names) == "AMCKXGX"
    assert atoms.resn2seq([]) == ""


def test_align_sequences_maps_matched_residues_only():
    # Two-residue deletion in the subject: the reference numbering skips 4-5.
    assert atoms.align_sequences("ACDEFGHIK", "ACDGHIK") == (
        "ACDGHIK",
        "ACDGHIK",
        [1, 2, 3, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7],
    )
    # A substitution stays aligned and keeps both numberings.
    assert atoms.align_sequences("ACDEFGHIK", "ACDQFGHIK") == (
        "ACDEFGHIK",
        "ACDQFGHIK",
        list(range(1, 10)),
        list(range(1, 10)),
    )
    # Unknown residues are aligned but never reported as matched.
    assert atoms.align_sequences("ACDEF", "ACXEF") == (
        "ACEF",
        "ACEF",
        [1, 2, 4, 5],
        [1, 2, 4, 5],
    )
    # Caller-provided numbering is passed through.
    assert atoms.align_sequences("AC", "AC", [7, 8], [30, 31]) == (
        "AC",
        "AC",
        [7, 8],
        [30, 31],
    )


def test_residue_index_mapping_mask(cif_atom_array):
    chain_a = cif_atom_array[cif_atom_array.chain_id == "A"]
    resolved = atoms.resn2seq(struc.get_residues(chain_a)[1])
    assert len(resolved) == 359

    # A reference with one extra tryptophan: only that position is unresolved.
    reference = resolved[:10] + "W" + resolved[10:]
    mask = atoms.get_residue_index_mapping_mask({"A": reference}, cif_atom_array)["A"]
    assert mask.shape == (360,)
    assert np.flatnonzero(mask == 0).tolist() == [10]

    # A reference that starts five residues in is fully resolved.
    mask = atoms.get_residue_index_mapping_mask({"A": resolved[5:]}, cif_atom_array)
    assert mask["A"].shape == (354,)
    assert mask["A"].all()


def test_write_cif_roundtrip(cif_atom_array, tmp_path):
    chain_a = cif_atom_array[cif_atom_array.chain_id == "A"]
    output = tmp_path / "nested" / "chain_a.cif"

    atoms.write_cif(chain_a, output)
    reread = atoms.atom_array_from_cif_file(output, use_author_fields=False)

    assert output.read_text().startswith("data_chain_a")
    assert reread.shape == chain_a.shape
    assert reread.res_name.tolist() == chain_a.res_name.tolist()
    with pytest.raises(ValueError, match="must end in .cif"):
        atoms.write_cif(chain_a, tmp_path / "chain_a.pdb")


def test_remove_all_hs():
    # explicit bond stereo - from PDB: 5j1x
    mol = Chem.MolFromSmiles("[H]/N=C(/N)NCCC[C@H](NC(=O)OC(C)(C)C)C(=O)O")
    mol = Chem.RemoveAllHs(mol, sanitize=False)
    assert mol.GetNumAtoms() == mol.GetNumHeavyAtoms()
    # hydrogen isotopes - from PDB: 1tuj
    mol2 = Chem.MolFromSmiles("[2H]C([2H])(C(=O)[O-])C([2H])([2H])[Si](C)(C)C")
    mol2 = Chem.RemoveAllHs(mol2, sanitize=False)
    assert mol2.GetNumAtoms() == mol2.GetNumHeavyAtoms()
    # more strange explicit Hs
    mol3 = Chem.MolFromSmiles(
        "[H]/N=C(\\N)c1ccc(O)c(C=NCCN=Cc2cc(/C(N)=N\\[H])ccc2O)c1"
    )
    mol3 = Chem.RemoveAllHs(mol3, sanitize=False)
    assert mol3.GetNumAtoms() == mol3.GetNumHeavyAtoms()


def test_generate_input_conformer_easy():
    # ligand_rdkit_canonical_smiles 4v2y__1__1.A__1.E
    thal_smiles = "O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1"
    mol = Chem.MolFromSmiles(thal_smiles)
    mol = generate_input_conformer(mol, addHs=True, minimize_maxIters=100)
    # check that it has a conformer
    assert mol.GetNumConformers() > 0
    # check that Z coords are not all zero (only for 2D mols)
    assert sum(abs(mol.GetConformer().GetPositions()[:, 2])) != 0


def test_generate_input_conformer_hard():
    # ligand_rdkit_canonical_smiles from system_id "102m__1__1.A__1.C"
    hard_smiles = "C=CC1=C(C)C2=Cc3c(C)c(CCC(=O)O)c4n3[Fe]35<-N6=C(C=c7c(C=C)c(C)c(n73)=CC1=N->52)C(C)=C(CCC(=O)O)C6=C4"
    mol = Chem.MolFromSmiles(hard_smiles)
    mol = generate_input_conformer(mol)
    # check that it has a conformer
    assert mol.GetNumConformers() > 0
    # check that Z coords are not all zero (only for 2D mols)
    assert sum(abs(mol.GetConformer().GetPositions()[:, 2])) != 0


def structure_from_fixture_archive(
    read_plinder_mount, tmp_path, archive_code, system_id, ligand_smiles
):
    with ZipFile(read_plinder_mount / "systems" / f"{archive_code}.zip") as archive:
        receptor = archive.extract(f"{system_id}/receptor.cif", tmp_path)
        ligand = archive.extract(
            f"{system_id}/ligand_files/{system_id.split('__')[-1]}.sdf",
            tmp_path,
        )
    ligand_id = system_id.split("__")[-1]
    return Structure(
        id=system_id,
        protein_path=receptor,
        ligand_sdfs={ligand_id: ligand},
        ligand_smiles={ligand_id: ligand_smiles},
    )


def test_structure_with_symmetry_in_ligand(read_plinder_mount, tmp_path):
    # structure with symmetry in ligand
    holo_struct = structure_from_fixture_archive(
        read_plinder_mount,
        tmp_path,
        "v2",
        "4v2y__1__1.A__1.E",
        "O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1",
    )
    tag = holo_struct.ligand_chain_ordered[0]
    holo_struct.input_ligand_templates[tag]
    holo_struct.ligand_template2resolved_atom_order_stacks[tag]
    # will get two different atom stacks (two matches for the re-ordering)
    assert np.shape(holo_struct.ligand_template2resolved_atom_order_stacks[tag]) == (
        2,
        2,
        19,
    )


def test_structure_partially_resolved_ligand(read_plinder_mount, tmp_path):
    # structure with partially resolved ligand that can be matched piecewise
    holo_struct = structure_from_fixture_archive(
        read_plinder_mount,
        tmp_path,
        "ng",
        "1ngx__1__1.A_1.B__1.E",
        "COCCO[C@@H](C)CO[C@H](C)CO[C@H](C)COC(C)CO[C@@H](C)CO[C@@H](C)"
        "CO[C@H](C)CO[C@H](C)COC[C@H](C)N",
    )
    tag = holo_struct.ligand_chain_ordered[0]
    holo_struct.input_ligand_templates[tag]
    holo_struct.ligand_template2resolved_atom_order_stacks[tag]
    # will get three different atom stacks (three partial matches for the re-ordering)
    assert np.shape(holo_struct.ligand_template2resolved_atom_order_stacks[tag]) == (
        2,
        3,
        28,
    )
