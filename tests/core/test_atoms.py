# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import biotite.structure as struc
import numpy as np
import pytest
from biotite.structure.atoms import AtomArray
from plinder.core.index.system import PlinderSystem
from plinder.core.structure import vendored as atoms
from plinder.core.structure.atoms import (
    atom_array_from_cif_file,
    generate_input_conformer,
)
from plinder.core.structure.models import BackboneDefinition


def test_pdb_loader(pdb_5a7w_hydrogen):
    arr = atoms.atom_array_from_pdb_file(pdb_5a7w_hydrogen)

    assert isinstance(arr, AtomArray)
    assert arr.shape == (5623,)


def test_atom_masks(cif_atom_array):
    print(cif_atom_array)
    arr = cif_atom_array.copy()
    mask = atoms.backbone_mask(arr, BackboneDefinition("dockq"))
    assert mask.sum() == 1794

    assert mask.shape == (3256,)
    assert set(arr[mask].atom_name) == set(atoms.DOCKQ_BACKBONE_ATOMS)

    assert atoms.apply_mask(arr, mask).shape == arr[mask].shape

    assert set(atoms.filter_atoms(arr, calpha_only=True).atom_name) == {"CA"}


@pytest.mark.parametrize(
    "backbone_only, calpha_only, expected_mask",
    [
        (True, True, 1794),
        (True, False, 1794),
        (False, True, 360),
    ],
)
def test_get_backbone_atom_masks(
    backbone_only, calpha_only, expected_mask, cif_atom_array
):
    arr = cif_atom_array.copy()
    arr_mask, stack_mask = atoms.get_backbone_atom_masks(
        arr, struc.stack([arr]), backbone_only=backbone_only, calpha_only=calpha_only
    )
    assert arr_mask.shape[0] == 3256
    assert arr_mask.shape == stack_mask.shape
    assert arr_mask.sum() == stack_mask.sum() == expected_mask


def test_resn2seq(cif_atom_array):
    assert atoms.resn2seq(cif_atom_array.res_name[0:5]) == "TTTTT"
    structure, numbering, resn = atoms._get_structure_and_res_info(cif_atom_array)
    assert isinstance(structure, AtomArray)
    assert set(numbering) == set(cif_atom_array.res_id)
    assert atoms.resn2seq(resn[0:2]) == "TT"


def test_get_seq_alignments(read_plinder_mount):
    pdb = PlinderSystem(system_id="19hc__1__1.A_1.B__1.K_1.M_1.N").receptor_pdb
    a = atoms.atom_array_from_pdb_file(pdb)
    b = atoms.atom_array_from_pdb_file(pdb)
    a_numbering, a_resn = struc.get_residues(a)
    b_numbering, b_resn = struc.get_residues(b)
    a_seq = atoms.resn2seq(a_resn).strip("X")
    b_seq = atoms.resn2seq(b_resn).strip("X")
    ident = atoms.get_seq_identity(a_seq, b_seq)
    assert isinstance(ident, float)
    assert ident == pytest.approx(1.0)

    alns = atoms.get_seq_alignments(a_seq, b_seq)
    mismatches, matches = atoms.calc_num_mismatches(alns)
    assert mismatches == 0
    assert matches == len(a_seq)

    a_seq_aln, b_seq_aln, a_numbering, b_numbering = atoms.align_sequences(a_seq, b_seq)
    assert a_seq_aln == b_seq_aln == a_seq == b_seq
    assert a_numbering == list(range(1, len(a_seq) + 1))
    assert b_numbering == list(range(1, len(b_seq) + 1))


def test_buried_sasa(read_plinder_mount):
    pdb = PlinderSystem(system_id="19hc__1__1.A_1.B__1.K_1.M_1.N").receptor_pdb
    arr = atoms.atom_array_from_pdb_file(pdb)
    a = arr[arr.chain_id == "A"]
    b = arr[arr.chain_id == "B"]
    dsasa = atoms.get_buried_sasa(a, b)
    assert isinstance(dsasa, int)
    assert dsasa == 2520


def test_atom_array_from_cif_file(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)
    assert isinstance(arr, AtomArray)


def test_generate_input_conformer_easy():
    from rdkit import Chem

    # ligand_rdkit_canonical_smiles 4v2y__1__1.A__1.E
    thal_smiles = "O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1"
    mol = Chem.MolFromSmiles(thal_smiles)
    mol = generate_input_conformer(mol, addHs=True, minimize_maxIters=100)
    # check that it has a conformer
    assert mol.GetNumConformers() > 0
    # check that Z coords are not all zero (only for 2D mols)
    assert sum(abs(mol.GetConformer().GetPositions()[:, 2])) != 0


def test_generate_input_conformer_hard():
    from rdkit import Chem

    # ligand_rdkit_canonical_smiles from system_id "102m__1__1.A__1.C"
    hard_smiles = "C=CC1=C(C)C2=Cc3c(C)c(CCC(=O)O)c4n3[Fe]35<-N6=C(C=c7c(C=C)c(C)c(n73)=CC1=N->52)C(C)=C(CCC(=O)O)C6=C4"
    mol = Chem.MolFromSmiles(hard_smiles)
    mol = generate_input_conformer(mol)
    # check that it has a conformer
    assert mol.GetNumConformers() > 0
    # check that Z coords are not all zero (only for 2D mols)
    assert sum(abs(mol.GetConformer().GetPositions()[:, 2])) != 0


def test_structure_with_symmetry_in_ligand(read_plinder_mount):
    # structure with symmetry in ligand
    holo_struct = PlinderSystem(system_id="4v2y__1__1.A__1.E").holo_structure
    tag = holo_struct.ligand_chain_ordered[0]
    holo_struct.input_ligand_templates[tag]
    holo_struct.ligand_template2resolved_atom_order_stacks[tag]
    # will get two different atom stacks (two matches for the re-ordering)
    assert np.shape(holo_struct.ligand_template2resolved_atom_order_stacks[tag]) == (
        2,
        2,
        19,
    )


def test_structure_partially_resolved_ligand(read_plinder_mount):
    # structure with partially resolved ligand that can be matched piecewise
    holo_struct = PlinderSystem(system_id="1ngx__1__1.A_1.B__1.E").holo_structure
    tag = holo_struct.ligand_chain_ordered[0]
    holo_struct.input_ligand_templates[tag]
    holo_struct.ligand_template2resolved_atom_order_stacks[tag]
    # will get three different atom stacks (three partial matches for the re-ordering)
    assert np.shape(holo_struct.ligand_template2resolved_atom_order_stacks[tag]) == (
        2,
        3,
        28,
    )
