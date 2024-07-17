# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from biotite.structure import AtomArray
from rdkit import Chem

from plinder.data.structure.atoms import (
    add_hydrogen_to_lig_array,
    atom_array_from_cif_file,
    atom_array_to_rdkit_mol,
    atom_vdw_radius,
    backbone_mask,
    filter_atoms,
    protonate_protein,
)


def test_atom_array_from_cif_file(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)
    assert isinstance(arr, AtomArray)
    assert atom_vdw_radius(arr[0]) == 1.64
    assert sum(backbone_mask(arr, "dockq")) > 0
    assert len(filter_atoms(arr)) > 0


def test_atom_array_to_rdkit_mol(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)

    lig_arr = arr[(arr.res_name == "ATP") & (arr.res_id == 380) & (arr.chain_id == "A")]
    assert isinstance(
        atom_array_to_rdkit_mol(
            lig_arr,
            smiles="C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
            chain_mapping={"A1": "A"},
            add_h=False,
        ),
        Chem.rdchem.Mol,
    )


def test_protonate_protein(protein_pdb_block):
    prot = protonate_protein(protein_pdb_block)
    print(prot[0])
    assert prot[0].split("\n")[-3].split()[-2] == "H"


def test_add_hydrogen_to_lig_array(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)

    lig_arr = arr[(arr.res_name == "ATP") & (arr.res_id == 380) & (arr.chain_id == "A")]
    lig_h = Chem.MolToPDBBlock(
        add_hydrogen_to_lig_array(
            lig_arr,
            smiles="C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
            chain_mapping={"A1": "A"},
            add_h=True,
        )
    )
    lig_h = lig_h.split("\n")[-32].split()[-1]
    assert lig_h == "H"
