# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from biotite.structure import AtomArray
from rdkit import Chem

from plinder.data.structure.atoms import (
    atom_array_from_cif_file,
    atom_array_to_rdkit_mol,
)


def test_atom_array_from_cif_file(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)
    assert isinstance(arr, AtomArray)


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
