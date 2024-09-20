# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from biotite.structure import AtomArray
from plinder.core.structure.atoms import (
    atom_array_from_cif_file,
    generate_input_conformer,
)


# TODO: likely not needed
def test_atom_array_from_cif_file(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)
    assert isinstance(arr, AtomArray)


def test_generate_input_conformer(cif_1qz5_unzipped):
    from rdkit import Chem

    # ligand_rdkit_canonical_smiles from system_id "102m__1__1.A__1.C"
    hard_smiles = "C=CC1=C(C)C2=Cc3c(C)c(CCC(=O)O)c4n3[Fe]35<-N6=C(C=c7c(C=C)c(C)c(n73)=CC1=N->52)C(C)=C(CCC(=O)O)C6=C4"
    mol = Chem.MolFromSmiles(hard_smiles)
    mol, _ = generate_input_conformer(mol)
    assert mol.GetNumConformers() > 0
