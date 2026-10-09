# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from plinder.data.annotations.ligand_utils import _check_stereo_vs_template


@pytest.mark.parametrize("different_chain", [True, False])
@pytest.mark.parametrize("mirror_first", [True, False])
def test_stereo_checks_each_connected_residue(different_chain, mirror_first):
    template = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)CO"))
    assert AllChem.EmbedMolecule(template, randomSeed=42) == 0
    template = Chem.RemoveHs(template)
    for atom in template.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{atom.GetIdx()}")
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)

    first, second = Chem.Mol(template), Chem.Mol(template)
    for atom in second.GetAtoms():
        if different_chain:
            atom.GetPDBResidueInfo().SetChainId("B")
        else:
            atom.GetPDBResidueInfo().SetInsertionCode("A")
    if mirror_first:
        conformer = first.GetConformer()
        positions = conformer.GetPositions()
        positions[:, 0] *= -1
        for index, position in enumerate(positions):
            conformer.SetAtomPosition(index, position)

    translation = (
        first.GetConformer().GetPositions()[4]
        + [1.45, 0, 0]
        - second.GetConformer().GetPositions()[4]
    )
    for index, position in enumerate(second.GetConformer().GetPositions()):
        second.GetConformer().SetAtomPosition(index, position + translation)

    combined = Chem.RWMol(Chem.CombineMols(first, second))
    # Connect the terminal hydroxyls without changing either stereocenter.
    combined.AddBond(4, first.GetNumAtoms() + 4, Chem.BondType.SINGLE)
    molecule = combined.GetMol()
    Chem.SanitizeMol(molecule)
    assert len(Chem.GetMolFrags(molecule)) == 1
    assert _check_stereo_vs_template(molecule, {"LIG": template}) is not mirror_first
