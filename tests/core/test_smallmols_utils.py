# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest
from rdkit import Chem


@pytest.mark.parametrize(
    ["smiles", "num_problems"],
    [
        ["CC(=O)OCCN(C)(C)C", 0],  # AtomValenceException
        ["c1ccnc1", 0],  # KekulizeException
    ],
)
def test_valence_issue_handling(smiles, num_problems):
    from plinder.core.structure.smallmols_utils import fix_valency_issues

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = fix_valency_issues(mol)
    problems = Chem.DetectChemistryProblems(mol)
    assert len(problems) == num_problems


@pytest.mark.parametrize(
    ["smiles", "num_charged_atoms"],
    [
        ["CC(=O)OCCN(C)(C)C", 0],
        ["[O-]C(=[O])CC[NH+](C)(C)", 0],
        ["[O-]C(=[OH+])CC[NH+](C)(C)", 0],
        ["OC(=[O])CC[N+](C)(C)(C)", 1],
        ["[O-]C(=[O])CC[N+](C)(C)(C)", 2],
        ["[O-]C(=[O])C.C[N+](C)(C)(C)", 2],
    ],
)
def test_uncharge_mol(smiles, num_charged_atoms):
    from plinder.core.structure.smallmols_utils import uncharge_mol

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = uncharge_mol(mol)
    assert (
        sum([at.GetFormalCharge() != 0 for at in mol.GetAtoms()]) == num_charged_atoms
    )


@pytest.mark.parametrize(
    ["smiles", "inchikey", "remove_stereo"],
    [
        ["CC/C=C/Cl", "DUDKKPVINWLFBI-ONEGZZNKSA-N", False],
        ["CC/C=C\\Cl", "DUDKKPVINWLFBI-ARJAWSKDSA-N", False],
        ["CC/C=C/Cl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", True],
        ["CC/C=C\\Cl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", True],
        ["CCC=CCl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", False],
        ["CCC=CCl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", True],
        ["C[C@@](F)(Cl)CBr", "REKDFINPOZVXJS-VKHMYHEASA-N", False],
        ["C[C@](F)(Cl)CBr", "REKDFINPOZVXJS-GSVOUGTGSA-N", False],
        ["C[C@](F)(Cl)CBr", "REKDFINPOZVXJS-UHFFFAOYSA-N", True],
        ["C[C@@](F)(Cl)CBr", "REKDFINPOZVXJS-UHFFFAOYSA-N", True],
        ["CC(F)(Cl)CBr", "REKDFINPOZVXJS-UHFFFAOYSA-N", True],
    ],
)
def test_inchikey(smiles, inchikey, remove_stereo):
    from plinder.core.structure.smallmols_similarity import smiles2inchikey

    assert inchikey == smiles2inchikey(smiles, remove_stereo=remove_stereo)


def test_matched_templates():
    from plinder.core.structure.smallmols_utils import (
        get_matched_template_v2,
        mol_assigned_bond_orders_by_template,
    )

    mol1 = Chem.MolFromSmiles("FC(Cl)(Br)C.CNCC1CCCCC1.CCC(OC)O")
    template = Chem.MolFromSmiles("F[C@@](Br)(Cl)CCCNCc1cc(C(=O)N/C=C/C(OC)=O)ccc1")
    matched_template = get_matched_template_v2(template, mol1)
    fixed_mol = mol_assigned_bond_orders_by_template(matched_template, mol1)
    fixed_mol_SMILES = Chem.CanonSmiles(Chem.MolToSmiles(fixed_mol))
    assert fixed_mol_SMILES.count("=") >= 2
    assert fixed_mol_SMILES == "C=CC(=O)OC.CC(F)(Cl)Br.CNCc1ccccc1"
