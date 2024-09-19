# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from rdkit import Chem


def test_valence_issue_handling():
    from plinder.data.utils.annotations.rdkit_utils import fix_valency_issues

    # AtomValenceException
    mol = Chem.MolFromSmiles("CC(=O)OCCN(C)(C)C", sanitize=False)
    mol = fix_valency_issues(mol)
    problems = Chem.DetectChemistryProblems(mol)
    assert len(problems) == 0

    # KekulizeException
    mol = Chem.MolFromSmiles("c1ccnc1", sanitize=False)
    mol = fix_valency_issues(mol)
    problems = Chem.DetectChemistryProblems(mol)
    assert len(problems) == 0


def test_matched_templates():
    from plinder.data.utils.annotations.rdkit_utils import (
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
