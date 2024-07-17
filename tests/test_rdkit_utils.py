# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

def test_valence_issue_handling():
    from rdkit import Chem
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
