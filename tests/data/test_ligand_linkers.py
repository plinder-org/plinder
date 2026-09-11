# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import pytest
from plinder.data.annotations import ligand_utils
from rdkit import Chem


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [("CC", 0), ("CCC", 1), ("CCCC", 2), ("C1CCCCC1", 0), ("CC(C)C", 0), ("CCOCC", 1)],
)
def test_linker_length_excludes_chain_ends_branches_and_rings(smiles, expected):
    molecule = Chem.MolFromSmiles(smiles)
    assert (
        ligand_utils.get_len_of_longest_linear_hydrocarbon_linker(molecule) == expected
    )


@pytest.mark.parametrize(
    ("carbons", "limit", "expected"), [(6, 3, 3), (5, 3, 3), (4, 3, 2), (55, 50, 50)]
)
def test_linker_search_limit_is_a_lower_bound(carbons, limit, expected):
    molecule = Chem.MolFromSmiles("C" * carbons)
    assert (
        ligand_utils.get_len_of_longest_linear_hydrocarbon_linker(
            molecule, max_count=limit
        )
        == expected
    )


@pytest.mark.parametrize("limit", [0, -1])
def test_linker_search_requires_a_positive_limit(limit):
    with pytest.raises(ValueError, match="max_count must be positive"):
        ligand_utils.get_len_of_longest_linear_hydrocarbon_linker(
            Chem.MolFromSmiles("CCC"), max_count=limit
        )


@pytest.mark.parametrize("name", [None, "example ligand"])
def test_linker_failure_is_unknown_and_safe_for_unnamed_molecules(
    monkeypatch, caplog, name
):
    molecule = Chem.MolFromSmiles("CCCCCC")
    if name is not None:
        molecule.SetProp("_Name", name)

    def fail(*_args, **_kwargs):
        raise RuntimeError("ring calculation failed")

    monkeypatch.setattr(Chem, "SanitizeMol", fail)
    assert ligand_utils.get_len_of_longest_linear_hydrocarbon_linker(molecule) is None
    assert "ring calculation failed" in caplog.text
    assert (name or "unnamed molecule") in caplog.text
    assert ligand_utils.is_excluded_mol("CCCCCC") is True


@pytest.mark.parametrize(
    ("carbons", "cutoff", "excluded"),
    [(14, 12, False), (15, 12, True), (62, 60, False), (63, 60, True)],
)
def test_artifact_linker_rule_honors_the_configured_cutoff(carbons, cutoff, excluded):
    assert (
        ligand_utils.is_excluded_mol(
            "C" * carbons, max_linear_hydrocarbon_linker=cutoff
        )
        is excluded
    )
