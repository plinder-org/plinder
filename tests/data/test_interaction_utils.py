# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import time

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from plinder.data.utils.annotations import cif_utils, interaction_utils
from rdkit import Chem


def test_peppr_tautomer_cache_requires_matching_atom_order(monkeypatch) -> None:
    molecule = Chem.MolFromSmiles("c1ncc[nH]1")
    assert molecule is not None
    calls = []

    def enumerate_tautomers(input_molecule: Chem.Mol) -> list[Chem.Mol]:
        calls.append(input_molecule)
        return [Chem.Mol(input_molecule)]

    monkeypatch.setattr(
        interaction_utils,
        "_ORIGINAL_GET_INTERCHANGEABLE_TAUTOMERS",
        enumerate_tautomers,
    )
    interaction_utils._TAUTOMER_CACHE.clear()
    try:
        first = interaction_utils._cached_get_interchangeable_tautomers(molecule)
        second = interaction_utils._cached_get_interchangeable_tautomers(
            Chem.Mol(molecule)
        )
        reordered = interaction_utils._cached_get_interchangeable_tautomers(
            Chem.RenumberAtoms(
                molecule,
                list(reversed(range(molecule.GetNumAtoms()))),
            )
        )
    finally:
        interaction_utils._TAUTOMER_CACHE.clear()

    assert len(calls) == 2
    assert Chem.MolToSmiles(first[0]) == Chem.MolToSmiles(second[0])
    assert Chem.MolToSmiles(first[0]) == Chem.MolToSmiles(reordered[0])
    assert first[0] is not second[0]


def _large_charged_molecule() -> Chem.Mol:
    molecule = Chem.MolFromSmiles("[O-]C(=O)" + "C" * 70)
    assert molecule is not None
    return molecule


def test_resonance_guard_preserves_fast_exact_result(monkeypatch) -> None:
    molecule = _large_charged_molecule()
    expected = (
        np.arange(molecule.GetNumAtoms()) % 2 == 0,
        np.arange(molecule.GetNumAtoms()) % 3 == 0,
        np.arange(molecule.GetNumAtoms()),
    )
    monkeypatch.setattr(
        interaction_utils,
        "_ORIGINAL_FIND_RESONANCE_CHARGES",
        lambda _molecule: expected,
    )
    interaction_utils._RESONANCE_CACHE.clear()
    try:
        actual = interaction_utils._bounded_find_resonance_charges(molecule)
    finally:
        interaction_utils._RESONANCE_CACHE.clear()

    assert all(np.array_equal(left, right) for left, right in zip(actual, expected))


def test_resonance_guard_falls_back_only_after_timeout(monkeypatch) -> None:
    molecule = _large_charged_molecule()

    def slow_exact_result(_molecule: Chem.Mol):
        time.sleep(1)
        raise AssertionError("worker should be terminated first")

    monkeypatch.setattr(
        interaction_utils,
        "_ORIGINAL_FIND_RESONANCE_CHARGES",
        slow_exact_result,
    )
    monkeypatch.setattr(interaction_utils, "_RESONANCE_TIMEOUT_SECONDS", 0.01)
    interaction_utils._RESONANCE_CACHE.clear()
    try:
        positive, negative, groups = interaction_utils._bounded_find_resonance_charges(
            molecule
        )
    finally:
        interaction_utils._RESONANCE_CACHE.clear()

    formal_charges = np.array([atom.GetFormalCharge() for atom in molecule.GetAtoms()])
    assert np.array_equal(positive, formal_charges > 0)
    assert np.array_equal(negative, formal_charges < 0)
    assert np.array_equal(groups, np.arange(molecule.GetNumAtoms()))


def test_symmetry_contacts_skip_placeholder_unit_cell(monkeypatch) -> None:
    asu = struc.AtomArray(1)
    asu.coord[:] = 0
    asu.element[:] = "C"
    asu.res_name[:] = "ALA"
    unit_cell = asu + asu
    unit_cell.box = np.eye(3)
    monkeypatch.setattr(
        cif_utils,
        "get_structure_with_altloc",
        lambda *_args, **_kwargs: asu,
    )
    monkeypatch.setattr(
        cif_utils,
        "get_unit_cell_with_altloc",
        lambda *_args, **_kwargs: unit_cell,
    )

    contacts = interaction_utils.get_symmetry_mate_contacts(pdbx.CIFFile())

    assert contacts == {}
