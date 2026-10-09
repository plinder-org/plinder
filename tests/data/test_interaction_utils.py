# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import time

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
import pytest
from biotite.structure import info
from rdkit import Chem

from plinder.data.annotations import cif_utils, interaction_utils


@pytest.mark.parametrize(
    ("receptor_name", "ligand_name", "contact_atom", "positive"),
    [
        ("ASP", "NH4", "OD1", False),
        ("ASP", "NH4", "OD2", False),
        ("LYS", "ACT", "NZ", True),
    ],
)
def test_salt_bridge_polarity(receptor_name, ligand_name, contact_atom, positive):
    receptor = info.residue(receptor_name)
    receptor = receptor[receptor.element != "H"]
    receptor.chain_id[:] = "R"
    receptor.res_id[:] = 10
    ligand = info.residue(ligand_name)
    ligand = ligand[ligand.element != "H"]
    ligand.chain_id[:] = "L"
    ligand.coord += (
        receptor.coord[receptor.atom_name == contact_atom][0]
        - ligand.coord[0]
        + [0, 0, 3]
    )
    # Binding-site indices must be translated back to the full receptor.
    distant = receptor.copy()
    distant.res_id[:] = 1
    distant.coord += 100
    receptor = distant + receptor

    interactions, _, failed = interaction_utils.run_peppr_interactions(
        receptor, ligand, struc.AtomArray(0), struc.AtomArray(0), "L"
    )

    assert "salt_bridges" not in failed
    bridges = [s for s in interactions["R"][10] if s.startswith("type:salt_bridges")]
    assert bridges == [f"type:salt_bridges__protispos:{positive}"]


@pytest.mark.parametrize("elements", [["NA", "MG"], ["MG", "NA"]])
def test_metal_bridge_indices_refer_to_original_array(elements):
    receptor = info.residue("ASP")
    receptor = receptor[receptor.element != "H"]
    receptor.chain_id[:] = "R"
    ligand = info.residue("ACT")
    ligand = ligand[ligand.element != "H"]
    ligand.chain_id[:] = "L"
    center = receptor.coord[receptor.atom_name == "OD2"][0] + [0, 0, 2]
    ligand.coord += center + [0, 0, 2] - ligand.coord[1]
    metals = struc.AtomArray(2)
    metals.res_name = np.array(elements)
    metals.element = np.array(elements)
    metals.coord[:] = [100, 100, 100]
    metals.coord[elements.index("MG")] = center

    bridges = interaction_utils.find_metal_bridges(receptor, ligand, metals)
    assert bridges
    assert all(
        metal_indices.tolist() == [elements.index("MG")]
        for _, _, metal_indices in bridges
    )

    interactions, _, failed = interaction_utils.run_peppr_interactions(
        receptor, ligand, struc.AtomArray(0), metals, "L"
    )
    assert "metal_complexes" not in failed
    labels = [
        label
        for labels in interactions["R"].values()
        for label in labels
        if label.startswith("type:metal_complexes")
    ]
    assert labels
    assert set(labels) == {"type:metal_complexes__metal_type:MG"}


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
