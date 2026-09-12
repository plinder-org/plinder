# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import biotite.structure as struc
import numpy as np

from plinder.core.structure.atoms import atom_array_from_cif_file
from plinder.core.structure.ccd_template import (
    UNRESOLVED_ANNOTATION,
    add_missing_atoms,
    ccd_component_template,
)
from plinder.core.structure.structure import Structure

CIF_6M92 = "xx/pdb_00006m92/pdb_00006m92_xyz-enrich.cif.gz"


def _lysine_30(atoms):
    return atoms[(atoms.chain_id == "A") & (atoms.res_id == 30)]


def test_add_missing_atoms_completes_side_chains_with_nan_and_mask(test_dir):
    atoms = atom_array_from_cif_file(test_dir / CIF_6M92, use_author_fields=False)
    atoms = atoms[atoms.res_name != "HOH"]
    before = _lysine_30(atoms)
    assert before.atom_name.tolist() == ["N", "CA", "C", "O", "CB"]

    completed = add_missing_atoms(atoms)

    after = _lysine_30(completed)
    assert after.atom_name.tolist() == [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "CE",
        "NZ",
    ]
    assert after.element.tolist()[5:] == ["C", "C", "C", "N"]
    assert np.isnan(after.coord[5:]).all() and not np.isnan(after.coord[:5]).any()
    assert (
        after.get_annotation(UNRESOLVED_ANNOTATION).tolist() == [False] * 5 + [True] * 4
    )
    assert after.res_name.tolist() == ["LYS"] * 9 and set(after.chain_id) == {"A"}

    # deposited atoms are untouched, in order, and residues stay contiguous
    mask = completed.get_annotation(UNRESOLVED_ANNOTATION)
    kept = completed[~mask]
    assert kept.atom_name.tolist() == atoms.atom_name.tolist()
    np.testing.assert_array_equal(kept.coord, atoms.coord)
    assert struc.get_residue_count(completed) == struc.get_residue_count(atoms)
    # leaving atoms (OXT) and hydrogens are never added
    added = completed[mask]
    assert (
        "OXT" not in set(added.atom_name)
        and not np.isin(added.element, ["H", "D"]).any()
    )
    # the input is not modified
    assert UNRESOLVED_ANNOTATION not in atoms.get_annotation_categories()


def test_add_missing_atoms_extends_bonds_and_remaps_existing_ones(test_dir):
    atoms = atom_array_from_cif_file(test_dir / CIF_6M92, use_author_fields=False)
    atoms = atoms[atoms.res_name != "HOH"]
    assert atoms.bonds is not None

    completed = add_missing_atoms(atoms)

    assert completed.bonds is not None
    mask = completed.get_annotation(UNRESOLVED_ANNOTATION)
    # every deposited bond survives, as a bond between the same two atoms
    old_pairs = {
        (
            atoms.chain_id[i],
            atoms.res_id[i],
            atoms.atom_name[i],
            atoms.chain_id[j],
            atoms.res_id[j],
            atoms.atom_name[j],
        )
        for i, j, _ in atoms.bonds.as_array()
    }
    new_pairs = {
        (
            completed.chain_id[i],
            completed.res_id[i],
            completed.atom_name[i],
            completed.chain_id[j],
            completed.res_id[j],
            completed.atom_name[j],
        )
        for i, j, _ in completed.bonds.as_array()
    }
    assert old_pairs <= new_pairs
    # the completed lysine side chain is chained CB-CG-CD-CE-NZ
    lysine = _lysine_30(completed)
    names = lysine.atom_name.tolist()
    bonded = {frozenset((names[i], names[j])) for i, j, _ in lysine.bonds.as_array()}
    assert {
        frozenset(p) for p in [("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")]
    } <= bonded
    assert mask.sum() > 0


def test_add_missing_atoms_is_conservative():
    atoms = struc.AtomArray(6)
    atoms.coord = np.zeros((6, 3))
    atoms.chain_id = np.array(["A"] * 6)
    atoms.res_id = np.array([1, 1, 1, 2, 2, 2])
    atoms.res_name = np.array(["ALA"] * 3 + ["NOTACCD"] * 3)
    atoms.atom_name = np.array(
        [
            "N",
            "CA",
            "CX",
            "C1",
            "C2",
            "C3",
        ]
    )  # CX is not an ALA atom
    atoms.element = np.array(["N", "C", "C", "C", "C", "C"])
    atoms.hetero = np.array([False] * 3 + [True] * 3)

    completed = add_missing_atoms(atoms)

    # off-template naming and unknown components are left exactly as they were
    assert completed.array_length() == 6
    assert not completed.get_annotation(UNRESOLVED_ANNOTATION).any()
    assert ccd_component_template("NOTACCD") is None
    template = ccd_component_template("ALA")
    assert (
        template.heavy == ("N", "CA", "C", "O", "CB", "OXT")
        and "OXT" in template.leaving
    )


def test_completed_atoms_match_the_wwpdb_unobserved_records(test_dir):
    from plinder.data.annotations.cif_utils import _iter_category_rows, read_mmcif_file

    block = list(read_mmcif_file(test_dir / CIF_6M92).values())[0]
    recorded = {
        (row["label_asym_id"], int(row["label_seq_id"]), row["label_atom_id"])
        for row in _iter_category_rows(
            block,
            "pdbx_unobs_or_zero_occ_atoms",
            [
                "polymer_flag",
                "occupancy_flag",
                "label_asym_id",
                "label_seq_id",
                "label_atom_id",
            ],
        )
        if row["polymer_flag"] == "Y" and row["occupancy_flag"] == "1"
    }
    assert len(recorded) == 136

    from plinder.data.annotations.cif_utils import get_label_asym_sequences

    sequences = get_label_asym_sequences(block)
    structure = Structure(
        id="6m92",
        protein_path=test_dir / CIF_6M92,
        protein_sequence=sequences,
        complete_missing_atoms=True,
    )
    atoms = structure.protein_atom_array
    mask = atoms.get_annotation(UNRESOLVED_ANNOTATION)
    completed = {
        (str(c), int(r), str(n))
        for c, r, n in zip(
            atoms.chain_id[mask], atoms.res_id[mask], atoms.atom_name[mask]
        )
    }

    # two independent sources agree: wwPDB's records and the CCD-template diff
    assert completed == recorded
    per_chain = structure.protein_unresolved_atom_mask
    assert [int(m.sum()) for m in per_chain] == [
        len([r for r in recorded if r[0] == chain])
        for chain in structure.protein_chain_ordered
    ]
    assert all(len(m) == len(c) for m, c in zip(per_chain, structure.protein_coords))
    assert all(
        coordinates.dtype == atoms.coord.dtype
        for coordinates in structure.protein_coords + structure.protein_calpha_coords
    )
    default = Structure(
        id="6m92", protein_path=test_dir / CIF_6M92, protein_sequence=sequences
    )
    assert not any(m.any() for m in default.protein_unresolved_atom_mask)
    assert default.protein_atom_array.array_length() == atoms.array_length() - len(
        recorded
    )
