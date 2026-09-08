# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from plinder.core.structure.ccd_template import (
    ccd_heavy_atom_names,
    unresolved_atoms_from_template,
)
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import read_mmcif_file
from plinder.data.annotations.protein_utils import get_unobserved_atoms

CIF_8PN3 = "xx/pdb_00008pn3/pdb_00008pn3_xyz-enrich.cif.gz"
CIF_6M92 = "xx/pdb_00006m92/pdb_00006m92_xyz-enrich.cif.gz"
CIF_19HC = "xx/pdb_000019hc/pdb_000019hc_xyz-enrich.cif.gz"
CIF_5LWX = "xx/pdb_00005lwx/pdb_00005lwx_xyz-enrich.cif.gz"


def test_ccd_template_diff_ignores_leaving_atoms_and_unknown_components():
    heavy, leaving = ccd_heavy_atom_names("A2G")
    assert (
        "O1" in leaving and "C8" in heavy and not any(n.startswith("H") for n in heavy)
    )
    resolved = heavy - {"O1", "C8"}  # glycosidic O1 leaves on linkage; C8 truly missing

    assert unresolved_atoms_from_template("A2G", resolved) == ["C8"]
    assert unresolved_atoms_from_template("ALA", ["N", "CA", "C", "O", "CB"]) == []
    assert unresolved_atoms_from_template("ALA", ["N", "CA", "C", "O"]) == ["CB"]
    assert unresolved_atoms_from_template("NOTACCDCODE", ["C1"]) is None


def test_unobserved_atom_records_are_indexed_by_residue_and_chain(test_dir):
    by_residue, by_chain = get_unobserved_atoms(
        list(read_mmcif_file(test_dir / CIF_8PN3).values())[0]
    )
    assert by_residue["C"][13] == ["CA", "C", "O", "CB"]
    assert by_chain["C"] == [("ALA", "13", name) for name in ("CA", "C", "O", "CB")]

    _, by_chain = get_unobserved_atoms(
        list(read_mmcif_file(test_dir / CIF_19HC).values())[0]
    )
    assert by_chain["R"] == [("HEM", "302", "CGD"), ("HEM", "302", "O2D")]


def test_pocket_residues_carry_unobserved_atoms(test_dir):
    entry = Entry.from_cif_file(test_dir / CIF_6M92, include_interfaces=False)

    assert entry.chains["A"].residues[30].unresolved_atom_names == [
        "CG",
        "CD",
        "CE",
        "NZ",
    ]
    assert entry.chains["A"].residues[31].unresolved_atom_names == []


def test_ligand_unresolved_atoms_from_records_and_pocket_export(test_dir):
    entry = Entry.from_cif_file(test_dir / CIF_19HC, include_interfaces=False)
    ligands = {
        ligand.asym_id: ligand
        for system in entry.systems.values()
        for ligand in system.ligands
    }

    heme = ligands["R"]
    assert heme.ccd_code == "HEM"
    assert heme.unresolved_atoms == ["R:HEM:302:CGD", "R:HEM:302:O2D"]
    assert heme.num_unresolved_heavy_atoms == 2
    complete = [ligand for ligand in ligands.values() if not ligand.unresolved_atoms]
    assert complete and all(l.num_unresolved_heavy_atoms == 0 for l in complete)

    row = heme.format(entry.chains)
    assert row["ligand_unresolved_atoms"] == heme.unresolved_atoms
    for item in row["ligand_pocket_unresolved_atoms"]:
        instance_chain, residue_number, atom_name = item.split("_")
        residue = entry.chains[instance_chain.split(".")[-1]].residues[
            int(residue_number)
        ]
        assert atom_name in residue.unresolved_atom_names
        assert int(residue_number) in heme.neighboring_residues[instance_chain]


def test_glycan_leaving_atoms_are_not_reported_as_unresolved(test_dir):
    # TODO: 5lwx glycans are linked, so dropping O1 is right; an unlinked sugar
    # with a missing O1 should be reported once the exclusion is link-aware.
    block = list(read_mmcif_file(test_dir / CIF_5LWX).values())[0]
    _, by_chain = get_unobserved_atoms(block)
    assert any(atom == "O1" for rows in by_chain.values() for _, _, atom in rows)

    entry = Entry.from_cif_file(test_dir / CIF_5LWX, include_interfaces=False)
    sugars = [
        ligand
        for system in entry.systems.values()
        for ligand in system.ligands
        if {"NAG", "BMA", "MAN"} & set(ligand.ccd_code.split("-"))
    ]
    assert sugars
    assert not any(
        item.endswith(":O1") for ligand in sugars for item in ligand.unresolved_atoms
    )
