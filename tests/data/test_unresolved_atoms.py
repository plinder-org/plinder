# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pytest
from plinder.core.structure.ccd_template import (
    ccd_heavy_atom_names,
    unresolved_atoms_from_template,
)
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import read_mmcif_file
from plinder.data.annotations.protein_utils import UnobservedAtom, get_unobserved_atoms

CIF_8PN3 = "xx/pdb_00008pn3/pdb_00008pn3_xyz-enrich.cif.gz"
CIF_6M92 = "xx/pdb_00006m92/pdb_00006m92_xyz-enrich.cif.gz"
CIF_19HC = "xx/pdb_000019hc/pdb_000019hc_xyz-enrich.cif.gz"
CIF_5LWX = "xx/pdb_00005lwx/pdb_00005lwx_xyz-enrich.cif.gz"
CIF_7GJ7 = "xx/pdb_00007gj7/pdb_00007gj7_xyz-enrich.cif.gz"
CIF_6LU7 = "xx/pdb_00006lu7/pdb_00006lu7_xyz-enrich.cif.gz"


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
    # a leaving atom departs only when the residue is linked through its anchor
    nag_heavy, _ = ccd_heavy_atom_names("NAG")
    assert unresolved_atoms_from_template("NAG", nag_heavy - {"O1"}, {"C1"}) == []
    assert unresolved_atoms_from_template("NAG", nag_heavy - {"O1"}, set()) == ["O1"]


def test_unobserved_atom_records_are_indexed_by_residue_and_chain(test_dir):
    by_residue, by_chain = get_unobserved_atoms(
        list(read_mmcif_file(test_dir / CIF_8PN3).values())[0]
    )
    assert by_residue["C"][13] == ["CA", "C", "O", "CB"]
    assert by_chain["C"] == [
        UnobservedAtom("ALA", "13", "13", name, False)
        for name in ("CA", "C", "O", "CB")
    ]

    _, by_chain = get_unobserved_atoms(
        list(read_mmcif_file(test_dir / CIF_19HC).values())[0]
    )
    assert by_chain["R"] == [
        UnobservedAtom("HEM", "302", ".", "CGD", True),
        UnobservedAtom("HEM", "302", ".", "O2D", True),
    ]


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
    assert heme.unresolved_atoms == ["302:HEM:R:.:CGD", "302:HEM:R:.:O2D"]
    assert heme.num_unresolved_heavy_atoms == 2
    # same residue address as covalent links
    assert "47:CYS:A:47:SG__301:HEM:F:.:CAB" in ligands["F"].covalent_linkages
    complete = [ligand for ligand in ligands.values() if not ligand.unresolved_atoms]
    assert complete and all(l.num_unresolved_heavy_atoms == 0 for l in complete)

    row = heme.format(entry.chains)
    assert row["ligand_unresolved_atoms"] == heme.unresolved_atoms
    assert row["ligand_pocket_unresolved_atoms"]
    for item in row["ligand_pocket_unresolved_atoms"]:
        auth_seq, comp_id, asym_id, label_seq, atom_name = item.split(":")
        residue = entry.chains[asym_id].residues[int(label_seq)]
        assert (residue.auth_number, residue.name) == (auth_seq, comp_id)
        assert atom_name in residue.unresolved_atom_names
        assert any(
            int(label_seq) in numbers and chain.endswith(f".{asym_id}")
            for chain, numbers in heme.neighboring_residues.items()
        )


def test_leaving_atoms_are_missing_only_on_unlinked_residues(test_dir):
    entry = Entry.from_cif_file(test_dir / CIF_5LWX, include_interfaces=False)
    ligands = {
        ligand.instance_chain: ligand
        for system in entry.systems.values()
        for ligand in system.ligands
    }
    # NAG-NAG on Asn216: C1 is linked, so the departed O1 is not missing
    assert ligands["1.C"].is_covalent and ligands["1.C"].unresolved_atoms == []
    # free NAG-NAG: the reducing-end O1 is a real gap, the second sugar's O1 departed
    assert ligands["1.B"].unresolved_atoms == ["1:NAG:B:.:O1"]
    # free single NAG
    assert ligands["1.H"].unresolved_atoms == ["605:NAG:H:.:O1"]


@pytest.mark.parametrize(
    "cif, instance_chain",
    [(CIF_6LU7, "1.B"), (CIF_7GJ7, "1.E"), (CIF_19HC, "1.F")],
)
def test_covalent_ligands_have_no_missing_atoms(test_dir, cif, instance_chain):
    # 6lu7 N3 (peptide-like, OXT departed on every peptide bond) and 7gj7 Q0I take
    # the CCD template path; 19hc heme F has wwPDB records
    entry = Entry.from_cif_file(test_dir / cif, include_interfaces=False)
    ligand = next(
        ligand
        for system in entry.systems.values()
        for ligand in system.ligands
        if ligand.instance_chain == instance_chain
    )
    assert ligand.is_covalent
    assert ligand.unresolved_atoms == []
    assert ligand.num_unresolved_heavy_atoms == 0
