# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import numpy as np
import pytest

from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import build_biounit, read_mmcif_file
from plinder.data.annotations.contact_areas import (
    chain_pair_contact_areas,
    partner_contact_areas,
    tessellation_atom_mask,
)
from plinder.data.annotations.interface_utils import (
    INTERFACE_ANNOTATION_SCHEMA,
    protein_interfaces_to_table,
)

CIF_1QZ5 = "xx/pdb_00001qz5/pdb_00001qz5_xyz-enrich.cif.gz"
CIF_7CMA = "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"


def _water_chains(atoms) -> set[str]:
    return {
        str(chain_id)
        for chain_id in np.unique(atoms.chain_id)
        if np.all(atoms.res_name[atoms.chain_id == chain_id] == "HOH")
    }


def test_chain_pair_contact_areas_match_reference_and_exclude_solvent(test_dir):
    from voronotalt import biotite_interface  # noqa: F401  registers from_biotite_atoms
    from voronotalt import voronotalt_python as vlt

    biounit = build_biounit(read_mmcif_file(test_dir / CIF_1QZ5), "1")
    water_chains = _water_chains(biounit)
    assert water_chains, "fixture must carry a water chain to prove its exclusion"
    assert not tessellation_atom_mask(biounit)[
        np.isin(biounit.chain_id, sorted(water_chains))
    ].any()

    areas = chain_pair_contact_areas(biounit)

    assert areas
    assert all(chain_a < chain_b for chain_a, chain_b in areas)
    assert all(area > 0 for area in areas.values())
    assert not any(water_chains & set(pair) for pair in areas)
    # Voronota-LT's own biotite shim (waters and hydrogens dropped, default probe)
    # is the reference for the vectorised ball construction.
    reference = vlt.MolecularRadicalTessellation.from_biotite_atoms(biounit)
    expected = {
        tuple(sorted((summary.ID1_chain, summary.ID2_chain))): summary.area
        for summary in reference.inter_chain_contact_summaries
    }
    assert areas.keys() == expected.keys()
    for pair, area in areas.items():
        assert area == pytest.approx(expected[pair], rel=1e-9)


def _three_chain_toy(distance_ab: float, occluder: bool):
    """Two carbon atoms in chains A and B ``distance_ab`` apart on the x axis,
    optionally with a chain-C carbon at their midpoint; no ins_code annotation."""
    import biotite.structure as struc

    count = 3 if occluder else 2
    atoms = struc.AtomArray(count)
    atoms.coord = np.zeros((count, 3))
    atoms.coord[1, 0] = distance_ab
    if occluder:
        atoms.coord[2, 0] = distance_ab / 2.0
    atoms.chain_id = np.array(["A", "B", "C"][:count])
    atoms.res_id = np.arange(1, count + 1)
    atoms.res_name = np.array(["ALA"] * count)
    atoms.atom_name = np.array(["CA"] * count)
    atoms.element = np.array(["C"] * count)
    atoms.hetero = np.zeros(count, dtype=bool)
    if "ins_code" in atoms.get_annotation_categories():
        atoms.del_annotation("ins_code")
    return atoms


def test_chain_pair_contact_areas_follow_contact_geometry():
    # Two carbons 4 Å apart overlap once each is grown by the 1.4 Å probe.
    touching = chain_pair_contact_areas(_three_chain_toy(4.0, occluder=False))
    assert list(touching) == [("A", "B")]
    assert touching[("A", "B")] > 0

    # Beyond the summed probe-expanded radii there is no contact at all.
    assert chain_pair_contact_areas(_three_chain_toy(12.0, occluder=False)) == {}

    # A third chain sitting between them takes over part of the shared plane.
    occluded = chain_pair_contact_areas(_three_chain_toy(4.0, occluder=True))
    assert occluded.get(("A", "B"), 0.0) < touching[("A", "B")]
    assert occluded[("A", "C")] > 0 and occluded[("B", "C")] > 0
    assert all(chain_a < chain_b for chain_a, chain_b in occluded)

    # Hydrogens are ignored, so adding one changes nothing.
    with_h = _three_chain_toy(4.0, occluder=True)
    hydrogen = with_h[2:3].copy()
    hydrogen.element = np.array(["H"])
    hydrogen.atom_name = np.array(["HA"])
    hydrogen.coord = hydrogen.coord + 0.7
    assert chain_pair_contact_areas(with_h + hydrogen) == pytest.approx(occluded)


def test_partner_contact_areas_sums_over_member_chains_and_skips_internal_pairs():
    pair_areas = {
        ("1.A", "1.C"): 10.0,
        ("1.C", "1.D"): 5.0,  # internal to the two-chain ligand
        ("1.A", "1.D"): 2.5,
        ("1.B", "1.C"): 1.0,
        ("1.A", "1.B"): 7.0,  # receptor-receptor, not a partner pair
    }

    assert partner_contact_areas(pair_areas, {"1.C", "1.D"}) == {
        "1.A": 12.5,
        "1.B": 1.0,
    }
    assert partner_contact_areas(pair_areas, {"1.Z"}) == {}


def test_ligand_contact_area_counts_receptor_chains_only(test_dir):
    """ATP in 1qz5 touches the receptor chain and a calcium ion: both are reported
    as partners, but only the receptor counts towards ``contact_area``."""
    entry = Entry.from_cif_file(test_dir / CIF_1QZ5, include_interfaces=False)

    assert entry.failed_contact_area_biounit_ids == []
    [atp] = [
        ligand
        for system in entry.systems.values()
        for ligand in system.ligands
        if ligand.ccd_code == "ATP"
    ]
    partners = atp.chain_contact_areas
    ion_partners = [c for c in partners if c.split(".")[-1] in entry.ligand_like_chains]
    receptor_partners = [c for c in partners if c not in ion_partners]
    assert ion_partners and receptor_partners
    assert all(area > 0 for area in partners.values())
    assert atp.contact_area == pytest.approx(
        sum(partners[chain_id] for chain_id in receptor_partners)
    )
    assert atp.contact_area > 100.0
    assert atp.contact_area < sum(partners.values())

    row = atp.format(entry.chains)
    assert row["ligand_contact_area"] == atp.contact_area
    assert row["ligand_contact_area_chains"] == sorted(partners)
    assert row["ligand_contact_area_values"] == [partners[c] for c in sorted(partners)]


def test_protein_interface_reports_its_contact_area(test_dir):
    entry = Entry.from_custom_cif_file(
        pdb_id="custom_7cma",
        cif_file=test_dir / CIF_7CMA,
        structure_mode="as_is",
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )
    [interface] = entry.interfaces

    assert interface.system_id == "custom_7cma__1__1.A--1.B"
    assert interface.contact_area is not None and interface.contact_area > 100.0
    assert interface.to_row()["interface_contact_area"] == interface.contact_area
    table = protein_interfaces_to_table(entry.interfaces)
    assert table.schema.equals(INTERFACE_ANNOTATION_SCHEMA, check_metadata=False)
    assert table.column("interface_contact_area").to_pylist() == [
        interface.contact_area
    ]


def test_contact_area_size_guard_records_the_skipped_assembly(test_dir):
    entry = Entry.from_cif_file(
        test_dir / CIF_1QZ5, include_interfaces=False, tessellation_atom_limit=10
    )

    assert entry.failed_contact_area_biounit_ids == ["1"]
    ligands = [ligand for system in entry.systems.values() for ligand in system.ligands]
    assert ligands
    for ligand in ligands:
        assert ligand.contact_area is None
        assert ligand.chain_contact_areas == {}
        row = ligand.format(entry.chains)
        assert row["ligand_contact_area"] is None
        assert row["ligand_contact_area_chains"] == []
        assert row["ligand_contact_area_values"] == []


def test_failed_tessellation_records_the_assembly_and_leaves_areas_null(
    test_dir, monkeypatch
):
    from plinder.data.annotations import aggregate_annotations

    def explode(_atoms):
        raise RuntimeError("voronota exploded")

    monkeypatch.setattr(aggregate_annotations, "chain_pair_contact_areas", explode)
    entry = Entry.from_custom_cif_file(
        pdb_id="custom_7cma",
        cif_file=test_dir / CIF_7CMA,
        structure_mode="as_is",
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )

    assert entry.failed_contact_area_biounit_ids == ["1"]
    [interface] = entry.interfaces
    assert interface.contact_area is None
    assert interface.to_row()["interface_contact_area"] is None


def test_protein_only_ingest_never_tessellates(test_dir, monkeypatch):
    from plinder.data.annotations import aggregate_annotations

    def explode(_atoms):
        raise AssertionError("tessellation must not run without consumers")

    monkeypatch.setattr(aggregate_annotations, "chain_pair_contact_areas", explode)
    entry = Entry.from_cif_file(
        test_dir / CIF_1QZ5,
        include_ligands=False,
        include_interfaces=False,
        protein_only=True,
    )

    assert entry.chains
    assert entry.failed_contact_area_biounit_ids == []
