# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pytest

from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import (
    find_ph_mentions,
    get_entry_ph_range,
    get_ligand_of_interest,
    parse_ph_range,
    read_mmcif_file,
)


def _block(test_dir, relative_path):
    return list(read_mmcif_file(test_dir / relative_path).values())[0]


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("7.5", (7.5, 7.5)),
        ("7.0-8.0", (7.0, 8.0)),
        ("pH8.0", (8.0, 8.0)),
        ("pH 6.5 to 7.5", (6.5, 7.5)),
        ("6.5 - 7.5", (6.5, 7.5)),
        ("7.5/8.0", (7.5, 8.0)),
        ("?", None),
        (".", None),
        ("", None),
        (None, None),
        ("15", None),
        ("0.1 M HEPES pH 7.0", None),
    ],
)
def test_parse_ph_range(text, expected):
    assert parse_ph_range(text) == expected


def test_find_ph_mentions_in_free_text():
    assert find_ph_mentions("0.1 M HEPES, 18-26% PEG 3350, 1-6% Tacsimate pH 7.0") == (
        7.0,
        7.0,
    )
    assert find_ph_mentions("crystals at pH 6.5-7.5 and PH8") == (6.5, 8.0)
    assert find_ph_mentions("20% PEG 3350, 0.2 M NaCl") is None
    assert find_ph_mentions("pH 25") is None
    assert find_ph_mentions(None) is None


def test_entry_ph_range_precedence(test_dir):
    # numeric pH column
    assert get_entry_ph_range(
        _block(test_dir, "xx/pdb_00001qz5/pdb_00001qz5_xyz-enrich.cif.gz")
    ) == (5.5, 5.5)
    # pH only inside pdbx_details ("... Tacsimate pH 7.0")
    assert get_entry_ph_range(
        _block(test_dir, "interfaces/bd/pdb_00007bdu/pdb_00007bdu_xyz-enrich.cif.gz")
    ) == (7.0, 7.0)
    # cryo-EM entry without crystal growth rows
    assert get_entry_ph_range(
        _block(test_dir, "interfaces/km/pdb_00007kmx/pdb_00007kmx_xyz-enrich.cif.gz")
    ) == (None, None)


def test_ligand_of_interest_flags(test_dir):
    assert get_ligand_of_interest(
        _block(test_dir, "xx/pdb_00008pn3/pdb_00008pn3_xyz-enrich.cif.gz")
    ) == (True, frozenset({"A2G"}))
    assert get_ligand_of_interest(
        _block(test_dir, "xx/pdb_00001atp/pdb_00001atp_xyz-enrich.cif.gz")
    ) == (None, None)
    assert get_ligand_of_interest(
        _block(test_dir, "interfaces/bd/pdb_00007bdu/pdb_00007bdu_xyz-enrich.cif.gz")
    ) == (False, None)


def test_entry_exports_ph_range_and_ligand_of_interest(test_dir):
    entry = Entry.from_cif_file(
        test_dir / "xx/pdb_00008pn3/pdb_00008pn3_xyz-enrich.cif.gz",
        include_interfaces=False,
    )
    ligands = [ligand for system in entry.systems.values() for ligand in system.ligands]

    assert (entry.pH_min, entry.pH_max) == (6.5, 6.5)
    assert entry.has_ligand_of_interest is True
    # the A2G glycans are covalently linked into glycopeptide ligands
    assert any("A2G" in ligand.ccd_code.split("-") for ligand in ligands)
    assert any(ligand.ccd_code == "GOL" for ligand in ligands)
    for ligand in ligands:
        assert ligand.is_subject_of_investigation is (
            "A2G" in ligand.ccd_code.split("-")
        )
    row = entry.format()
    assert (row["entry_pH_min"], row["entry_pH_max"]) == (6.5, 6.5)
    assert row["entry_has_ligand_of_interest"] is True
    assert ligands[0].format(entry.chains)["ligand_is_subject_of_investigation"] in {
        True,
        False,
    }

    legacy = Entry.from_cif_file(
        test_dir / "xx/pdb_00001atp/pdb_00001atp_xyz-enrich.cif.gz",
        include_interfaces=False,
    )
    assert legacy.has_ligand_of_interest is None
    assert all(
        ligand.is_subject_of_investigation is None
        for system in legacy.systems.values()
        for ligand in system.ligands
    )
