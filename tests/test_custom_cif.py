# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Tests for custom CIF processing with missing bond orders.

These tests verify that:
1. Missing _chem_comp_bond in CIF files is detected and raises an error
2. Known CCD compounds (ATP, etc.) are skipped — no SMILES needed
3. Bond orders can be assigned from SMILES and written to the CIF
4. The enriched CIF can be read back with correct bond information
"""

from __future__ import annotations

import shutil
from pathlib import Path

import biotite.structure.io.pdbx as pdbx
import pytest
from plinder.data.utils.annotations.biotite_utils import (
    MissingBondOrderError,
    assign_bond_orders_from_smiles,
    check_cif_bond_orders,
    get_unknown_ligand_ids,
)

TEST_DATA = Path(__file__).parent / "test_data" / "custom_cif"
BOLTZ_CIF = TEST_DATA / "boltz_8c3u_input_model_0.cif"
LIGAND_SMILES = "Cc1ccc2c(c1)NC(=O)C2(c3cc(ccc3O)c4ccc(cc4C(=O)O)C(=O)O)c5c[nH]nc5"


@pytest.fixture
def boltz_cif(tmp_path):
    """Copy the Boltz CIF to a temp dir so tests can modify it."""
    dst = tmp_path / "boltz_model.cif"
    shutil.copy(BOLTZ_CIF, dst)
    return dst


def test_boltz_cif_has_no_bond_orders(boltz_cif):
    """Boltz output CIF should have no _chem_comp_bond category."""
    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    assert "chem_comp_bond" not in block


def test_unknown_ligand_ids_detects_lig(boltz_cif):
    """LIG is not in CCD, so it should be flagged as unknown."""
    unknown = get_unknown_ligand_ids(boltz_cif)
    assert "LIG" in unknown


def test_known_compounds_not_flagged(boltz_cif):
    """Known CCD compounds like ATP should not be flagged as unknown."""
    # Inject a fake ATP HETATM into the CIF to verify it gets skipped
    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    atom_site = block["atom_site"]

    # Read all columns and append one ATP row
    columns = {}
    for col_name in atom_site.keys():
        arr = list(atom_site[col_name].as_array())
        # Copy the last row and modify it
        arr.append(arr[-1])
        columns[col_name] = arr

    # Set the last row to be ATP
    n = len(columns["group_PDB"]) - 1
    columns["group_PDB"][n] = "HETATM"
    columns["label_comp_id"][n] = "ATP"

    block["atom_site"] = pdbx.CIFCategory(columns)
    modified = boltz_cif.parent / "with_atp.cif"
    f.write(str(modified))

    unknown = get_unknown_ligand_ids(modified)
    assert "ATP" not in unknown, "ATP is a known CCD compound, should not be flagged"
    assert "LIG" in unknown, "LIG should still be flagged"


def test_check_cif_bond_orders_raises_on_unknown(boltz_cif):
    """check_cif_bond_orders should raise for unknown ligands without bonds."""
    with pytest.raises(MissingBondOrderError, match="unknown ligands"):
        check_cif_bond_orders(boltz_cif)


def test_assign_bond_orders_from_smiles(boltz_cif):
    """Assigning bond orders from SMILES should write _chem_comp_bond."""
    output = boltz_cif.parent / "enriched.cif"
    assign_bond_orders_from_smiles(
        boltz_cif,
        ligand_smiles={"LIG": LIGAND_SMILES},
        output_path=output,
    )

    f = pdbx.CIFFile.read(str(output))
    block = list(f.values())[0]
    assert "chem_comp_bond" in block

    bond_cat = block["chem_comp_bond"]
    comp_ids = bond_cat["comp_id"].as_array()
    orders = set(bond_cat["value_order"].as_array())

    assert all(c == "LIG" for c in comp_ids)
    assert len(comp_ids) > 0
    assert "SING" in orders or "AROM" in orders
    assert "DOUB" in orders or "AROM" in orders


def test_check_passes_after_enrichment(boltz_cif):
    """After enrichment, check_cif_bond_orders should not raise."""
    assign_bond_orders_from_smiles(
        boltz_cif,
        ligand_smiles={"LIG": LIGAND_SMILES},
    )
    check_cif_bond_orders(boltz_cif)


def test_assign_skips_known_compounds(boltz_cif):
    """Providing SMILES for a known compound should be silently skipped."""
    assign_bond_orders_from_smiles(
        boltz_cif,
        ligand_smiles={
            "LIG": LIGAND_SMILES,
            "ATP": "dummy_will_be_skipped",  # ATP is known, won't be processed
        },
    )
    check_cif_bond_orders(boltz_cif)


def test_assign_missing_smiles_raises(boltz_cif):
    """Not providing SMILES for an unknown ligand should raise."""
    with pytest.raises(MissingBondOrderError, match="need SMILES"):
        assign_bond_orders_from_smiles(
            boltz_cif,
            ligand_smiles={},  # LIG is unknown but no SMILES given
        )


def test_assign_invalid_smiles_raises(boltz_cif):
    """Invalid SMILES should raise ValueError."""
    with pytest.raises(ValueError, match="Invalid SMILES"):
        assign_bond_orders_from_smiles(
            boltz_cif,
            ligand_smiles={"LIG": "not_a_smiles!!!"},
        )


# ---------------------------------------------------------------------------
# Integration tests: Entry.from_custom_cif_file
# ---------------------------------------------------------------------------


def test_from_custom_cif_raises_without_smiles(boltz_cif):
    """from_custom_cif_file should raise when unknown ligands lack SMILES."""
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    with pytest.raises(MissingBondOrderError):
        Entry.from_custom_cif_file(
            pdb_id="8c3u",
            cif_file=boltz_cif,
        )


def test_from_custom_cif_with_smiles(boltz_cif):
    """from_custom_cif_file should succeed when SMILES are provided."""
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    entry = Entry.from_custom_cif_file(
        pdb_id="8c3u",
        cif_file=boltz_cif,
        ligand_smiles_dict={"LIG": LIGAND_SMILES},
    )
    assert entry.pdb_id == "8c3u"
    assert len(entry.systems) > 0, "Should detect at least one system"

    # Verify the CIF was enriched in-place
    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    assert "chem_comp_bond" in block
