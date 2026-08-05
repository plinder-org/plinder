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

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
import pytest
import yaml
from plinder.data.annotations.cif_utils import (
    MissingBondOrderError,
    apply_struct_conn_bonds,
    assign_bond_orders_from_smiles,
    build_biounit,
    check_cif_bond_orders,
    get_legacy_chain_instance_mapping,
    get_structure_with_altloc,
    get_unknown_ligand_ids,
)

CUSTOM_CIF_DIR = Path(__file__).parent / "test_data" / "custom_cif"
BOLTZ_CIF = CUSTOM_CIF_DIR / "boltz_8c3u_input_model_0.cif"
BOLTZ_INPUT_YAML = CUSTOM_CIF_DIR / "boltz_8c3u_input.yaml"


def _load_boltz_ligand_smiles() -> str:
    """Parse the ligand SMILES from the Boltz input YAML (single source of truth)."""
    config = yaml.safe_load(BOLTZ_INPUT_YAML.read_text())
    for seq in config["sequences"]:
        if "ligand" in seq:
            return seq["ligand"]["smiles"]
    raise ValueError(f"No ligand SMILES found in {BOLTZ_INPUT_YAML}")


LIGAND_SMILES = _load_boltz_ligand_smiles()


def test_deposited_first_altloc_accepts_numeric_ids():
    atoms = struc.AtomArray(5)
    atoms.coord = np.arange(15).reshape(5, 3)
    atoms.chain_id = np.array(["A"] * 5)
    atoms.res_id = np.array([1] * 5)
    atoms.res_name = np.array(["LIG"] * 5)
    atoms.atom_name = np.array(["C", "X", "X", "Y", "Y"])
    atoms.element = np.array(["C"] * 5)
    atoms.hetero = np.ones(5, dtype=bool)
    atoms.set_annotation("occupancy", np.array([1.0, 0.486, 0.514, 0.486, 0.514]))
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, atoms)
    block = cif_file.block
    original_ids = [".", "1", "2", "1", "2"]
    block["atom_site"]["label_alt_id"] = pdbx.CIFColumn(original_ids)

    filtered = get_structure_with_altloc(
        cif_file,
        model=1,
        use_author_fields=False,
        extra_fields=["occupancy"],
    )

    assert filtered.atom_name.tolist() == ["C", "X", "Y"]
    assert filtered.occupancy.tolist() == [1.0, 0.486, 0.486]
    assert filtered.selected_altloc_id.tolist() == ["1", "1", "1"]
    assert block["atom_site"]["label_alt_id"].as_array(str).tolist() == original_ids


def test_selected_altloc_survives_hydrogen_removal():
    atoms = struc.AtomArray(3)
    atoms.coord = np.arange(9).reshape(3, 3)
    atoms.chain_id = np.array(["A"] * 3)
    atoms.res_id = np.array([1] * 3)
    atoms.res_name = np.array(["LIG"] * 3)
    atoms.atom_name = np.array(["C", "H1", "H1"])
    atoms.element = np.array(["C", "H", "H"])
    atoms.hetero = np.ones(3, dtype=bool)
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, atoms)
    cif_file.block["atom_site"]["label_alt_id"] = pdbx.CIFColumn([".", "A", "B"])

    filtered = get_structure_with_altloc(cif_file)
    heavy_atoms = filtered[filtered.element != "H"]

    assert filtered.selected_altloc_id.tolist() == ["A", "A"]
    assert heavy_atoms.selected_altloc_id.tolist() == ["A"]


def test_build_biounit_normalizes_altlocs_and_uses_deposited_first(monkeypatch):
    atoms = struc.AtomArray(2)
    atoms.coord = np.zeros((2, 3))
    atoms.chain_id = np.array(["A", "A"])
    atoms.res_id = np.array([1, 1])
    atoms.res_name = np.array(["LIG", "LIG"])
    atoms.atom_name = np.array(["C1", "C1"])
    atoms.element = np.array(["C", "C"])
    atoms.hetero = np.ones(2, dtype=bool)
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, atoms)
    original_ids = ["1", "2"]
    cif_file.block["atom_site"]["label_alt_id"] = pdbx.CIFColumn(original_ids)

    def fake_get_assembly(cif, **kwargs):
        assert kwargs["altloc"] == "first"
        assert cif.block["atom_site"]["label_alt_id"].as_array(str).tolist() == [
            "A",
            "B",
        ]
        assembly = struc.AtomArray(1)
        assembly.coord = np.zeros((1, 3))
        assembly.chain_id = np.array(["A"])
        assembly.res_id = np.array([1])
        assembly.res_name = np.array(["LIG"])
        assembly.atom_name = np.array(["C1"])
        assembly.element = np.array(["C"])
        assembly.set_annotation("sym_id", np.array([0]))
        assembly.set_annotation("label_asym_id", np.array(["A"]))
        assembly.bonds = struc.BondList(1)
        return assembly

    monkeypatch.setattr(pdbx, "get_assembly", fake_get_assembly)

    assembly = build_biounit(cif_file, "1")

    assert assembly.chain_id.tolist() == ["1.A"]
    assert assembly.legacy_chain_id.tolist() == ["1.A"]
    assert (
        cif_file.block["atom_site"]["label_alt_id"].as_array(str).tolist()
        == original_ids
    )


def test_legacy_chain_instance_mapping_uses_global_operation_order():
    block = pdbx.CIFBlock()
    block["pdbx_struct_assembly_gen"] = pdbx.CIFCategory(
        {
            "assembly_id": ["1", "1", "1", "2"],
            "oper_expression": ["1", "(2-3)", "(4,5)(6-7)", "1"],
            "asym_id_list": ["A,C", "B,C", "A", "A"],
        }
    )

    assert get_legacy_chain_instance_mapping(block, "1") == {
        "1.A": "1.A",
        "2.A": "4.A",
        "3.A": "5.A",
        "4.A": "6.A",
        "5.A": "7.A",
        "1.B": "2.B",
        "2.B": "3.B",
        "1.C": "1.C",
        "2.C": "2.C",
        "3.C": "3.C",
    }


def test_apply_struct_conn_bonds_indexes_partners_per_assembly_instance():
    atoms = struc.AtomArray(8)
    atoms.chain_id = np.array(["1.A", "1.A", "1.B", "1.B", "2.A", "2.A", "2.B", "2.B"])
    atoms.set_annotation(
        "label_asym_id", np.array(["A", "A", "B", "B", "A", "A", "B", "B"])
    )
    atoms.set_annotation("sym_id", np.array([0, 0, 0, 0, 1, 1, 1, 1]))
    atoms.res_id = np.array([1, 1, 2, 2, 1, 1, 2, 2])
    atoms.atom_name = np.array(["C1", "X", "N1", "Y", "C1", "X", "N1", "Y"])
    atoms.bonds = struc.BondList(len(atoms))
    atoms.bonds.add_bond(0, 2, struc.BondType.DOUBLE)

    block = pdbx.CIFBlock()
    block["struct_conn"] = pdbx.CIFCategory(
        {
            "conn_type_id": ["covale"],
            "ptnr1_label_asym_id": ["A"],
            "ptnr1_label_seq_id": ["1"],
            "ptnr1_label_atom_id": ["C1"],
            "ptnr1_label_comp_id": ["L1"],
            "ptnr2_label_asym_id": ["B"],
            "ptnr2_label_seq_id": ["2"],
            "ptnr2_label_atom_id": ["N1"],
            "ptnr2_label_comp_id": ["L2"],
            "ptnr1_auth_seq_id": ["1"],
            "ptnr2_auth_seq_id": ["2"],
        }
    )

    apply_struct_conn_bonds(atoms, block)

    neighbors0, bond_types0 = atoms.bonds.get_bonds(0)
    assert list(neighbors0) == [2]
    assert list(bond_types0) == [struc.BondType.DOUBLE]
    neighbors4, bond_types4 = atoms.bonds.get_bonds(4)
    assert list(neighbors4) == [6]
    assert list(bond_types4) == [struc.BondType.SINGLE]


@pytest.fixture
def boltz_cif(tmp_path):
    """Copy Boltz CIF to temp dir — tests modify it in-place."""
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
    # Inject fake ATP HETATMs into the CIF (enough to match CCD atom count)
    import biotite.structure.info as info
    from plinder.core.structure.atoms import is_hydrogen_isotope

    atp_ref = info.residue("ATP")
    atp_heavy = atp_ref[~is_hydrogen_isotope(atp_ref.element)]

    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    atom_site = block["atom_site"]

    columns = {}
    for col_name in atom_site.keys():
        columns[col_name] = list(atom_site[col_name].as_array())

    # Add ATP atoms with correct CCD atom names
    for i in range(len(atp_heavy)):
        for col_name in columns:
            columns[col_name].append(columns[col_name][-1])
        n = len(columns["group_PDB"]) - 1
        columns["group_PDB"][n] = "HETATM"
        columns["label_comp_id"][n] = "ATP"
        columns["label_atom_id"][n] = atp_heavy.atom_name[i]
        if "type_symbol" in columns:
            columns["type_symbol"][n] = atp_heavy.element[i]

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


def test_assign_atom_count_mismatch_raises(boltz_cif):
    """Default positional path should raise when heavy-atom count differs."""
    # Truncated SMILES — fewer atoms than the CIF ligand
    short_smiles = "CC"
    with pytest.raises(ValueError, match="Atom count mismatch"):
        assign_bond_orders_from_smiles(
            boltz_cif,
            ligand_smiles={"LIG": short_smiles},
        )


def test_assign_element_mismatch_raises(boltz_cif):
    """Default positional path should raise when elements don't match.

    Same heavy-atom count as LIG but with a different first-atom element
    (N instead of C) to force a position-0 element mismatch.
    """
    # LIG has 35 heavy atoms starting with C (methyl group). Build a
    # SMILES with the same count but starting with N to trigger a
    # position-0 element mismatch.
    lig_atom_count = 35
    mismatched_smiles = "N" + "C" * (lig_atom_count - 1)
    with pytest.raises(ValueError, match="Element mismatch.*position 0"):
        assign_bond_orders_from_smiles(
            boltz_cif,
            ligand_smiles={"LIG": mismatched_smiles},
        )


def test_assign_force_substructure_match_succeeds(boltz_cif):
    """Opt-in substructure match path should still work end-to-end."""
    assign_bond_orders_from_smiles(
        boltz_cif,
        ligand_smiles={"LIG": LIGAND_SMILES},
        force_substructure_match=True,
    )
    check_cif_bond_orders(boltz_cif)


def test_assign_rejects_divergent_atom_naming(boltz_cif, tmp_path):
    """Multi-instance custom comp_ids with divergent atom names must raise.

    Duplicate the LIG residue and rename atom 0 of the second instance —
    mmCIF ``_chem_comp_bond`` keys by comp_id, so biotite would silently
    fail to apply bonds to the second instance. We must raise a clear
    error rather than emit chemically wrong/incomplete bonds.
    """
    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    atom_site = block["atom_site"]
    columns = {col: list(atom_site[col].as_array()) for col in atom_site.keys()}
    lig_indices = [i for i, c in enumerate(columns["label_comp_id"]) if c == "LIG"]
    next_atom_id = max(int(x) for x in columns["id"]) + 1 if "id" in columns else None

    new_indices = []
    for src in lig_indices:
        for col_name in columns:
            columns[col_name].append(columns[col_name][src])
        n = len(columns["label_comp_id"]) - 1
        columns["label_asym_id"][n] = "C"
        if "auth_asym_id" in columns:
            columns["auth_asym_id"][n] = "C"
        if next_atom_id is not None:
            columns["id"][n] = str(next_atom_id)
            next_atom_id += 1
        new_indices.append(n)
    # Rename one atom in the duplicate instance to break alignment
    columns["label_atom_id"][new_indices[0]] = "X_RENAMED"

    block["atom_site"] = pdbx.CIFCategory(columns)
    divergent = tmp_path / "divergent.cif"
    f.write(str(divergent))

    with pytest.raises(ValueError, match="disagree on heavy-atom"):
        assign_bond_orders_from_smiles(divergent, ligand_smiles={"LIG": LIGAND_SMILES})


def test_assign_handles_multi_instance_comp_id(boltz_cif, tmp_path):
    """Multi-instance custom comp_ids must enrich without atom-count mismatch.

    Duplicate the LIG residue in the CIF so the file has 2 instances of
    the same comp_id. The positional path used to fail with an atom-count
    mismatch (2 * 35 != 35); now ``enrich_cif_with_smiles_bonds`` picks
    one representative instance and writes a single ``_chem_comp_bond``
    entry that biotite applies to all copies.
    """
    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    atom_site = block["atom_site"]

    columns = {col: list(atom_site[col].as_array()) for col in atom_site.keys()}
    lig_indices = [i for i, c in enumerate(columns["label_comp_id"]) if c == "LIG"]
    assert lig_indices, "test setup expects original LIG atoms"

    # Duplicate every LIG row, change the chain to 'C' to mark the second
    # instance as a distinct copy (same comp_id, different chain/res_id).
    next_atom_id = max(int(x) for x in columns["id"]) + 1 if "id" in columns else None
    for src in lig_indices:
        for col_name in columns:
            columns[col_name].append(columns[col_name][src])
        n = len(columns["label_comp_id"]) - 1
        columns["label_asym_id"][n] = "C"
        if "auth_asym_id" in columns:
            columns["auth_asym_id"][n] = "C"
        if next_atom_id is not None:
            columns["id"][n] = str(next_atom_id)
            next_atom_id += 1

    block["atom_site"] = pdbx.CIFCategory(columns)
    duplicated = tmp_path / "duplicated.cif"
    f.write(str(duplicated))

    # Should NOT raise atom count mismatch
    output = tmp_path / "enriched.cif"
    assign_bond_orders_from_smiles(
        duplicated, ligand_smiles={"LIG": LIGAND_SMILES}, output_path=output
    )
    block = list(pdbx.CIFFile.read(str(output)).values())[0]
    bond_cat = block["chem_comp_bond"]
    lig_bonds = sum(1 for c in bond_cat["comp_id"].as_array() if c == "LIG")
    # Bonds defined exactly once for the comp_id, regardless of N copies
    from rdkit import Chem

    template = Chem.MolFromSmiles(LIGAND_SMILES)
    expected_bonds = Chem.RemoveHs(template, sanitize=False).GetNumBonds()
    assert (
        lig_bonds == expected_bonds
    ), f"Expected {expected_bonds} LIG bonds (one per template bond), got {lig_bonds}"


# ---------------------------------------------------------------------------
# Integration tests: Entry.from_custom_cif_file
# ---------------------------------------------------------------------------


def test_custom_cif_modes_distinguish_assembled_from_deposited_pdb(test_dir):
    from plinder.data.annotations.aggregate_annotations import Entry

    cif = test_dir / "interfaces/cm/pdb_00007cm8/pdb_00007cm8_xyz-enrich.cif.gz"
    assembled = Entry.from_custom_cif_file(
        pdb_id="custom_7cm8",
        cif_file=cif,
        structure_mode="as_is",
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )
    deposited = Entry.from_custom_cif_file(
        pdb_id=None,
        cif_file=cif,
        structure_mode="pdb",
        assembly_ids=["1"],
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )

    assert assembled.pdb_id == "custom_7cm8"
    assert assembled.interfaces == []
    assert deposited.pdb_id == "7cm8"
    assert [interface.system_id for interface in deposited.interfaces] == [
        "7cm8__1__1.A--2.A"
    ]
    assert set(deposited.biounit_chain_ids) == {"1"}


def test_custom_as_is_mode_detects_interface_in_supplied_coordinates(test_dir):
    from plinder.data.annotations.aggregate_annotations import Entry

    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    entry = Entry.from_custom_cif_file(
        pdb_id="custom_7cma",
        cif_file=cif,
        structure_mode="as_is",
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )

    assert [interface.system_id for interface in entry.interfaces] == [
        "custom_7cma__1__1.A--1.B"
    ]


def test_custom_interface_only_mode_preserves_saved_ligands(test_dir, tmp_path):
    from plinder.data.annotations.aggregate_annotations import Entry

    ligand_file = tmp_path / "custom_7cma/ligand_files/A.sdf"
    ligand_file.parent.mkdir(parents=True)
    ligand_file.write_text("existing canonical ligand")
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"

    Entry.from_custom_cif_file(
        pdb_id="custom_7cma",
        cif_file=cif,
        structure_mode="as_is",
        save_folder=tmp_path,
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )

    assert ligand_file.read_text() == "existing canonical ligand"


def test_custom_pdb_output_does_not_infer_data_dir(
    cif_6i41, tmp_path, monkeypatch
):
    from plinder.data.annotations import aggregate_annotations as agg
    from plinder.data.annotations.aggregate_annotations import Entry

    def unexpected_reference_lookup(_data_dir):
        pytest.fail("save_folder must not be interpreted as data_dir")

    monkeypatch.setattr(agg, "get_artifact_codes", unexpected_reference_lookup)
    Entry.from_custom_cif_file(
        pdb_id=None,
        cif_file=cif_6i41,
        structure_mode="pdb",
        assembly_ids=["1"],
        save_folder=tmp_path / "standalone-output",
        include_ligands=True,
        include_interfaces=False,
    )

    assert list((tmp_path / "standalone-output/6i41/ligand_files").glob("*.sdf"))


def test_custom_cif_rejects_legacy_pdb_path():
    from plinder.data.annotations.aggregate_annotations import Entry

    with pytest.raises(ValueError, match="must be mmCIF"):
        Entry.from_custom_cif_file(
            pdb_id="custom",
            cif_file=Path("model.pdb"),
            include_ligands=False,
            include_interfaces=True,
        )


def test_custom_cif_rejects_assembly_selection_in_as_is_mode(boltz_cif):
    from plinder.data.annotations.aggregate_annotations import Entry

    with pytest.raises(ValueError, match="only valid in pdb mode"):
        Entry.from_custom_cif_file(
            pdb_id="custom",
            cif_file=boltz_cif,
            assembly_ids=["1"],
        )


def test_custom_pdb_mode_uses_deposited_entry_id(boltz_cif):
    from plinder.data.annotations.aggregate_annotations import Entry

    with pytest.raises(ValueError, match="pdb_id must be None"):
        Entry.from_custom_cif_file(
            pdb_id="alias",
            cif_file=boltz_cif,
            structure_mode="pdb",
        )


def test_selected_assembly_ids_validates_requested_subset():
    from plinder.data.annotations.aggregate_annotations import _selected_assembly_ids

    assert _selected_assembly_ids(["1", "2"], ["2", "2", "1"]) == ["2", "1"]
    assert _selected_assembly_ids(["1", "2"], "2") == ["2"]
    with pytest.raises(ValueError, match="absent from the mmCIF"):
        _selected_assembly_ids(["1"], ["2"])
    with pytest.raises(ValueError, match="must not be empty"):
        _selected_assembly_ids(["1"], [])


def test_from_custom_cif_warns_on_multi_model(boltz_cif, tmp_path, monkeypatch):
    """Multi-model CIFs (NMR ensembles, multi-sample) warn and use model 1."""
    from plinder.data.annotations import aggregate_annotations as agg
    from plinder.data.annotations.aggregate_annotations import Entry

    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    atom_site = block["atom_site"]
    columns = {col: list(atom_site[col].as_array()) for col in atom_site.keys()}

    # If the input has no model-num column, add one with all "1"s first.
    if "pdbx_PDB_model_num" not in columns:
        columns["pdbx_PDB_model_num"] = ["1"] * len(columns["label_comp_id"])

    # Duplicate every atom under model "2" to create a 2-model CIF.
    n_orig = len(columns["label_comp_id"])
    for i in range(n_orig):
        for col_name in columns:
            columns[col_name].append(columns[col_name][i])
        columns["pdbx_PDB_model_num"][n_orig + i] = "2"

    block["atom_site"] = pdbx.CIFCategory(columns)
    multi = tmp_path / "two_models.cif"
    f.write(str(multi))

    # plinder's setup_logger sets propagate=False, so caplog can't see
    # records via the root logger. Capture LOG.warning calls directly.
    warnings: list[str] = []
    real_warning = agg.LOG.warning
    monkeypatch.setattr(
        agg.LOG,
        "warning",
        lambda msg, *a, **kw: warnings.append(str(msg)) or real_warning(msg, *a, **kw),
    )

    entry = Entry.from_custom_cif_file(
        pdb_id="8c3u",
        cif_file=multi,
        ligand_smiles_dict={"LIG": LIGAND_SMILES},
    )

    # A warning was emitted naming the model count
    assert any(
        "2 models" in w for w in warnings
    ), f"Expected warning about 2 models, got: {warnings}"
    # Parsing succeeded using model 1 — entry has the same systems as
    # the single-model run.
    single_entry = Entry.from_custom_cif_file(
        pdb_id="8c3u",
        cif_file=boltz_cif,
        ligand_smiles_dict={"LIG": LIGAND_SMILES},
    )
    assert sorted(entry.systems.keys()) == sorted(single_entry.systems.keys())


def test_from_custom_cif_raises_without_smiles(boltz_cif):
    """from_custom_cif_file should raise when unknown ligands lack SMILES."""
    from plinder.data.annotations.aggregate_annotations import Entry

    with pytest.raises(MissingBondOrderError):
        Entry.from_custom_cif_file(
            pdb_id="8c3u",
            cif_file=boltz_cif,
        )


def test_from_custom_cif_with_smiles(boltz_cif):
    """from_custom_cif_file should succeed when SMILES are provided.

    The input CIF must not be mutated on disk — bond-order enrichment
    happens on an in-memory copy.
    """
    from plinder.data.annotations.aggregate_annotations import Entry

    before_bytes = boltz_cif.read_bytes()

    entry = Entry.from_custom_cif_file(
        pdb_id="8c3u",
        cif_file=boltz_cif,
        ligand_smiles_dict={"LIG": LIGAND_SMILES},
    )
    assert entry.pdb_id == "8c3u"
    assert len(entry.systems) > 0, "Should detect at least one system"

    # Input file on disk must be byte-identical — no side effects
    assert (
        boltz_cif.read_bytes() == before_bytes
    ), "from_custom_cif_file should not mutate the input CIF on disk"
    # And the original CIF should still have no _chem_comp_bond (unknown LIG)
    f = pdbx.CIFFile.read(str(boltz_cif))
    block = list(f.values())[0]
    assert "chem_comp_bond" not in block


def test_from_custom_cif_user_smiles_takes_precedence(boltz_cif):
    """User-supplied SMILES wins over the CCD placeholder for custom residues.

    biotite ships a generic placeholder for the CCD code ``LIG`` — if we
    used it, the ligand's ``smiles`` field would be wrong AND the stereo
    check would silently pass any 3D conformer. This test asserts:
      1. ``lig.smiles`` equals the canonical form of the user SMILES
         (not the CCD placeholder).
      2. With correct stereo, ``resolved_stereo_matches_template`` is True.
      3. With inverted stereo, it flips to False — proving the check
         actually uses the user-provided template.
    """
    import shutil

    from plinder.data.annotations.aggregate_annotations import Entry
    from plinder.data.annotations.ligand_utils import _get_ccd_smiles
    from rdkit import Chem

    # Sanity: the biotite CCD placeholder for "LIG" is a different molecule
    placeholder = _get_ccd_smiles("LIG")
    canonical_user = Chem.MolToSmiles(Chem.MolFromSmiles(LIGAND_SMILES))
    assert (
        placeholder is not None and placeholder != canonical_user
    ), "Expected the biotite LIG placeholder to differ from the user SMILES"

    assert "[C@@]" in LIGAND_SMILES, "YAML SMILES must have the stereo center"
    inverted = LIGAND_SMILES.replace("[C@@]", "[C@]")
    assert LIGAND_SMILES != inverted

    for expected_stereo, smi in [(False, inverted), (True, LIGAND_SMILES)]:
        copy = boltz_cif.parent / f"copy_{expected_stereo}.cif"
        shutil.copy(boltz_cif, copy)
        entry = Entry.from_custom_cif_file(
            pdb_id="8c3u",
            cif_file=copy,
            ligand_smiles_dict={"LIG": smi},
        )
        ligs = [
            l
            for sys in entry.systems.values()
            for l in sys.ligands
            if l.ccd_code == "LIG"
        ]
        assert ligs, "LIG ligand not found in systems"
        for lig in ligs:
            expected_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            assert (
                lig.smiles == expected_canonical
            ), f"lig.smiles should match user SMILES, got {lig.smiles}"
            assert (
                lig.smiles != placeholder
            ), "lig.smiles fell back to CCD placeholder — user SMILES did not win"
            assert lig.resolved_stereo_matches_template is expected_stereo, (
                f"expected stereo_matches={expected_stereo} for "
                f"{'correct' if expected_stereo else 'inverted'} SMILES, "
                f"got {lig.resolved_stereo_matches_template}"
            )


def test_from_custom_cif_save_fixed_roundtrip(boltz_cif, tmp_path):
    """Full round-trip: bad input -> fix -> save -> reload should pass validation.

    1. Input CIF has no _chem_comp_bond (fails check_cif_bond_orders).
    2. from_custom_cif_file with save_fixed_cif writes the enriched CIF.
    3. Reloading the saved file:
       - has _chem_comp_bond
       - passes check_cif_bond_orders
       - produces an equivalent Entry without needing SMILES again
    """
    from plinder.data.annotations.aggregate_annotations import Entry

    # 1. Input is bad — confirm it fails validation
    with pytest.raises(MissingBondOrderError):
        check_cif_bond_orders(boltz_cif)

    fixed_cif = tmp_path / "fixed.cif"
    input_bytes_before = boltz_cif.read_bytes()

    # 2. Fix + save
    entry1 = Entry.from_custom_cif_file(
        pdb_id="8c3u",
        cif_file=boltz_cif,
        ligand_smiles_dict={"LIG": LIGAND_SMILES},
        save_fixed_cif=fixed_cif,
    )
    assert fixed_cif.is_file(), "save_fixed_cif target should be written"
    assert (
        boltz_cif.read_bytes() == input_bytes_before
    ), "Input CIF must remain untouched"

    # 3. Reload the saved fixed CIF and confirm it's self-sufficient
    block = list(pdbx.CIFFile.read(str(fixed_cif)).values())[0]
    assert "chem_comp_bond" in block
    check_cif_bond_orders(fixed_cif)  # must not raise

    entry2 = Entry.from_custom_cif_file(
        pdb_id="8c3u",
        cif_file=fixed_cif,  # no ligand_smiles_dict needed — already enriched
    )
    assert sorted(entry1.systems.keys()) == sorted(
        entry2.systems.keys()
    ), "Systems from the round-tripped fixed CIF must match the original run"


def test_save_fixed_cif_refuses_to_overwrite_input(boltz_cif):
    """save_fixed_cif pointing at the input path must raise, not overwrite."""
    from plinder.data.annotations.aggregate_annotations import Entry

    with pytest.raises(ValueError, match="must not point at the input"):
        Entry.from_custom_cif_file(
            pdb_id="8c3u",
            cif_file=boltz_cif,
            ligand_smiles_dict={"LIG": LIGAND_SMILES},
            save_fixed_cif=boltz_cif,
        )


def test_save_fixed_cif_refuses_to_overwrite_existing(boltz_cif, tmp_path):
    """save_fixed_cif pointing at an existing file must raise, not overwrite."""
    from plinder.data.annotations.aggregate_annotations import Entry

    existing = tmp_path / "existing.cif"
    existing.write_text("DO NOT OVERWRITE ME")

    with pytest.raises(FileExistsError):
        Entry.from_custom_cif_file(
            pdb_id="8c3u",
            cif_file=boltz_cif,
            ligand_smiles_dict={"LIG": LIGAND_SMILES},
            save_fixed_cif=existing,
        )
    assert existing.read_text() == "DO NOT OVERWRITE ME"


# ---------------------------------------------------------------------------
# atoms_to_rdkit_mol unit tests
# ---------------------------------------------------------------------------


def test_atoms_to_rdkit_mol_error():
    """atoms_to_rdkit_mol raises ValueError on empty input."""
    import biotite.structure as struc
    from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

    empty = struc.AtomArray(0)
    try:
        atoms_to_rdkit_mol(empty)
        assert False, "Should have raised ValueError"
    except (ValueError, Exception):
        pass
