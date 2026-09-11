# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
import pytest
from biotite.structure.io import pdbx

from plinder.core import PlinderInterface, PlinderRelease
from plinder.data.annotations.cif_utils import read_mmcif_file


@pytest.fixture
def interface_release(tmp_path) -> PlinderRelease:
    index = tmp_path / "index"
    index.mkdir()
    system_id = "2y4i__1__1.A--1.B"
    pd.DataFrame({
        "entry_pdb_id": ["2y4i"],
        "system_id": [system_id],
        "system_biounit_id": ["1"],
        "interface_chain_1": ["1.A"],
        "interface_chain_2": ["1.B"],
        "interface_chain_1_residue_numbers": [[22, 23]],
        "interface_chain_1_residue_indices": [[0, 1]],
        "interface_chain_2_residue_numbers": [[39]],
        "interface_chain_2_residue_indices": [[0]],
        "interface_num_contact_residue_pairs": [2],
    }).to_parquet(index / "interface_annotation_table.parquet", index=False)
    pd.DataFrame({
        "entry_pdb_id": ["2y4i", "2y4i"],
        "chain_asym_id": ["A", "B"],
        "chain_sequence": ["AC", "DEF"],
    }).to_parquet(index / "entry_chains.parquet", index=False)
    return PlinderRelease(data_dir=tmp_path)


def test_interface_loads_annotation_sequences_and_residue_views(
    interface_release, cif_2y4i, tmp_path
):
    interface = PlinderInterface(
        system_id="2y4i__1__1.A--1.B",
        release=interface_release,
        source_mmcif=cif_2y4i,
        reconstruction_dir=tmp_path / "reconstructed",
    )

    assert interface.entry_pdb_id == "2y4i"
    assert interface.biounit_id == "1"
    assert interface.chains == ("1.A", "1.B")
    assert interface.sequences == {"1.A": "AC", "1.B": "DEF"}
    assert set(interface.chain_structures) == {"1.A", "1.B"}
    assert all(len(atoms) > 0 for atoms in interface.chain_structures.values())
    assert set(interface.interface_residue_masks) == {"1.A", "1.B"}
    assert len(interface.interface_structure) > 0
    assert set(interface.interface_structure.chain_id) == {"1.A", "1.B"}

    output = interface.reconstruct()
    assert output == tmp_path / "reconstructed" / "interface.cif"
    cif_file = read_mmcif_file(output)
    block = cif_file[interface.system_id]
    assert {
        "entry",
        "entity",
        "entity_poly",
        "entity_poly_seq",
        "chem_comp",
        "struct_asym",
        "atom_site",
    }.issubset(block)
    atoms = pdbx.get_structure(
        cif_file,
        model=1,
        use_author_fields=False,
        include_bonds=True,
    )
    assert len(set(atoms.chain_id.astype(str))) == 2
    assert set(block["atom_site"]["label_asym_id"].as_array(str)).issubset(
        block["struct_asym"]["id"].as_array(str)
    )
    assert set(block["atom_site"]["label_entity_id"].as_array(str)).issubset(
        block["entity"]["id"].as_array(str)
    )
    assert interface.interface_cif == output
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        interface.reconstruct()


def test_interface_uses_explicit_release_for_default_reconstruction_dir(
    interface_release, cif_2y4i
):
    interface = PlinderInterface(
        system_id="2y4i__1__1.A--1.B",
        release=interface_release,
        source_mmcif=cif_2y4i,
    )

    assert interface.reconstruction_dir == (
        interface_release.data_dir / "reconstructed_interfaces" / "2y4i__1__1.A--1.B"
    )


def test_interface_rejects_source_with_different_residue_numbering(
    interface_release, cif_2y4i, tmp_path
):
    interface = PlinderInterface(
        system_id="2y4i__1__1.A--1.B",
        release=interface_release,
        source_mmcif=cif_2y4i,
        reconstruction_dir=tmp_path / "reconstructed",
    )
    annotation = interface.annotation.copy()
    annotation["interface_chain_1_residue_numbers"] = [999, 23]
    interface._annotation = annotation

    with pytest.raises(ValueError, match="do not match the release annotation"):
        _ = interface.interface_residue_masks


def test_interface_reports_an_unknown_id(interface_release):
    interface = PlinderInterface(
        system_id="2y4i__1__1.A--1.Z",
        release=interface_release,
    )
    with pytest.raises(ValueError, match="is not in the release"):
        _ = interface.annotation
