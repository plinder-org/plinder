# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pyarrow as pa
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import read_mmcif_file
from plinder.data.annotations.protein_utils import get_modified_residues
from plinder.data.pipeline.collate import ENTRY_CHAIN_SCHEMA

CIF_6M92 = "xx/pdb_00006m92/pdb_00006m92_xyz-enrich.cif.gz"


def test_modified_residues_come_from_seqres_with_deposited_parent(test_dir):
    block = list(read_mmcif_file(test_dir / CIF_6M92).values())[0]

    assert get_modified_residues(block) == {"C": ["33:SEP:C:18>SER (modified residue)"]}


def test_modified_residues_fall_back_to_the_ccd_parent(test_dir):
    block = list(read_mmcif_file(test_dir / CIF_6M92).values())[0]
    assert "pdbx_struct_mod_residue" in block
    del block["pdbx_struct_mod_residue"]

    assert get_modified_residues(block) == {"C": ["33:SEP:C:18>SER"]}


def test_entry_chain_table_carries_modifications_and_noncanonical_sequence(test_dir):
    entry = Entry.from_cif_file(test_dir / CIF_6M92, include_interfaces=False)

    assert entry.chains["C"].modified_residues == ["33:SEP:C:18>SER (modified residue)"]
    assert entry.chains["A"].modified_residues == []
    assert entry.chain_to_seqres["C"][17] == "S"
    assert "(SEP)" in entry.chain_to_seqres_noncanonical["C"]
    assert "(" not in entry.chain_to_seqres_noncanonical["A"]

    frame = entry.chains_to_df()
    table = pa.Table.from_pandas(frame, schema=ENTRY_CHAIN_SCHEMA, preserve_index=False)
    rows = {row["chain_asym_id"]: row for row in table.to_pylist()}
    assert rows["C"]["chain_modified_residues"] == [
        "33:SEP:C:18>SER (modified residue)"
    ]
    assert (
        rows["C"]["chain_sequence_noncanonical"]
        == entry.chain_to_seqres_noncanonical["C"]
    )
    assert rows["A"]["chain_modified_residues"] == []
