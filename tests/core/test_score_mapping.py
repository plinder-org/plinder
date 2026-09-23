# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path

import pandas as pd

from plinder.core.scores.mapping import (
    expand_residue_positions,
    pack_residue_identities,
    unpack_residue_identities,
)


def test_residue_identity_bitset_roundtrip() -> None:
    values = [True, False, True, True, False, False, True, False, True]

    packed = pack_residue_identities(values)

    assert packed == bytes([0b01001101, 0b00000001])
    assert unpack_residue_identities(packed, len(values)) == values


def test_expand_chain_relative_residue_positions(tmp_path: Path) -> None:
    lookup = tmp_path / "alignment_chain_lookup.parquet"
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "chain_asym_id": ["A", "B"],
            "selected_residue_numbers": [[10, 30], [101, 205]],
        }
    ).to_parquet(lookup, index=False)
    compact = pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "target_entry": ["2def"],
            "query_chain_mapped": ["A"],
            "target_chain_mapped": ["B"],
            "query_selected_residue_positions": [[2, 1]],
            "target_selected_residue_positions": [[0, 2]],
        }
    )

    expanded = expand_residue_positions(compact, chain_lookup=lookup)

    assert expanded.loc[0, "query_selected_residue_numbers"] == [30, 10]
    assert expanded.loc[0, "target_selected_residue_numbers"] == [-1, 205]
