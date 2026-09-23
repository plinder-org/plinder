# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Compact residue mappings shared by search and score reconstruction."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd


def pack_residue_identities(values: Iterable[bool]) -> bytes:
    """Pack residue identity flags least-significant bit first."""
    result = bytearray()
    for index, value in enumerate(values):
        byte_index = index // 8
        if byte_index == len(result):
            result.append(0)
        if value:
            result[byte_index] |= 1 << (index % 8)
    return bytes(result)


def unpack_residue_identities(value: object, count: int) -> list[bool]:
    """Expand ``count`` flags from a packed identity bitset."""
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError("residue identity bitset must be bytes-like")
    packed = bytes(value)
    if len(packed) * 8 < count:
        raise ValueError(
            f"residue identity bitset has {len(packed) * 8} bits for {count} residues"
        )
    return [bool(packed[index // 8] & (1 << (index % 8))) for index in range(count)]


def expand_residue_positions(
    alignments: pd.DataFrame,
    *,
    chain_lookup: Path,
) -> pd.DataFrame:
    """Expand compact chain-relative positions to label-sequence numbers."""
    if alignments.empty:
        result = alignments.copy()
        result["query_selected_residue_numbers"] = pd.Series(dtype="object")
        result["target_selected_residue_numbers"] = pd.Series(dtype="object")
        return result.drop(
            columns=[
                "query_selected_residue_positions",
                "target_selected_residue_positions",
            ],
            errors="ignore",
        )

    entry_ids = sorted(
        set(alignments["query_entry"].astype(str))
        | set(alignments["target_entry"].astype(str))
    )
    lookup = pd.read_parquet(
        chain_lookup,
        columns=["entry_pdb_id", "chain_asym_id", "selected_residue_numbers"],
        filters=[("entry_pdb_id", "in", entry_ids)],
    )
    numbers = {
        (str(row.entry_pdb_id), str(row.chain_asym_id)): [
            int(value) for value in row.selected_residue_numbers
        ]
        for row in lookup.itertuples(index=False)
    }

    query_numbers: list[list[int]] = []
    target_numbers: list[list[int]] = []
    for row in alignments.itertuples(index=False):
        query = numbers[(str(row.query_entry), str(row.query_chain_mapped))]
        target = numbers.get((str(row.target_entry), str(row.target_chain_mapped)), [])
        query_numbers.append(
            [
                query[int(position) - 1]
                for position in row.query_selected_residue_positions
            ]
        )
        target_numbers.append(
            [
                -1 if int(position) == 0 else target[int(position) - 1]
                for position in row.target_selected_residue_positions
            ]
        )

    result = alignments.copy()
    result["query_selected_residue_numbers"] = query_numbers
    result["target_selected_residue_numbers"] = target_numbers
    return result.drop(
        columns=[
            "query_selected_residue_positions",
            "target_selected_residue_positions",
        ]
    )
