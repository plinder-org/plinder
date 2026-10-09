# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Compact residue mappings shared by search and score reconstruction."""

from __future__ import annotations

import re
from bisect import bisect_left
from collections.abc import Iterable
from pathlib import Path

import pandas as pd


def cigar_alignment_sql(path: Path, chain_lookup: Path) -> str:
    """DuckDB relation recovering sparse residue positions without a Python UDF."""
    source = str(path).replace("'", "''")
    lookup = str(chain_lookup).replace("'", "''")
    return f"""WITH native AS (
        SELECT a AS original,
            CASE WHEN a.source='mmseqs' THEN q.selected_residue_numbers
                ELSE list_transform(q.selected_residue_indices, i -> i+1) END AS qnative,
            coalesce(CASE WHEN a.source='mmseqs' THEN t.selected_residue_numbers
                ELSE list_transform(t.selected_residue_indices, i -> i+1) END, []) AS tnative,
            regexp_extract_all(a.cigar, '([0-9]+)([MID=X])', 1)::BIGINT[] AS lengths,
            regexp_extract_all(a.cigar, '([0-9]+)([MID=X])', 2) AS operations
        FROM read_parquet('{source}') a
        JOIN read_parquet('{lookup}') q
            ON a.query_entry=q.entry_pdb_id AND a.query_chain_mapped=q.chain_asym_id
        LEFT JOIN read_parquet('{lookup}') t
            ON a.target_entry=t.entry_pdb_id AND a.target_chain_mapped=t.chain_asym_id
    ), consumed AS (
        SELECT *,
            list_transform(list_zip(lengths, operations), r -> CASE WHEN r[2]='D' THEN 0 ELSE r[1] END) AS qlengths,
            list_transform(list_zip(lengths, operations), r -> CASE WHEN r[2]='I' THEN 0 ELSE r[1] END) AS tlengths
        FROM native
    ), starts AS (
        SELECT *,
            list_transform(range(1, len(lengths)+1), i -> original.query_start + coalesce(list_sum(list_slice(qlengths, 1, i-1)), 0)) AS qstarts,
            list_transform(range(1, len(lengths)+1), i -> original.target_start + coalesce(list_sum(list_slice(tlengths, 1, i-1)), 0)) AS tstarts
        FROM consumed
    ), matched AS (
        SELECT *, list_transform(qnative, n -> list_position(
            list_transform(range(1, len(lengths)+1), j ->
                n >= qstarts[j] AND n < qstarts[j]+qlengths[j]
                AND operations[j] IN ('M', '=', 'X')), true)) AS matching_runs
        FROM starts
    ), selected AS (
        SELECT *, list_filter(list_grade_up(qnative), i -> matching_runs[i] IS NOT NULL) AS qpositions
        FROM matched
    ) SELECT original.*,
        qpositions AS query_selected_residue_positions,
        list_transform(qpositions, i -> coalesce(list_position(tnative,
            tstarts[matching_runs[i]] + qnative[i] - qstarts[matching_runs[i]]), 0)) AS target_selected_residue_positions
    FROM selected"""


def decode_cigar_residue_positions(
    alignments: pd.DataFrame, *, chain_lookup: Path
) -> pd.DataFrame:
    """Recover selected residue pairs from CIGARs and the chain lookup.

    Only selected query residues are visited; target position zero denotes a
    residue outside the target's selected pocket/interface residues.
    """
    entries = sorted(set(alignments.query_entry) | set(alignments.target_entry))
    lookup = (
        pd.read_parquet(
            chain_lookup,
            columns=[
                "entry_pdb_id",
                "chain_asym_id",
                "selected_residue_numbers",
                "selected_residue_indices",
            ],
            filters=[("entry_pdb_id", "in", entries)],
        )
        if entries
        else pd.DataFrame()
    )
    maps: dict[tuple[str, str, str], tuple[list[int], dict[int, int]]] = {}
    for row in lookup.itertuples(index=False):
        for backend, numbers in (
            ("mmseqs", [int(n) for n in row.selected_residue_numbers]),
            ("foldseek", [int(i) + 1 for i in row.selected_residue_indices]),
        ):
            maps[(row.entry_pdb_id, row.chain_asym_id, backend)] = (
                sorted(numbers),
                {n: i for i, n in enumerate(numbers, start=1)},
            )
    queries, targets = [], []
    for row in alignments.itertuples(index=False):
        selected, qmap = maps[(row.query_entry, row.query_chain_mapped, row.source)]
        _, tmap = maps.get(
            (row.target_entry, row.target_chain_mapped, row.source), ([], {})
        )
        q, t = int(row.query_start), int(row.target_start)
        query, target = [], []
        for length_text, operation in re.findall(r"(\d+)([MID=X])", row.cigar):
            length = int(length_text)
            if operation in {"M", "=", "X"}:
                lo, hi = bisect_left(selected, q), bisect_left(selected, q + length)
                for number in selected[lo:hi]:
                    query.append(qmap[number])
                    target.append(tmap.get(t + number - q, 0))
            if operation in {"M", "I", "=", "X"}:
                q += length
            if operation in {"M", "D", "=", "X"}:
                t += length
        queries.append(query)
        targets.append(target)
    result = alignments.copy()
    result["query_selected_residue_positions"] = queries
    result["target_selected_residue_positions"] = targets
    return result


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
