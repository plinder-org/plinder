# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Load full chain alignments from optional CIGAR shards."""

from __future__ import annotations

import re
from typing import Literal

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.utils.dec import timeit

POSITION_COLUMNS = [
    "alignment_index",
    "alignment_column",
    "operation",
    "query_seqres_position",
    "target_seqres_position",
    "query_coordinate_index",
    "target_coordinate_index",
]
_CIGAR_RUN = re.compile(r"([1-9][0-9]*)([MID=X])")


def expand_alignment_cigar(
    cigar: str,
    *,
    query_start: int,
    target_start: int,
    source: Literal["foldseek", "mmseqs"],
    include_gaps: bool = True,
    alignment_index: int = 0,
) -> pd.DataFrame:
    """Expand one CIGAR into backend-native residue positions.

    MMseqs positions are one-based SEQRES positions. Foldseek positions are
    zero-based indices over resolved coordinate residues.
    """
    if source not in {"foldseek", "mmseqs"}:
        raise ValueError("source must be 'foldseek' or 'mmseqs'")
    if query_start < 1 or target_start < 1:
        raise ValueError("query_start and target_start must be positive")
    runs = _CIGAR_RUN.findall(cigar)
    if (
        not runs
        or "".join(f"{length}{operation}" for length, operation in runs) != cigar
    ):
        raise ValueError(f"invalid alignment CIGAR: {cigar!r}")

    query_position = query_start
    target_position = target_start
    records: list[dict[str, object]] = []
    alignment_column = 0
    for length_text, operation in runs:
        query_consumed = operation in {"M", "I", "=", "X"}
        target_consumed = operation in {"M", "D", "=", "X"}
        for _ in range(int(length_text)):
            query_value = query_position if query_consumed else None
            target_value = target_position if target_consumed else None
            if include_gaps or (query_value is not None and target_value is not None):
                records.append(
                    {
                        "alignment_index": alignment_index,
                        "alignment_column": alignment_column,
                        "operation": operation,
                        "query_seqres_position": (
                            query_value if source == "mmseqs" else None
                        ),
                        "target_seqres_position": (
                            target_value if source == "mmseqs" else None
                        ),
                        "query_coordinate_index": (
                            query_value - 1
                            if source == "foldseek" and query_value is not None
                            else None
                        ),
                        "target_coordinate_index": (
                            target_value - 1
                            if source == "foldseek" and target_value is not None
                            else None
                        ),
                    }
                )
            query_position += int(query_consumed)
            target_position += int(target_consumed)
            alignment_column += 1
    result = pd.DataFrame.from_records(records, columns=POSITION_COLUMNS)
    for column in POSITION_COLUMNS[3:]:
        result[column] = result[column].astype("Int64")
    return result


@timeit
def map_chain_alignment(
    query_entry: str,
    query_chain: str,
    target_entry: str,
    target_chain: str,
    *,
    source: Literal["foldseek", "mmseqs"],
    search_db: str = "holo",
    include_gaps: bool = True,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Load and expand full alignments for a pair of label-asym chains."""
    selected_release = release or PlinderRelease()
    shard = selected_release.fetch(
        "alignment_cigar_shard",
        search_db=search_db,
        alignment_type=source,
        shard=query_entry[-3:-1],
    )
    alignments = pd.read_parquet(
        shard,
        columns=["query_start", "target_start", "cigar"],
        filters=[
            ("query_entry", "==", query_entry),
            ("query_chain_mapped", "==", query_chain),
            ("target_entry", "==", target_entry),
            ("target_chain_mapped", "==", target_chain),
        ],
    )
    frames = [
        expand_alignment_cigar(
            str(row.cigar),
            query_start=int(row.query_start),
            target_start=int(row.target_start),
            source=source,
            include_gaps=include_gaps,
            alignment_index=index,
        )
        for index, row in enumerate(alignments.itertuples(index=False))
    ]
    if not frames:
        return pd.DataFrame(columns=POSITION_COLUMNS)
    return pd.concat(frames, ignore_index=True)


__all__ = ["expand_alignment_cigar", "map_chain_alignment"]
