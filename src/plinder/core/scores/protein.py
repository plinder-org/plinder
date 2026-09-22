# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Query compact chain-level protein similarities."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.scores.mapping import unpack_residue_identities
from plinder.core.scores.query import Filters, read_score_table
from plinder.core.utils.dec import timeit

OVERLAP_COLUMNS = [
    "query_entry",
    "query_chain",
    "target_entry",
    "target_chain",
    "source",
    "kind",
    "overlapping_residues",
    "identical_overlapping_residues",
    "query_total_residues",
    "target_total_residues",
    "query_overlap_fraction",
    "target_overlap_fraction",
]


@timeit
def query_protein_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query chain-level Foldseek and MMseqs similarities in integer percent."""
    dataset = (release or PlinderRelease()).fetch("protein_similarity_scores")
    return read_score_table(dataset, columns=columns, filters=filters)


def _selected_residue_positions(
    lookup: Path,
    *,
    entry: str,
    chain: str,
) -> dict[int, int]:
    rows = pd.read_parquet(
        lookup,
        columns=["selected_residue_numbers"],
        filters=[("entry_pdb_id", "==", entry), ("chain_asym_id", "==", chain)],
    )
    if len(rows) != 1:
        raise KeyError(f"unknown protein chain {entry}_{chain}")
    return {
        int(number): position
        for position, number in enumerate(
            rows.iloc[0]["selected_residue_numbers"], start=1
        )
    }


def _pocket_residue_numbers(
    release: PlinderRelease,
    *,
    entry: str,
    chain: str,
) -> set[int]:
    rows = pd.read_parquet(
        release.fetch("ligand_pocket_residues"),
        columns=["residue_label_seq_id"],
        filters=[("entry_pdb_id", "==", entry), ("chain_asym_id", "==", chain)],
    )
    return set(rows["residue_label_seq_id"].astype(int))


def _interface_residue_numbers(
    release: PlinderRelease,
    *,
    entry: str,
    chain: str,
) -> set[int]:
    rows = pd.read_parquet(
        release.fetch("interface_annotations"),
        columns=[
            "interface_chain_1",
            "interface_chain_1_residue_numbers",
            "interface_chain_2",
            "interface_chain_2_residue_numbers",
        ],
        filters=[("entry_pdb_id", "==", entry)],
    )
    numbers: set[int] = set()
    for row in rows.itertuples(index=False):
        for instance_chain, residues in (
            (row.interface_chain_1, row.interface_chain_1_residue_numbers),
            (row.interface_chain_2, row.interface_chain_2_residue_numbers),
        ):
            if str(instance_chain).split(".", maxsplit=1)[-1] == chain:
                numbers.update(int(number) for number in residues)
    return numbers


def _residue_positions(
    release: PlinderRelease,
    lookup: Path,
    *,
    entry: str,
    chain: str,
    kind: Literal["pocket", "interface"],
) -> set[int]:
    selected = _selected_residue_positions(lookup, entry=entry, chain=chain)
    if kind == "pocket":
        numbers = _pocket_residue_numbers(release, entry=entry, chain=chain)
    else:
        numbers = _interface_residue_numbers(release, entry=entry, chain=chain)
    return {selected[number] for number in numbers}


def _overlap_fraction(overlap: int, total: int) -> float:
    return overlap / total if total else 0.0


@timeit
def query_chain_overlap(
    query_entry: str,
    query_chain: str,
    target_entry: str,
    target_chain: str,
    *,
    kind: Literal["pocket", "interface"],
    source: Literal["foldseek", "mmseqs"] | None = None,
    search_db: str = "holo",
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Count aligned pocket or interface residues between two protein chains.

    Chain IDs are label asym IDs. Residues are collected across all ligand
    pockets or all protein interfaces involving each chain. One row is returned
    per available search backend.
    """
    if kind not in {"pocket", "interface"}:
        raise ValueError("kind must be 'pocket' or 'interface'")
    selected_release = release or PlinderRelease()
    lookup = selected_release.fetch("alignment_chain_lookup")
    query_positions = _residue_positions(
        selected_release,
        lookup,
        entry=query_entry,
        chain=query_chain,
        kind=kind,
    )
    target_positions = _residue_positions(
        selected_release,
        lookup,
        entry=target_entry,
        chain=target_chain,
        kind=kind,
    )
    query_total = len(query_positions)
    target_total = len(target_positions)

    sources = (source,) if source is not None else ("foldseek", "mmseqs")
    records: list[dict[str, object]] = []
    for alignment_type in sources:
        try:
            shard = selected_release.fetch(
                "alignment_shard",
                search_db=search_db,
                alignment_type=alignment_type,
                shard=query_entry[-3:-1],
            )
        except FileNotFoundError:
            if source is not None:
                raise
            continue
        alignments = pd.read_parquet(
            shard,
            columns=[
                "query_selected_residue_positions",
                "target_selected_residue_positions",
                "selected_residue_identity_bits",
            ],
            filters=[
                ("query_entry", "==", query_entry),
                ("query_chain_mapped", "==", query_chain),
                ("target_entry", "==", target_entry),
                ("target_chain_mapped", "==", target_chain),
            ],
        )
        for row in alignments.itertuples(index=False):
            aligned_query_positions = [
                int(position) for position in row.query_selected_residue_positions
            ]
            aligned_target_positions = [
                int(position) for position in row.target_selected_residue_positions
            ]
            identities = unpack_residue_identities(
                row.selected_residue_identity_bits,
                len(aligned_target_positions),
            )
            overlap_mask = [
                query_position in query_positions
                and target_position in target_positions
                for query_position, target_position in zip(
                    aligned_query_positions,
                    aligned_target_positions,
                    strict=True,
                )
            ]
            overlap = sum(overlap_mask)
            records.append(
                {
                    "query_entry": query_entry,
                    "query_chain": query_chain,
                    "target_entry": target_entry,
                    "target_chain": target_chain,
                    "source": alignment_type,
                    "kind": kind,
                    "overlapping_residues": overlap,
                    "identical_overlapping_residues": sum(
                        identity and overlaps
                        for identity, overlaps in zip(
                            identities, overlap_mask, strict=True
                        )
                    ),
                    "query_total_residues": query_total,
                    "target_total_residues": target_total,
                    "query_overlap_fraction": _overlap_fraction(overlap, query_total),
                    "target_overlap_fraction": _overlap_fraction(overlap, target_total),
                }
            )
    return pd.DataFrame.from_records(records, columns=OVERLAP_COLUMNS)


__all__ = ["query_chain_overlap", "query_protein_similarity"]
