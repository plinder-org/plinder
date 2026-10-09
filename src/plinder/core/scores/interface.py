# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Query complete protein-interface similarities."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.scores.entries import InterfaceView
from plinder.core.scores.mapping import unpack_residue_identities
from plinder.core.scores.query import Filters, read_score_table
from plinder.core.utils.dec import timeit

INTERFACE_SIDE_COLUMNS = [
    "custom_structure_id",
    "custom_chain",
    "plinder_interface_id",
    "plinder_chain",
    "plinder_partner_chain",
    "source",
    "interface_side_qcov",
    "interface_side_fident",
    "plinder_residue_numbers",
    "custom_residue_numbers",
]


def calculate_interface_side_scores(
    alignments: pd.DataFrame,
    *,
    interfaces: Mapping[str, InterfaceView],
) -> pd.DataFrame:
    """Map known interface sides onto custom chains, without requiring a partner.

    Alignments run from the release chain to the custom chain. Coverage and
    identity are integer percentages of *all* residues in the known interface
    side, so gaps and unaligned residues contribute zero. Each backend retains
    its best single alignment per interface side and custom chain; alignments
    are never combined to inflate coverage.
    """
    sides = defaultdict(list)
    for interface in interfaces.values():
        for chain, partner, residues in (
            (
                interface.chain_1,
                interface.chain_2,
                interface.chain_1_residue_number_to_index,
            ),
            (
                interface.chain_2,
                interface.chain_1,
                interface.chain_2_residue_number_to_index,
            ),
        ):
            if residues:
                sides[(interface.pdb_id, chain.split(".", 1)[-1])].append(
                    (interface.id, chain, partner, frozenset(residues))
                )
    best = {}
    for row in alignments.itertuples(index=False):
        matching = sides.get((str(row.query_entry), str(row.query_chain_mapped)), [])
        if not matching:
            continue
        release_numbers = list(row.query_selected_residue_numbers)
        custom_numbers = list(row.target_selected_residue_numbers)
        identities = unpack_residue_identities(
            row.selected_residue_identity_bits, len(release_numbers)
        )
        pairs = {
            int(release_number): (int(custom_number), identical)
            for release_number, custom_number, identical in zip(
                release_numbers, custom_numbers, identities, strict=True
            )
            if int(custom_number) > 0
        }
        for interface_id, chain, partner, residues in matching:
            numbers = sorted(residues.intersection(pairs))
            if not numbers:
                continue
            identical_count = sum(pairs[number][1] for number in numbers)
            key = (
                str(row.target_entry),
                str(row.target_chain_mapped),
                interface_id,
                chain,
                str(row.source),
            )
            rank = (len(numbers), identical_count)
            record = dict(zip(INTERFACE_SIDE_COLUMNS[:4], key[:4], strict=True))
            record.update(
                plinder_partner_chain=partner,
                source=key[4],
                interface_side_qcov=int(100 * len(numbers) / len(residues) + 0.5),
                interface_side_fident=int(100 * identical_count / len(residues) + 0.5),
                plinder_residue_numbers=numbers,
                custom_residue_numbers=[pairs[number][0] for number in numbers],
            )
            if key not in best or rank > best[key][0]:
                best[key] = (rank, record)
    result = pd.DataFrame(
        [record for _, record in best.values()], columns=INTERFACE_SIDE_COLUMNS
    )
    return result.sort_values(
        [
            "interface_side_qcov",
            "interface_side_fident",
            *INTERFACE_SIDE_COLUMNS[:4],
            "source",
        ],
        ascending=[False, False, True, True, True, True, True],
        ignore_index=True,
    )


@timeit
def query_interface_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query complete directed protein-interface similarities."""
    dataset = (release or PlinderRelease()).fetch("interface_similarity_scores")
    return read_score_table(dataset, columns=columns, filters=filters)


@timeit
def query_half_interface_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query directed coverage between individual interface sides."""
    dataset = (release or PlinderRelease()).fetch("interface_half_similarity_scores")
    return read_score_table(dataset, columns=columns, filters=filters)


__all__ = ["query_interface_similarity", "query_half_interface_similarity"]
