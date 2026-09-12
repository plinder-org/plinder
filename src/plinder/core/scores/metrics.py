# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Similarity metric definitions shared by scoring and clustering."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping

PROTEIN_SCORE_NAMES = (
    "protein_lddt",
    "protein_lddt_qcov",
    "protein_qcov",
    "protein_fident",
    "protein_fident_qcov",
    "protein_seqsim",
    "protein_seqsim_qcov",
)

LIGAND_SCORE_NAMES = (
    "pocket_qcov",
    "pocket_fident",
    "pocket_fident_qcov",
    "pli_fident",
    "pli_qcov",
    "pli_unique_qcov",
    "shape",
    "color",
    "sucos_shape",
    "sucos_shape_pocket_qcov",
)

SCORE_NAMES = PROTEIN_SCORE_NAMES + LIGAND_SCORE_NAMES

GATED_LIGAND_DIAGNOSTIC_METRICS = frozenset({"shape", "color", "sucos_shape"})
NON_CLUSTERING_LIGAND_METRICS = GATED_LIGAND_DIAGNOSTIC_METRICS | {
    "pocket_fident",
    "pocket_fident_qcov",
    "pli_fident",
}

FOLDSEEK_ONLY_PROTEIN_SCORE_NAMES = frozenset({"protein_lddt", "protein_lddt_qcov"})

PROTEIN_CLUSTER_METRICS = tuple(
    f"{score_name}{suffix}"
    for score_name in PROTEIN_SCORE_NAMES
    if score_name not in FOLDSEEK_ONLY_PROTEIN_SCORE_NAMES
    for suffix in ("_max", "_weighted_max", "_weighted_sum")
)

# Shape, color, raw SuCOS, and pocket fident remain queryable scores but are not
# release clustering metrics. The 3D diagnostics are evaluated only after the
# positive-pocket gate.
LIGAND_CLUSTER_METRICS = tuple(
    metric
    for metric in LIGAND_SCORE_NAMES
    if metric not in NON_CLUSTERING_LIGAND_METRICS
)

# Fingerprint-derived ligand-similarity metrics. Each is computed on the shared
# unique-SMILES node universe (``ligand_smiles_id``) but from its own fingerprint
# and similarity measure: ECFP4/1024 Tanimoto and MHFP6/2048 estimated Jaccard.
CHEMICAL_CLUSTER_METRICS = (
    "tanimoto_similarity_ecfp4_1024",
    "jaccard_similarity_mhfp6_2048",
)

# Each chemical metric also publishes one convenience cluster column from its
# 90-percent reciprocal set cover, plus the distinct-PDB count of that cluster
# under ``<column>_num_pdb_ids``.
CHEMICAL_CLUSTER_SUMMARY_THRESHOLD = 90
CHEMICAL_CLUSTER_SUMMARY_COLUMNS = {
    "tanimoto_similarity_ecfp4_1024": "ligand_tanimoto_ecfp4_1024_90_cluster",
    "jaccard_similarity_mhfp6_2048": "ligand_jaccard_mhfp6_2048_90_cluster",
}

DEFAULT_CLUSTER_METRICS = LIGAND_CLUSTER_METRICS + CHEMICAL_CLUSTER_METRICS


def maximum_weight_bipartite_assignment(
    query_nodes: Iterable[str],
    target_nodes: Iterable[str],
    primary_weights: Mapping[tuple[str, str], float],
    *,
    secondary_weights: Mapping[tuple[str, str], float] | None = None,
) -> list[tuple[str, str]]:
    """Maximize integer-stepped primary weights, then secondary weights.

    Zero-weight dummy nodes allow either side to remain unmatched. Secondary
    weights are normalized so that their complete assignment can never
    outweigh one primary-weight unit. Returned pairs are deterministic and
    ordered by their contribution to the objective.
    """
    query = sorted(set(map(str, query_nodes)))
    target = sorted(set(map(str, target_nodes)))
    size = max(len(query), len(target))
    if size == 0:
        return []

    secondary_weights = secondary_weights or {}
    finite_secondary = [
        float(value)
        for value in secondary_weights.values()
        if math.isfinite(float(value)) and float(value) > 0
    ]
    secondary_scale = max([1.0, *finite_secondary])
    primary_multiplier = float(size + 1)
    weights = [[0.0] * size for _ in range(size)]
    components: dict[tuple[int, int], tuple[float, float]] = {}
    maximum = 0.0
    for query_index, query_node in enumerate(query):
        for target_index, target_node in enumerate(target):
            pair = (query_node, target_node)
            primary = float(primary_weights.get(pair, 0.0))
            secondary = float(secondary_weights.get(pair, 0.0))
            if not math.isfinite(primary) or primary < 0:
                primary = 0.0
            if not math.isfinite(secondary) or secondary < 0:
                secondary = 0.0
            secondary /= secondary_scale
            weight = primary * primary_multiplier + secondary
            weights[query_index][target_index] = weight
            components[(query_index, target_index)] = (primary, secondary)
            maximum = max(maximum, weight)
    if maximum == 0:
        return []

    # Shortest-augmenting-path Hungarian algorithm for a square minimum-cost
    # problem. Converting max weights to costs preserves the optimum.
    costs = [[maximum - weight for weight in row] for row in weights]
    row_potential = [0.0] * (size + 1)
    column_potential = [0.0] * (size + 1)
    column_match = [0] * (size + 1)
    predecessor = [0] * (size + 1)
    for row in range(1, size + 1):
        column_match[0] = row
        minimum = [math.inf] * (size + 1)
        used = [False] * (size + 1)
        column = 0
        while True:
            used[column] = True
            matched_row = column_match[column]
            delta = math.inf
            next_column = 0
            for candidate in range(1, size + 1):
                if used[candidate]:
                    continue
                reduced_cost = (
                    costs[matched_row - 1][candidate - 1]
                    - row_potential[matched_row]
                    - column_potential[candidate]
                )
                if reduced_cost < minimum[candidate]:
                    minimum[candidate] = reduced_cost
                    predecessor[candidate] = column
                if minimum[candidate] < delta:
                    delta = minimum[candidate]
                    next_column = candidate
            for candidate in range(size + 1):
                if used[candidate]:
                    row_potential[column_match[candidate]] += delta
                    column_potential[candidate] -= delta
                else:
                    minimum[candidate] -= delta
            column = next_column
            if column_match[column] == 0:
                break
        while True:
            previous = predecessor[column]
            column_match[column] = column_match[previous]
            column = previous
            if column == 0:
                break

    row_to_column = [-1] * size
    for column in range(1, size + 1):
        if column_match[column] > 0:
            row_to_column[column_match[column] - 1] = column - 1
    weighted_pairs = [
        (
            *components[(row, column)],
            query[row],
            target[column],
        )
        for row, column in enumerate(row_to_column[: len(query)])
        if 0 <= column < len(target) and weights[row][column] > 0
    ]
    weighted_pairs.sort(key=lambda pair: (-pair[0], -pair[1], pair[2], pair[3]))
    return [
        (query_node, target_node) for _, _, query_node, target_node in weighted_pairs
    ]


def is_chemical_cluster_metric(metric: str) -> bool:
    """Return whether a metric is a fingerprint-derived ligand-similarity metric.

    Chemical metrics cluster the unique-SMILES node universe and are symmetric,
    which sets them apart from the pocket/protein score metrics throughout the
    clustering pipeline.
    """
    return metric in CHEMICAL_CLUSTER_METRICS


def is_ligand_level_metric(metric: str) -> bool:
    """Return whether graph nodes for this score are individual ligands."""
    return (
        metric in LIGAND_SCORE_NAMES
        or metric in PROTEIN_CLUSTER_METRICS
        or metric in CHEMICAL_CLUSTER_METRICS
        or metric.startswith(("pocket_", "pli_"))
    )
