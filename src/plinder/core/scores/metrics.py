# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Similarity metric definitions shared by scoring and clustering."""

from __future__ import annotations

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
    "pli_qcov",
    "pli_unique_qcov",
    "shape",
    "color",
    "sucos_shape",
    "sucos_shape_pocket_qcov",
)

SCORE_NAMES = PROTEIN_SCORE_NAMES + LIGAND_SCORE_NAMES

GATED_LIGAND_DIAGNOSTIC_METRICS = frozenset({"shape", "color", "sucos_shape"})

PROTEIN_CLUSTER_METRICS = tuple(
    f"{score_name}{suffix}"
    for score_name in PROTEIN_SCORE_NAMES
    for suffix in ("_max", "_weighted_max", "_weighted_sum")
)

# Shape, color, and raw SuCOS remain available as diagnostic scores. They are
# evaluated only after the positive-pocket gate, so missing values do not mean
# dissimilar ligands and must not be interpreted as absent clustering edges.
LIGAND_CLUSTER_METRICS = tuple(
    metric
    for metric in LIGAND_SCORE_NAMES
    if metric not in GATED_LIGAND_DIAGNOSTIC_METRICS
)

CHEMICAL_CLUSTER_METRICS = ("tanimoto_similarity_ecfp4_1024",)

DEFAULT_CLUSTER_METRICS = (
    PROTEIN_CLUSTER_METRICS + LIGAND_CLUSTER_METRICS + CHEMICAL_CLUSTER_METRICS
)


def is_ligand_level_metric(metric: str) -> bool:
    """Return whether graph nodes for this score are individual ligands."""
    return metric in LIGAND_SCORE_NAMES or metric.startswith(("pocket_", "pli_"))
