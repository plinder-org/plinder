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
NON_CLUSTERING_LIGAND_METRICS = GATED_LIGAND_DIAGNOSTIC_METRICS | {
    "pocket_fident",
    "pocket_fident_qcov",
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

CHEMICAL_CLUSTER_METRICS = ("tanimoto_similarity_ecfp4_1024",)

DEFAULT_CLUSTER_METRICS = LIGAND_CLUSTER_METRICS + CHEMICAL_CLUSTER_METRICS


def is_ligand_level_metric(metric: str) -> bool:
    """Return whether graph nodes for this score are individual ligands."""
    return (
        metric in LIGAND_SCORE_NAMES
        or metric in PROTEIN_CLUSTER_METRICS
        or metric in CHEMICAL_CLUSTER_METRICS
        or metric.startswith(("pocket_", "pli_"))
    )
