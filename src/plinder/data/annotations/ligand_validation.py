# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Validation of a ligand and its own pocket, independent of system grouping."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

from plinder.data.annotations.get_ligand_validation import (
    ResidueListValidation,
    ResidueValidation,
    ResidueValidationThresholds,
)

if TYPE_CHECKING:
    from plinder.data.annotations.protein_utils import Chain


def validation_columns(prefix: str) -> list[str]:
    """Return all summary fields, including fields unavailable for an entry."""
    fields = [
        key
        for key in ResidueListValidation.model_fields
        if key not in {"thresholds", "percent_outliers"}
    ]
    fields.extend(
        f"percent_outliers_{kind}"
        for kind in ("geometry", "density", "chirality", "clashes")
    )
    return [f"{prefix}_validation_{key}" for key in fields]


def format_selection_validation(
    selection: Mapping[str, Sequence[int]],
    chains: Mapping[str, Chain],
    *,
    prefix: str,
    thresholds: ResidueValidationThresholds | None = None,
) -> dict[str, float | int | None]:
    """Aggregate selected residues, never other ligands sharing a system.

    A multi-chain covalent ligand includes all its members. Duplicate residue
    numbers within a chain are removed; distinct assembly instances remain
    distinct. Missing residue validation remains in the coverage denominator.
    """
    residues: list[ResidueValidation | None] = []
    selected_thresholds = thresholds
    for instance_chain, numbers in sorted(selection.items()):
        chain = chains[instance_chain.split(".", maxsplit=1)[-1]]
        if selected_thresholds is None and chain.validation is not None:
            selected_thresholds = chain.validation.thresholds
        residues.extend(
            chain.residues[number].validation for number in sorted(set(numbers))
        )
    summary = ResidueListValidation.from_residues(
        residues, selected_thresholds or ResidueValidationThresholds()
    )
    if summary is None:
        return dict.fromkeys(validation_columns(prefix))
    return {f"{prefix}_{key}": value for key, value in summary.format().items()}
