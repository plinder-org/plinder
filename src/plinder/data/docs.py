# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd

from plinder.data import column_descriptions
from plinder.data.utils.annotations.aggregate_annotations import Entry, System
from plinder.data.utils.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
)
from plinder.data.utils.annotations.ligand_utils import Ligand
from plinder.data.utils.annotations.protein_utils import Chain

TSV_DIR = Path(column_descriptions.__file__).parent
VALIDATION_TYPES = [
    "system_pocket",
    "system_ligand",
    "ligand_interacting_ligand_chains",
    "ligand_neighboring_ligand_chains",
    "ligand_protein_chains",
    "system_ligand_chains",
    "system_protein_chains",
]
VALIDATION_OUTLIER_KEYS = [
    "chirality",
    "clashes",
    "density",
    "geometry",
]
MAPPING_NAMES = [
    "CATH",
    "ECOD",
    "ECOD_t_name",
    "Pfam",
    "SCOP2",
    "SCOP2B",
    "PANTHER",
    "UniProt",
    "kinase_name",
]
CHAIN_TYPES = [
    "system_protein_chains",
    "system_ligand_chains",
    "ligand_interacting_ligand_chains",
    "ligand_neighboring_ligand_chains",
    "ligand_protein_chains",
]


def get_cluster_column_descriptions(
    plindex: pd.DataFrame
) -> list[tuple[str, str | None, str | None]]:
    rows: list[tuple[str, str | None, str | None]] = []
    component_columns = [c for c in plindex.columns if c.endswith("__component")]
    for column in component_columns:
        metric, threshold, directed, cluster = column.split("__")
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {directed} {cluster} built from {metric} metric with {threshold} threshold",
            )
        )
    community_columns = [c for c in plindex.columns if c.endswith("__community")]
    for column in community_columns:
        metric, threshold, cluster = column.split("__")
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {cluster} built from {metric} metric with {threshold} threshold",
            )
        )
    return rows


def get_all_column_descriptions(
    *,
    plindex: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if plindex is not None:
        make_column_descriptions(plindex=plindex)
    dfs = []
    for tsv in TSV_DIR.glob("*.tsv"):
        dfs.append(pd.read_csv(tsv, sep="\t"))
    return pd.concat(dfs).reset_index(drop=True)


def make_column_descriptions(*, plindex: pd.DataFrame) -> None:
    output_dir = TSV_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    Entry.document_properties_to_tsv(prefix="entry", filename=output_dir / "entry.tsv")
    EntryValidation.document_properties_to_tsv(
        prefix="entry_validation", filename=output_dir / "entry_validation.tsv"
    )
    System.document_properties_to_tsv(
        prefix="system", filename=output_dir / "system.tsv"
    )
    for validation_type in VALIDATION_TYPES:
        ResidueListValidation.document_properties_to_tsv(
            prefix=f"{validation_type}_validation",
            filename=output_dir / f"{validation_type}_validation.tsv",
            nested=validation_type.endswith("_chains"),
        )
        with open(output_dir / f"{validation_type}_validation.tsv", "a") as f:
            typ = "list[float]" if validation_type.endswith("_chains") else "float"
            for key in VALIDATION_OUTLIER_KEYS:
                name = f"{validation_type}_validation_percent_outliers_{key}"
                f.write(f"{name}\t{typ}\tPercent outliers for {key}\n")
    for chain_type in CHAIN_TYPES:
        Chain.document_properties_to_tsv(
            prefix=chain_type,
            filename=output_dir / f"{chain_type}.tsv",
            nested=True,
        )
    with open(output_dir / "system_pocket.tsv", "w") as f:
        f.write("Name\tType\tDescription\n")
        for key in MAPPING_NAMES:
            name = f"system_pocket_{key}"
            f.write(f"{name}\tstr\t{key} domain for the pocket\n")
    Ligand.document_properties_to_tsv(
        prefix="ligand", filename=output_dir / "ligands.tsv"
    )
    with open(output_dir / "similarity_clusters.tsv", "w") as f:
        f.write("Name\tType\tDescription\n")
        rows = get_cluster_column_descriptions(plindex)
        for row in rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
