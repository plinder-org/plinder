# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd

from plinder.data import column_descriptions
from plinder.data.annotations.aggregate_annotations import Entry, System
from plinder.data.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
)
from plinder.data.annotations.ligand_utils import Ligand
from plinder.data.annotations.protein_utils import Chain

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
    "Pfam",
    "SCOP2",
    "SCOP2B",
    "UniProt",
]
CHAIN_TYPES = [
    "system_protein_chains",
    "system_ligand_chains",
    "ligand_interacting_ligand_chains",
    "ligand_neighboring_ligand_chains",
    "ligand_protein_chains",
]
DERIVED_LIGAND_COLUMNS = [
    (
        "ligand_is_3d_score_able",
        "bool | None",
        "Whether the canonical ASU SDF loads and supports finite shape, color, and SuCOS self-scoring",
    ),
    (
        "ligand_smiles_id",
        "int | None",
        "Integer node ID assigned to this exact canonical SMILES for ligand similarity scoring",
    ),
    (
        "ligand_max_cofactor_similarity",
        "float | None",
        "Maximum ECFP4/1024 Tanimoto similarity (percent) to any CCD structure in the cofactor list",
    ),
    (
        "ligand_most_similar_cofactor",
        "str | None",
        "CCD code of the cofactor with maximum ECFP4/1024 Tanimoto similarity",
    ),
    (
        "ligand_is_cofactor_like",
        "bool | None",
        "Whether maximum similarity to a listed CCD cofactor is at least 90 percent",
    ),
    (
        "ligand_tanimoto_ecfp4_1024_90_cluster",
        "str | None",
        "Connected-component ID among unique ligand SMILES using Tanimoto similarity of at least 90 percent",
    ),
    (
        "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids",
        "int | None",
        "Number of distinct PDB entries containing this ligand or another ligand in its 90-percent Tanimoto component",
    ),
]


def get_cluster_column_descriptions(
    plindex: pd.DataFrame,
) -> list[tuple[str, str | None, str | None]]:
    rows: list[tuple[str, str | None, str | None]] = []
    component_columns = [c for c in plindex.columns if c.endswith("component")]
    for column in component_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        ligand_level = parts[2] == "ligand"
        half_interface = parts[-1].startswith("chain_")
        direction = (
            "reciprocal-minimum"
            if ligand_level or metric.startswith("interface_")
            else parts[-2]
            if parts[-2] in {"weak", "strong"}
            else "directed"
        )
        cluster = "component"
        level = (
            "ligand-level "
            if ligand_level
            else f"{parts[-1].removesuffix('_component').replace('_', ' ')} "
            if half_interface
            else ""
        )
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {level}{direction} {cluster} built from "
                f"{metric} metric with {threshold} threshold",
            )
        )
    community_columns = [c for c in plindex.columns if c.endswith("community")]
    for column in community_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        ligand_level = parts[2] == "ligand"
        half_interface = parts[-1].startswith("chain_")
        cluster = "community"
        level = (
            "ligand-level "
            if ligand_level
            else f"{parts[-1].removesuffix('_community').replace('_', ' ')} "
            if half_interface
            else ""
        )
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {level}greedy centroid {cluster} built from "
                f"reciprocal-minimum {metric} with {threshold} threshold; each "
                "member meets the threshold in both directions to its centroid",
            )
        )
    directed_cover_columns = [
        c for c in plindex.columns if c.endswith("directed_set_cover")
    ]
    for column in directed_cover_columns:
        parts = column.split("__")
        metric, threshold = parts[:2]
        ligand_level = parts[2] == "ligand"
        half_interface = parts[-1].startswith("chain_")
        level = (
            "ligand-level "
            if ligand_level
            else f"{parts[-1].removesuffix('_directed_set_cover').replace('_', ' ')} "
            if half_interface
            else ""
        )
        rows.append(
            (
                column,
                "str",
                f"Cluster ID for {level}directed set cover built from "
                f"{metric} with {threshold} threshold; each member's "
                "query-to-centroid score meets the threshold",
            )
        )
    column_order = {column: index for index, column in enumerate(plindex.columns)}
    return sorted(rows, key=lambda row: column_order[row[0]])


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
    with (output_dir / "ligands.tsv").open("a") as f:
        for name, typ, description in DERIVED_LIGAND_COLUMNS:
            f.write(f"{name}\t{typ}\t{description}\n")
    with open(output_dir / "similarity_clusters.tsv", "w") as f:
        f.write("Name\tType\tDescription\n")
        rows = get_cluster_column_descriptions(plindex)
        for row in rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
