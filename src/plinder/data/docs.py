# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd

from plinder.data.utils.annotations.aggregate_annotations import Entry, System
from plinder.data.utils.annotations.protein_utils import Chain
from plinder.data.utils.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
)
from plinder.data.utils.annotations.ligand_utils import Ligand


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
    "BIRD",
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


def get_cluster_column_descriptions(plindex: pd.DataFrame) -> list[tuple[str, str | None, str | None]]:
    rows = []
    component_columns = [c for c in plindex.columns if c.endswith("__component")]
    for column in component_columns:
        metric, threshold, directed, cluster = column.split("__")
        rows.append(
            (
                column,
                str,
                f"Cluster ID for {directed} {cluster} built from {metric} metric with {threshold} threshold",
            )
        )
    community_columns = [c for c in plindex.columns if c.endswith("__community")]
    for column in community_columns:
        metric, threshold, cluster = column.split("__")
        rows.append(
            (
                column,
                str,
                f"Cluster ID for {cluster} built from {metric} metric with {threshold} threshold",
            )
        )
    return rows


def get_all_column_descriptions(*, plindex: pd.DataFrame) -> list[tuple[str, str | None, str | None]]:
    descriptions = []
    descriptions.extend(list(Entry.document_properties("entry")))
    descriptions.extend(list(EntryValidation.document_properties("entry_validation")))
    descriptions.extend(list(System.document_properties("system")))
    for validation_type in VALIDATION_TYPES:
        descriptions.extend(
            list(ResidueListValidation.document_properties(f"{validation_type}_validation"))
        )
        descriptions.extend([
            (
                f"{validation_type}_validation_percent_outliers_{key}",
                float,
                f"Percent outliers for {key}",
            )
            for key in VALIDATION_OUTLIER_KEYS
        ])
    for chain_type in CHAIN_TYPES:
        descriptions.extend(list(Chain.document_properties(chain_type)))
        descriptions.extend([
            (
                f"{chain_type}_{key}",
                dict[str, tuple[str, str]],
                f"Domains and ranges for {key}",
            )
            for key in MAPPING_NAMES
        ])
    descriptions.extend([
        (
            f"system_pocket_{key}",
            dict[str, tuple[str, str]],
            f"Domains and ranges for {key}",
        )
        for key in MAPPING_NAMES
    ])
    descriptions.extend(list(Ligand.document_properties("ligand")))
    descriptions.extend(get_cluster_column_descriptions(plindex))

    hack_path = Path.cwd() / "column_descriptions" 
    for tsv in ["extra.tsv", "posebusters_checks.tsv", "qc.tsv"]:
        if not (hack_path / tsv).is_file():
            raise ValueError(f"missing {hack_path / tsv}")
        df = pd.read_csv(hack_path / tsv, sep="\t")
        df["Type"] = df["Type"].apply(eval)
        descriptions.extend(df.itertuples(index=False))
    return descriptions


def make_column_descriptions(*, plindex: pd.DataFrame, output_dir: Path) -> None:
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
        )
        with open(output_dir / f"{validation_type}_validation.tsv", "a") as f:
            for key in VALIDATION_OUTLIER_KEYS:
                name = f"{validation_type}_validation_percent_outliers_{key}"
                f.write(f"{name}\tfloat\tPercent outliers for {key}\n")
    for chain_type in CHAIN_TYPES:
        Chain.document_properties_to_tsv(
            prefix=chain_type,
            filename=output_dir / f"{chain_type}.tsv",
        )
        with open(output_dir / f"{chain_type}.tsv", "a") as f:
            for key in MAPPING_NAMES:
                name = f"{chain_type}_{key}"
                f.write(
                    f"{name}\tdict[str, tuple[str, str]]\tDomains and ranges for {key}\n"
                )
    with open(output_dir / "system_pocket.tsv", "w") as f:
        f.write("Name\tType\tDescription\n")
        for key in MAPPING_NAMES:
            name = f"system_pocket_{key}"
            f.write(
                f"{name}\tdict[str, tuple[str, str]]\tDomains and ranges for {key}\n"
            )
    Ligand.document_properties_to_tsv(
        prefix="ligand", filename=output_dir / "ligands.tsv"
    )
    with open(output_dir / "similarity_clusters.tsv", "w") as f:
        f.write("Name\tType\tDescription\n")
        rows = get_cluster_column_descriptions(plindex)
        for row in rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")


if __name__ == "__main__":
    from plinder.data.pipeline import config

    import types

    cfg = config.get_config()
    plindex = pd.read_parquet(Path(cfg.data.plinder_dir) / cfg.data.index / cfg.data.index_file)
    # plindex = pd.read_parquet(Path("column_descriptions/h8.parquet"))
    count = 0
    aliases = 0
    columns = set()
    for row in get_all_column_descriptions(plindex=plindex):
        print(row[0], row[1], row[2], type(row[1]))
        if isinstance(row[1], types.GenericAlias):
            aliases += 1
        count += 1
        columns.add(row[0])
    print(f"total {count} columns, aliases={aliases}")
    print("in plindex not in columns", len(set(plindex.columns).difference(columns)))
    print("in plindex not in columns", set(plindex.columns).difference(columns))
    print("in columns not in plindex", len(columns.difference(plindex.columns)))
