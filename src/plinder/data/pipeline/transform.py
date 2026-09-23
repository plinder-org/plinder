# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import re
from pathlib import Path

import numpy as np
import pandas as pd

AFFINITY_RECORD_COLUMNS = [
    "pdbid_ligid",
    "pdb_id",
    "ligand_het_id",
    "reactant_set_id",
    "monomer_id",
    "curation_source",
    "article_doi",
    "bindingdb_entry_doi",
    "source_row",
    "target_sequences",
    "target_sequence",
    "endpoint",
    "raw_value",
    "relation",
    "value_nm",
    "pchembl",
]


def transform_bindingdb_measurements(
    *, raw_affinity_path: Path, chunksize: int = 100_000
) -> pd.DataFrame:
    """Keep each Ki/Kd measurement linked to its source row and target.

    PDB/CCD cross-references are candidate matches, not evidence that an
    assay measured a particular deposited complex. ``relation`` applies to
    the logarithmic pKi/pKd value, so a raw ``Ki > 100 nM`` is ``pKi < 7``.
    """
    header = pd.read_csv(raw_affinity_path, sep="\t", nrows=0).columns.tolist()
    sequence_columns = [
        column
        for column in header
        if re.fullmatch(r"BindingDB Target Chain\s+Sequence(?: \d+)?", column)
    ]
    if not sequence_columns:
        raise ValueError("BindingDB TSV has no target-chain sequence columns")
    required = [
        "Ligand HET ID in PDB",
        "PDB ID(s) for Ligand-Target Complex",
        "Ki (nM)",
        "Kd (nM)",
    ]
    missing = sorted(set(required).difference(header))
    if missing:
        raise ValueError(f"BindingDB TSV is missing columns: {missing}")
    reactant_column = "BindingDB Reactant_set_id"
    monomer_column = "BindingDB MonomerID"
    source_columns = {
        "curation_source": "Curation/DataSource",
        "article_doi": "Article DOI",
        "bindingdb_entry_doi": "BindingDB Entry DOI",
    }
    columns = required + sequence_columns
    if reactant_column in header:
        columns.append(reactant_column)
    if monomer_column in header:
        columns.append(monomer_column)
    columns.extend(column for column in source_columns.values() if column in header)

    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        raw_affinity_path,
        sep="\t",
        usecols=columns,
        dtype="string",
        chunksize=chunksize,
    ):
        chunk = chunk.loc[
            chunk[required[0]].notna() & chunk[required[1]].notna()
        ].copy()
        if chunk.empty:
            continue
        chunk["source_row"] = chunk.index
        chunk["target_sequences"] = chunk[sequence_columns].apply(
            lambda row: [
                value.strip().upper()
                for value in row
                if isinstance(value, str) and value.strip()
            ],
            axis=1,
        )
        chunk["target_sequence"] = chunk["target_sequences"].map(
            lambda values: values[0] if len(values) == 1 else None
        )
        chunk["reactant_set_id"] = (
            chunk[reactant_column] if reactant_column in chunk else pd.NA
        )
        chunk["monomer_id"] = (
            chunk[monomer_column] if monomer_column in chunk else pd.NA
        )
        for output, source in source_columns.items():
            chunk[output] = chunk[source] if source in chunk else pd.NA
        for endpoint in ("Ki", "Kd"):
            source = f"{endpoint} (nM)"
            records = chunk.loc[chunk[source].notna()].copy()
            records["raw_value"] = records[source].str.strip()
            records = records.loc[records["raw_value"].ne("")].copy()
            if records.empty:
                continue
            parsed = records["raw_value"].str.extract(
                r"^(?P<bound><=|>=|<|>|=)?\s*"
                r"(?P<value>(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)$"
            )
            records["value_nm"] = pd.to_numeric(parsed["value"], errors="coerce")
            records.loc[records["value_nm"] <= 0, "value_nm"] = np.nan
            records["relation"] = (
                parsed["bound"]
                .fillna("=")
                .map({"=": "=", "<": ">", "<=": ">=", ">": "<", ">=": "<="})
            )
            records.loc[records["value_nm"].isna(), "relation"] = None
            with np.errstate(divide="ignore", invalid="ignore"):
                records["pchembl"] = 9 - np.log10(records["value_nm"])
            records["endpoint"] = endpoint
            records["pdb_id"] = records[required[1]].str.split(r"[,;]")
            records = records.explode("pdb_id")
            records["pdb_id"] = records["pdb_id"].str.strip().str.upper()
            records = records.loc[records["pdb_id"].str.fullmatch(r"[A-Z0-9]{4}")]
            records["ligand_het_id"] = records[required[0]].str.strip().str.upper()
            records["pdbid_ligid"] = records["pdb_id"] + "_" + records["ligand_het_id"]
            parts.append(records[AFFINITY_RECORD_COLUMNS])

    if not parts:
        return pd.DataFrame(columns=AFFINITY_RECORD_COLUMNS)
    return pd.concat(parts, ignore_index=True).drop_duplicates(
        ["source_row", "pdb_id", "ligand_het_id", "endpoint"]
    )


def strict_bindingdb_candidates(records: pd.DataFrame) -> pd.DataFrame:
    """Summarize only comparable, uncensored measurements per target/endpoint."""
    columns = ["pdbid_ligid", "target_sequence", "endpoint", "pchembl", "count"]
    if records.empty:
        return pd.DataFrame(columns=columns)
    # Multi-chain or missing-target records cannot be checked against one
    # receptor sequence; do not let them disappear before the ambiguity test.
    unusable_keys = records.loc[
        records["target_sequence"].isna(), "pdbid_ligid"
    ].unique()
    eligible = records.loc[~records["pdbid_ligid"].isin(unusable_keys)].copy()
    if eligible.empty:
        return pd.DataFrame(columns=columns)
    eligible["valid"] = (
        eligible["relation"].eq("=") & eligible["pchembl"].notna()
    ).fillna(False)
    grouped = eligible.groupby(["pdbid_ligid", "target_sequence"], sort=False).agg(
        endpoint=("endpoint", "first"),
        endpoint_count=("endpoint", "nunique"),
        all_valid=("valid", "all"),
        pchembl=("pchembl", "median"),
        count=("pchembl", "size"),
    )
    grouped = grouped.loc[grouped["all_valid"] & grouped["endpoint_count"].eq(1)]
    return grouped.reset_index()[columns]
