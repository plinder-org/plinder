# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd

from plinder.core.index.query import query_table
from plinder.core.release import PlinderRelease
from plinder.core.scores.query import Filter, Filters, read_score_table
from plinder.core.utils.dec import timeit


@timeit
def query_ligand_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
) -> pd.DataFrame:
    """
    Query the ligand similarity database
    and return the results.

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str | set[str]]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame
        The ligand similarity results.
    """
    dataset = PlinderRelease().fetch("ligand_scores")
    return read_score_table(
        dataset,
        columns=columns,
        filters=filters,
    )


@timeit
def map_cross_similarity(
    df: pd.DataFrame, target_ligands: set[str], metric: str
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            {
                "system_id": pd.Series(dtype="object"),
                metric: pd.Series(dtype="float64"),
            }
        )
    target_is_requested = df["target_ligand_id"].isin(target_ligands)
    df["updated_query_ligand_id"] = df["target_ligand_id"].where(
        target_is_requested, df["query_ligand_id"]
    )
    idx = df.groupby("updated_query_ligand_id")[metric].idxmax()
    df = df.loc[idx]

    ligand_ids = set(df["query_ligand_id"].astype(int))
    ligand_occurrences = query_table(
        "annotation",
        columns=["system_id", "ligand_smiles_id"],
        filters=[("ligand_smiles_id", "in", ligand_ids)],
    )
    id_column = "ligand_smiles_id"
    ligand_to_system: dict[int, set[str]] = {}
    for ligand_id, group in ligand_occurrences.groupby(id_column):
        ligand_to_system[int(ligand_id)] = set(group["system_id"])
    df["query_system"] = df["query_ligand_id"].map(ligand_to_system)
    return (
        df.explode("query_system")
        .rename(
            columns={
                "query_system": "system_id",
            }
        )
        .drop_duplicates("system_id")[["system_id", metric]]
        .reset_index(drop=True)
    )


@timeit
def cross_similarity(
    *,
    query_ligands: set[str],
    target_ligands: set[str],
    metric: str | None = None,
) -> pd.DataFrame:
    """
    Query the ligand similarity database for
    a cross similarity between a set of query
    and target ligands.

    Parameters
    ----------
    query_ligands : set[str]
        the set of query ligands
    target_ligands : set[str]
        the set of target ligands

    Returns
    -------
    df : pd.DataFrame
        the cross similarity results
    """
    dataset = PlinderRelease().fetch("ligand_scores")
    if metric is None:
        metric = "tanimoto_similarity_ecfp4_1024"
    filters: list[list[Filter]] = [
        [
            ("query_ligand_id", "in", query_ligands),
            ("target_ligand_id", "in", target_ligands),
        ],
        [
            ("query_ligand_id", "in", target_ligands),
            ("target_ligand_id", "in", query_ligands),
        ],
    ]
    columns = ["query_ligand_id", "target_ligand_id", metric]
    similarities = read_score_table(
        dataset,
        columns=columns,
        filters=filters,
    )
    return map_cross_similarity(similarities, target_ligands, metric)
