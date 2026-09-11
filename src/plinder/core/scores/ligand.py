# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from typing import cast

import pandas as pd
from duckdb import sql

from plinder.core.release import PlinderRelease
from plinder.core.scores.query import FILTER, FILTERS, make_query
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import TANIMOTO_SCORE_SCHEMA

LOG = setup_logger(__name__)


@timeit
def query_ligand_similarity(
    *,
    columns: list[str] | None = None,
    filters: FILTERS = None,
) -> pd.DataFrame | None:
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
    df : pd.DataFrame | None
        the protein similarity results
    """
    dataset = PlinderRelease().fetch("ligand_scores")
    schema = TANIMOTO_SCORE_SCHEMA
    metric = schema.names[-1]
    query = make_query(
        schema=schema,
        dataset=dataset,
        filters=filters,
        columns=columns,
    )
    if query is None:
        LOG.warning(f"try minimally passing filters=[('{metric}', '>', 50)]")
        return None

    return sql(query).to_df()


@timeit
def map_cross_similarity(
    df: pd.DataFrame, target_ligands: set[str], metric: str
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame({
            "system_id": pd.Series(dtype="object"),
            metric: pd.Series(dtype="float64"),
        })
    updated_query_ligands = []
    for q, t in zip(df["query_ligand_id"], df["target_ligand_id"]):
        if t in target_ligands:
            updated_query_ligands.append(t)
        else:
            updated_query_ligands.append(q)
    df["updated_query_ligand_id"] = updated_query_ligands
    idx = df.groupby("updated_query_ligand_id")[metric].idxmax()
    df = df.loc[idx]

    from plinder.core.scores.index import query_index

    ligand_ids = set(df["query_ligand_id"].astype(int))
    ligand_occurrences = query_index(
        columns=["system_id", "ligand_smiles_id"],
        filters=[
            FILTER((
                "ligand_smiles_id",
                "in",
                cast(set[str], ligand_ids),
            ))
        ],
    )
    id_column = "ligand_smiles_id"
    ligand_to_system: dict[int, set[str]] = {}
    for ligand_id, group in ligand_occurrences.groupby(id_column):
        ligand_to_system[int(ligand_id)] = set(group["system_id"])
    df["query_system"] = df["query_ligand_id"].map(ligand_to_system)
    return (
        df
        .explode("query_system")
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
    schema = TANIMOTO_SCORE_SCHEMA
    if metric is None:
        metric = "tanimoto_similarity_ecfp4_1024"
    filters = [
        [
            FILTER(("query_ligand_id", "in", query_ligands)),
            FILTER(("target_ligand_id", "in", target_ligands)),
        ],
        [
            FILTER(("query_ligand_id", "in", target_ligands)),
            FILTER(("target_ligand_id", "in", query_ligands)),
        ],
    ]
    columns = ["query_ligand_id", "target_ligand_id", metric]
    query = make_query(
        schema=schema,
        dataset=dataset,
        columns=columns,
        filters=filters,
    )
    assert query is not None
    return map_cross_similarity(sql(query).to_df(), target_ligands, metric)
