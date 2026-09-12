# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from typing import cast

import pandas as pd

from plinder.core.index.query import query_table
from plinder.core.scores.query import Filter, Filters, read_score_table
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit


@timeit
def query_protein_similarity(
    *,
    search_db: str,
    columns: list[str] | None = None,
    filters: Filters = None,
) -> pd.DataFrame:
    """
    Query the protein similarity database for
    a given search_db and return the results.

    Parameters
    ----------
    search_db : str
        the name of the search database
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str | set[str]]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame
        The protein similarity results.
    """
    if search_db not in ["apo", "holo", "pred"]:
        raise ValueError(f"search_db={search_db} not in ['apo', 'holo', 'pred']")
    if filters and not isinstance(filters[0], list):
        filters = [
            condition
            for condition in cast(list[Filter], filters)
            if condition[0] != "search_db"
        ]
    cfg = get_config()
    dataset = cpl.get_plinder_path(rel=f"{cfg.data.scores}/search_db={search_db}")
    return read_score_table(
        dataset,
        columns=columns,
        filters=filters,
    )


@timeit
def map_cross_similarity(
    df: pd.DataFrame, target_systems: set[str], metric: str
) -> pd.DataFrame:
    target_is_requested = df["target_system"].isin(target_systems)
    df["updated_query_system"] = df["target_system"].where(
        target_is_requested, df["query_system"]
    )
    df["updated_target_system"] = df["query_system"].where(
        target_is_requested, df["target_system"]
    )
    idx = df.groupby("updated_query_system")["similarity"].idxmax()
    return df.loc[idx][
        ["updated_query_system", "updated_target_system", "similarity"]
    ].rename(
        columns={
            "updated_query_system": "query_system",
            "updated_target_system": "target_system",
            "similarity": metric,
        }
    )


@timeit
def cross_similarity(
    *,
    query_systems: set[str],
    target_systems: set[str],
    metric: str,
) -> pd.DataFrame:
    cfg = get_config()
    dataset = cpl.get_plinder_path(rel=f"{cfg.data.scores}/search_db=holo")
    filters: list[list[Filter]] = [
        [
            ("metric", "==", metric),
            ("query_system", "in", query_systems),
            ("target_system", "in", target_systems),
        ],
        [
            ("metric", "==", metric),
            ("query_system", "in", target_systems),
            ("target_system", "in", query_systems),
        ],
    ]
    columns = ["query_system", "target_system", "similarity"]
    similarities = read_score_table(
        dataset,
        columns=columns,
        filters=filters,
    )
    return map_cross_similarity(similarities, target_systems, metric)


@timeit
def multi_query_protein_similarity(
    *,
    system_id: str,
    search_db: str,
    filter_criteria: dict[str, int],
) -> pd.DataFrame:
    """
    Searches the protein similarity database for systems satisfying ALL filter criteria

    Parameters
    ----------
    system_id : str
        the system_id to search for
    search_db : str
        the search database to search in
    filter_criteria : dict[str, int]
        Metric and threshold pairs. For example::

            {
                "pocket_fident": 100,
                "protein_fident_weighted_sum": 95,
                "protein_fident_qcov_weighted_sum": 80,
                "pocket_lddt": 20,
                "protein_lddt_weighted_sum": 20,
            }

    Returns
    -------
    df : pd.DataFrame
        the protein similarity results across all metrics in filter_criteria
    """
    empty_df = pd.DataFrame(
        columns=["query_system", "target_system"] + list(filter_criteria.keys())
    )
    if search_db == "holo":
        target_systems_df = query_table("annotation", columns=["system_id"])
        target_systems = set(target_systems_df["system_id"])
        if not target_systems:
            return empty_df
    filters: list[list[Filter]] = []
    for metric, threshold in filter_criteria.items():
        conditions: list[Filter] = [
            ("metric", "==", metric),
            ("similarity", ">=", threshold),
            ("query_system", "==", system_id),
        ]
        if search_db == "holo":
            conditions.append(("target_system", "in", target_systems))
        filters.append(conditions)
    links = query_protein_similarity(
        search_db=search_db,
        columns=["query_system", "target_system", "metric", "similarity"],
        filters=filters,
    )
    if links.empty:
        return empty_df
    links = links.iloc[
        links.groupby(["query_system", "target_system", "metric"], observed=True)[
            "similarity"
        ].idxmax()
    ]
    links = links.pivot(
        index=["query_system", "target_system"],
        columns="metric",
        values="similarity",
    ).reset_index()
    if not filter_criteria.keys() <= set(links.columns):
        return empty_df
    keep = pd.Series(True, index=links.index)
    for metric, threshold in filter_criteria.items():
        keep &= links[metric] >= threshold
    return links.loc[keep]
