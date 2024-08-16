# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import FILTER, FILTERS, make_query
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import PROTEIN_SIMILARITY_SCHEMA

LOG = setup_logger(__name__)


@timeit
def query_protein_similarity(
    *,
    search_db: str,
    columns: list[str] | None = None,
    filters: FILTERS = None,
) -> pd.DataFrame | None:
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
    df : pd.DataFrame | None
        the protein similarity results
    """
    if search_db not in ["apo", "holo", "pred"]:
        raise ValueError(f"search_db={search_db} not in ['apo', 'holo', 'pred']")
    cfg = get_config()
    if isinstance(filters, list):
        if len(filters) and not isinstance(filters[0], list):
            for i, (col, _, _) in enumerate(filters):
                if col == "search_db":
                    filters.pop(i)
                    break
    dataset = cpl.get_plinder_path(rel=f"{cfg.data.scores}/search_db={search_db}")
    query = make_query(
        schema=PROTEIN_SIMILARITY_SCHEMA,
        dataset=dataset,
        filters=filters,
        columns=columns,
    )
    if query is None:
        LOG.warning("try minimally passing filters=[('similarity', '>', 90)]")
        return None
    return sql(query).to_df()


@timeit
def map_cross_similarity(
    df: pd.DataFrame, target_systems: set[str], metric: str
) -> pd.DataFrame:
    updated_query_systems = []
    updated_target_systems = []
    for q, t in zip(df["query_system"], df["target_system"]):
        if t in target_systems:
            updated_query_systems.append(t)
            updated_target_systems.append(q)
        else:
            updated_query_systems.append(q)
            updated_target_systems.append(t)
    df["updated_query_system"] = updated_query_systems
    df["updated_target_system"] = updated_target_systems
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
    filters: list[list[FILTER]] = [
        [
            FILTER(("metric", "==", metric)),
            FILTER(("query_system", "in", query_systems)),
            FILTER(("target_system", "in", target_systems)),
        ],
        [
            FILTER(("metric", "==", metric)),
            FILTER(("query_system", "in", target_systems)),
            FILTER(("target_system", "in", query_systems)),
        ],
    ]
    columns = ["query_system", "target_system", "similarity"]
    query = make_query(
        schema=PROTEIN_SIMILARITY_SCHEMA,
        dataset=dataset,
        columns=columns,
        filters=filters,
    )
    assert query is not None
    return map_cross_similarity(sql(query).to_df(), target_systems, metric)
