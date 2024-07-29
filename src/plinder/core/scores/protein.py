# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import ensure_dataset, make_query
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
    filters: list[tuple[str, str, str]] | None = None,
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
    filters : list[tuple[str, str, str]]
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
        for i, (col, _, _) in enumerate(filters):
            if col == "search_db":
                filters = filters[:i] + filters[i + 1 :]
                break
    dataset = ensure_dataset(rel=f"{cfg.data.scores}/search_db={search_db}")
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
def cross_similarity(
    *,
    query_systems: list[str],
    target_systems: list[str],
    metric: str,
) -> pd.DataFrame:
    cfg = get_config()
    dataset = ensure_dataset(rel=f"{cfg.data.scores}/search_db=holo")
    filters = [
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
    query = make_query(
        schema=PROTEIN_SIMILARITY_SCHEMA,
        dataset=dataset,
        columns=columns,
        filters=filters,
    )
    assert query is not None
    return sql(query).to_df()
