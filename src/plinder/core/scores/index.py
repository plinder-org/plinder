# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import ensure_dataset, make_query_no_schema
from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def query_index(
    *,
    columns: list[str] | None = None,
    filters: list[tuple[str, str, str]] | None = None,
) -> pd.DataFrame:
    """
    Query the index database.

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame | None
        the index results
    """
    cfg = get_config()
    dataset = ensure_dataset(rel=cfg.data.index)
    if columns is None:
        columns = ["system_id", "entry_pdb_id"]
    query = make_query_no_schema(
        dataset=dataset,
        columns=columns,
        filters=filters,
        filename=cfg.data.index_file,
        allow_no_filters=True,
    )
    assert query is not None
    return sql(query).to_df()
