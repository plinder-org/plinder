from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import ensure_dataset, make_query_no_schema
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def query_index(
    *,
    columns: list[str] | None = None,
    filters: list[tuple[str, str, str]] | None = None,
) -> pd.DataFrame | None:
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
    raise NotImplementedError(
        "duckdb reads ((1988116, 487) vs (1748019, 487)) pyarrow :grimacing:"
    )
    dataset = ensure_dataset(rel="index")
    query = make_query_no_schema(
        dataset=dataset,
        columns=columns,
        filters=filters,
        nested=True,
        allow_no_filters=True,
    )
    if query is None:
        LOG.warning("try minimally passing filters=[('system_type', '==', 'holo')]")
        return None

    return sql(query).to_df()
