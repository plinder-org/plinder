# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import FILTERS, make_query
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


@timeit
def query_links(
    *,
    columns: list[str] | None = None,
    filters: FILTERS = None,
) -> pd.DataFrame:
    """
    Query the linked systems dataset

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame
        the linked systems results
    """
    cfg = get_config()
    dataset = cpl.get_plinder_path(rel=cfg.data.links)
    query = make_query(
        dataset=dataset,
        filters=filters,
        columns=columns or ["*"],
        allow_no_filters=True,
        include_filename=True,
    )
    assert query is not None
    df = sql(query).to_df()
    df["kind"] = df["filename"].apply(lambda x: Path(x).stem.split("_links")[0])
    return df
