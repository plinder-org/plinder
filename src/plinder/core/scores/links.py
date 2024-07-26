# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations
from pathlib import Path

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import ensure_dataset, make_query
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import STRUCTURE_LINK_SCHEMA

LOG = setup_logger(__name__)


@timeit
def query_links(
    *,
    columns: list[str] | None = None,
    filters: list[tuple[str, str, str]] | None = None,
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
    dataset = ensure_dataset(rel=cfg.data.linked_systems)
    query = make_query(
        dataset=dataset,
        filters=filters,
        columns=columns or ["*"],
        schema=STRUCTURE_LINK_SCHEMA,
        allow_no_filters=True,
        include_filename=True,
    )
    assert query is not None
    df = sql(query).to_df()
    df["filename"] = df["filename"].apply(lambda x: Path(x).stem)
    return df
