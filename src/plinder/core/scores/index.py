# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import FILTERS, make_query
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def query_index(
    *,
    columns: list[str] | None = None,
    splits: list[str] | None = None,
    filters: FILTERS = None,
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
    dataset = cpl.get_plinder_path(rel=f"{cfg.data.index}/{cfg.data.index_file}")
    if columns is None:
        columns = ["system_id", "entry_pdb_id"]
    if "system_id" not in columns and "*" not in columns:
        columns = ["system_id"] + columns
    query = make_query(
        dataset=dataset,
        columns=columns,
        filters=filters,
        allow_no_filters=True,
    )
    assert query is not None
    df = sql(query).to_df()
    if splits is None:
        splits = ["train", "val"]
    split = cpl.get_plinder_path(rel=f"{cfg.data.splits}/{cfg.data.split_file}")
    split_df = pd.read_parquet(split)
    split_dict = dict(zip(split_df["system_id"], split_df["split"]))
    df["split"] = df["system_id"].map(lambda x: split_dict.get(x, "unassigned"))
    if "*" not in splits:
        df = df[df["split"].isin(splits)].reset_index(drop=True)
    return df
