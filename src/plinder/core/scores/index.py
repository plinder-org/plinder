# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from typing import cast

import pandas as pd

from plinder.core.index.query import Filters, query_table
from plinder.core.scores.query import FILTERS
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def _filter_columns(filters: FILTERS) -> set[str]:
    return {
        condition[0]
        for item in filters or []
        for condition in (item if isinstance(item, list) else [item])
    }


def query_index(
    *,
    columns: list[str] | None = None,
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
    if columns is None:
        columns = ["system_id", "entry_pdb_id"]
    if "system_id" not in columns and "*" not in columns:
        columns = ["system_id"] + columns
    requested_columns = set(columns) | _filter_columns(filters)
    needs_entry_metadata = "*" in requested_columns or any(
        column.startswith("entry_") and column != "entry_pdb_id"
        for column in requested_columns
    )
    return query_table(
        "annotation",
        columns=columns,
        filters=cast(Filters, filters),
        joins=["entry_metadata"] if needs_entry_metadata else None,
    )
