# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import FILTERS, make_query
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import CLUSTER_SCHEMA

LOG = setup_logger(__name__)


@timeit
def query_clusters(
    *,
    columns: list[str] | None = None,
    filters: FILTERS = None,
) -> pd.DataFrame | None:
    """
    Query the cluster database.

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame | None
        the cluster results
    """

    cfg = get_config()
    dataset = cpl.get_plinder_path(rel=cfg.data.clusters)
    query = make_query(
        schema=CLUSTER_SCHEMA,
        dataset=dataset,
        filters=filters,
        columns=columns,
        nested=True,
        allow_no_filters=True,
    )
    if query is None:
        LOG.warning("try minimally passing filters=[('metric', '==', 'pli_qcov')]")
        return None
    return sql(query).to_df()
