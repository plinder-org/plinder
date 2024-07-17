from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import ensure_dataset, make_query
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import TANIMOTO_SCORE_SCHEMA

LOG = setup_logger(__name__)


@timeit
def query_ligand_similarity(
    *,
    columns: list[str] | None = None,
    filters: list[tuple[str, str, str]] | None = None,
) -> pd.DataFrame | None:
    """
    Query the ligand similarity database
    and return the results.

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame | None
        the protein similarity results
    """
    cfg = get_config()
    dataset = ensure_dataset(rel=cfg.data.ligand_scores)
    query = make_query(
        schema=TANIMOTO_SCORE_SCHEMA,
        dataset=dataset,
        filters=filters,
        columns=columns,
    )
    if query is None:
        LOG.warning(
            "try minimally passing filters=[('tanimoto_similarity_max', '>', 0.5)]"
        )
        return None

    return sql(query).to_df()
