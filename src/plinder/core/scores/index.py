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
    # START patch-1
    # TODO-1: remove this patch after binding_affinity is fixed
    if "system_has_binding_affinity" in columns or "ligand_binding_affinity" in columns:
        raise ValueError(
            "columns containing binding_affinity have been removed until bugfix"
            "see: https://github.com/plinder-org/plinder/issues/94"
        )
    # END patch-1
    query = make_query(
        dataset=dataset,
        columns=columns,
        filters=filters,
        allow_no_filters=True,
    )
    assert query is not None
    df = sql(query).to_df()
    # START patch-2
    # TODO-2: remove this patch after entry_release_date is fixed
    if "entry_release_date" in df.columns:
        from importlib import resources

        df_fixed_time = pd.read_csv(
            resources.files("plinder") / "data/utils/annotations/static_files/dates.csv"
        )[["entry_release_date", "entry_pdb_id"]]
        if "entry_pdb_id" not in df.columns:
            # hacky fix - assuming standard pdb names - to be removed
            df["entry_pdb_id"] = df.system_id.apply(lambda x: x[:4])
        df = df.drop("entry_release_date", axis=1).merge(
            df_fixed_time, on="entry_pdb_id"
        )
    # END patch-2
    if splits is None:
        splits = ["train", "val"]
    split = cpl.get_plinder_path(rel=f"{cfg.data.splits}/{cfg.data.split_file}")
    split_df = pd.read_parquet(split)
    split_dict = dict(zip(split_df["system_id"], split_df["split"]))
    df["split"] = df["system_id"].map(lambda x: split_dict.get(x, "unassigned"))
    if "*" not in splits:
        df = df[df["split"].isin(splits)].reset_index(drop=True)
    return df
