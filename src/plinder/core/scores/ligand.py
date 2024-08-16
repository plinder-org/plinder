# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import pandas as pd
from duckdb import sql

from plinder.core.scores.query import FILTER, FILTERS, make_query
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import TANIMOTO_SCORE_SCHEMA

LOG = setup_logger(__name__)


@timeit
def query_ligand_similarity(
    *,
    columns: list[str] | None = None,
    filters: FILTERS = None,
) -> pd.DataFrame | None:
    """
    Query the ligand similarity database
    and return the results.

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str | set[str]]]
        the filters to apply

    Returns
    -------
    df : pd.DataFrame | None
        the protein similarity results
    """
    cfg = get_config()
    dataset = cpl.get_plinder_path(rel=cfg.data.ligand_scores)
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


@timeit
def map_cross_similarity(
    df: pd.DataFrame, target_ligands: set[str], metric: str
) -> pd.DataFrame:
    updated_query_ligands = []
    for q, t in zip(df["query_ligand_id"], df["target_ligand_id"]):
        if t in target_ligands:
            updated_query_ligands.append(t)
        else:
            updated_query_ligands.append(q)
    df["updated_query_ligand_id"] = updated_query_ligands
    idx = df.groupby("updated_query_ligand_id")["tanimoto_similarity_max"].idxmax()
    df = df.loc[idx]

    cfg = get_config()
    dataset = cpl.get_plinder_path(
        rel=f"{cfg.data.fingerprints}/{cfg.data.fingerprint_file}"
    )
    ligands_per_system = pd.read_parquet(dataset)
    ligand_to_system: dict[int, set[str]] = {}
    for ligand_id, group in ligands_per_system.groupby("number_id_by_inchikeys"):
        ligand_to_system[int(ligand_id)] = set(group["system_id"])
    df["query_system"] = df["query_ligand_id"].map(ligand_to_system)
    return (
        df.explode("query_system")
        .rename(
            columns={
                "query_system": "system_id",
            }
        )
        .drop_duplicates("system_id")[["system_id", metric]]
        .reset_index(drop=True)
    )


@timeit
def cross_similarity(
    *,
    query_ligands: set[str],
    target_ligands: set[str],
    metric: str = "tanimoto_similarity_max",
) -> pd.DataFrame:
    """
    Query the ligand similarity database for
    a cross similarity between a set of query
    and target ligands.

    Parameters
    ----------
    query_ligands : set[str]
        the set of query ligands
    target_ligands : set[str]
        the set of target ligands

    Returns
    -------
    df : pd.DataFrame
        the cross similarity results
    """
    cfg = get_config()
    dataset = cpl.get_plinder_path(rel=cfg.data.ligand_scores)
    filters = [
        [
            FILTER(("query_ligand_id", "in", query_ligands)),
            FILTER(("target_ligand_id", "in", target_ligands)),
        ],
        [
            FILTER(("query_ligand_id", "in", target_ligands)),
            FILTER(("target_ligand_id", "in", query_ligands)),
        ],
    ]
    columns = ["query_ligand_id", "target_ligand_id", "tanimoto_similarity_max"]
    query = make_query(
        schema=TANIMOTO_SCORE_SCHEMA,
        dataset=dataset,
        columns=columns,
        filters=filters,
    )
    assert query is not None
    return map_cross_similarity(sql(query).to_df(), target_ligands, metric)
