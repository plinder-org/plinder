# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from plinder.core.index.query import query_table
from plinder.core.release import PlinderRelease
from plinder.core.scores.query import Filter, Filters, read_score_table
from plinder.core.utils.dec import timeit


def _ligand_ids(values: Iterable[int | str]) -> set[int]:
    if isinstance(values, (str, bytes)):
        raise TypeError("ligand IDs must be provided as a collection")
    try:
        return {int(value) for value in values}
    except (TypeError, ValueError) as exc:
        raise ValueError("ligand IDs must be integers") from exc


@timeit
def query_ligand_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query complete directed ligand-pair similarities.

    Parameters
    ----------
    columns : list[str], default=None
        the columns to return
    filters : list[tuple[str, str, str | set[str]]]
        the filters to apply
    release : PlinderRelease | None
        Explicit local release, or the configured release when omitted.

    Returns
    -------
    df : pd.DataFrame
        The ligand similarity results.
    """
    dataset = (release or PlinderRelease()).fetch("ligand_similarity_scores")
    return read_score_table(
        dataset,
        columns=columns,
        filters=filters,
    )


@timeit
def query_ligand_chemical_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query ECFP4 Tanimoto similarities between canonical ligand SMILES."""
    dataset = (release or PlinderRelease()).fetch("ligand_scores")
    return read_score_table(dataset, columns=columns, filters=filters)


@timeit
def _map_cross_chemical_similarity(
    df: pd.DataFrame,
    target_ligands: set[int],
    metric: str,
    release: PlinderRelease | None,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            {
                "system_id": pd.Series(dtype="object"),
                metric: pd.Series(dtype="float64"),
            }
        )
    target_is_requested = df["target_ligand_id"].isin(target_ligands)
    df["updated_query_ligand_id"] = df["target_ligand_id"].where(
        target_is_requested, df["query_ligand_id"]
    )
    idx = df.groupby("updated_query_ligand_id")[metric].idxmax()
    df = df.loc[idx]

    ligand_ids = set(df["query_ligand_id"].astype(int))
    ligand_occurrences = query_table(
        "annotation",
        columns=["system_id", "ligand_smiles_id"],
        filters=[("ligand_smiles_id", "in", ligand_ids)],
        release=release,
    )
    id_column = "ligand_smiles_id"
    ligand_to_system: dict[int, set[str]] = {}
    for ligand_id, group in ligand_occurrences.groupby(id_column):
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
def cross_ligand_chemical_similarity(
    *,
    query_ligands: Iterable[int | str],
    target_ligands: Iterable[int | str],
    metric: str | None = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """
    Query the ligand similarity database for
    a cross similarity between a set of query
    and target ligands.

    Parameters
    ----------
    query_ligands : Iterable[int | str]
        The query ligand IDs.
    target_ligands : Iterable[int | str]
        The target ligand IDs.

    Returns
    -------
    df : pd.DataFrame
        the cross similarity results
    """
    query_ids = _ligand_ids(query_ligands)
    target_ids = _ligand_ids(target_ligands)
    dataset = (release or PlinderRelease()).fetch("ligand_scores")
    if metric is None:
        metric = "tanimoto_similarity_ecfp4_1024"
    filters: list[list[Filter]] = [
        [
            ("query_ligand_id", "in", query_ids),
            ("target_ligand_id", "in", target_ids),
        ],
        [
            ("query_ligand_id", "in", target_ids),
            ("target_ligand_id", "in", query_ids),
        ],
    ]
    columns = ["query_ligand_id", "target_ligand_id", metric]
    similarities = read_score_table(
        dataset,
        columns=columns,
        filters=filters,
    )
    return _map_cross_chemical_similarity(
        similarities,
        target_ids,
        metric,
        release,
    )
