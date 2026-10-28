# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Read filtered similarity-score Parquet datasets."""

from __future__ import annotations

from pathlib import Path
from typing import Any, TypeAlias

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.scores.store import read_partitioned_scores

Filter: TypeAlias = tuple[str, str, Any]
Filters: TypeAlias = list[Filter] | list[list[Filter]] | None


def read_score_table(
    dataset: Path,
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
) -> pd.DataFrame:
    """Read selected score rows without allowing an unbounded scan."""
    if not filters:
        raise ValueError("at least one score filter is required")
    return pd.read_parquet(dataset, columns=columns, filters=filters)


def query_score_table(
    artifact: str,
    *,
    release: PlinderRelease,
    columns: list[str] | None = None,
    filters: Filters = None,
) -> pd.DataFrame:
    """Route filtered reads before downloading, with legacy-release fallback."""
    if not filters or any(not branch for branch in filters):
        raise ValueError("at least one score filter is required")
    result = read_partitioned_scores(
        release, artifact, columns=columns, filters=filters
    )
    if result is not None:
        return result
    return read_score_table(release.fetch(artifact), columns=columns, filters=filters)


__all__ = ["Filter", "Filters", "read_score_table"]
