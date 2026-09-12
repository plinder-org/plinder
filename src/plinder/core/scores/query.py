# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Read filtered similarity-score Parquet datasets."""

from __future__ import annotations

from pathlib import Path
from typing import Any, TypeAlias

import pandas as pd

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


__all__ = ["Filter", "Filters", "read_score_table"]
