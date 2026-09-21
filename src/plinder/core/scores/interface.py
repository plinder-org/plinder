# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Query complete protein-interface similarities."""

from __future__ import annotations

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.scores.query import Filters, read_score_table
from plinder.core.utils.dec import timeit


@timeit
def query_interface_similarity(
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query complete directed protein-interface similarities."""
    dataset = (release or PlinderRelease()).fetch("interface_similarity_scores")
    return read_score_table(dataset, columns=columns, filters=filters)


__all__ = ["query_interface_similarity"]
