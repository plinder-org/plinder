# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
The plinder.core.scores subpackage provides a consistent API for querying
the various parquet collections in the PLINDER dataset. The preferred
parquet reader engine is duckdb, but much of the code previously used
pandas and pyarrow directly. The internal query API supports converting
the same pyarrow query filters used in pd.read_parquet into raw SQL for
duckdb to execute.
"""
from .clusters import query_clusters
from .index import query_index
from .ligand import cross_similarity as cross_ligand_similarity
from .ligand import query_ligand_similarity
from .links import query_links
from .protein import (
    cross_similarity as cross_protein_similarity,
)
from .protein import (
    query_protein_similarity,
)

__all__ = [
    "query_ligand_similarity",
    "cross_ligand_similarity",
    "query_protein_similarity",
    "cross_protein_similarity",
    "query_clusters",
    "query_links",
    "query_index",
]
