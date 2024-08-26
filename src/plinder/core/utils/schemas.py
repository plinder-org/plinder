# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pyarrow as pa

PROTEIN_SIMILARITY_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("target_system", pa.string()),
        ("protein_mapping", pa.string()),
        ("mapping", pa.string()),
        ("protein_mapper", pa.dictionary(pa.int8(), pa.string())),
        ("source", pa.dictionary(pa.int8(), pa.string(), ordered=True)),
        ("metric", pa.dictionary(pa.int8(), pa.string(), ordered=True)),
        ("similarity", pa.int8()),
    ]
)


NETWORKX_CLUSTER_SCHEMA = pa.schema(
    [
        ("system_id", pa.dictionary(pa.int32(), pa.string())),
        ("component", pa.dictionary(pa.int32(), pa.string())),
        ("community", pa.dictionary(pa.int32(), pa.string())),
    ]
)


GRAPHTOOL_CLUSTER_SCHEMA = pa.schema(
    [
        ("system_id", pa.string()),
        ("component", pa.dictionary(pa.int32(), pa.string())),
        ("metric", pa.dictionary(pa.int32(), pa.string())),
        ("directed", pa.dictionary(pa.int32(), pa.string())),
        ("threshold", pa.int8()),
    ]
)


CLUSTER_SCHEMA = pa.schema(
    [
        ("system_id", pa.string()),
        ("label", pa.string()),
        ("metric", pa.string()),
        ("cluster", pa.string()),
        ("directed", pa.bool_()),
        ("threshold", pa.int8()),
    ]
)


TANIMOTO_SCORE_SCHEMA = pa.schema(
    [
        pa.field("query_ligand_id", pa.int32()),
        pa.field("target_ligand_id", pa.int32()),
        pa.field("tanimoto_similarity_max", pa.int8()),
    ]
)


CLUSTER_DATASET_SCHEMA = pa.schema(
    [
        ("metric", pa.string()),
        ("directed", pa.bool_()),
        ("threshold", pa.int8()),
        ("system_id", pa.string()),
        ("component", pa.string()),
    ]
)


SPLIT_DATASET_SCHEMA = pa.schema(
    [
        ("system_id", pa.string()),
        ("split", pa.string()),
        ("cluster", pa.string()),
        ("cluster_for_val_split", pa.string()),
    ]
)


# subject to criteria used in save_linked_structures.py
# TODO: this schema is now out of date since addition of
#       scores.json contents but it now contains >50 columns
STRUCTURE_LINK_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("target_system", pa.string()),
        ("protein_qcov_weighted_sum", pa.float32()),
        ("protein_fident_weighted_sum", pa.float32()),
        ("pocket_fident", pa.float32()),
        ("target_id", pa.string()),
        ("sort_score", pa.float32()),
    ]
)
