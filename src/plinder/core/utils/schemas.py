# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pyarrow as pa

MAPPED_ALIGNMENT_REQUIRED_COLUMNS = frozenset(
    {
        "query_entry",
        "target_entry",
        "query_chain_mapped",
        "target_chain_mapped",
        "source",
        "qcov",
        "fident",
        "seqsim",
        "query_selected_residue_numbers",
        "target_selected_residue_numbers",
        "selected_residue_identity",
    }
)


def mapped_alignment_schema(*, alignment_type: str) -> pa.Schema:
    """Return the compact release schema, including typed empty shards."""
    fields = [
        pa.field("query_entry", pa.string()),
        pa.field("target_entry", pa.string()),
        pa.field("query_chain_mapped", pa.string()),
        pa.field("target_chain_mapped", pa.string()),
        pa.field("source", pa.string()),
        pa.field("qcov", pa.float64()),
        pa.field("fident", pa.float64()),
        pa.field("seqsim", pa.float64()),
        pa.field("query_selected_residue_numbers", pa.list_(pa.int32())),
        pa.field("target_selected_residue_numbers", pa.list_(pa.int32())),
        pa.field("selected_residue_identity", pa.binary()),
    ]
    if alignment_type == "foldseek":
        fields.append(pa.field("lddt", pa.float64()))
    elif alignment_type != "mmseqs":
        raise ValueError(f"unknown alignment type: {alignment_type}")
    return pa.schema(fields)


def mapped_alignment_schema_is_current(
    columns: set[str], *, alignment_type: str
) -> bool:
    """Return whether a mapped alignment has the compact V3 schema."""
    required = MAPPED_ALIGNMENT_REQUIRED_COLUMNS
    if alignment_type == "foldseek":
        required = required | {"lddt"}
    elif alignment_type != "mmseqs":
        raise ValueError(f"unknown alignment type: {alignment_type}")
    return required.issubset(columns)


PROTEIN_SIMILARITY_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("query_ligand_id", pa.string()),
        ("target_system", pa.string()),
        ("target_ligand_id", pa.string()),
        ("protein_mapping", pa.string()),
        ("mapping", pa.string()),
        ("protein_mapper", pa.dictionary(pa.int8(), pa.string())),
        ("source", pa.dictionary(pa.int8(), pa.string(), ordered=True)),
        ("metric", pa.dictionary(pa.int8(), pa.string(), ordered=True)),
        ("similarity", pa.int8()),
    ]
)

INTERFACE_SIMILARITY_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("target_system", pa.string()),
        ("mapping", pa.string()),
        ("source", pa.dictionary(pa.int8(), pa.string(), ordered=True)),
        ("metric", pa.dictionary(pa.int8(), pa.string(), ordered=True)),
        ("iface1_qcov", pa.float32()),
        ("iface2_qcov", pa.float32()),
        ("similarity", pa.int8()),
    ]
)

INTERFACE_SCORE_SHARD_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("target_system", pa.string()),
        ("mapping", pa.string()),
        ("source", pa.string()),
        ("metric", pa.string()),
        ("iface1_qcov", pa.float32()),
        ("iface2_qcov", pa.float32()),
        ("similarity", pa.int8()),
    ]
)

INTERFACE_QCOV_EXPORT_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("target_system", pa.string()),
        ("iface1_qcov", pa.float32()),
        ("iface2_qcov", pa.float32()),
        ("similarity", pa.int8()),
    ]
)

LIGAND_3D_CANDIDATE_SCHEMA = pa.schema(
    [
        ("query_system", pa.string()),
        ("query_ligand_id", pa.string()),
        ("query_entry", pa.string()),
        ("query_ligand_asym_id", pa.string()),
        ("target_system", pa.string()),
        ("target_ligand_id", pa.string()),
        ("target_entry", pa.string()),
        ("target_ligand_asym_id", pa.string()),
        ("protein_mapping", pa.string()),
        ("protein_mapper", pa.string()),
        ("pocket_qcov", pa.float64()),
    ]
)

LIGAND_3D_PAIR_CANDIDATE_SCHEMA = pa.schema(
    [
        ("query_entry", pa.string()),
        ("query_ligand_asym_id", pa.string()),
        ("target_entry", pa.string()),
        ("target_ligand_asym_id", pa.string()),
        ("pocket_qcov", pa.float64()),
    ]
)

LIGAND_3D_SCORE_SCHEMA = pa.schema(
    [
        ("query_entry", pa.string()),
        ("query_ligand_asym_id", pa.string()),
        ("target_entry", pa.string()),
        ("target_ligand_asym_id", pa.string()),
        ("shape", pa.float64()),
        ("color", pa.float64()),
        ("sucos_shape", pa.float64()),
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


LIGAND_CLUSTER_SCHEMA = pa.schema(
    [
        ("ligand_id", pa.string()),
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
        pa.field("tanimoto_similarity_ecfp4_1024", pa.float32()),
    ]
)

LEGACY_TANIMOTO_SCORE_SCHEMA = pa.schema(
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
