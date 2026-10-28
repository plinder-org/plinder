# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Manifest-directed downloads for partitioned score exports.

Only positive equality/membership predicates are used for routing. All original
predicates are still applied by Arrow, including every branch of an OR filter.
"""

from __future__ import annotations

import base64
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from plinder.core.release import PlinderRelease
from plinder.core.utils import cpl

# Query and target aliases that share the same entry (or numeric-ID) bucket.
SCORE_KEYS: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {
    "ligand_similarity_scores": (
        ("query_system", "query_ligand_id"),
        ("target_system", "target_ligand_id"),
    ),
    "interface_similarity_scores": (("query_system",), ("target_system",)),
    "interface_half_similarity_scores": (
        ("query_half_interface_id",),
        ("target_half_interface_id",),
    ),
    "protein_similarity_scores": (("query_entry",), ("target_entry",)),
    "monomer_similarity_scores": (("query_entry",), ("target_entry",)),
    "ligand_scores": (("query_ligand_id",), ("target_ligand_id",)),
}


def score_bucket(value: Any, *, numeric: bool = False) -> str:
    """Stable entry-prefix buckets; numeric ligand IDs use ranges of 1024."""
    if numeric:
        return f"n{int(value) // 1024}"
    text = str(value)
    if re.match(r"^[a-zA-Z0-9]{4}(?:__|$)", text):
        return "p" + text[1:3].lower()
    return "h" + hashlib.sha256(text.split("__", 1)[0].encode()).hexdigest()[:2]


def _file(root: Path, relative: str, release: PlinderRelease) -> Path:
    path = Path(relative)
    if (
        path.is_absolute()
        or "\\" in relative
        or any(part in {"", ".", ".."} for part in relative.split("/"))
    ):
        raise ValueError(f"unsafe score-store path: {relative!r}")
    local = root / path
    local.resolve().relative_to(root.resolve())
    return local


def _files(root: Path, relatives: list[str], release: PlinderRelease) -> list[Path]:
    paths = [_file(root, relative, release) for relative in relatives]
    if release.data_dir is None:
        cpl.download_paths(paths=paths)
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"missing score-store file: {path}")
    return paths


def _to_pandas(table: pa.Table, partition_dictionaries: dict[str, str]) -> pd.DataFrame:
    """Preserve source Hive categories, including those absent from the rows."""
    result = table.to_pandas()
    for name, encoded in partition_dictionaries.items():
        if name not in result:
            continue
        dictionary = pa.ipc.open_stream(base64.b64decode(encoded)).read_all()
        result[name] = pd.Categorical(
            result[name],
            categories=dictionary.column(0).to_pandas().rename(None),
            ordered=table.schema.field(name).type.ordered,
        )
    return result


def read_partitioned_scores(
    release: PlinderRelease,
    artifact: str,
    *,
    columns: list[str] | None,
    filters: Any,
) -> pd.DataFrame | None:
    """Return None only when this release has no partitioned export yet."""
    try:
        manifest_path = release.fetch("score_store_manifest", artifact=artifact)
    except FileNotFoundError:
        return None
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("version") != 1 or manifest.get("artifact") != artifact:
        raise ValueError(f"unsupported score-store manifest: {manifest_path}")
    schema = pa.ipc.read_schema(pa.BufferReader(base64.b64decode(manifest["schema"])))
    partition_dictionaries = manifest.get("partition_dictionaries", {})
    # Bind filters/columns against the schema even if no shard matches.
    expression = pq.filters_to_expression(filters)
    empty = ds.dataset(pa.Table.from_batches([], schema=schema))
    empty_result = _to_pandas(
        empty.to_table(columns=columns, filter=expression), partition_dictionaries
    )
    query_keys, target_keys = SCORE_KEYS[artifact]
    numeric = artifact == "ligand_scores"
    data: dict[str, list[str]] = manifest["data"]
    targets: dict[str, list[str]] = manifest["targets"]
    branches = [filters] if isinstance(filters[0], tuple) else filters
    selected: set[str] = set()
    index_cache: dict[str, pd.DataFrame] = {}
    for branch in branches:
        buckets = set(data)
        for column, operator, value in branch:
            if column not in query_keys or operator not in {"=", "==", "in"}:
                continue
            values = value if operator == "in" else [value]
            buckets &= {score_bucket(v, numeric=numeric) for v in values}
        for column, operator, value in branch:
            if not buckets or len(buckets) == 1:
                break  # The original Arrow filter handles target predicates.
            if column not in target_keys or operator not in {"=", "==", "in"}:
                continue
            values = value if operator == "in" else [value]
            strings = {str(int(v)) if numeric else str(v) for v in values}
            matching: set[str] = set()
            for bucket in {score_bucket(v, numeric=numeric) for v in values}:
                if bucket not in index_cache:
                    paths = _files(
                        manifest_path.parent, targets.get(bucket, []), release
                    )
                    index_cache[bucket] = (
                        pd.read_parquet(paths)
                        if paths
                        else pd.DataFrame(columns=["column", "value", "shard"])
                    )
                index = index_cache[bucket]
                matching.update(
                    index.loc[
                        (index["column"] == column) & index["value"].isin(strings),
                        "shard",
                    ]
                )
            buckets &= matching
        selected.update(buckets)
    relatives = [path for bucket in sorted(selected) for path in data[bucket]]
    paths = _files(manifest_path.parent, relatives, release)
    if not paths:
        return empty_result
    return _to_pandas(
        ds.dataset(paths, format="parquet", schema=schema).to_table(
            columns=columns, filter=expression
        ),
        partition_dictionaries,
    )
