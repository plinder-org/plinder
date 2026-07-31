# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Plan, shard, validate, and install the local V3 annotation index."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterable, Sequence, cast
from uuid import uuid4

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.data.annotations.interface_utils import (
    INTERFACE_ANNOTATION_SCHEMA,
    min_interface_residues_from_schema,
)
from plinder.data.pipeline.ingest import completed_entry_metrics

COLLATION_VERSION = 2
STAGING_RELATIVE = Path("index/.staging/v3_collation")
MANIFEST_NAME = "entries.parquet"
PLAN_NAME = "plan.json"
FINAL_MARKER_NAME = "collation.json"
REPAIR_REQUIRED_STATUS = "requires_downstream_repair"
RETIRED_ENRICHMENT_MARKERS = ("ecod", "panther", "kinase")
SYSTEM_LIGAND_FLAGS = (
    "lipinski",
    "cofactor",
    "fragment",
    "monosaccharide",
    "oligosaccharide",
    "mononucleotide",
    "oligonucleotide",
    "monopeptide",
    "oligopeptide",
    "artifact",
    "other",
    "covalent",
    "invalid",
    "ion",
)
AGGREGATED_COLUMNS = (
    "biounit_num_ligands",
    "biounit_num_unique_ccd_codes",
    "biounit_num_proper_ligands",
    *(f"system_ligand_has_{name}" for name in SYSTEM_LIGAND_FLAGS),
    "system_protein_chains_total_length",
    "system_unique_ccd_codes",
    "system_proper_unique_ccd_codes",
    "ligand_is_3d_score_able",
)

MANIFEST_SCHEMA = pa.schema(
    [
        ("pdb_id", pa.string()),
        ("code", pa.string()),
        ("annotation_path", pa.string()),
        ("chain_path", pa.string()),
        ("biounit_chain_path", pa.string()),
        ("metadata_path", pa.string()),
        ("interface_path", pa.string()),
        ("source_path", pa.string()),
        ("ligand_path", pa.string()),
        ("annotation_size", pa.int64()),
        ("annotation_mtime_ns", pa.int64()),
        ("chain_size", pa.int64()),
        ("chain_mtime_ns", pa.int64()),
        ("biounit_chain_size", pa.int64()),
        ("biounit_chain_mtime_ns", pa.int64()),
        ("metadata_size", pa.int64()),
        ("metadata_mtime_ns", pa.int64()),
        ("interface_size", pa.int64()),
        ("interface_mtime_ns", pa.int64()),
        ("source_size", pa.int64()),
        ("source_mtime_ns", pa.int64()),
        ("ligand_size", pa.int64()),
        ("ligand_mtime_ns", pa.int64()),
        ("interface_min_residues", pa.int64()),
    ]
)

ENTRY_CHAIN_SCHEMA = pa.schema(
    [
        ("entry_pdb_id", pa.string()),
        ("chain_asym_id", pa.string()),
        ("chain_auth_id", pa.string()),
        ("chain_entity_id", pa.string()),
        ("chain_type", pa.string()),
        ("chain_receptor_type", pa.string()),
        ("chain_length", pa.int64()),
        ("chain_num_unresolved_residues", pa.int64()),
        ("chain_is_holo", pa.bool_()),
        ("chain_is_ligand_like", pa.bool_()),
        ("chain_uniprot_ids", pa.list_(pa.string())),
    ]
)

BIOUNIT_CHAIN_SCHEMA = pa.schema(
    [
        ("entry_pdb_id", pa.string()),
        ("biounit_id", pa.string()),
        ("chain_instance", pa.string()),
        ("chain_asym_id", pa.string()),
        ("chain_role", pa.string()),
    ]
)

ENTRY_SOURCE_SCHEMA = pa.schema(
    [
        ("entry_pdb_id", pa.string()),
        ("source_mmcif_major_revision", pa.int64()),
        ("source_mmcif_minor_revision", pa.int64()),
    ]
)

SIDECAR_SCHEMAS = {
    "entry_chains": ENTRY_CHAIN_SCHEMA,
    "entry_biounit_chains": BIOUNIT_CHAIN_SCHEMA,
    "interfaces": INTERFACE_ANNOTATION_SCHEMA,
    "entry_sources": ENTRY_SOURCE_SCHEMA,
}


def staging_dir(data_dir: Path) -> Path:
    """Return the private collation workspace under the release directory."""
    return data_dir.resolve() / STAGING_RELATIVE


def manifest_path(data_dir: Path) -> Path:
    return staging_dir(data_dir) / MANIFEST_NAME


def plan_path(data_dir: Path) -> Path:
    return staging_dir(data_dir) / PLAN_NAME


def _temporary_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.{os.getpid()}.{uuid4().hex}.tmp")


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_table_atomic(
    table: pa.Table,
    path: Path,
    *,
    row_group_size: int = 100_000,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        pq.write_table(
            table,
            temporary,
            compression="zstd",
            row_group_size=row_group_size,
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _stat_values(path: Path, prefix: str) -> dict[str, int]:
    stat = path.stat()
    return {
        f"{prefix}_size": stat.st_size,
        f"{prefix}_mtime_ns": stat.st_mtime_ns,
    }


def _row_signature(rows: Sequence[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in sorted(rows, key=lambda item: str(item["pdb_id"])):
        digest.update(
            json.dumps(row, sort_keys=True, separators=(",", ":")).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _entry_manifest_row(data_dir: Path, entry_dir: Path) -> dict[str, Any]:
    pdb_id = entry_dir.name.lower()
    if re.fullmatch(r"[0-9][a-z0-9]{3}", pdb_id) is None:
        raise ValueError(f"invalid raw-entry directory: {entry_dir}")
    code = pdb_id[1:3]
    if entry_dir.parent.name.lower() != code:
        raise ValueError(f"raw-entry shard does not match PDB ID: {entry_dir}")
    annotation_path = entry_dir.parent / f"{pdb_id}.parquet"
    ligand_path = data_dir / "ligands" / f"{pdb_id}.parquet"
    paths = {
        "annotation": annotation_path,
        "chain": entry_dir / "entry_chains.parquet",
        "biounit_chain": entry_dir / "entry_biounit_chains.parquet",
        "metadata": entry_dir / "entry_metadata.parquet",
        "interface": entry_dir / "interfaces.parquet",
        "source": entry_dir / "entry_source.parquet",
        "ligand": ligand_path,
    }
    optional = {"annotation", "ligand"}
    if annotation_path.is_file() != ligand_path.is_file():
        raise FileNotFoundError(
            f"incomplete V3 entry {pdb_id}; ligand outputs disagree: "
            f"annotation={annotation_path.is_file()}, ligand={ligand_path.is_file()}"
        )
    path_stats: dict[str, os.stat_result] = {}
    missing: list[Path] = []
    for name, path in paths.items():
        try:
            path_stat = path.stat()
        except FileNotFoundError:
            if name not in optional:
                missing.append(path)
            continue
        if not stat.S_ISREG(path_stat.st_mode):
            if name not in optional:
                missing.append(path)
            continue
        path_stats[name] = path_stat
    if missing:
        formatted = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"incomplete V3 entry {pdb_id}; missing: {formatted}")
    if annotation_path.is_file():
        if pq.ParquetFile(annotation_path).metadata.num_rows < 1:
            raise ValueError(f"ligand annotation is empty for V3 entry {pdb_id}")
    else:
        if pq.ParquetFile(paths["interface"]).metadata.num_rows < 1:
            raise ValueError(
                f"materialized V3 entry {pdb_id} has no ligand or interface rows"
            )
        if completed_entry_metrics(data_dir, pdb_id) is None:
            raise ValueError(
                f"interface-only V3 entry {pdb_id} has no successful ingest marker"
            )
    row: dict[str, Any] = {
        "pdb_id": pdb_id,
        "code": code,
        "interface_min_residues": min_interface_residues_from_schema(
            pq.read_schema(paths["interface"])
        ),
        **{f"{name}_path": str(path.resolve()) for name, path in paths.items()},
    }
    for name in paths:
        row[f"{name}_size"] = None
        row[f"{name}_mtime_ns"] = None
    for name, path_stat in path_stats.items():
        row[f"{name}_size"] = path_stat.st_size
        row[f"{name}_mtime_ns"] = path_stat.st_mtime_ns
    return row


def plan_collation(data_dir: Path, *, threads: int = 1) -> dict[str, Any]:
    """Inventory every materialized V3 entry and atomically publish a plan."""
    if threads < 1:
        raise ValueError("planning threads must be positive")
    data_dir = data_dir.resolve()
    raw_entries = data_dir / "raw_entries"
    if not raw_entries.is_dir():
        raise FileNotFoundError(f"missing raw-entry dataset: {raw_entries}")
    output = manifest_path(data_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(output)
    writer: Any | None = None
    seen: set[str] = set()
    code_signatures: dict[str, str] = {}
    code_counts: dict[str, int] = {}
    interface_min_residues: set[int] = set()
    total_entries = 0
    try:
        writer = pq.ParquetWriter(temporary, MANIFEST_SCHEMA, compression="zstd")
        code_dirs = sorted(
            path
            for path in raw_entries.iterdir()
            if path.is_dir() and re.fullmatch(r"[a-z0-9]{2}", path.name.lower())
        )
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for code_dir in code_dirs:
                entry_dirs = sorted(
                    path
                    for path in code_dir.iterdir()
                    if path.is_dir()
                    and re.fullmatch(r"[0-9][a-z0-9]{3}", path.name.lower())
                )
                rows = list(
                    executor.map(
                        lambda path: _entry_manifest_row(data_dir, path),
                        entry_dirs,
                    )
                )
                if not rows:
                    continue
                duplicates = sorted(
                    str(row["pdb_id"]) for row in rows if str(row["pdb_id"]) in seen
                )
                if duplicates:
                    raise ValueError(f"duplicate raw-entry annotations: {duplicates}")
                seen.update(str(row["pdb_id"]) for row in rows)
                code = code_dir.name.lower()
                code_signatures[code] = _row_signature(rows)
                code_counts[code] = len(rows)
                interface_min_residues.update(
                    int(row["interface_min_residues"]) for row in rows
                )
                total_entries += len(rows)
                writer.write_table(pa.Table.from_pylist(rows, schema=MANIFEST_SCHEMA))
        if total_entries == 0:
            raise ValueError(f"no materialized V3 entries found in {raw_entries}")
        if len(interface_min_residues) != 1:
            raise ValueError(
                "mixed interface.min_interface_residues values in ingest outputs: "
                f"{sorted(interface_min_residues)}"
            )
        writer.close()
        writer = None
        temporary.replace(output)
    finally:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)

    summary: dict[str, Any] = {
        "version": COLLATION_VERSION,
        "status": "complete",
        "data_dir": str(data_dir),
        "manifest": str(output),
        "entry_count": total_entries,
        "code_count": len(code_signatures),
        "codes": sorted(code_signatures),
        "code_entry_counts": code_counts,
        "code_signatures": code_signatures,
        "interface_min_residues": interface_min_residues.pop(),
    }
    _write_json_atomic(plan_path(data_dir), summary)
    return summary


def load_plan(data_dir: Path) -> dict[str, Any]:
    """Load a completed plan and verify its manifest is present."""
    path = plan_path(data_dir)
    if not path.is_file():
        raise FileNotFoundError(f"missing collation plan: {path}")
    plan = cast(dict[str, Any], json.loads(path.read_text()))
    if plan.get("status") != "complete" or plan.get("version") != COLLATION_VERSION:
        raise ValueError(f"invalid collation plan: {path}")
    if int(plan.get("interface_min_residues", 0)) < 1:
        raise ValueError(f"invalid interface threshold in collation plan: {path}")
    if not manifest_path(data_dir).is_file():
        raise FileNotFoundError(
            f"missing collation manifest: {manifest_path(data_dir)}"
        )
    return plan


def _load_manifest_rows(data_dir: Path, code: str) -> list[dict[str, Any]]:
    normalized = code.lower()
    if re.fullmatch(r"[a-z0-9]{2}", normalized) is None:
        raise ValueError(f"invalid two-character code: {code!r}")
    table = pq.read_table(
        manifest_path(data_dir),
        filters=[("code", "=", normalized)],
    )
    return cast(list[dict[str, Any]], table.to_pylist())


def _verify_manifest_row(row: dict[str, Any]) -> None:
    for prefix in (
        "annotation",
        "chain",
        "biounit_chain",
        "metadata",
        "interface",
        "source",
        "ligand",
    ):
        path = Path(str(row[f"{prefix}_path"]))
        if row[f"{prefix}_size"] is None:
            if path.exists():
                raise RuntimeError(
                    "collation input appeared after planning; rerun the plan "
                    f"stage: {path}"
                )
            continue
        try:
            values = _stat_values(path, prefix)
        except FileNotFoundError as exc:
            raise RuntimeError(
                f"collation input disappeared after planning: {path}"
            ) from exc
        if any(values[key] != row[key] for key in values):
            raise RuntimeError(
                "collation input changed after planning; rerun the plan stage: "
                f"{path}"
            )


def _verify_manifest_inputs(
    rows: Sequence[dict[str, Any]], *, threads: int = 1
) -> None:
    if threads < 1:
        raise ValueError("verification threads must be positive")
    if threads == 1:
        for row in rows:
            _verify_manifest_row(row)
        return
    batch_size = threads * 64
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for start in range(0, len(rows), batch_size):
            for _ in executor.map(
                _verify_manifest_row, rows[start : start + batch_size]
            ):
                pass


def _configure_duckdb(
    connection: duckdb.DuckDBPyConnection,
    *,
    threads: int,
    memory_limit: str,
    scratch_dir: Path | None,
) -> None:
    if threads < 1:
        raise ValueError("threads must be positive")
    connection.execute(f"SET threads = {threads:d}")
    connection.execute("SET preserve_insertion_order = false")
    escaped_memory = memory_limit.replace("'", "''")
    connection.execute(f"SET memory_limit = '{escaped_memory}'")
    if scratch_dir is not None:
        scratch_dir.mkdir(parents=True, exist_ok=True)
        escaped_scratch = str(scratch_dir.resolve()).replace("'", "''")
        connection.execute(f"SET temp_directory = '{escaped_scratch}'")


def _fetch_scalar(connection: duckdb.DuckDBPyConnection, query: str) -> Any:
    """Execute a query that must return exactly one scalar row."""
    row = connection.execute(query).fetchone()
    if row is None:
        raise RuntimeError(f"query unexpectedly returned no rows: {query}")
    return row[0]


def _quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _copy_query(
    connection: duckdb.DuckDBPyConnection,
    query: str,
    path: Path,
    *,
    row_group_size: int,
) -> None:
    escaped = str(path.resolve()).replace("'", "''")
    connection.execute(
        f"COPY ({query}) TO '{escaped}' "
        f"(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE {row_group_size:d})"
    )


def _copy_query_atomic(
    connection: duckdb.DuckDBPyConnection,
    query: str,
    path: Path,
    *,
    row_group_size: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        _copy_query(connection, query, temporary, row_group_size=row_group_size)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _relation_columns(
    connection: duckdb.DuckDBPyConnection, relation: str
) -> list[str]:
    rows = connection.execute(f"PRAGMA table_info('{relation}')").fetchall()
    return [str(row[1]) for row in rows]


def _require_columns(
    columns: Iterable[str], required: Iterable[str], label: str
) -> None:
    missing = sorted(set(required).difference(columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def _build_annotation_view(
    connection: duckdb.DuckDBPyConnection,
    annotation_paths: Sequence[str],
    ligand_paths: Sequence[str],
    *,
    empty: bool = False,
) -> None:
    annotation = connection.read_parquet(list(annotation_paths), union_by_name=True)
    ligand = connection.read_parquet(list(ligand_paths), union_by_name=True)
    if empty:
        annotation = annotation.filter("false")
        ligand = ligand.filter("false")
    annotation.create_view("raw_annotation", replace=True)
    ligand.create_view("raw_ligand_ability", replace=True)
    raw_columns = _relation_columns(connection, "raw_annotation")
    required = {
        "entry_pdb_id",
        "system_biounit_id",
        "system_id",
        "system_type",
        "system_protein_chains_length",
        "ligand_id",
        "ligand_unique_ccd_code",
        "ligand_is_proper",
        *(f"ligand_is_{name}" for name in SYSTEM_LIGAND_FLAGS),
    }
    _require_columns(raw_columns, required, "raw annotation")
    retained = [
        column
        for column in raw_columns
        if column not in AGGREGATED_COLUMNS
        and not any(
            marker in column.casefold() for marker in RETIRED_ENRICHMENT_MARKERS
        )
    ]
    selected = ", ".join(_quote_identifier(column) for column in retained)
    connection.execute(
        f"CREATE OR REPLACE TEMP VIEW base_annotation AS SELECT {selected} "
        "FROM raw_annotation"
    )
    ability_columns = _relation_columns(connection, "raw_ligand_ability")
    _require_columns(
        ability_columns,
        {"ligand_id", "ligand_is_3d_score_able"},
        "ligand scoreability",
    )
    conflicts = connection.execute(
        "SELECT ligand_id FROM raw_ligand_ability GROUP BY ligand_id "
        "HAVING count(DISTINCT ligand_is_3d_score_able) > 1 LIMIT 10"
    ).fetchall()
    if conflicts:
        raise ValueError(f"conflicting 3D scoreability values: {conflicts}")
    connection.execute(
        "CREATE OR REPLACE TEMP VIEW ligand_ability AS "
        "SELECT ligand_id, bool_or(ligand_is_3d_score_able) "
        "AS ligand_is_3d_score_able FROM raw_ligand_ability GROUP BY ligand_id"
    )
    missing_abilities = _fetch_scalar(
        connection,
        "SELECT count(*) FROM base_annotation AS b "
        "LEFT JOIN ligand_ability AS a USING (ligand_id) "
        "WHERE b.system_type = 'holo' AND coalesce(b.ligand_is_proper, false) "
        "AND a.ligand_id IS NULL",
    )
    if missing_abilities:
        raise ValueError(
            f"{missing_abilities} proper holo rows lack 3D-scoreability annotations"
        )

    connection.execute(
        "CREATE OR REPLACE TEMP VIEW biounit_aggregate AS "
        "SELECT entry_pdb_id, system_biounit_id, "
        "count(system_id)::BIGINT AS biounit_num_ligands, "
        "count(DISTINCT ligand_unique_ccd_code)::BIGINT "
        "AS biounit_num_unique_ccd_codes, "
        "sum(CASE WHEN coalesce(ligand_is_proper, false) THEN 1 ELSE 0 END)::BIGINT "
        "AS biounit_num_proper_ligands "
        "FROM base_annotation GROUP BY entry_pdb_id, system_biounit_id"
    )
    flag_sql = ", ".join(
        f"bool_or(ligand_is_{name}) AS system_ligand_has_{name}"
        for name in SYSTEM_LIGAND_FLAGS
    )
    connection.execute(
        "CREATE OR REPLACE TEMP VIEW system_aggregate AS "
        f"SELECT system_id, {flag_sql}, "
        "string_agg(DISTINCT ligand_unique_ccd_code, '-' "
        "ORDER BY ligand_unique_ccd_code) AS system_unique_ccd_codes "
        "FROM base_annotation GROUP BY system_id"
    )
    connection.execute(
        "CREATE OR REPLACE TEMP VIEW proper_system_aggregate AS "
        "SELECT system_id, string_agg(DISTINCT ligand_unique_ccd_code, '-' "
        "ORDER BY ligand_unique_ccd_code) AS system_proper_unique_ccd_codes "
        "FROM base_annotation WHERE coalesce(ligand_is_proper, false) "
        "GROUP BY system_id"
    )
    raw_select = ", ".join(f"b.{_quote_identifier(column)}" for column in retained)
    flag_select = ", ".join(
        f"s.system_ligand_has_{name}" for name in SYSTEM_LIGAND_FLAGS
    )
    connection.execute(
        "CREATE OR REPLACE TEMP VIEW collated_annotation AS SELECT "
        f"{raw_select}, u.biounit_num_ligands, "
        "u.biounit_num_unique_ccd_codes, u.biounit_num_proper_ligands, "
        f"{flag_select}, list_sum(b.system_protein_chains_length)::BIGINT "
        "AS system_protein_chains_total_length, s.system_unique_ccd_codes, "
        "p.system_proper_unique_ccd_codes, "
        "CASE WHEN b.system_type = 'holo' AND coalesce(b.ligand_is_proper, false) "
        "THEN a.ligand_is_3d_score_able ELSE false END "
        "AS ligand_is_3d_score_able FROM base_annotation AS b "
        "JOIN biounit_aggregate AS u USING (entry_pdb_id, system_biounit_id) "
        "JOIN system_aggregate AS s USING (system_id) "
        "LEFT JOIN proper_system_aggregate AS p USING (system_id) "
        "LEFT JOIN ligand_ability AS a USING (ligand_id)"
    )


def _normalize_table(path: Path, schema: pa.Schema) -> pa.Table:
    table = pq.read_table(path)
    if table.num_rows and not set(schema.names).issubset(table.column_names):
        missing = sorted(set(schema.names).difference(table.column_names))
        raise ValueError(f"{path} is missing sidecar columns: {missing}")
    arrays = []
    for field in schema:
        if field.name in table.column_names:
            column = table[field.name].combine_chunks().cast(field.type, safe=False)
        else:
            column = pa.nulls(table.num_rows, type=field.type)
        arrays.append(column)
    return pa.Table.from_arrays(arrays, schema=schema)


def _collate_sidecars(
    paths: Sequence[Path],
    output: Path,
    schema: pa.Schema,
    *,
    sort_columns: Sequence[str],
    row_group_size: int,
) -> None:
    tables = [_normalize_table(path, schema) for path in paths]
    table = (
        pa.concat_tables(tables) if tables else pa.Table.from_batches([], schema=schema)
    )
    if table.num_rows:
        table = table.sort_by([(column, "ascending") for column in sort_columns])
    _write_table_atomic(table, output, row_group_size=row_group_size)


def _collate_entry_metadata(
    paths: Sequence[Path],
    output: Path,
    *,
    row_group_size: int,
) -> None:
    """Collate the evolving entry schema by column name."""
    tables = [pq.read_table(path) for path in paths]
    if not tables:
        raise ValueError("cannot collate an empty entry-metadata shard")
    table = pa.concat_tables(tables, promote_options="default")
    if "entry_pdb_id" not in table.column_names:
        raise ValueError("entry metadata is missing entry_pdb_id")
    table = table.sort_by([("entry_pdb_id", "ascending")])
    _write_table_atomic(table, output, row_group_size=row_group_size)


def _shard_paths(data_dir: Path, code: str) -> dict[str, Path]:
    root = staging_dir(data_dir)
    return {
        "annotation": root / "annotations" / f"{code}.parquet",
        "entry_chains": root / "entry_chains" / f"{code}.parquet",
        "entry_biounit_chains": root / "entry_biounit_chains" / f"{code}.parquet",
        "entry_metadata": root / "entry_metadata" / f"{code}.parquet",
        "interfaces": root / "interfaces" / f"{code}.parquet",
        "entry_sources": root / "entry_sources" / f"{code}.parquet",
        "metrics": root / "shards" / f"{code}.json",
    }


def _completed_shard(
    paths: dict[str, Path], *, signature: str
) -> dict[str, Any] | None:
    metrics_path = paths["metrics"]
    if not metrics_path.is_file():
        return None
    try:
        metrics = cast(dict[str, Any], json.loads(metrics_path.read_text()))
    except (OSError, json.JSONDecodeError):
        return None
    outputs = [path for name, path in paths.items() if name != "metrics"]
    if (
        metrics.get("status") == "complete"
        and metrics.get("signature") == signature
        and all(path.is_file() for path in outputs)
    ):
        return metrics
    return None


def collate_shard(
    data_dir: Path,
    code: str,
    *,
    force: bool = False,
    threads: int = 1,
    memory_limit: str = "8GB",
    scratch_dir: Path | None = None,
    row_group_size: int = 100_000,
) -> dict[str, Any]:
    """Collate one two-character shard with bounded memory and atomic writes."""
    data_dir = data_dir.resolve()
    plan = load_plan(data_dir)
    normalized = code.lower()
    if normalized not in plan["code_signatures"]:
        raise KeyError(f"two-character code is not in collation plan: {normalized}")
    rows = _load_manifest_rows(data_dir, normalized)
    signature = _row_signature(rows)
    if signature != plan["code_signatures"][normalized]:
        raise RuntimeError(f"manifest signature mismatch for shard {normalized}")
    _verify_manifest_inputs(rows, threads=threads)
    paths = _shard_paths(data_dir, normalized)
    if not force and (completed := _completed_shard(paths, signature=signature)):
        return completed
    metrics: dict[str, Any] = {
        "version": COLLATION_VERSION,
        "status": "running",
        "code": normalized,
        "signature": signature,
        "entry_count": len(rows),
        "outputs": {
            name: str(path) for name, path in paths.items() if name != "metrics"
        },
    }
    try:
        connection = duckdb.connect()
        try:
            shard_scratch = (
                scratch_dir / normalized if scratch_dir is not None else None
            )
            _configure_duckdb(
                connection,
                threads=threads,
                memory_limit=memory_limit,
                scratch_dir=shard_scratch,
            )
            ligand_rows = [row for row in rows if row["annotation_size"] is not None]
            if ligand_rows:
                schema_rows = ligand_rows
                empty_annotation = False
            else:
                schema_rows = [
                    row
                    for row in cast(
                        list[dict[str, Any]],
                        pq.read_table(manifest_path(data_dir)).to_pylist(),
                    )
                    if row["annotation_size"] is not None
                ][:1]
                if not schema_rows:
                    raise ValueError(
                        "interface collation requires at least one ligand-bearing "
                        "entry to define the ligand annotation schema"
                    )
                empty_annotation = True
            _build_annotation_view(
                connection,
                [str(row["annotation_path"]) for row in schema_rows],
                [str(row["ligand_path"]) for row in schema_rows],
                empty=empty_annotation,
            )
            _copy_query_atomic(
                connection,
                "SELECT * FROM collated_annotation "
                "ORDER BY entry_pdb_id, system_id, ligand_id",
                paths["annotation"],
                row_group_size=row_group_size,
            )
        finally:
            connection.close()
        _collate_sidecars(
            [Path(str(row["chain_path"])) for row in rows],
            paths["entry_chains"],
            ENTRY_CHAIN_SCHEMA,
            sort_columns=("entry_pdb_id", "chain_asym_id"),
            row_group_size=row_group_size,
        )
        _collate_entry_metadata(
            [Path(str(row["metadata_path"])) for row in rows],
            paths["entry_metadata"],
            row_group_size=row_group_size,
        )
        _collate_sidecars(
            [Path(str(row["interface_path"])) for row in rows],
            paths["interfaces"],
            INTERFACE_ANNOTATION_SCHEMA,
            sort_columns=("entry_pdb_id", "system_id"),
            row_group_size=row_group_size,
        )
        _collate_sidecars(
            [Path(str(row["biounit_chain_path"])) for row in rows],
            paths["entry_biounit_chains"],
            BIOUNIT_CHAIN_SCHEMA,
            sort_columns=("entry_pdb_id", "biounit_id", "chain_instance"),
            row_group_size=row_group_size,
        )
        _collate_sidecars(
            [Path(str(row["source_path"])) for row in rows],
            paths["entry_sources"],
            ENTRY_SOURCE_SCHEMA,
            sort_columns=("entry_pdb_id",),
            row_group_size=row_group_size,
        )
        metrics["counts"] = {
            name: pq.ParquetFile(path).metadata.num_rows
            for name, path in paths.items()
            if name != "metrics"
        }
        metrics["status"] = "complete"
    except BaseException as exc:
        metrics["status"] = "failed"
        metrics["error"] = repr(exc)
        metrics["traceback"] = traceback.format_exc()
        raise
    finally:
        _write_json_atomic(paths["metrics"], metrics)
    return metrics


def _load_completed_shards(
    data_dir: Path, plan: dict[str, Any], *, verify_threads: int
) -> tuple[dict[str, list[Path]], dict[str, int]]:
    manifest_rows = cast(
        list[dict[str, Any]],
        pq.read_table(manifest_path(data_dir)).to_pylist(),
    )
    rows_by_code: dict[str, list[dict[str, Any]]] = {}
    for row in manifest_rows:
        rows_by_code.setdefault(str(row["code"]), []).append(row)
    planned_codes = [str(code) for code in plan["codes"]]
    if set(rows_by_code) != set(planned_codes):
        raise RuntimeError("collation manifest codes differ from the completed plan")
    for code in planned_codes:
        rows = rows_by_code[code]
        if len(rows) != int(plan["code_entry_counts"][code]):
            raise RuntimeError(f"collation manifest row count changed for shard {code}")
        signature = _row_signature(rows)
        if signature != str(plan["code_signatures"][code]):
            raise RuntimeError(f"collation manifest signature changed for shard {code}")
    _verify_manifest_inputs(manifest_rows, threads=verify_threads)

    files: dict[str, list[Path]] = {
        "annotation": [],
        "entry_chains": [],
        "entry_biounit_chains": [],
        "entry_metadata": [],
        "interfaces": [],
        "entry_sources": [],
    }
    expected_counts = {name: 0 for name in files}
    for code in planned_codes:
        rows = rows_by_code[code]
        signature = _row_signature(rows)
        paths = _shard_paths(data_dir, code)
        metrics = _completed_shard(paths, signature=signature)
        if metrics is None:
            raise RuntimeError(f"collation shard is incomplete or stale: {code}")
        for name in files:
            files[name].append(paths[name])
            expected_counts[name] += int(metrics["counts"][name])
    return files, expected_counts


def _install_final_tables_fail_closed(
    temporary_paths: dict[str, Path],
    final_paths: dict[str, Path],
    marker_path: Path,
) -> None:
    """Install validated tables while making partial generations unreadable."""
    marker_path.unlink(missing_ok=True)
    final_paths["annotation"].unlink(missing_ok=True)
    try:
        for name in (
            "entry_chains",
            "entry_biounit_chains",
            "entry_metadata",
            "interfaces",
            "entry_sources",
        ):
            temporary_paths[name].replace(final_paths[name])
        temporary_paths["annotation"].replace(final_paths["annotation"])
    except BaseException:
        final_paths["annotation"].unlink(missing_ok=True)
        raise


def _validate_final_tables(
    paths: dict[str, Path],
    *,
    expected_counts: dict[str, int],
    min_interface_residues: int,
    threads: int,
    memory_limit: str,
    scratch_dir: Path | None,
) -> dict[str, Any]:
    if min_interface_residues < 1:
        raise ValueError("minimum interface residues must be positive")
    connection = duckdb.connect()
    try:
        _configure_duckdb(
            connection,
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
        )
        for name, path in paths.items():
            connection.read_parquet(str(path)).create_view(name, replace=True)
        actual_counts = {
            name: int(_fetch_scalar(connection, f"SELECT count(*) FROM {name}"))
            for name in paths
        }
        if actual_counts != expected_counts:
            raise ValueError(
                f"final collation row counts differ: {actual_counts} != {expected_counts}"
            )
        duplicate_queries = {
            "annotation": (
                "SELECT count(*) FROM (SELECT system_id, ligand_id FROM annotation "
                "GROUP BY system_id, ligand_id HAVING count(*) > 1)"
            ),
            "entry_chains": (
                "SELECT count(*) FROM (SELECT entry_pdb_id, chain_asym_id "
                "FROM entry_chains GROUP BY entry_pdb_id, chain_asym_id "
                "HAVING count(*) > 1)"
            ),
            "entry_biounit_chains": (
                "SELECT count(*) FROM (SELECT entry_pdb_id, biounit_id, "
                "chain_instance FROM entry_biounit_chains GROUP BY entry_pdb_id, "
                "biounit_id, chain_instance HAVING count(*) > 1)"
            ),
            "entry_metadata": (
                "SELECT count(*) FROM (SELECT entry_pdb_id FROM entry_metadata "
                "GROUP BY entry_pdb_id HAVING count(*) > 1)"
            ),
            "interfaces": (
                "SELECT count(*) FROM (SELECT system_id FROM interfaces "
                "GROUP BY system_id HAVING count(*) > 1)"
            ),
            "entry_sources": (
                "SELECT count(*) FROM (SELECT entry_pdb_id FROM entry_sources "
                "GROUP BY entry_pdb_id HAVING count(*) > 1)"
            ),
        }
        duplicates = {
            name: int(_fetch_scalar(connection, query))
            for name, query in duplicate_queries.items()
        }
        if any(duplicates.values()):
            raise ValueError(f"duplicate keys in final collation: {duplicates}")
        invalid_chain_metadata = {
            "nonpositive_lengths": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM entry_chains "
                    "WHERE chain_length IS NULL OR chain_length <= 0",
                )
            ),
            "negative_unresolved": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM entry_chains WHERE "
                    "chain_num_unresolved_residues IS NULL OR "
                    "chain_num_unresolved_residues < 0",
                )
            ),
            "unresolved_exceeds_length": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM entry_chains WHERE "
                    "chain_num_unresolved_residues > chain_length",
                )
            ),
        }
        if any(invalid_chain_metadata.values()):
            raise ValueError(
                "invalid entry-chain sequence metadata: " f"{invalid_chain_metadata}"
            )
        invalid_interfaces = {
            "system_id": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM interfaces WHERE system_id != "
                    "entry_pdb_id || '__' || system_biounit_id || '__' || "
                    "interface_chain_1 || '--' || interface_chain_2",
                )
            ),
            "same_chain": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM interfaces WHERE "
                    "interface_chain_1 >= interface_chain_2",
                )
            ),
            "short_side": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM interfaces WHERE "
                    "len(interface_chain_1_residue_numbers) < "
                    f"{min_interface_residues} OR "
                    "len(interface_chain_2_residue_numbers) < "
                    f"{min_interface_residues}",
                )
            ),
            "mapping_length": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(*) FROM interfaces WHERE "
                    "len(interface_chain_1_residue_numbers) != "
                    "len(interface_chain_1_residue_indices) OR "
                    "len(interface_chain_2_residue_numbers) != "
                    "len(interface_chain_2_residue_indices)",
                )
            ),
        }
        if any(invalid_interfaces.values()):
            raise ValueError(f"invalid protein interfaces: {invalid_interfaces}")
        annotation_columns = _relation_columns(connection, "annotation")
        retired = sorted(
            column
            for column in annotation_columns
            if any(marker in column.casefold() for marker in RETIRED_ENRICHMENT_MARKERS)
        )
        if retired:
            raise ValueError(f"retired enrichment columns remain: {retired}")
        _require_columns(
            annotation_columns,
            AGGREGATED_COLUMNS,
            "collated annotation",
        )
        orphan_sources = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM (SELECT DISTINCT a.entry_pdb_id FROM annotation a "
                "ANTI JOIN entry_sources s USING (entry_pdb_id))",
            )
        )
        orphan_chains = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM (SELECT DISTINCT a.entry_pdb_id FROM annotation a "
                "ANTI JOIN entry_chains c USING (entry_pdb_id))",
            )
        )
        orphan_biounits = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM (SELECT DISTINCT entry_pdb_id, "
                "system_biounit_id AS biounit_id FROM annotation) a ANTI JOIN "
                "(SELECT DISTINCT entry_pdb_id, biounit_id FROM entry_biounit_chains) b "
                "USING (entry_pdb_id, biounit_id)",
            )
        )
        orphan_metadata = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM (SELECT entry_pdb_id FROM entry_sources "
                "ANTI JOIN entry_metadata USING (entry_pdb_id))",
            )
        )
        orphan_interface_instances = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM ("
                "SELECT i.system_id FROM interfaces i LEFT JOIN "
                "entry_biounit_chains b1 ON b1.entry_pdb_id = i.entry_pdb_id "
                "AND b1.biounit_id = i.system_biounit_id "
                "AND b1.chain_instance = i.interface_chain_1 LEFT JOIN "
                "entry_biounit_chains b2 ON b2.entry_pdb_id = i.entry_pdb_id "
                "AND b2.biounit_id = i.system_biounit_id "
                "AND b2.chain_instance = i.interface_chain_2 "
                "WHERE b1.chain_instance IS NULL OR b2.chain_instance IS NULL)",
            )
        )
        orphan_interface_chains = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM ("
                "SELECT i.system_id FROM interfaces i LEFT JOIN entry_chains c1 "
                "ON c1.entry_pdb_id = i.entry_pdb_id AND c1.chain_asym_id = "
                "split_part(i.interface_chain_1, '.', 2) LEFT JOIN entry_chains c2 "
                "ON c2.entry_pdb_id = i.entry_pdb_id AND c2.chain_asym_id = "
                "split_part(i.interface_chain_2, '.', 2) WHERE "
                "c1.chain_asym_id IS NULL OR c2.chain_asym_id IS NULL OR "
                "c1.chain_receptor_type != 'protein' OR "
                "c2.chain_receptor_type != 'protein')",
            )
        )
        if (
            orphan_sources
            or orphan_chains
            or orphan_biounits
            or orphan_metadata
            or orphan_interface_instances
            or orphan_interface_chains
        ):
            raise ValueError(
                "referential-integrity failures: "
                f"sources={orphan_sources}, chains={orphan_chains}, "
                f"biounits={orphan_biounits}, metadata={orphan_metadata}, "
                f"interface_instances={orphan_interface_instances}, "
                f"interface_chains={orphan_interface_chains}"
            )
        missing_abilities = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM annotation WHERE system_type = 'holo' "
                "AND coalesce(ligand_is_proper, false) "
                "AND ligand_is_3d_score_able IS NULL",
            )
        )
        invalid_nonproper_abilities = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM annotation WHERE NOT (system_type = 'holo' "
                "AND coalesce(ligand_is_proper, false)) "
                "AND coalesce(ligand_is_3d_score_able, false)",
            )
        )
        invalid_systems = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM (SELECT system_id FROM annotation "
                "GROUP BY system_id HAVING bool_and(coalesce(ligand_is_ion, false) "
                "OR coalesce(ligand_is_artifact, false)))",
            )
        )
        mismatched_biounits = int(
            _fetch_scalar(
                connection,
                "SELECT count(*) FROM annotation WHERE "
                "split_part(ligand_id, '__', 2) IS DISTINCT FROM system_biounit_id "
                "OR split_part(system_id, '__', 2) IS DISTINCT FROM system_biounit_id",
            )
        )
        if (
            missing_abilities
            or invalid_nonproper_abilities
            or invalid_systems
            or mismatched_biounits
        ):
            raise ValueError(
                "annotation validation failed: "
                f"missing_scoreability={missing_abilities}, "
                f"nonproper_scoreability={invalid_nonproper_abilities}, "
                f"all_ion_or_artifact_systems={invalid_systems}, "
                f"mismatched_biounits={mismatched_biounits}"
            )
        return {
            "row_counts": actual_counts,
            "duplicate_key_counts": duplicates,
            "entry_count": int(
                _fetch_scalar(
                    connection,
                    "SELECT count(DISTINCT entry_pdb_id) FROM entry_metadata",
                )
            ),
            "system_count": int(
                _fetch_scalar(
                    connection, "SELECT count(DISTINCT system_id) FROM annotation"
                )
            ),
            "ligand_count": int(
                _fetch_scalar(
                    connection, "SELECT count(DISTINCT ligand_id) FROM annotation"
                )
            ),
            "interface_count": int(
                _fetch_scalar(connection, "SELECT count(*) FROM interfaces")
            ),
        }
    finally:
        connection.close()


def finalize_collation(
    data_dir: Path,
    *,
    threads: int = 4,
    memory_limit: str = "16GB",
    scratch_dir: Path | None = None,
    row_group_size: int = 100_000,
) -> dict[str, Any]:
    """Merge completed shards, validate them, and install final index files."""
    data_dir = data_dir.resolve()
    plan = load_plan(data_dir)
    shard_files, expected_counts = _load_completed_shards(
        data_dir, plan, verify_threads=threads
    )
    final_paths = {
        "annotation": data_dir / "index" / "annotation_table.parquet",
        "entry_chains": data_dir / "index" / "entry_chains.parquet",
        "entry_biounit_chains": data_dir / "index" / "entry_biounit_chains.parquet",
        "entry_metadata": data_dir / "index" / "entry_metadata.parquet",
        "interfaces": data_dir / "index" / "interface_annotation_table.parquet",
        "entry_sources": data_dir / "index" / "entry_sources.parquet",
    }
    temporary_paths = {
        name: _temporary_path(path) for name, path in final_paths.items()
    }
    for path in final_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect()
    try:
        _configure_duckdb(
            connection,
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
        )
        sort_columns = {
            "annotation": "entry_pdb_id, system_id, ligand_id",
            "entry_chains": "entry_pdb_id, chain_asym_id",
            "entry_biounit_chains": "entry_pdb_id, biounit_id, chain_instance",
            "entry_metadata": "entry_pdb_id",
            "interfaces": "entry_pdb_id, system_id",
            "entry_sources": "entry_pdb_id",
        }
        for name, paths in shard_files.items():
            connection.read_parquet(
                [str(path) for path in paths], union_by_name=True
            ).create_view(f"sharded_{name}", replace=True)
            _copy_query(
                connection,
                f"SELECT * FROM sharded_{name} ORDER BY {sort_columns[name]}",
                temporary_paths[name],
                row_group_size=row_group_size,
            )
    finally:
        connection.close()
    try:
        validation = _validate_final_tables(
            temporary_paths,
            expected_counts=expected_counts,
            min_interface_residues=int(plan["interface_min_residues"]),
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
        )
        _install_final_tables_fail_closed(
            temporary_paths,
            final_paths,
            data_dir / "index" / FINAL_MARKER_NAME,
        )
    finally:
        for path in temporary_paths.values():
            path.unlink(missing_ok=True)
    report: dict[str, Any] = {
        "version": COLLATION_VERSION,
        "status": "complete",
        "plan": str(plan_path(data_dir)),
        "manifest": str(manifest_path(data_dir)),
        "code_count": len(plan["codes"]),
        "interface_min_residues": int(plan["interface_min_residues"]),
        "outputs": {name: str(path) for name, path in final_paths.items()},
        **validation,
    }
    _write_json_atomic(data_dir / "index" / FINAL_MARKER_NAME, report)
    return report


def repair_collation(
    data_dir: Path,
    pdb_ids: Sequence[str],
    *,
    threads: int = 4,
    memory_limit: str = "16GB",
    scratch_dir: Path | None = None,
    row_group_size: int = 100_000,
) -> dict[str, Any]:
    """Atomically replace selected entries in an existing collated index.

    This is intended for corrections that require re-annotating a bounded set of
    entries.  Unselected rows are read from the installed index, so release-only
    columns are preserved without rebuilding every raw entry.
    """
    data_dir = data_dir.resolve()
    selected = sorted({str(pdb_id).lower() for pdb_id in pdb_ids})
    invalid = [
        pdb_id
        for pdb_id in selected
        if re.fullmatch(r"[0-9][a-z0-9]{3}", pdb_id) is None
    ]
    if invalid:
        raise ValueError(f"invalid repair PDB IDs: {invalid[:10]}")
    if not selected:
        raise ValueError("no PDB IDs selected for index repair")

    final_paths = {
        "annotation": data_dir / "index" / "annotation_table.parquet",
        "entry_chains": data_dir / "index" / "entry_chains.parquet",
        "entry_biounit_chains": data_dir / "index" / "entry_biounit_chains.parquet",
        "entry_metadata": data_dir / "index" / "entry_metadata.parquet",
        "interfaces": data_dir / "index" / "interface_annotation_table.parquet",
        "entry_sources": data_dir / "index" / "entry_sources.parquet",
    }
    missing_final = [path for path in final_paths.values() if not path.is_file()]
    if missing_final:
        raise FileNotFoundError(f"missing installed index tables: {missing_final}")

    marker_path = data_dir / "index" / FINAL_MARKER_NAME
    try:
        installed_marker = cast(dict[str, Any], json.loads(marker_path.read_text()))
        installed_interface_min_residues = int(
            installed_marker["interface_min_residues"]
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(
            "targeted repair requires a collated interface threshold"
        ) from exc

    rows = [
        _entry_manifest_row(
            data_dir,
            data_dir / "raw_entries" / pdb_id[1:3] / pdb_id,
        )
        for pdb_id in selected
    ]
    replacement_interface_thresholds = {
        int(row["interface_min_residues"]) for row in rows
    }
    if replacement_interface_thresholds != {installed_interface_min_residues}:
        raise ValueError(
            "targeted repair interface.min_interface_residues differs from the "
            f"installed release: replacement={sorted(replacement_interface_thresholds)}, "
            f"installed={installed_interface_min_residues}"
        )
    _verify_manifest_inputs(rows, threads=threads)
    temporary_paths = {
        name: _temporary_path(final_paths[name])
        for name in ("annotation", "entry_metadata", "interfaces")
    }
    replacement_paths = {
        name: _temporary_path(path.with_name(f"repair-{path.name}"))
        for name, path in final_paths.items()
        if name != "annotation"
    }

    try:
        for name, schema in SIDECAR_SCHEMAS.items():
            source_key = {
                "entry_chains": "chain_path",
                "entry_biounit_chains": "biounit_chain_path",
                "interfaces": "interface_path",
                "entry_sources": "source_path",
            }[name]
            sort_columns = {
                "entry_chains": ("entry_pdb_id", "chain_asym_id"),
                "entry_biounit_chains": (
                    "entry_pdb_id",
                    "biounit_id",
                    "chain_instance",
                ),
                "interfaces": ("entry_pdb_id", "system_id"),
                "entry_sources": ("entry_pdb_id",),
            }[name]
            _collate_sidecars(
                [Path(str(row[source_key])) for row in rows],
                replacement_paths[name],
                schema,
                sort_columns=sort_columns,
                row_group_size=row_group_size,
            )
        _collate_entry_metadata(
            [Path(str(row["metadata_path"])) for row in rows],
            replacement_paths["entry_metadata"],
            row_group_size=row_group_size,
        )

        connection = duckdb.connect()
        try:
            _configure_duckdb(
                connection,
                threads=threads,
                memory_limit=memory_limit,
                scratch_dir=scratch_dir,
            )
            connection.register(
                "repaired_entries", pa.table({"entry_pdb_id": selected})
            )
            ligand_rows = [row for row in rows if row["annotation_size"] is not None]
            if ligand_rows:
                _build_annotation_view(
                    connection,
                    [str(row["annotation_path"]) for row in ligand_rows],
                    [str(row["ligand_path"]) for row in ligand_rows],
                )
            connection.read_parquet(str(final_paths["annotation"])).create_view(
                "installed_annotation", replace=True
            )
            repaired_annotation = (
                "UNION ALL BY NAME SELECT * FROM collated_annotation"
                if ligand_rows
                else ""
            )
            _copy_query(
                connection,
                "SELECT * FROM (SELECT installed.* FROM installed_annotation "
                "AS installed ANTI JOIN repaired_entries USING (entry_pdb_id) "
                f"{repaired_annotation}) "
                "ORDER BY entry_pdb_id, system_id, ligand_id",
                temporary_paths["annotation"],
                row_group_size=row_group_size,
            )
            for name in (
                "entry_chains",
                "entry_biounit_chains",
                "entry_sources",
            ):
                connection.read_parquet(str(final_paths[name])).create_view(
                    f"installed_{name}", replace=True
                )
                connection.read_parquet(str(replacement_paths[name])).create_view(
                    f"replacement_{name}", replace=True
                )
                differences = int(
                    _fetch_scalar(
                        connection,
                        "SELECT count(*) FROM ("
                        f"(SELECT installed.* FROM installed_{name} AS installed "
                        "INNER JOIN repaired_entries USING (entry_pdb_id) "
                        f"EXCEPT ALL SELECT * FROM replacement_{name}) "
                        "UNION ALL "
                        f"(SELECT * FROM replacement_{name} EXCEPT ALL "
                        f"SELECT installed.* FROM installed_{name} AS installed "
                        "INNER JOIN repaired_entries USING (entry_pdb_id)))",
                    )
                )
                if differences:
                    raise ValueError(
                        f"targeted repair changed {name}; rebuild the protein "
                        "scoring inputs and mapped alignments instead"
                    )
            for name, order_by in (
                ("entry_metadata", "entry_pdb_id"),
                ("interfaces", "entry_pdb_id, system_id"),
            ):
                connection.read_parquet(str(final_paths[name])).create_view(
                    f"installed_{name}", replace=True
                )
                connection.read_parquet(str(replacement_paths[name])).create_view(
                    f"replacement_{name}", replace=True
                )
                _copy_query(
                    connection,
                    f"SELECT * FROM (SELECT installed.* FROM installed_{name} "
                    "AS installed ANTI JOIN repaired_entries USING (entry_pdb_id) "
                    f"UNION ALL BY NAME SELECT * FROM replacement_{name}) "
                    f"ORDER BY {order_by}",
                    temporary_paths[name],
                    row_group_size=row_group_size,
                )
        finally:
            connection.close()

        validation_paths = {**final_paths, **temporary_paths}
        expected_counts = {
            name: pq.ParquetFile(path).metadata.num_rows
            for name, path in validation_paths.items()
        }
        validation = _validate_final_tables(
            validation_paths,
            expected_counts=expected_counts,
            min_interface_residues=installed_interface_min_residues,
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
        )
        marker = data_dir / "index" / FINAL_MARKER_NAME
        marker.unlink(missing_ok=True)
        final_paths["annotation"].unlink(missing_ok=True)
        temporary_paths["entry_metadata"].replace(final_paths["entry_metadata"])
        temporary_paths["interfaces"].replace(final_paths["interfaces"])
        temporary_paths["annotation"].replace(final_paths["annotation"])
        # The nonredundant index and every release-only ligand enrichment were
        # derived from the pre-repair annotation.  Do not leave a stale
        # nonredundant table looking publishable while scores, fingerprints,
        # clusters, and final index enrichment are being repaired.
        (data_dir / "index" / "annotation_table_nonredundant.parquet").unlink(
            missing_ok=True
        )
    finally:
        for path in [*temporary_paths.values(), *replacement_paths.values()]:
            path.unlink(missing_ok=True)

    report: dict[str, Any] = {
        "version": COLLATION_VERSION,
        "status": REPAIR_REQUIRED_STATUS,
        "mode": "targeted_repair",
        "repaired_entry_count": len(selected),
        "interface_min_residues": installed_interface_min_residues,
        "repaired_entry_digest": hashlib.sha256(
            ("\n".join(selected) + "\n").encode("utf-8")
        ).hexdigest(),
        "outputs": {name: str(path) for name, path in final_paths.items()},
        "required_downstream_artifacts": [
            "ligand_similarity",
            "scores",
            "ligand_clusters",
            "ligand_sampling",
            "interface_scores",
            "interface_clusters",
            "interface_sampling",
            "annotation_table_nonredundant",
        ],
        **validation,
    }
    _write_json_atomic(data_dir / "index" / FINAL_MARKER_NAME, report)
    return report


def finalize_repair_marker(data_dir: Path) -> dict[str, Any] | None:
    """Mark a targeted collation repair complete after final index enrichment."""
    marker_path = data_dir / "index" / FINAL_MARKER_NAME
    if not marker_path.is_file():
        return None
    report = cast(dict[str, Any], json.loads(marker_path.read_text()))
    if (
        report.get("mode") != "targeted_repair"
        or report.get("status") != REPAIR_REQUIRED_STATUS
    ):
        return None
    nonredundant = data_dir / "index" / "annotation_table_nonredundant.parquet"
    if not nonredundant.is_file():
        raise FileNotFoundError(
            "targeted repair cannot be finalized before rebuilding "
            "annotation_table_nonredundant.parquet"
        )
    report["status"] = "complete"
    report["downstream_repair_complete"] = True
    _write_json_atomic(marker_path, report)
    return report


def run_collation(
    data_dir: Path,
    *,
    threads: int = 1,
    memory_limit: str = "8GB",
    scratch_dir: Path | None = None,
    row_group_size: int = 100_000,
    force: bool = False,
) -> dict[str, Any]:
    """Run all collation phases locally; intended for tests and small releases."""
    plan = plan_collation(data_dir, threads=threads)
    for code in plan["codes"]:
        collate_shard(
            data_dir,
            code,
            force=force,
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
            row_group_size=row_group_size,
        )
    return finalize_collation(
        data_dir,
        threads=threads,
        memory_limit=memory_limit,
        scratch_dir=scratch_dir,
        row_group_size=row_group_size,
    )


def planned_code_batch(
    data_dir: Path, *, batch_index: int, batch_size: int
) -> list[str]:
    """Return one fixed-size slice of the planned two-character codes."""
    if batch_index < 0:
        raise ValueError("batch_index must be non-negative")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    codes = [str(code) for code in load_plan(data_dir)["codes"]]
    start = batch_index * batch_size
    return codes[start : start + batch_size]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "shard", "finalize", "repair", "run"):
        command = commands.add_parser(name)
        command.add_argument("data_dir", type=Path)
        if name == "shard":
            command.add_argument("codes", nargs="*")
            command.add_argument("--batch-index", type=int)
            command.add_argument("--batch-size", type=int)
            command.add_argument("--pdb-manifest", type=Path)
        if name == "repair":
            command.add_argument("--pdb-manifest", type=Path, required=True)
        if name in {"plan", "shard", "finalize", "repair", "run"}:
            command.add_argument("--threads", type=int, default=1)
        if name in {"shard", "finalize", "repair", "run"}:
            command.add_argument("--memory-limit", default="8GB")
            command.add_argument("--scratch-dir", type=Path)
            command.add_argument("--row-group-size", type=int, default=100_000)
        if name in {"shard", "run"}:
            command.add_argument("--force", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "plan":
        result = plan_collation(args.data_dir, threads=args.threads)
    elif args.command == "shard":
        codes = args.codes
        if args.pdb_manifest is not None:
            if codes:
                raise ValueError("pass a PDB manifest or explicit codes, not both")
            from plinder.data.pipeline.ingest import load_manifest

            codes = sorted({pdb_id[1:3] for pdb_id in load_manifest(args.pdb_manifest)})
            if args.batch_index is not None or args.batch_size is not None:
                if args.batch_index is None or args.batch_size is None:
                    raise ValueError(
                        "both --batch-index and --batch-size are required for a slice"
                    )
                start = args.batch_index * args.batch_size
                codes = codes[start : start + args.batch_size]
                args.batch_index = None
                args.batch_size = None
        if args.batch_index is not None or args.batch_size is not None:
            if codes or args.batch_index is None or args.batch_size is None:
                raise ValueError(
                    "pass either explicit codes or both --batch-index and --batch-size"
                )
            codes = planned_code_batch(
                args.data_dir,
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
        if not codes:
            raise ValueError("no collation shard codes selected")
        result = {
            "shards": [
                collate_shard(
                    args.data_dir,
                    code,
                    force=args.force,
                    threads=args.threads,
                    memory_limit=args.memory_limit,
                    scratch_dir=args.scratch_dir,
                    row_group_size=args.row_group_size,
                )
                for code in codes
            ]
        }
    elif args.command == "finalize":
        result = finalize_collation(
            args.data_dir,
            threads=args.threads,
            memory_limit=args.memory_limit,
            scratch_dir=args.scratch_dir,
            row_group_size=args.row_group_size,
        )
    elif args.command == "repair":
        from plinder.data.pipeline.ingest import load_manifest

        result = repair_collation(
            args.data_dir,
            load_manifest(args.pdb_manifest),
            threads=args.threads,
            memory_limit=args.memory_limit,
            scratch_dir=args.scratch_dir,
            row_group_size=args.row_group_size,
        )
    else:
        result = run_collation(
            args.data_dir,
            threads=args.threads,
            memory_limit=args.memory_limit,
            scratch_dir=args.scratch_dir,
            row_group_size=args.row_group_size,
            force=args.force,
        )
    printable = result
    if args.command == "plan":
        printable = {
            key: value
            for key, value in result.items()
            if key not in {"code_entry_counts", "code_signatures", "codes"}
        }
    print(json.dumps(printable, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
