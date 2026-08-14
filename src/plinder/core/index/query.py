# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Query published PLINDER tables, with optional grain-preserving joins."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any, TypeAlias

import duckdb
import pandas as pd
import pyarrow.dataset as ds

from plinder.core.release import RELEASE_TABLES, PlinderRelease

Filter: TypeAlias = tuple[str, str, Any]
Filters: TypeAlias = list[Filter | list[Filter]] | None
JoinKeys: TypeAlias = tuple[tuple[str, str], ...]

DISABLED_ANNOTATION_COLUMNS = frozenset(
    {"system_has_binding_affinity", "ligand_binding_affinity"}
)


def _is_repeated_entry_column(column: str) -> bool:
    return column.startswith("entry_") and column != "entry_pdb_id"


# Every relationship points from a base table to a table whose join columns are
# unique.  These joins therefore preserve the base table's documented row grain.
TABLE_JOINS: dict[str, dict[str, JoinKeys]] = {
    "annotation": {
        "entry_metadata": (("entry_pdb_id", "entry_pdb_id"),),
        "entry_sources": (("entry_pdb_id", "entry_pdb_id"),),
        "ligand_pocket_membership": (("ligand_id", "ligand_id"),),
    },
    "entry_chains": {
        "alignment_chain_lookup": (
            ("entry_pdb_id", "entry_pdb_id"),
            ("chain_asym_id", "chain_asym_id"),
        ),
        "entry_metadata": (("entry_pdb_id", "entry_pdb_id"),),
        "entry_sources": (("entry_pdb_id", "entry_pdb_id"),),
    },
    "entry_biounit_chains": {
        "entry_chains": (
            ("entry_pdb_id", "entry_pdb_id"),
            ("chain_asym_id", "chain_asym_id"),
        ),
        "entry_metadata": (("entry_pdb_id", "entry_pdb_id"),),
        "entry_sources": (("entry_pdb_id", "entry_pdb_id"),),
    },
    "entry_metadata": {
        "entry_sources": (("entry_pdb_id", "entry_pdb_id"),),
    },
    "interface_annotations": {
        "entry_metadata": (("entry_pdb_id", "entry_pdb_id"),),
        "entry_sources": (("entry_pdb_id", "entry_pdb_id"),),
        "interface_membership": (("system_id", "system_id"),),
    },
    "alignment_chain_lookup": {
        "entry_chains": (
            ("entry_pdb_id", "entry_pdb_id"),
            ("chain_asym_id", "chain_asym_id"),
        ),
    },
    "ligand_pocket_membership": {
        "annotation": (("ligand_id", "ligand_id"),),
    },
    "interface_membership": {
        "interface_annotations": (("system_id", "system_id"),),
    },
    "linked_apo_structures": {
        "entry_metadata": (("source_entry_id", "entry_pdb_id"),),
        "entry_sources": (("source_entry_id", "entry_pdb_id"),),
    },
}


def _quote_identifier(value: str) -> str:
    return f'"{value.replace(chr(34), chr(34) * 2)}"'


def _parquet_sql(path: Path) -> str:
    source = path.resolve()
    if source.is_dir():
        source = source / "**" / "*.parquet"
    escaped = source.as_posix().replace("'", "''")
    return f"read_parquet('{escaped}', union_by_name = true)"


def _schema_names(path: Path) -> list[str]:
    if path.is_dir():
        dataset = ds.dataset(
            path,
            format="parquet",
            partitioning="hive",
            exclude_invalid_files=True,
        )
    else:
        dataset = ds.dataset(path, format="parquet", exclude_invalid_files=True)
    return list(dataset.schema.names)


def _condition_sql(
    condition: Filter,
    *,
    columns: dict[str, str],
    parameters: list[Any],
) -> str:
    if len(condition) != 3:
        raise ValueError(
            f"filters must be (column, operator, value): got {condition!r}"
        )
    column, operator, value = condition
    if column not in columns:
        raise ValueError(
            f"filter column {column!r} is unavailable; choose from: "
            f"{', '.join(columns)}"
        )
    field = columns[column]
    operator = operator.casefold().strip()
    operator = {"==": "=", "not in": "not in"}.get(operator, operator)
    if value is None:
        if operator in {"=", "is"}:
            return f"{field} IS NULL"
        if operator in {"!=", "<>", "is not"}:
            return f"{field} IS NOT NULL"
        raise ValueError(f"operator {operator!r} cannot compare with None")
    if operator in {"in", "not in"}:
        if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
            raise TypeError(f"operator {operator!r} requires a non-string iterable")
        values = list(value)
        if not values:
            return "FALSE" if operator == "in" else "TRUE"
        parameters.extend(values)
        placeholders = ", ".join("?" for _ in values)
        return f"{field} {operator.upper()} ({placeholders})"
    allowed = {"=", "!=", "<>", "<", "<=", ">", ">="}
    if operator not in allowed:
        raise ValueError(
            f"unsupported filter operator {operator!r}; choose from: "
            f"{', '.join(sorted(allowed | {'==', 'in', 'not in'}))}"
        )
    parameters.append(value)
    return f"{field} {operator} ?"


def _filters_sql(
    filters: Filters,
    *,
    columns: dict[str, str],
) -> tuple[str, list[Any]]:
    parameters: list[Any] = []
    clauses: list[str] = []
    for item in filters or []:
        if isinstance(item, list):
            if not item:
                raise ValueError("OR filter groups must not be empty")
            group = [
                _condition_sql(
                    condition,
                    columns=columns,
                    parameters=parameters,
                )
                for condition in item
            ]
            clauses.append(f"({' OR '.join(group)})")
        else:
            clauses.append(_condition_sql(item, columns=columns, parameters=parameters))
    if not clauses:
        return "", parameters
    return f"\nWHERE {' AND '.join(clauses)}", parameters


def _filter_columns(filters: Filters) -> set[str]:
    """Return every column referenced by a filter expression."""
    return {
        condition[0]
        for item in filters or []
        for condition in (item if isinstance(item, list) else [item])
        if len(condition) == 3
    }


def _column_sql(alias: str, column: str) -> str:
    return f"{_quote_identifier(alias)}.{_quote_identifier(column)}"


def query_table(
    table_name: str = "annotation",
    *,
    columns: list[str] | None = None,
    filters: Filters = None,
    joins: list[str] | None = None,
    release: PlinderRelease | None = None,
) -> pd.DataFrame:
    """Query one release table with optional grain-preserving sidecar joins.

    Joined tables must be listed explicitly.  Their registered join columns are
    unique, so joining never creates extra base rows. If a sidecar repeats a
    non-key column from the base table, its non-null values take precedence and
    unmatched base values are retained.

    Parameters
    ----------
    table_name : str
        Name in :data:`plinder.core.RELEASE_TABLES`.
    columns : list[str] | None
        Output columns. ``None`` returns every available column.
    filters : list | None
        Conditions combined with AND. A nested list is an OR group.
    joins : list[str] | None
        Registered sidecars to left-join to the base table.
    release : PlinderRelease | None
        Explicit local release, or the configured release when omitted.
    """
    if table_name not in RELEASE_TABLES:
        choices = ", ".join(sorted(RELEASE_TABLES))
        raise KeyError(f"unknown release table {table_name!r}; choose from: {choices}")
    release = release or PlinderRelease()
    selected_joins = list(dict.fromkeys(joins or []))
    allowed_joins = TABLE_JOINS.get(table_name, {})
    invalid_joins = [name for name in selected_joins if name not in allowed_joins]
    if invalid_joins:
        choices = ", ".join(sorted(allowed_joins)) or "none"
        raise ValueError(
            f"table {table_name!r} cannot join {invalid_joins}; "
            f"grain-preserving joins: {choices}"
        )

    requested_filter_columns = _filter_columns(filters)
    requested_columns = set(columns or [])
    requested_columns.discard("*")
    disabled_requested = DISABLED_ANNOTATION_COLUMNS.intersection(
        requested_columns | requested_filter_columns
    )
    annotation_is_available = (
        table_name == "annotation" or "annotation" in selected_joins
    )
    if annotation_is_available and disabled_requested:
        raise ValueError(
            "binding_affinity columns are disabled in the current dataset: "
            f"{sorted(disabled_requested)}"
        )

    table_names = [table_name, *selected_joins]
    aliases = {name: f"t{index}" for index, name in enumerate(table_names)}
    paths = {
        name: release.fetch(str(RELEASE_TABLES[name]["artifact"]))
        for name in table_names
    }
    schemas = {name: _schema_names(path) for name, path in paths.items()}

    output_order = [
        name
        for name in schemas[table_name]
        if table_name != "annotation"
        or (
            name not in DISABLED_ANNOTATION_COLUMNS
            and not _is_repeated_entry_column(name)
        )
    ]
    output_columns = {
        name: _column_sql(aliases[table_name], name) for name in output_order
    }
    joined_column_owner: dict[str, str] = {}
    join_sql: list[str] = []
    for join_name in selected_joins:
        keys = allowed_joins[join_name]
        base_missing = [
            column for column, _ in keys if column not in schemas[table_name]
        ]
        join_missing = [
            column for _, column in keys if column not in schemas[join_name]
        ]
        if base_missing or join_missing:
            raise ValueError(
                f"join {table_name!r} -> {join_name!r} is unavailable; "
                f"missing base columns {base_missing}, sidecar columns {join_missing}"
            )
        side_key_columns = {side_column for _, side_column in keys}
        visible_columns = list(schemas[join_name])
        if join_name == "annotation":
            visible_columns = [
                name
                for name in visible_columns
                if name not in DISABLED_ANNOTATION_COLUMNS
                and not _is_repeated_entry_column(name)
            ]
        for name in visible_columns:
            if name in side_key_columns:
                continue
            if name in joined_column_owner:
                previous = joined_column_owner[name]
                raise ValueError(
                    f"joined column {name!r} is provided by both {previous!r} "
                    f"and {join_name!r}"
                )
            joined_column_owner[name] = join_name
            if name not in output_columns:
                output_order.append(name)
                output_columns[name] = _column_sql(aliases[join_name], name)
            else:
                output_columns[name] = (
                    f"COALESCE({_column_sql(aliases[join_name], name)}, "
                    f"{output_columns[name]})"
                )
        conditions = " AND ".join(
            f"{_quote_identifier(aliases[table_name])}.{_quote_identifier(base)} = "
            f"{_quote_identifier(aliases[join_name])}.{_quote_identifier(side)}"
            for base, side in keys
        )
        join_sql.append(
            f"LEFT JOIN {_parquet_sql(paths[join_name])} "
            f"AS {_quote_identifier(aliases[join_name])} ON {conditions}"
        )

    requested = output_order if columns is None or columns == ["*"] else columns
    if not requested:
        raise ValueError("columns must not be empty")
    duplicates = sorted({name for name in requested if requested.count(name) > 1})
    if duplicates:
        raise ValueError(f"duplicate output columns: {duplicates}")
    missing = [name for name in requested if name not in output_columns]
    if missing:
        raise ValueError(
            f"columns {missing} are unavailable; choose from: "
            f"{', '.join(output_order)}"
        )
    select_sql = ", ".join(
        f"{output_columns[name]} AS {_quote_identifier(name)}" for name in requested
    )
    where_sql, parameters = _filters_sql(filters, columns=output_columns)
    query = (
        f"SELECT {select_sql}\n"
        f"FROM {_parquet_sql(paths[table_name])} "
        f"AS {_quote_identifier(aliases[table_name])}\n"
        f"{' '.join(join_sql)}{where_sql}"
    )
    connection = duckdb.connect()
    try:
        return connection.execute(query, parameters).fetch_df()
    finally:
        connection.close()


__all__ = ["TABLE_JOINS", "query_table"]
