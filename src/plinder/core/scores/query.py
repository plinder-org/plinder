# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import List, NewType, Set, Tuple, Union

import numpy as np
import pyarrow as pa

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)
SQL_OP_MAP = {"==": "="}
FILTER = NewType("FILTER", Tuple[str, str, Union[str, Set[str]]])
FILTERS = Union[List[List[FILTER]], List[FILTER], None]


def _handle_condition_by_schema(
    schema: pa.Schema, col: str, val: str | set[str]
) -> str:
    """
    Deal with quoting hell for the varied conditions
    and types in the schema.

    Parameters
    ----------
    schema : pa.Schema
        the schema to use
    col : str
        the field name to get the type
    val : str
        the value provided in the filter

    Returns
    -------
    val : str
        the value to use in the query
    """
    try:
        dt = schema.types[schema.names.index(col)].to_pandas_dtype()
        if np.issubdtype(dt, np.number):
            val = val
        else:
            val = f"'{val}'"
    except NotImplementedError:
        val = f"'{val}'"
    if isinstance(val, (set, list, tuple)):
        val = f"({str(val)[1:-1]})"
    return val


def _handle_condition_by_type(val: str | set[str]) -> str:
    """
    Deal with quoting hell for the varied types passed
    in the filters.

    Parameters
    ----------
    val : str
        the value provided in the filter

    Returns
    -------
    val : str
        the value to use in the query
    """
    if isinstance(val, str):
        val = f"'{val}'"
    elif isinstance(val, (list, set)):
        val = f"({str(val)[1:-1]})"
    else:
        LOG.debug(f"_handle_condition_by_type: type(val)={type(val)} not mutated")
    return val


def _handle_inner_filter(
    filter: FILTER,
    schema: pa.Schema | None,
) -> str:
    """
    Convert the tuple format of a filter into a valid SQL substring

    Parameters
    ----------
    filter : tuple[str, str, str | set[str]]
        the filter to apply
    schema : pa.Schema | None
        the schema to use

    Returns
    -------
    str
        the formatted SQL filter
    """
    if len(filter) != 3:
        raise ValueError(f"filters must be (column, operator, value): got {filter}")
    col, op, val = filter
    if op == "in":
        val = _handle_condition_by_type(val)
    else:
        if schema is not None:
            if col not in schema.names:
                raise ValueError(f"column={col} not in schema={schema.names}")
            val = _handle_condition_by_schema(schema, col, val)
        else:
            val = _handle_condition_by_type(val)
    return f"{col} {SQL_OP_MAP.get(op, op)} {val}"


def _handle_filters(
    filters: FILTERS,
    allow_no_filters: bool,
    schema: pa.Schema | None,
) -> list[str | list[str]]:
    """
    Format the filters for the query

    Parameters
    ----------
    filters : list[tuple[str, str, str | set[str]]], default=None
        the filters to apply
    allow_no_filters : bool
        if True, allow no filters
    schema : pa.Schema, default=None
        the schema to validate the filters against

    Returns
    -------
    filters : list[tuple[str, str, str | set[str]]]
        the formatted filters
    """
    missing = False
    if filters is None:
        missing = True
    elif not len(filters):
        missing = True
    if missing and not allow_no_filters:
        raise ValueError("no filters provided, aborting query generation!")
    wheres: list[str | list[str]] = []
    if filters is not None:
        for filter in filters:
            if isinstance(filter, list):
                inner_wheres: list[str] = []
                for inner_filter in filter:
                    inner_wheres.append(_handle_inner_filter(inner_filter, schema))
                wheres.append(inner_wheres)
            else:
                wheres.append(_handle_inner_filter(filter, schema))
    return wheres


def _handle_columns(
    columns: list[str] | None,
    include_filename: bool,
    schema: pa.Schema | None,
) -> list[str]:
    """
    Format the columns for the query

    Parameters
    ----------
    columns : list[str], default=None
        the columns to select
    include_filename : bool
        if True, include the filename in the columns
    schema : pa.Schema, default=None
        the schema to validate the columns against
    """
    cols: list[str]
    if schema is not None:
        cols = schema.names if columns is None or not len(columns) else columns
        if cols != ["*"]:
            for col in cols:
                if col not in schema.names:
                    if include_filename and col == "filename":
                        continue
                    raise ValueError(f"column={col} not in schema={schema.names}")
    else:
        cols = columns if columns is not None else ["*"]
    return cols


def _handle_source(
    dataset: Path,
    nested: bool,
    include_filename: bool,
    filename: str | None,
) -> str:
    """
    Format the source of data for the query

    Parameters
    ----------
    dataset : Path
        the dataset to query
    nested : bool
        if True, select all parquet files in the dataset
    include_filename : bool
        if True, include the filename in the result
    filename : str, default=None
        a specific filename to read from if present

    Returns
    -------
    content : str
        the formatted source
    """
    if dataset.is_file():
        content = f"'{dataset.as_posix()}'"
    else:
        content = f"'{dataset}/*.parquet'"
    if nested:
        content = f"'{dataset}/**/*.parquet'"
    if filename is not None:
        content = f"'{dataset}/{filename}'"
    if include_filename:
        content = f"read_parquet({content}, filename = true)"
    return content


def _format_query(
    columns: list[str],
    filters: list[str | list[str]],
    source: str,
) -> str:
    """
    Pretty format the SQL query for readability of the SQL
    syntax when debugging.

    Parameters
    ----------
    columns : list[str]
        the columns to select
    filters : list[list[str]] | list[str]
        the filters to apply
    source : str
        the source of the data

    Returns
    -------
    qry : str
        the formatted SQL query
    """
    margin = "\n    "
    select = margin.join(["", f",{margin}".join(columns)])
    ands = []
    union = "AND"
    for ors in filters:
        if isinstance(ors, list):
            inner = f"{margin}    "
            union = "OR"
            ands.append("".join([f"({inner}", f"{inner}AND ".join(ors), f"{margin})"]))
        else:
            ands.append(f"{ors}")
    where = margin.join(["", f"{margin}{union} {margin}".join(ands)])
    qry = f"""\
SELECT {select}
FROM
    {source}"""
    if len(filters):
        qry += f"\nWHERE {where}"
    if ";" in qry:
        raise ValueError(f"query={qry} contains a semicolon!")
    LOG.debug("\n" + qry + ";")
    return qry


def make_query(
    dataset: Path,
    schema: pa.Schema | None = None,
    columns: list[str] | None = None,
    filters: FILTERS = None,
    nested: bool = False,
    allow_no_filters: bool = False,
    include_filename: bool = False,
    filename: str | None = None,
) -> str | None:
    """
    Convert the kwargs commonly passed to pd.read_parquet
    into a duckdb SQL query string. Fill in the blanks and
    raise an error if no filters are provided.

    Parameters
    ----------
    schema : pa.Schema
        the pyarrow schema of the dataset
    dataset : Path
        the path to the dataset
    columns : list[str], default=None
        the columns to select
    filters : list[tuple[str, str, str]], default=None
        the filters to apply
    nested : bool, default=False
        if True, select all parquet files in the dataset
    allow_no_filters : bool, default=False
        if True, allow no filters
    include_filename : bool, default=False
        if True, include the filename in the result
    filename : str, default=None
        a specific filename to read from if present

    Returns
    -------
    query : str | None
        the duckdb SQL query string
    """
    wheres = _handle_filters(filters, allow_no_filters, schema)
    cols = _handle_columns(columns, include_filename, schema)
    content = _handle_source(dataset, nested, include_filename, filename)
    return f"{_format_query(cols, wheres, content)};"
