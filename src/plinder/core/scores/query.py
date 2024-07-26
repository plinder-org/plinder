# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import numpy as np
import pyarrow as pa

from plinder.core.utils import cpl
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)
SQL_OP_MAP = {"==": "="}


def ensure_dataset(rel: str) -> Path:
    """
    Check if the dataset exists and download it if not.

    Parameters
    ----------
    rel : str
        the relative path to the dataset

    Returns
    -------
    root : Path
        the local path to the dataset
    """
    root = cpl.get_plinder_path(rel=rel)
    local = Path(root.fspath)
    LOG.debug(f"ensure_dataset: root={root}, local={local}")
    if not local.exists():
        LOG.info(f"downloading {rel}, stand by...")
        cpl.download_many(rel=rel)
    return local


def _handle_condition_by_schema(schema: pa.Schema, col: str, val: str) -> str:
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
    return val


def _handle_condition_by_type(val: str) -> str:
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
    return val


def make_query(
    schema: pa.Schema,
    dataset: Path,
    columns: list[str] | None = None,
    filters: list[tuple[str, str, str]] | None = None,
    nested: bool = False,
    allow_no_filters: bool = False,
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

    Returns
    -------
    query : str | None
        the duckdb SQL query string
    """
    if filters is None or not len(filters):
        if not allow_no_filters:
            LOG.error("no filters provided, aborting query generation!")
            return None
        else:
            filters = []
    wheres = []
    for filter in filters:
        if len(filter) != 3:
            raise ValueError(f"filters must be (column, operator, value): got {filter}")
        col, op, val = filter
        if col not in schema.names:
            raise ValueError(f"column={col} not in schema={schema.names}")
        val = _handle_condition_by_schema(schema, col, val)
        wheres.append(f"{col} {SQL_OP_MAP.get(op, op)} {val}")
    cs = schema.names if columns is None or not len(columns) else columns
    for col in cs:
        if col not in schema.names:
            raise ValueError(f"column={col} not in schema={schema.names}")
    margin = "\n            "
    select = margin.join(["", f",{margin}".join(cs)])
    where = margin.join(["", f" AND {margin}".join(wheres)])
    glob = f"{margin}'{dataset}/*.parquet'"
    if nested:
        glob = f"{margin}'{dataset}/**/*.parquet'"
    qry = dedent(
        f"""\
        SELECT {select}
        FROM {glob}
        """
    )
    if not allow_no_filters:
        qry += f"\nWHERE {where}"
    if ";" in qry:
        raise ValueError(f"query={qry} contains a semicolon!")
    LOG.debug("\n" + qry + ";")
    return f"{qry};"


def make_query_no_schema(
    dataset: Path,
    columns: list[str] | None = None,
    filters: list[tuple[str, str, str]] | None = None,
    nested: bool = False,
    allow_no_filters: bool = False,
    filename: str | None = None,
) -> str:
    """
    Query a parquet dataset without a provided schema

    Parameters
    ----------
    dataset : Path
        the dataset to query
    columns : list[str], default=None
        the columns to select
    filters : list[tuple[str, str, str]], default=None
        the filters to apply
    nested : bool, default=False
        if True, select all parquet files in the dataset
    allow_no_filters : bool, default=False
        if True, allow no filters
    filename : str, default=None
        a specific filename to read from if present

    Returns
    -------
    str
        the duckdb SQL query string
    """
    margin = "\n            "
    wheres = []
    if filters is not None:
        for filter in filters:
            if len(filter) != 3:
                raise ValueError(
                    f"filters must be (column, operator, value): got {filter}"
                )
            col, op, val = filter
            val = _handle_condition_by_type(val)
            wheres.append(f"{col} {SQL_OP_MAP.get(op, op)} {val}")
    cs = columns if columns is not None else ["*"]
    select = margin.join(["", f",{margin}".join(cs)])
    where = margin.join(["", f" AND {margin}".join(wheres)])
    glob = f"{margin}'{dataset}/*.parquet'"
    if nested:
        glob = f"{margin}'{dataset}/**/*.parquet'"
    if filename is not None:
        glob = f"{margin}'{dataset}/{filename}'"
    qry = dedent(
        f"""\
        SELECT {select}
        FROM {glob}
        """
    )
    if not allow_no_filters:
        qry += f"\nWHERE {where}"
    if ";" in qry:
        raise ValueError(f"query={qry} contains a semicolon!")
    LOG.debug("\n" + qry + ";")
    return f"{qry};"
