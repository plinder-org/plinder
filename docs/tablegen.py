from __future__ import annotations

from glob import glob
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from google.cloud import storage
from itables import to_html_datatable
from sphinx.util import logging

ROWS_PER_PAGE = 10
# The columns that should be displayed first in the table
PRIMARY_COLUMNS = [
    "system_id",
]
CACHE_FILE = "index.parquet"

CHECKMARK = '<svg class="marks" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 448 512"><!--!Font Awesome Free 6.6.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free Copyright 2024 Fonticons, Inc.--><path fill="#5c7ec0" d="M438.6 105.4c12.5 12.5 12.5 32.8 0 45.3l-256 256c-12.5 12.5-32.8 12.5-45.3 0l-128-128c-12.5-12.5-12.5-32.8 0-45.3s32.8-12.5 45.3 0L160 338.7 393.4 105.4c12.5-12.5 32.8-12.5 45.3 0z"/></svg>'
CROSSMARK = '<svg class="marks" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 384 512"><!--!Font Awesome Free 6.6.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free Copyright 2024 Fonticons, Inc.--><path fill="#cccccc" d="M342.6 150.6c12.5-12.5 12.5-32.8 0-45.3s-32.8-12.5-45.3 0L192 210.7 86.6 105.4c-12.5-12.5-32.8-12.5-45.3 0s-12.5 32.8 0 45.3L146.7 256 41.4 361.4c-12.5 12.5-12.5 32.8 0 45.3s32.8 12.5 45.3 0L192 301.3 297.4 406.6c12.5 12.5 32.8 12.5 45.3 0s12.5-32.8 0-45.3L237.3 256 342.6 150.6z"/></svg>'


logger = logging.getLogger(__name__)


def generate_table(description_dir: Path, output_html_path: Path) -> None:
    """
    Generate a HTML table containing the column descriptions for
    `annotation_table.parquet`.

    The table is generated from TSV files.

    Parameters
    ----------
    description_dir : Path
        Directory containing the TSV files with the column descriptions.
    output_html_path : Path
        Path to the output HTML file.
        If the file already exists, it this function does nothing.

    Notes
    -----
    The functions caches the output file:
    If the file already exists, it does nothing.
    This skips the time consuming part of downloading and reading the annotation table.
    """
    column_descriptions = []
    for description_file in glob(str(description_dir / "*.tsv")):
        with open(description_file, "r") as f:
            column_descriptions_in_file = pd.read_csv(f, sep="\t")
        column_descriptions.append(column_descriptions_in_file)
    column_descriptions = pd.concat(column_descriptions, ignore_index=True)
#     # TODO: Remove as soon as wrong column names are fixed
#     column_descriptions = column_descriptions[
#         ~column_descriptions["Name"].str.contains("Kinase")
#     ]

    annotation_table = _get_annotation_table("2024-06", "v2", Path(CACHE_FILE))

    is_mandatory = np.zeros(column_descriptions.shape[0], dtype=bool)
    examples = [None] * column_descriptions.shape[0]
    for i, (column_name, data_type) in enumerate(
        zip(column_descriptions["Name"], column_descriptions["Type"])
    ):
        try:
            column = annotation_table[column_name]
        except KeyError:
            logger.warning(
                f"Column '{column_name}' is in column descriptions, "
                "but not found in annotation table."
            )
            continue
        is_value = _is_value(column, data_type)
        # Columns are considered mandatory if they contain a value in all rows
        is_mandatory[i] = is_value.all()
        # Get the first non-empty value as an example
        # whitelist ligand_rdkit_validation because we're only sampling the first 1000 rows
        if not is_value.any():
            if column_name in ["ligand_rdkit_validation"]:
                examples[i] = repr(None)
            else:
                logger.warning(f"Column '{column_name}' has no values.")
                examples[i] = repr(None)
        else:
            examples[i] = repr(column[is_value].iloc[0])
    column_descriptions["Mandatory"] = is_mandatory
    column_descriptions["Example"] = examples

    # Put important columns first
    column_descriptions = _prepend(column_descriptions, PRIMARY_COLUMNS)

    column_descriptions["Name"] = _to_monospace(column_descriptions["Name"])
    column_descriptions["Type"] = _to_monospace(column_descriptions["Type"])
    column_descriptions["Mandatory"] = _to_marks(column_descriptions["Mandatory"])
    column_descriptions["Example"] = _to_example(_to_monospace(
        column_descriptions["Example"])
    )  # fmt: skip

    html = to_html_datatable(
        column_descriptions,
        display_logo_when_loading=False,
        # Only one value, as the user can't change the number of rows
        lengthMenu=[ROWS_PER_PAGE],
        layout={
            "topStart": "search",
            "topEnd": "paging",
            "bottomStart": None,
            "bottomEnd": None,
        },
        # Use 'table' class from pydata theme to inherit the style of the usual tables
        classes="table display compact",
        # Without setting width to 100%,
        # the table would get a small horizontal scrollbar
        style="width:100%;overflow-wrap:break-word",
        autoWidth=False,
    )
    with open(output_html_path, "w") as f:
        # Inline HTML in Markdown does not allow empty lines
        f.write(_remove_empty_lines(html))


def _remove_empty_lines(string: str) -> str:
    """
    Remove empty lines from a multiline string.

    Parameters
    ----------
    html : str
        HTML string.

    Returns
    -------
    str
        HTML string without empty lines.
    """
    return "\n".join(line for line in string.splitlines() if len(line.strip()) != 0)


def _to_monospace(values: pd.Series) -> pd.Series:
    """
    Convert a series of strings to monospace.

    Parameters
    ----------
    values : pd.Series
        Series of strings.

    Returns
    -------
    pd.Series
        Series of strings in monospace.
    """
    return values.map(lambda x: f"<code>{x}</code>")


def _to_marks(values: pd.Series) -> pd.Series:
    """
    Convert a series of booleans to check/cross marks.

    Parameters
    ----------
    values : pd.Series
        Series of booleans.

    Returns
    -------
    pd.Series
        Series of check/cross marks.
    """
    return values.map(lambda x: CHECKMARK if x else CROSSMARK)


def _to_example(values: pd.Series) -> pd.Series:
    """
    Put the given series of text into a `div` with the `.example` class.

    Parameters
    ----------
    values : pd.Series
        Series of strings.

    Returns
    -------
    pd.Series
        Series of strings in a `div` with the `.example` class.
    """
    return values.map(lambda x: f'<div class="example">{x}</div>')


def _get_annotation_table(
    release: str, iteration: str, cache_path: Path
) -> pd.DataFrame:
    """
    Read the annotation table from the Parquet file.

    Parameters
    ----------
    release : str
        Time stamp of the last RCSB sync.
    iteration : str
        Iterative development within a release, i.e. the version
    cache_path : str

    Returns
    -------
    pd.DataFrame
        Annotation table.
    """
    if not cache_path.exists():
        storage.Client.create_anonymous_client().get_bucket("plinder").blob(
            blob_name=f"{release}/{iteration}/index/annotation_table.parquet"
        ).download_to_filename(cache_path)
    # only load the first 1000 rows to avoid memory issues
    df = duckdb.sql(f"select * from read_parquet('{cache_path.as_posix()}') limit 1000;").to_df()
    return df


def _is_value(column: pd.Series, data_type) -> pd.Series:
    """
    Get a mask for non-empty, non-NaN values.

    Parameters
    ----------
    column : pd.Series
        Column of the annotation table to check.
    data_type : str
        Data type of the column.
        This is not the `Series.dtype`, but the data type as a string to
        distinguish between different types of `object` columns.

    Returns
    -------
    mask : pd.Series
        Whether the elements in the series are actually values.
    """
    if data_type in ("str", "list[str]", "list[list[str]]"):
        # For strings and lists check if they are non-empty
        return ~pd.Series([obj is None or pd.isna(obj).all() if isinstance(obj, np.ndarray) else pd.isna(obj) or len(obj) == 0 for obj in column])
    return column.notna()


def _prepend(
    column_descriptions: pd.DataFrame, column_names: list[str]
) -> pd.DataFrame:
    """
    Put the specified column descriptions first in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to reorder.
    columns : list[str]
        Columns that should be displayed first.

    Returns
    -------
    pd.DataFrame
        DataFrame with the specified column descriptions first.
    """
    primary_indices = []
    for column_name in column_names:
        indices = np.where(column_descriptions["Name"] == column_name)[0]
        if len(indices) == 0:
            logger.warning(
                f"Primary column '{column_name}' not found in column descriptions."
            )
            continue
        primary_indices.append(indices[0])
    # Remove the primary columns from the original order...
    indices = np.arange(column_descriptions.shape[0])
    indices = indices[~np.isin(indices, primary_indices)]
    # ...and insert them at the beginning
    indices = np.concatenate([primary_indices, indices])
    column_descriptions = column_descriptions.iloc[indices]
    column_descriptions.reset_index(drop=True, inplace=True)
    return column_descriptions
