from __future__ import annotations

from glob import glob
from pathlib import Path

import pandas as pd
from itables import to_html_datatable

ROWS_PER_PAGE = 10

CHECKMARK = '<svg class="marks" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 448 512"><!--!Font Awesome Free 6.6.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free Copyright 2024 Fonticons, Inc.--><path fill="#da5da4" d="M438.6 105.4c12.5 12.5 12.5 32.8 0 45.3l-256 256c-12.5 12.5-32.8 12.5-45.3 0l-128-128c-12.5-12.5-12.5-32.8 0-45.3s32.8-12.5 45.3 0L160 338.7 393.4 105.4c12.5-12.5 32.8-12.5 45.3 0z"/></svg>'
CROSSMARK = '<svg class="marks" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 384 512"><!--!Font Awesome Free 6.6.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free Copyright 2024 Fonticons, Inc.--><path fill="#cccccc" d="M342.6 150.6c12.5-12.5 12.5-32.8 0-45.3s-32.8-12.5-45.3 0L192 210.7 86.6 105.4c-12.5-12.5-32.8-12.5-45.3 0s-12.5 32.8 0 45.3L146.7 256 41.4 361.4c-12.5 12.5-12.5 32.8 0 45.3s32.8 12.5 45.3 0L192 301.3 297.4 406.6c12.5 12.5 32.8 12.5 45.3 0s12.5-32.8 0-45.3L237.3 256 342.6 150.6z"/></svg>'


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
    """
    column_descriptions = []
    for description_file in glob(str(description_dir / "*.tsv")):
        with open(description_file, "r") as f:
            column_descriptions_in_file = pd.read_csv(f, sep="\t")
        column_descriptions.append(column_descriptions_in_file)
    column_descriptions = pd.concat(column_descriptions, ignore_index=True)

    column_descriptions["Mandatory"] = [True] * column_descriptions.shape[0]
    column_descriptions["Example"] = [
        '"aasdnkasdhjkashdkhaskjdhashkdhjaksdhkashjd"'
    ] * column_descriptions.shape[0]

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


"""
show(
    annotations,
    classes="cell-border",
    columns= [{ "width": '10%' }, None, None, None, None]
)
"""
