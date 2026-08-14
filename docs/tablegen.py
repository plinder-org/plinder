from __future__ import annotations

from html import escape
from pathlib import Path

import pandas as pd
from itables import to_html_datatable

ROWS_PER_PAGE = 10
DESCRIPTION_COLUMNS = ["Name", "Type", "Description"]


def generate_table(description_dir: Path, output_html_path: Path) -> None:
    """Render the checked-in release-table column descriptions as HTML."""
    frames: list[pd.DataFrame] = []
    for description_file in sorted((description_dir / "tables").glob("*.tsv")):
        frame = pd.read_csv(description_file, sep="\t")
        if frame.columns.tolist() != DESCRIPTION_COLUMNS:
            raise ValueError(
                f"{description_file} must contain columns {DESCRIPTION_COLUMNS}"
            )
        if frame["Name"].duplicated().any():
            raise ValueError(f"{description_file} contains duplicate column names")
        frame = frame.apply(
            lambda column: column.map(lambda value: escape(str(value), quote=True))
        )
        frame.insert(0, "Table", escape(description_file.stem, quote=True))
        frames.append(frame)
    if not frames:
        raise FileNotFoundError(
            f"no release-table descriptions found below {description_dir / 'tables'}"
        )

    descriptions = pd.concat(frames, ignore_index=True)
    for column in ("Table", "Name", "Type"):
        descriptions[column] = descriptions[column].map(
            lambda value: f"<code>{value}</code>"
        )

    html = to_html_datatable(
        descriptions,
        display_logo_when_loading=False,
        lengthMenu=[ROWS_PER_PAGE],
        layout={
            "topStart": "search",
            "topEnd": "paging",
            "bottomStart": None,
            "bottomEnd": None,
        },
        classes="table display compact",
        style="width:100%;overflow-wrap:break-word",
        autoWidth=False,
    )
    output_html_path.write_text(
        "\n".join(line for line in html.splitlines() if line.strip()),
        encoding="utf-8",
    )
