from __future__ import annotations

from html import escape
from pathlib import Path

import pandas as pd
from itables import to_html_datatable

ROWS_PER_PAGE = 10
DESCRIPTION_COLUMNS = ["Name", "Type", "Description"]


def _render_description_table(frame: pd.DataFrame) -> str:
    """Render one release table's column descriptions."""
    return to_html_datatable(
        frame,
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


def generate_table(description_dir: Path, output_html_path: Path) -> None:
    """Render one searchable HTML table per checked-in release table."""
    sections: list[str] = []
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
        table_name = escape(description_file.stem, quote=True)
        table_id = description_file.stem.replace("_", "-")
        table_html = _render_description_table(frame)
        sections.extend(
            [
                f'<section class="release-column-table" id="columns-{table_id}">',
                f"<h3>{table_name}</h3>",
                *(line for line in table_html.splitlines() if line.strip()),
                "</section>",
            ]
        )
    if not sections:
        raise FileNotFoundError(
            f"no release-table descriptions found below {description_dir / 'tables'}"
        )

    output_html_path.write_text(
        "\n".join(line for line in sections if line.strip()),
        encoding="utf-8",
    )
