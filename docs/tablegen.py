from glob import glob
from pathlib import Path
from xml.etree import ElementTree

import pandas as pd


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
    _write_table(
        column_descriptions,
        output_html_path,
        # 'table' class is also used by Sphinx for all other tables,
        # hence this class applies the same styling to this table as well.
        table_class="table",
        table_id="tableReference",
    )


def _write_table(
    column_descriptions: pd.DataFrame,
    output_html_path: Path,
    table_class: str | None = None,
    table_id: str | None = None,
) -> None:
    table = ElementTree.Element("table")
    if table_id is not None:
        table.set("id", table_id)
    if table_class is not None:
        table.set("class", table_class)

    thead = ElementTree.SubElement(table, "thead")
    tr = ElementTree.SubElement(thead, "tr")
    for column_name in column_descriptions.columns:
        th = ElementTree.SubElement(tr, "th")
        th.text = column_name

    tbody = ElementTree.SubElement(table, "tbody")
    for _, row in column_descriptions.iterrows():
        tr = ElementTree.SubElement(tbody, "tr")
        for column_name, value in row.items():
            td = ElementTree.SubElement(tr, "td")
            td.text = str(value)

    with open(output_html_path, "w") as f:
        f.write(ElementTree.tostring(table).decode())
