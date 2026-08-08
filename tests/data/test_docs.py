# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from plinder.data import docs


def test_ligand_cluster_column_descriptions():
    import pandas as pd

    columns = [
        "shape__50__strong__component",
        "shape__50__ligand__component",
        "color__70__community",
        "color__70__ligand__community",
        "pocket_qcov__50__ligand__directed_set_cover",
        "pocket_qcov__50__ligand__directed_set_cover__is_centroid",
        "pocket_qcov__50__ligand__directed_set_cover__coverage_count",
        "pocket_qcov__50__ligand__directed_set_cover__coverage_fraction",
        "interface_side_qcov__70__chain_1_component",
        "interface_side_qcov__70__chain_2_community",
        "interface_side_qcov__70__chain_1_directed_set_cover",
    ]

    rows = docs.get_cluster_column_descriptions(pd.DataFrame(columns=columns))

    assert [row[0] for row in rows] == columns
    descriptions = {name: description for name, _, description in rows}
    assert "ligand-level reciprocal-minimum component" in descriptions[columns[1]]
    assert "ligand-level greedy centroid community" in descriptions[columns[3]]
    assert "ligand-level directed set cover" in descriptions[columns[4]]
    assert "query-to-centroid score" in descriptions[columns[4]]
    assert "published centroid" in descriptions[columns[5]]
    assert "Number of directed-cover query nodes" in descriptions[columns[6]]
    assert "Fraction of its directed weak component" in descriptions[columns[7]]
    assert "chain 1 reciprocal-minimum component" in descriptions[columns[8]]
    assert "chain 2 greedy centroid community" in descriptions[columns[9]]
    assert "chain 1 directed set cover" in descriptions[columns[10]]


def test_description_markers_are_explicit_and_independent():
    from plinder.data.annotations.utils import (
        DocBaseModel,
        description_excluded_from_column_docs,
        description_excluded_from_flat_export,
    )
    from pydantic import Field

    class Example(DocBaseModel):
        visible: int = Field(description="Visible in both places")
        excluded: int = Field(description="[EXCLUDE] Internal implementation detail")
        custom_export: int = Field(
            description="[CUSTOM_EXPORT] Emitted by a custom formatter"
        )

    documented = {
        name: description
        for name, _, description in Example.document_properties(prefix="example")
    }
    assert documented == {
        "example_visible": "Visible in both places",
        "example_custom_export": "Emitted by a custom formatter",
    }
    assert description_excluded_from_column_docs(
        Example.model_fields["excluded"].description
    )
    assert description_excluded_from_flat_export(
        Example.model_fields["excluded"].description
    )
    assert description_excluded_from_flat_export(
        Example.model_fields["custom_export"].description
    )
    assert not description_excluded_from_column_docs(
        Example.model_fields["custom_export"].description
    )


def test_annotation_models_use_only_readable_description_markers():
    from plinder.data.annotations.aggregate_annotations import Entry, System
    from plinder.data.annotations.get_ligand_validation import ResidueListValidation
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain, Residue

    for model in (Entry, System, Ligand, Chain, Residue, ResidueListValidation):
        descriptions = model.get_descriptions_and_types()
        assert not any(
            str(description).lstrip().startswith("__")
            for description, _ in descriptions.values()
        )
        assert all(
            not str(description).lstrip().startswith("[")
            or str(description).lstrip().startswith(("[EXCLUDE]", "[CUSTOM_EXPORT]"))
            for description, _ in descriptions.values()
        )
        assert all(
            not description.startswith(("[EXCLUDE]", "[CUSTOM_EXPORT]"))
            for _, _, description in model.document_properties(prefix="test")
        )


def test_annotation_descriptions_follow_arrow_schema_order():
    import pyarrow as pa
    from plinder.data.annotations.aggregate_annotations import System
    from plinder.data.annotations.ligand_utils import Ligand

    schema = pa.schema(
        [
            ("ligand_id_legacy", pa.string()),
            ("ligand__members", pa.struct([("1.A", pa.list_(pa.int64()))])),
            ("ligand_member_asym_ids", pa.list_(pa.string())),
            ("system_id_legacy", pa.string()),
            ("pli_qcov__50__ligand__component", pa.string()),
            (
                "pli_qcov__50__ligand__directed_set_cover__is_centroid",
                pa.bool_(),
            ),
        ]
    )

    descriptions = docs.get_table_column_descriptions(
        table_name="annotation", schema=schema
    )

    assert descriptions["Name"].tolist() == schema.names
    assert descriptions["Type"].tolist() == [str(field.type) for field in schema]
    assert descriptions["Description"].str.len().gt(0).all()
    by_name = descriptions.set_index("Name")["Description"].to_dict()
    ligand_descriptions = {
        name: description
        for name, _, description in Ligand.document_properties("ligand")
    }
    system_descriptions = {
        name: description
        for name, _, description in System.document_properties("system")
    }
    assert by_name["ligand_id_legacy"] == ligand_descriptions["ligand_id_legacy"]
    assert by_name["ligand__members"] == ligand_descriptions["ligand__members"]
    assert (
        by_name["ligand_member_asym_ids"]
        == ligand_descriptions["ligand_member_asym_ids"]
    )
    assert by_name["system_id_legacy"] == system_descriptions["system_id_legacy"]


def test_table_descriptions_reject_missing_column_prose():
    import pyarrow as pa
    import pytest

    with pytest.raises(ValueError, match="undocumented_column"):
        docs.get_table_column_descriptions(
            table_name="annotation",
            schema=pa.schema([("undocumented_column", pa.string())]),
        )


def test_checked_in_descriptions_cover_every_table():
    from plinder.core.release import RELEASE_TABLES

    assert {path.stem for path in docs.TABLE_TSV_DIR.glob("*.tsv")} == set(
        RELEASE_TABLES
    )
    for table_name in RELEASE_TABLES:
        descriptions = docs.get_column_descriptions(table_name)
        assert not descriptions.empty
        assert list(descriptions.columns) == ["Name", "Type", "Description"]
        assert descriptions["Description"].notna().all()


def test_write_column_descriptions_uses_release_table_schemas(tmp_path, monkeypatch):
    import pyarrow as pa
    import pyarrow.parquet as pq

    release_dir = tmp_path / "release"
    table_path = release_dir / "index" / "entry_metadata.parquet"
    table_path.parent.mkdir(parents=True)
    pq.write_table(
        pa.table(
            {
                "entry_pdb_id": ["1abc"],
                "entry_source_taxonomy_ids": [[9606]],
            }
        ),
        table_path,
    )
    monkeypatch.setattr(
        docs,
        "RELEASE_TABLES",
        {
            "entry_metadata": {
                "artifact": "entry_metadata",
                "row_grain": "PDB entry",
                "primary_key": ("entry_pdb_id",),
            }
        },
    )
    output_dir = tmp_path / "descriptions"
    output_dir.mkdir()
    stale_path = output_dir / "stale.tsv"
    stale_path.write_text("Name\tType\tDescription\n")

    docs.write_column_descriptions(data_dir=release_dir, output_dir=output_dir)

    written = docs.get_column_descriptions("entry_metadata", description_dir=output_dir)
    assert written["Name"].tolist() == [
        "entry_pdb_id",
        "entry_source_taxonomy_ids",
    ]
    assert not stale_path.exists()
