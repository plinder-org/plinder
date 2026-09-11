# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import pytest
from plinder.data import docs


def _import_tablegen(monkeypatch):
    import importlib
    import sys
    from pathlib import Path
    from types import ModuleType

    itables = ModuleType("itables")
    itables.to_html_datatable = lambda frame, **_: frame.to_html(
        index=False, escape=False
    )
    monkeypatch.setitem(sys.modules, "itables", itables)
    monkeypatch.delitem(sys.modules, "tablegen", raising=False)
    repository = Path(__file__).resolve().parents[2]
    monkeypatch.syspath_prepend(str(repository / "docs"))
    return importlib.import_module("tablegen")


def test_tablegen_renders_checked_in_table_descriptions(tmp_path, monkeypatch):
    tablegen = _import_tablegen(monkeypatch)

    description_dir = tmp_path / "column_descriptions"
    table_dir = description_dir / "tables"
    table_dir.mkdir(parents=True)
    (table_dir / "alpha.tsv").write_text(
        "Name\tType\tDescription\n"
        "entry<id>\tlist<element: string>\tApply <operation> & keep metadata\n",
        encoding="utf-8",
    )
    (table_dir / "beta.tsv").write_text(
        "Name\tType\tDescription\n" "score\tdouble\tSimilarity score\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "table.html"

    tablegen.generate_table(description_dir, output_path)

    html = output_path.read_text(encoding="utf-8")
    assert '<section class="release-column-table" id="columns-alpha">' in html
    assert '<section class="release-column-table" id="columns-beta">' in html
    assert "<h3>alpha</h3>" in html
    assert "<h3>beta</h3>" in html
    assert "<code>" not in html
    assert ">entry&lt;id&gt;<" in html
    assert ">list&lt;element: string&gt;<" in html
    assert "Apply &lt;operation&gt; &amp; keep metadata" in html
    assert ">score<" in html
    assert "Similarity score" in html
    assert html.count("<table") == 2
    assert all(line.strip() for line in html.splitlines())


def test_tablegen_rejects_invalid_description_columns(tmp_path, monkeypatch):
    import pytest

    tablegen = _import_tablegen(monkeypatch)

    description_dir = tmp_path / "column_descriptions"
    table_dir = description_dir / "tables"
    table_dir.mkdir(parents=True)
    (table_dir / "broken.tsv").write_text(
        "Name\tDescription\nentry_id\tStable entry identifier\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="must contain columns"):
        tablegen.generate_table(description_dir, tmp_path / "table.html")


def test_dataset_preserves_notebook_link_target():
    from pathlib import Path

    dataset_doc = (
        Path(__file__).resolve().parents[2] / "docs" / "dataset.md"
    ).read_text()

    assert "(annotation-tables-index)=" in dataset_doc


def test_ligand_cluster_column_descriptions():
    import pandas as pd

    columns = [
        "tanimoto_similarity_ecfp4_1024__70__ligand__set_cover",
        "tanimoto_similarity_ecfp4_1024__70__ligand__set_cover__is_centroid",
        "pocket_qcov__50__ligand__directed_set_cover",
        "pocket_qcov__50__ligand__directed_set_cover__is_centroid",
        "pocket_qcov__50__ligand__directed_set_cover__coverage_count",
        "pocket_qcov__50__ligand__directed_set_cover__coverage_fraction",
        "interface_side_qcov__70__chain_1_directed_set_cover",
    ]

    rows = docs.get_cluster_column_descriptions(pd.DataFrame(columns=columns))

    assert [row[0] for row in rows] == columns
    descriptions = {name: description for name, _, description in rows}
    assert "ligand-level set cover" in descriptions[columns[0]]
    assert "direct threshold-qualified edge" in descriptions[columns[0]]
    assert "published centroid" in descriptions[columns[1]]
    assert "ligand-level directed set cover" in descriptions[columns[2]]
    assert "query-to-centroid score" in descriptions[columns[2]]
    assert "published centroid" in descriptions[columns[3]]
    assert "Number of directed-cover query nodes" in descriptions[columns[4]]
    assert "Fraction of its directed weak component" in descriptions[columns[5]]
    assert "chain 1 directed set cover" in descriptions[columns[6]]


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

    annotation_schema = pa.schema(
        [
            ("ligand_id_legacy", pa.string()),
            ("ligand__members", pa.struct([("1.A", pa.list_(pa.int64()))])),
            ("ligand_member_asym_ids", pa.list_(pa.string())),
            ("system_id_legacy", pa.string()),
        ]
    )
    cluster_schema = pa.schema(
        [
            ("ligand_id", pa.string()),
            (
                "tanimoto_similarity_ecfp4_1024__50__ligand__set_cover",
                pa.string(),
            ),
            (
                "pli_qcov__50__ligand__directed_set_cover__is_centroid",
                pa.bool_(),
            ),
        ]
    )

    descriptions = docs.get_table_column_descriptions(
        table_name="annotation", schema=annotation_schema
    )
    cluster_descriptions = docs.get_table_column_descriptions(
        table_name="ligand_clusters", schema=cluster_schema
    )

    assert descriptions["Name"].tolist() == annotation_schema.names
    assert descriptions["Type"].tolist() == [
        str(field.type) for field in annotation_schema
    ]
    assert descriptions["Description"].str.len().gt(0).all()
    assert cluster_descriptions["Name"].tolist() == cluster_schema.names
    assert cluster_descriptions["Description"].str.len().gt(0).all()
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


def test_annotation_descriptions_reject_repeated_entry_metadata():
    import pyarrow as pa
    import pytest

    with pytest.raises(ValueError, match="columns owned by 'entry_metadata'"):
        docs.get_table_column_descriptions(
            table_name="annotation",
            schema=pa.schema(
                [
                    ("entry_pdb_id", pa.string()),
                    ("entry_resolution", pa.float64()),
                ]
            ),
        )


@pytest.mark.parametrize(
    "column",
    [
        "system_id_no_biounit",
        "system_ligand_chains",
        "ligand_rdkit_canonical_smiles",
        "system_protein_chains_auth_id",
        "system_protein_chains_validation_average_rsr",
        "system_ligand_validation_average_rsr",
    ],
)
def test_annotation_descriptions_reject_moved_or_retired_columns(column: str):
    import pyarrow as pa

    with pytest.raises(ValueError, match="owned by sidecars or retired"):
        docs.get_table_column_descriptions(
            table_name="annotation",
            schema=pa.schema([("ligand_id", pa.string()), (column, pa.string())]),
        )


def test_system_validation_descriptions_reject_unrelated_columns():
    import pyarrow as pa

    with pytest.raises(ValueError, match="non-validation columns"):
        docs.get_table_column_descriptions(
            table_name="system_validation",
            schema=pa.schema(
                [("system_id", pa.string()), ("ligand_smiles", pa.string())]
            ),
        )


def test_table_descriptions_reject_retired_cover_modes():
    import pyarrow as pa
    import pytest

    invalid_schemas = [
        (
            "ligand_clusters",
            pa.schema(
                [
                    (
                        "tanimoto_similarity_ecfp4_1024__50__ligand__"
                        "directed_set_cover",
                        pa.string(),
                    )
                ]
            ),
        ),
        (
            "ligand_clusters",
            pa.schema([("pocket_qcov__50__ligand__community", pa.string())]),
        ),
        (
            "interface_clusters",
            pa.schema([("interface_qcov__50__component", pa.string())]),
        ),
        (
            "interface_clusters",
            pa.schema(
                [
                    (
                        "interface_qcov__50__directed_set_cover__is_centroid",
                        pa.string(),
                    )
                ]
            ),
        ),
        (
            "interface_clusters",
            pa.schema(
                [
                    (
                        "interface_qcov__50__directed_set_cover__coverage_count",
                        pa.int64(),
                    )
                ]
            ),
        ),
        (
            "interface_clusters",
            pa.schema(
                [
                    (
                        "interface_side_qcov__50__chain_1_directed_set_cover__"
                        "coverage_fraction",
                        pa.float64(),
                    )
                ]
            ),
        ),
    ]
    for table_name, schema in invalid_schemas:
        with pytest.raises(ValueError, match="not published by the current pipeline"):
            docs.get_table_column_descriptions(
                table_name=table_name,
                schema=schema,
            )


def test_table_descriptions_treat_every_chemical_metric_alike():
    import pyarrow as pa
    import pytest
    from plinder.core.scores.metrics import (
        CHEMICAL_CLUSTER_METRICS,
        CHEMICAL_CLUSTER_SUMMARY_COLUMNS,
    )

    assert set(CHEMICAL_CLUSTER_SUMMARY_COLUMNS) == set(CHEMICAL_CLUSTER_METRICS)
    fields = [("ligand_id", pa.string())]
    for metric in CHEMICAL_CLUSTER_METRICS:
        summary = CHEMICAL_CLUSTER_SUMMARY_COLUMNS[metric]
        fields.extend(
            [
                (f"{metric}__90__ligand__set_cover", pa.string()),
                (summary, pa.string()),
                (f"{summary}_num_pdb_ids", pa.int32()),
                (f"{metric}__90__ligand__set_cover__is_centroid", pa.bool_()),
            ]
        )

    descriptions = docs.get_table_column_descriptions(
        table_name="ligand_clusters", schema=pa.schema(fields)
    )

    assert descriptions["Name"].tolist() == [name for name, _ in fields]
    by_name = descriptions.set_index("Name")["Description"]
    assert (
        "MHFP6/2048" in by_name["jaccard_similarity_mhfp6_2048__90__ligand__set_cover"]
    )
    assert "MHFP6/2048" in by_name["ligand_jaccard_mhfp6_2048_90_cluster"]
    assert "ECFP4/1024" in by_name["ligand_tanimoto_ecfp4_1024_90_cluster"]
    for metric in CHEMICAL_CLUSTER_METRICS:
        with pytest.raises(ValueError, match="not published by the current pipeline"):
            docs.get_table_column_descriptions(
                table_name="ligand_clusters",
                schema=pa.schema(
                    [(f"{metric}__50__ligand__directed_set_cover", pa.string())]
                ),
            )
        with pytest.raises(ValueError):
            docs.get_table_column_descriptions(
                table_name="annotation",
                schema=pa.schema(
                    [(CHEMICAL_CLUSTER_SUMMARY_COLUMNS[metric], pa.string())]
                ),
            )


def test_table_descriptions_accept_published_interface_cover_columns():
    import pyarrow as pa

    names = [
        "interface_qcov__50__directed_set_cover",
        "interface_side_qcov__50__chain_1_directed_set_cover",
        "interface_side_qcov__50__chain_2_directed_set_cover",
    ]
    descriptions = docs.get_table_column_descriptions(
        table_name="interface_clusters",
        schema=pa.schema([(name, pa.string()) for name in names]),
    )

    assert descriptions["Name"].tolist() == names


@pytest.mark.parametrize(
    ("table_name", "column"),
    [
        ("entry_metadata", "entry_failed_assembly_ids"),
        ("annotation", "ligand_failed_interaction_types"),
    ],
)
def test_failure_descriptions_match_model_fields(table_name, column):
    import pyarrow as pa

    generated = docs.get_table_column_descriptions(
        table_name=table_name,
        schema=pa.schema([(column, pa.list_(pa.field("element", pa.string())))]),
    ).set_index("Name")
    checked_in = docs.get_column_descriptions(table_name).set_index("Name")
    assert checked_in.loc[column].to_dict() == generated.loc[column].to_dict()


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


def test_checked_in_cluster_descriptions_match_published_cover_modes():
    from plinder.core.scores.metrics import (
        DEFAULT_CLUSTER_METRICS,
        is_chemical_cluster_metric,
    )

    ligand_names = docs.get_column_descriptions("ligand_clusters")["Name"].tolist()
    metric_names = set(DEFAULT_CLUSTER_METRICS)
    ligand_cluster_names = [
        name for name in ligand_names if name.split("__", maxsplit=1)[0] in metric_names
    ]
    assert ligand_cluster_names
    assert not any(
        "__component" in name or "__community" in name for name in ligand_cluster_names
    )
    for name in ligand_cluster_names:
        if is_chemical_cluster_metric(name.split("__", maxsplit=1)[0]):
            assert "__ligand__set_cover" in name
            assert "__directed_set_cover" not in name
        else:
            assert "__ligand__directed_set_cover" in name

    interface_names = docs.get_column_descriptions("interface_clusters")[
        "Name"
    ].tolist()
    interface_cluster_names = [
        name
        for name in interface_names
        if name.startswith(("interface_qcov__", "interface_side_qcov__"))
    ]
    assert interface_cluster_names
    assert all("directed_set_cover" in name for name in interface_cluster_names)
    assert not any(
        "__component" in name or "__community" in name
        for name in interface_cluster_names
    )


def test_checked_in_annotation_keeps_only_the_entry_join_key():
    annotation_names = docs.get_column_descriptions("annotation")["Name"].tolist()
    entry_names = [name for name in annotation_names if name.startswith("entry_")]

    assert entry_names == ["entry_pdb_id"]


def test_linked_apo_descriptions_match_release_schema():
    from plinder.core.utils.schemas import STRUCTURE_LINK_SCHEMA

    descriptions = docs.get_column_descriptions("linked_apo_structures")

    assert descriptions["Name"].tolist() == STRUCTURE_LINK_SCHEMA.names
    generated = docs.get_table_column_descriptions(
        table_name="linked_apo_structures",
        schema=STRUCTURE_LINK_SCHEMA,
    )
    assert generated["Name"].tolist() == STRUCTURE_LINK_SCHEMA.names


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


def test_exported_column_descriptions_are_unique_within_each_model():
    """Each exported column's description is its only user-facing explanation,
    so two columns of one model must never share the same text. This catches
    copy-pasted docstrings such as a ``proper_*`` property repeating its
    unfiltered sibling, which the annotation table docs then publish verbatim."""
    from collections import defaultdict

    from plinder.data.annotations.aggregate_annotations import Entry, System
    from plinder.data.annotations.get_ligand_validation import (
        EntryValidation,
        ResidueListValidation,
    )
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain

    for model in (Entry, System, Ligand, Chain, ResidueListValidation, EntryValidation):
        columns_by_description: dict[str, list[str]] = defaultdict(list)
        for name, _, description in model.document_properties(model.__name__.lower()):
            columns_by_description[description].append(name)
        duplicated = {
            description: names
            for description, names in columns_by_description.items()
            if len(names) > 1
        }
        assert not duplicated, f"{model.__name__}: {duplicated}"
