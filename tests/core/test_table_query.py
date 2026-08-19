# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from plinder.core import PlinderRelease, query_table
from plinder.core.index.query import DISABLED_ANNOTATION_COLUMNS


def _write_table(
    release: PlinderRelease,
    artifact: str,
    data: dict[str, list[object]],
) -> None:
    path = release.path(artifact)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data).to_parquet(path, index=False)


@pytest.fixture
def local_release(tmp_path: Path) -> PlinderRelease:
    release = PlinderRelease(tmp_path)
    _write_table(
        release,
        "annotation_table",
        {
            "ligand_id": ["1abc__1__1.L", "2def__1__1.M", "3ghi__1__1.N"],
            "system_id": [
                "1abc__1__1.A__1.L",
                "2def__1__1.B__1.M",
                "3ghi__1__1.C__1.N",
            ],
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "entry_resolution": [99.0, 99.0, 99.0],
            "entry_release_date": ["1900-01-01"] * 3,
            "system_has_binding_affinity": [False, True, False],
            "ligand_binding_affinity": [None, "Kd=10nM", None],
        },
    )
    _write_table(
        release,
        "system_validation",
        {
            "system_id": [
                "1abc__1__1.A__1.L",
                "2def__1__1.B__1.M",
                "3ghi__1__1.C__1.N",
            ],
            "system_pocket_validation_average_rscc": [0.95, 0.8, None],
        },
    )
    _write_table(
        release,
        "entry_metadata",
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "entry_resolution": [1.5, 2.5, None],
            "entry_release_date": ["2020-01-01", "2021-01-01", "2022-01-01"],
        },
    )
    _write_table(
        release,
        "ligand_pocket_membership",
        {
            "ligand_id": ["1abc__1__1.L", "2def__1__1.M"],
            "system_id": ["1abc__1__1.A__1.L", "2def__1__1.B__1.M"],
            "pocket_cluster": ["cluster_a", "cluster_b"],
        },
    )
    _write_table(
        release,
        "linked_apo_structures",
        {
            "reference_system_id": ["1abc__1__1.A__1.L"],
            "linked_structure_id": ["2def_B"],
            "source_entry_id": ["2def"],
            "rank": [1],
        },
    )
    _write_table(
        release,
        "ligand_mmp_pairs",
        {
            "ligand_smiles_id_1": [0],
            "ligand_smiles_id_2": [1],
            "transformation": ["*O>>*N"],
            "shared_core_smiles": ["*C"],
            "ligand_1_shared_core_fraction": [1 / 3],
        },
    )
    return release


def test_query_table_infers_related_tables_without_changing_ligand_rows(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=[
            "ligand_id",
            "entry_resolution",
            "entry_release_date",
            "pocket_cluster",
        ],
        filters=[
            ("entry_resolution", "<=", 2.0),
            [
                ("pocket_cluster", "==", "cluster_a"),
                ("ligand_id", "==", "not a ligand"),
            ],
        ],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {
            "ligand_id": "1abc__1__1.L",
            "entry_resolution": 1.5,
            "entry_release_date": "2020-01-01",
            "pocket_cluster": "cluster_a",
        }
    ]


def test_requested_entry_metadata_is_joined_automatically(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=["ligand_id", "entry_resolution"],
        release=local_release,
    )

    assert result["entry_resolution"].tolist()[:2] == [1.5, 2.5]
    assert pd.isna(result["entry_resolution"].iloc[2])


def test_requested_system_validation_is_joined_automatically(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=["ligand_id", "system_pocket_validation_average_rscc"],
        filters=[("system_pocket_validation_average_rscc", ">=", 0.9)],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {
            "ligand_id": "1abc__1__1.L",
            "system_pocket_validation_average_rscc": 0.95,
        }
    ]


def test_sparse_sidecar_does_not_erase_base_identifiers(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=["ligand_id", "system_id", "pocket_cluster"],
        filters=[("system_id", "==", "3ghi__1__1.C__1.N")],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {
            "ligand_id": "3ghi__1__1.N",
            "system_id": "3ghi__1__1.C__1.N",
            "pocket_cluster": None,
        }
    ]


def test_annotation_release_dates_come_from_entry_metadata(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=["ligand_id", "entry_release_date"],
        filters=[("entry_release_date", ">=", "2022-01-01")],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {
            "ligand_id": "3ghi__1__1.N",
            "entry_release_date": "2022-01-01",
        }
    ]


def test_annotation_adds_entry_metadata_only_when_requested(
    local_release: PlinderRelease,
) -> None:
    result = query_table("annotation", release=local_release)

    assert "entry_release_date" not in result.columns
    assert "entry_resolution" not in result.columns
    requested = query_table(
        "annotation",
        columns=["ligand_id", "entry_resolution"],
        release=local_release,
    )
    assert requested["entry_resolution"].tolist()[:2] == [1.5, 2.5]
    assert pd.isna(requested["entry_resolution"].iloc[2])


def test_annotation_binding_affinity_columns_are_disabled(
    local_release: PlinderRelease,
) -> None:
    default = query_table("annotation", release=local_release)
    assert DISABLED_ANNOTATION_COLUMNS.isdisjoint(default.columns)

    with pytest.raises(ValueError, match="binding_affinity columns are disabled"):
        query_table(
            "annotation",
            columns=["ligand_binding_affinity"],
            release=local_release,
        )

    joined = query_table(
        "ligand_pocket_membership",
        joins=["annotation"],
        release=local_release,
    )
    assert DISABLED_ANNOTATION_COLUMNS.isdisjoint(joined.columns)

    with pytest.raises(ValueError, match="binding_affinity columns are disabled"):
        query_table(
            "ligand_pocket_membership",
            columns=["ligand_binding_affinity"],
            joins=["annotation"],
            release=local_release,
        )


def test_query_table_uses_bound_filter_parameters(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=["ligand_id"],
        filters=[("ligand_id", "==", "1abc__1__1.L' OR true --")],
        release=local_release,
    )

    assert result.empty


def test_query_linked_apo_with_source_entry_metadata(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "linked_apo_structures",
        columns=["linked_structure_id", "entry_resolution"],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {"linked_structure_id": "2def_B", "entry_resolution": 2.5}
    ]


def test_query_ligand_mmp_pairs(local_release: PlinderRelease) -> None:
    result = query_table(
        "ligand_mmp_pairs",
        columns=[
            "ligand_smiles_id_1",
            "ligand_smiles_id_2",
            "ligand_1_shared_core_fraction",
        ],
        filters=[("ligand_1_shared_core_fraction", ">=", 0.3)],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {
            "ligand_smiles_id_1": 0,
            "ligand_smiles_id_2": 1,
            "ligand_1_shared_core_fraction": pytest.approx(1 / 3),
        }
    ]


def test_query_table_rejects_unregistered_related_tables(
    local_release: PlinderRelease,
) -> None:
    with pytest.raises(ValueError, match="available related tables"):
        query_table(
            "entry_metadata",
            joins=["annotation"],
            release=local_release,
        )


def test_query_table_rejects_unavailable_columns(
    local_release: PlinderRelease,
) -> None:
    with pytest.raises(ValueError, match="columns .* are unavailable"):
        query_table(
            "annotation",
            columns=["not_a_column"],
            release=local_release,
        )


def test_query_table_requires_a_choice_for_ambiguous_related_columns(
    tmp_path: Path,
) -> None:
    release = PlinderRelease(tmp_path)
    _write_table(
        release,
        "annotation_table",
        {"ligand_id": ["1abc__1__1.L"], "entry_pdb_id": ["1abc"]},
    )
    _write_table(
        release,
        "entry_metadata",
        {
            "entry_pdb_id": ["1abc"],
            "metadata_only": [3],
            "shared_value": [1],
        },
    )
    _write_table(
        release,
        "entry_sources",
        {"entry_pdb_id": ["1abc"], "shared_value": [2]},
    )

    with pytest.raises(ValueError, match="available from multiple related tables"):
        query_table(
            "annotation",
            columns=["ligand_id", "shared_value"],
            release=release,
        )

    inferred = query_table(
        "annotation",
        columns=["ligand_id", "shared_value", "metadata_only"],
        release=release,
    )
    assert inferred.to_dict("records") == [
        {
            "ligand_id": "1abc__1__1.L",
            "shared_value": 1,
            "metadata_only": 3,
        }
    ]

    selected = query_table(
        "annotation",
        columns=["ligand_id", "shared_value"],
        joins=["entry_metadata"],
        release=release,
    )
    assert selected.to_dict("records") == [
        {"ligand_id": "1abc__1__1.L", "shared_value": 1}
    ]
