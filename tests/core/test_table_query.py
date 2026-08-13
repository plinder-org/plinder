# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from plinder.core import PlinderRelease, query_table


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
            "ligand_id": ["1abc__1__1.L", "2def__1__1.M", "3ghi__1__1.N"],
            "pocket_cluster": ["cluster_a", "cluster_b", "cluster_a"],
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
    return release


def test_query_table_filters_joined_fields_without_changing_ligand_grain(
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
        joins=["entry_metadata", "ligand_pocket_membership"],
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


def test_joined_sidecar_replaces_a_retired_repeated_column(
    local_release: PlinderRelease,
) -> None:
    result = query_table(
        "annotation",
        columns=["ligand_id", "entry_resolution"],
        joins=["entry_metadata"],
        release=local_release,
    )

    assert result["entry_resolution"].tolist()[:2] == [1.5, 2.5]
    assert pd.isna(result["entry_resolution"].iloc[2])


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
        joins=["entry_metadata"],
        release=local_release,
    )

    assert result.to_dict("records") == [
        {"linked_structure_id": "2def_B", "entry_resolution": 2.5}
    ]


def test_query_table_rejects_joins_that_would_change_row_grain(
    local_release: PlinderRelease,
) -> None:
    with pytest.raises(ValueError, match="grain-preserving joins"):
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
