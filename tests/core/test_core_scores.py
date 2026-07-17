# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
import pytest
from plinder.core import scores
from plinder.core.scores import index as index_module
from plinder.core.scores import ligand as ligand_module
from plinder.core.scores.protein import multi_query_protein_similarity


def test_query_index(read_plinder_mount):
    df = scores.query_index(columns=["system_id"], splits=["*"])
    assert len(df.index) == 57


@pytest.mark.parametrize(
    "system_id, correct_release_date",
    [
        ("6cex__1__1.D__1.M", "2018-04-04"),
        ("8grn__1__1.A__1.C", "2023-06-14"),
    ],
)
def test_entry_release_date(system_id, correct_release_date):
    df = scores.query_index(
        splits=["*"],
        columns=["entry_release_date"],
        filters=[("system_id", "==", system_id)],
    )
    assert df.iloc[0].entry_release_date == correct_release_date


def test_query_protein_similarity(read_plinder_mount):
    df = scores.query_protein_similarity(
        search_db="holo",
        filters=[
            ("metric", "==", "pocket_lddt"),
            ("similarity", ">=", 90),
        ],
    )
    assert df is not None
    assert len(df.index)


def test_query_protein_similarity_empty(read_plinder_mount):
    with pytest.raises(ValueError):
        scores.query_protein_similarity(
            search_db="holo",
            filters=[],
        )


def test_query_protein_similarity_removes_search_db(read_plinder_mount):
    df = scores.query_protein_similarity(
        search_db="holo",
        filters=[
            ("search_db", "==", "holo"),
            ("metric", "==", "pocket_lddt"),
            ("similarity", ">=", 90),
        ],
    )
    assert df is not None
    assert len(df.index)


def test_query_protein_similarity_raises(read_plinder_mount):
    with pytest.raises(ValueError):
        scores.query_protein_similarity(
            search_db="test",
            filters=[
                ("search_db", "==", "holo"),
                ("metric", "==", "pocket_lddt"),
                ("similarity", ">=", 90),
            ],
        )


def test_query_protein_cross_similarity(read_plinder_mount):
    df = scores.cross_protein_similarity(
        query_systems=["8t49__1__1.G__1.AB", "6cex__1__1.D__1.M"],
        target_systems=["4r2g__3__1.P__1.AB", "4ln4__1__1.F__1.U"],
        metric="pocket_lddt",
    )
    assert len(df.index)


def test_query_ligand_similarity(read_plinder_mount):
    df = scores.query_ligand_similarity(
        filters=[
            ("query_ligand_id", "<", "100"),
        ]
    )
    assert df is not None
    assert len(df.index)


def test_query_ligand_similarity_empty(read_plinder_mount):
    with pytest.raises(ValueError):
        scores.query_ligand_similarity(filters=[])


def test_query_ligand_cross_similarity(read_plinder_mount):
    df = scores.cross_ligand_similarity(
        query_ligands=[29, 51], target_ligands=[49918, 36689]
    )
    assert len(df.index)


def test_v3_ligand_cross_similarity_maps_nodes_through_index(monkeypatch):
    class DataConfig:
        plinder_iteration = "v3"

    class Config:
        data = DataConfig()

    calls = []

    def fake_query_index(*, columns, filters, splits):
        calls.append((columns, filters, splits))
        return pd.DataFrame(
            {
                "system_id": ["1aaa__1__1.A__1.X", "2bbb__1__1.B__1.Y"],
                "ligand_smiles_id": [0, 0],
            }
        )

    monkeypatch.setattr(ligand_module, "get_config", Config)
    monkeypatch.setattr(index_module, "query_index", fake_query_index)
    result = ligand_module.map_cross_similarity(
        pd.DataFrame(
            {
                "query_ligand_id": [0],
                "target_ligand_id": [1],
                "tanimoto_similarity_ecfp4_1024": [95.0],
            }
        ),
        target_ligands={1},
        metric="tanimoto_similarity_ecfp4_1024",
    )

    assert set(result["system_id"]) == {
        "1aaa__1__1.A__1.X",
        "2bbb__1__1.B__1.Y",
    }
    assert calls and calls[0][0] == ["system_id", "ligand_smiles_id"]


def test_v3_ligand_cross_similarity_returns_empty_without_querying_index(
    monkeypatch,
):
    class DataConfig:
        plinder_iteration = "v3"

    class Config:
        data = DataConfig()

    monkeypatch.setattr(ligand_module, "get_config", Config)
    monkeypatch.setattr(
        index_module,
        "query_index",
        lambda **_kwargs: pytest.fail("empty similarities must not query the index"),
    )

    result = ligand_module.map_cross_similarity(
        pd.DataFrame(
            columns=[
                "query_ligand_id",
                "target_ligand_id",
                "tanimoto_similarity_ecfp4_1024",
            ]
        ),
        target_ligands={1},
        metric="tanimoto_similarity_ecfp4_1024",
    )

    assert result.empty
    assert result.columns.tolist() == ["system_id", "tanimoto_similarity_ecfp4_1024"]


def test_query_links(read_plinder_mount):
    system_id = "4dd7__1__1.A__1.B"
    df = scores.query_links(filters=[("reference_system_id", "==", system_id)])
    assert len(df.index)


def test_query_links_columns(read_plinder_mount):
    system_id = "4dd7__1__1.A__1.B"
    df = scores.query_links(
        columns=["reference_system_id"],
        filters=[("reference_system_id", "==", system_id)],
    )
    assert len(df.index)
    assert "reference_system_id" in df.columns
    assert "kind" in df.columns


def test_multi_query_protein_similarity(read_plinder_mount):
    system_id = "8t49__1__1.G__1.AB"
    filter_criteria: dict[str, int] = {
        "protein_fident_qcov_weighted_sum": 0,
        "pocket_lddt": 90,
    }
    df = multi_query_protein_similarity(
        system_id=system_id,
        search_db="holo",
        filter_criteria=filter_criteria,
        splits=["*"],
    )
    assert len(df.index)
    assert all(k in df.columns for k in filter_criteria)
