# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
import pytest
from plinder.core import query_table, scores
from plinder.core.scores import ligand as ligand_module
from plinder.core.scores import protein as protein_module
from plinder.core.scores.protein import multi_query_protein_similarity


def test_custom_scoring_public_api():
    assert callable(scores.score_custom_cif_files)
    assert callable(scores.score_custom_sequence_file)
    assert scores.CustomProteinSearchConfig().max_seqs == 10_000


@pytest.fixture
def current_ligand_scores(read_plinder_mount, tmp_path, monkeypatch):
    source = read_plinder_mount / "ligand_scores" / "ligand_scores.parquet"
    scores_dir = tmp_path / "ligand_scores"
    scores_dir.mkdir()
    pd.read_parquet(source).rename(
        columns={
            "tanimoto_similarity_max": "tanimoto_similarity_ecfp4_1024",
        }
    ).to_parquet(scores_dir / "ligand_scores.parquet", index=False)

    original_fetch = ligand_module.PlinderRelease.fetch

    def fetch(release, name, **parameters):
        if name == "ligand_scores":
            return scores_dir
        return original_fetch(release, name, **parameters)

    monkeypatch.setattr(ligand_module.PlinderRelease, "fetch", fetch)
    return scores_dir


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_annotation_table():
    df = query_table("annotation", columns=["system_id"])
    assert len(df.index) == 57
    assert "split" not in df.columns


@pytest.mark.usefixtures("read_plinder_mount")
@pytest.mark.parametrize(
    "system_id, correct_release_date",
    [
        ("6cex__1__1.D__1.M", "2018-04-04"),
        ("8grn__1__1.A__1.C", "2023-06-14"),
    ],
)
def test_entry_release_date(system_id, correct_release_date):
    df = query_table(
        "annotation",
        columns=["entry_release_date"],
        filters=[("system_id", "==", system_id)],
    )
    assert df.iloc[0].entry_release_date == correct_release_date


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_protein_similarity():
    df = scores.query_protein_similarity(
        search_db="holo",
        filters=[
            ("metric", "==", "pocket_lddt"),
            ("similarity", ">=", 90),
        ],
    )
    assert df is not None
    assert len(df.index)


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_protein_similarity_empty():
    with pytest.raises(ValueError):
        scores.query_protein_similarity(
            search_db="holo",
            filters=[],
        )


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_protein_similarity_removes_search_db_without_mutating_filters():
    filters = [
        ("search_db", "==", "holo"),
        ("metric", "==", "pocket_lddt"),
        ("similarity", ">=", 90),
    ]
    df = scores.query_protein_similarity(search_db="holo", filters=filters)

    assert len(df.index)
    assert filters[0] == ("search_db", "==", "holo")


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_protein_similarity_raises():
    with pytest.raises(ValueError):
        scores.query_protein_similarity(
            search_db="test",
            filters=[
                ("search_db", "==", "holo"),
                ("metric", "==", "pocket_lddt"),
                ("similarity", ">=", 90),
            ],
        )


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_protein_cross_similarity():
    df = scores.cross_protein_similarity(
        query_systems=["8t49__1__1.G__1.AB", "6cex__1__1.D__1.M"],
        target_systems=["4r2g__3__1.P__1.AB", "4ln4__1__1.F__1.U"],
        metric="pocket_lddt",
    )
    assert len(df.index)


def test_query_ligand_similarity(current_ligand_scores):
    df = scores.query_ligand_similarity(
        filters=[
            ("query_ligand_id", "<", 100),
        ]
    )
    assert df is not None
    assert len(df.index)


@pytest.mark.usefixtures("read_plinder_mount")
def test_query_ligand_similarity_empty():
    with pytest.raises(ValueError):
        scores.query_ligand_similarity(filters=[])


def test_query_ligand_cross_similarity(current_ligand_scores, monkeypatch):
    monkeypatch.setattr(
        ligand_module,
        "map_cross_similarity",
        lambda frame, _target_ligands, _metric: frame,
    )
    df = scores.cross_ligand_similarity(
        query_ligands=["29", "51"], target_ligands=["49918", "36689"]
    )
    assert len(df.index)


def test_query_ligand_cross_similarity_rejects_non_numeric_ids(
    current_ligand_scores,
):
    with pytest.raises(ValueError, match="ligand IDs must be integers"):
        scores.cross_ligand_similarity(
            query_ligands=["not-an-id"], target_ligands=["49918"]
        )


def test_ligand_cross_similarity_maps_nodes_through_annotation(monkeypatch):
    calls = []

    def fake_query_table(table_name, *, columns, filters):
        calls.append((table_name, columns, filters))
        return pd.DataFrame(
            {
                "system_id": ["1aaa__1__1.A__1.X", "2bbb__1__1.B__1.Y"],
                "ligand_smiles_id": [0, 0],
            }
        )

    monkeypatch.setattr(ligand_module, "query_table", fake_query_table)
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
    assert calls and calls[0][:2] == (
        "annotation",
        ["system_id", "ligand_smiles_id"],
    )


def test_ligand_cross_similarity_returns_empty_without_querying_annotation(
    monkeypatch,
):
    monkeypatch.setattr(
        ligand_module,
        "query_table",
        lambda *_args, **_kwargs: pytest.fail(
            "empty similarities must not query annotation"
        ),
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


@pytest.mark.usefixtures("read_plinder_mount")
def test_multi_query_protein_similarity():
    system_id = "8t49__1__1.G__1.AB"
    filter_criteria: dict[str, int] = {
        "protein_fident_qcov_weighted_sum": 0,
        "pocket_lddt": 90,
    }
    df = multi_query_protein_similarity(
        system_id=system_id,
        search_db="holo",
        filter_criteria=filter_criteria,
    )
    assert len(df.index)
    assert all(k in df.columns for k in filter_criteria)


def test_multi_query_protein_similarity_requires_every_metric(monkeypatch):
    monkeypatch.setattr(
        protein_module,
        "query_table",
        lambda *_args, **_kwargs: pd.DataFrame({"system_id": ["target"]}),
    )
    monkeypatch.setattr(
        protein_module,
        "query_protein_similarity",
        lambda **_kwargs: pd.DataFrame(
            {
                "query_system": ["query"],
                "target_system": ["target"],
                "metric": ["pocket_lddt"],
                "similarity": [95],
            }
        ),
    )

    result = multi_query_protein_similarity(
        system_id="query",
        search_db="holo",
        filter_criteria={"pocket_lddt": 90, "pocket_fident": 50},
    )

    assert result.empty
