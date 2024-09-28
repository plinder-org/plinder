# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest
from plinder.core import scores
from plinder.core.scores.protein import multi_query_protein_similarity


def test_query_index(read_plinder_mount):
    df = scores.query_index(columns=["system_id"], splits=["*"])
    assert len(df.index) == 57


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
