import pytest

from plinder.core import scores


def test_query_protein_similarity(read_plinder_mount):

    df = scores.query_protein_similarity(
        search_db="holo",
        filters=[
            ("metric", "==", "pocket_lddt"),
            ("similarity", ">=", 90),
        ]
    )
    assert len(df.index)


def test_query_protein_similarity_empty(read_plinder_mount):
    df = scores.query_protein_similarity(
        search_db="holo",
        filters=[],
    )
    assert df is None


def test_query_protein_similarity_removes_search_db(read_plinder_mount):
    df = scores.query_protein_similarity(
        search_db="holo",
        filters=[
            ("search_db", "==", "holo"),
            ("metric", "==", "pocket_lddt"),
            ("similarity", ">=", 90),
        ],
    )
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


def test_query_ligand_similarity(read_plinder_mount):

    df = scores.query_ligand_similarity(
        filters=[
            ("query_ligand_id", "<", "100"),
        ]
    )
    assert len(df.index)


def test_query_ligand_similarity_empty(read_plinder_mount):

    df = scores.query_ligand_similarity(
        filters=[]
    )
    assert df is None
