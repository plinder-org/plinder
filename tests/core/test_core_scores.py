# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
import pytest

from plinder.core import PlinderRelease, query_table, scores
from plinder.core.scores import ligand as ligand_module


def test_custom_scoring_public_api():
    assert callable(scores.score_custom_cif_files)
    assert callable(scores.score_custom_sequence_file)
    assert scores.CustomProteinSearchConfig().max_seqs == 10_000


@pytest.fixture
def similarity_release(tmp_path):
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "chain_asym_id": ["A", "B"],
            "chain_auth_id": ["A", "B"],
            "selected_residue_numbers": [[10, 20, 30], [5, 6, 7, 8]],
            "selected_residue_indices": [[0, 1, 2], [0, 1, 2, 3]],
        }
    ).to_parquet(index / "alignment_chain_lookup.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "2def", "2def"],
            "chain_asym_id": ["A", "A", "B", "B"],
            "residue_label_seq_id": [10, 20, 6, 8],
        }
    ).to_parquet(index / "ligand_pocket_residues.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "interface_chain_1": ["1.A", "1.B"],
            "interface_chain_1_residue_numbers": [[20, 30], [8]],
            "interface_chain_2": ["1.C", "1.D"],
            "interface_chain_2_residue_numbers": [[40], [9]],
        }
    ).to_parquet(index / "interface_annotation_table.parquet", index=False)
    exports = tmp_path / "exports"
    exports.mkdir()
    pd.DataFrame(
        {
            "query_system": ["query", "query"],
            "query_ligand_id": ["query__1.A", "query__1.B"],
            "target_system": ["target", "other"],
            "target_ligand_id": ["target__1.X", "other__1.Y"],
            "pocket_qcov": [90, 40],
            "pocket_fident_qcov": [80, 20],
            "pli_qcov": [70, 10],
            "sucos_shape": [60, 30],
        }
    ).to_parquet(exports / "ligand_similarity_scores.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": ["query_interface", "query_interface"],
            "target_system": ["target_interface", "other_interface"],
            "iface1_qcov": [0.9, 0.4],
            "iface2_qcov": [0.8, 0.3],
            "similarity": [80, 30],
        }
    ).to_parquet(exports / "interface_similarity_scores.parquet", index=False)
    protein_scores = exports / "protein_similarity_scores/alignment_type=foldseek"
    protein_scores.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "target_entry": ["2def"],
            "query_chain_mapped": ["A"],
            "target_chain_mapped": ["B"],
            "source": ["foldseek"],
            "qcov": pd.Series([90], dtype="uint8"),
            "tcov": pd.Series([85], dtype="uint8"),
            "fident": pd.Series([70], dtype="uint8"),
            "seqsim": pd.Series([75], dtype="uint8"),
            "lddt": pd.Series([80], dtype="uint8"),
        }
    ).to_parquet(protein_scores / "shard=ab.parquet", index=False)
    alignment_shard = (
        tmp_path / "alignments/search_db=holo/alignment_type=mmseqs/shard=ab.parquet"
    )
    alignment_shard.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "target_entry": ["2def"],
            "query_chain_mapped": ["A"],
            "target_chain_mapped": ["B"],
            "query_selected_residue_positions": [[1, 2, 3]],
            "target_selected_residue_positions": [[2, 0, 4]],
            "selected_residue_identity_bits": [b"\x05"],
        }
    ).to_parquet(alignment_shard, index=False)
    ligand_scores = tmp_path / "ligand_scores"
    ligand_scores.mkdir()
    pd.DataFrame(
        {
            "query_ligand_id": [29, 51],
            "target_ligand_id": [49918, 36689],
            "tanimoto_similarity_ecfp4_1024": [95.0, 75.0],
        }
    ).to_parquet(ligand_scores / "scores.parquet", index=False)
    return PlinderRelease(tmp_path)


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


def test_query_ligand_similarity_reads_complete_export(similarity_release):
    result = scores.query_ligand_similarity(
        columns=["target_system", "pocket_qcov", "pli_qcov", "sucos_shape"],
        filters=[("query_system", "==", "query"), ("pocket_qcov", ">=", 80)],
        release=similarity_release,
    )

    assert result.to_dict("records") == [
        {
            "target_system": "target",
            "pocket_qcov": 90,
            "pli_qcov": 70,
            "sucos_shape": 60,
        }
    ]


def test_query_protein_similarity_reads_integer_scores(similarity_release):
    result = scores.query_protein_similarity(
        columns=["target_entry", "fident", "seqsim", "lddt"],
        filters=[("query_entry", "==", "1abc")],
        release=similarity_release,
    )

    assert result.to_dict("records") == [
        {"target_entry": "2def", "fident": 70, "seqsim": 75, "lddt": 80}
    ]


@pytest.mark.parametrize(
    ("kind", "expected"),
    [
        (
            "pocket",
            {
                "overlapping_residues": 1,
                "identical_overlapping_residues": 1,
                "query_total_residues": 2,
                "target_total_residues": 2,
                "query_overlap_fraction": pytest.approx(1 / 2),
                "target_overlap_fraction": pytest.approx(1 / 2),
            },
        ),
        (
            "interface",
            {
                "overlapping_residues": 1,
                "identical_overlapping_residues": 1,
                "query_total_residues": 2,
                "target_total_residues": 1,
                "query_overlap_fraction": pytest.approx(1 / 2),
                "target_overlap_fraction": pytest.approx(1.0),
            },
        ),
    ],
)
def test_query_chain_overlap_counts_selected_residues(
    similarity_release, kind, expected
):
    result = scores.query_chain_overlap(
        "1abc",
        "A",
        "2def",
        "B",
        kind=kind,
        source="mmseqs",
        release=similarity_release,
    )

    assert result.to_dict("records") == [
        {
            "query_entry": "1abc",
            "query_chain": "A",
            "target_entry": "2def",
            "target_chain": "B",
            "source": "mmseqs",
            "kind": kind,
            **expected,
        }
    ]


def test_query_interface_similarity_reads_complete_export(similarity_release):
    result = scores.query_interface_similarity(
        columns=["target_system", "iface1_qcov", "iface2_qcov", "similarity"],
        filters=[
            ("query_system", "==", "query_interface"),
            ("similarity", ">=", 70),
        ],
        release=similarity_release,
    )

    assert result.to_dict("records") == [
        {
            "target_system": "target_interface",
            "iface1_qcov": pytest.approx(0.9),
            "iface2_qcov": pytest.approx(0.8),
            "similarity": 80,
        }
    ]


@pytest.mark.parametrize(
    "query",
    [
        scores.query_protein_similarity,
        scores.query_ligand_similarity,
        scores.query_interface_similarity,
    ],
)
def test_complete_similarity_queries_require_filters(query, similarity_release):
    with pytest.raises(ValueError, match="at least one score filter"):
        query(filters=[], release=similarity_release)


def test_query_ligand_chemical_similarity(similarity_release):
    result = scores.query_ligand_chemical_similarity(
        filters=[("query_ligand_id", "==", 29)],
        release=similarity_release,
    )

    assert result["tanimoto_similarity_ecfp4_1024"].tolist() == [95.0]


def test_cross_ligand_chemical_similarity(similarity_release, monkeypatch):
    monkeypatch.setattr(
        ligand_module,
        "_map_cross_chemical_similarity",
        lambda frame, _target_ligands, _metric, _release: frame,
    )
    result = scores.cross_ligand_chemical_similarity(
        query_ligands=["29", "51"],
        target_ligands=["49918", "36689"],
        release=similarity_release,
    )

    assert len(result.index) == 2


def test_cross_ligand_chemical_similarity_rejects_non_numeric_ids(
    similarity_release,
):
    with pytest.raises(ValueError, match="ligand IDs must be integers"):
        scores.cross_ligand_chemical_similarity(
            query_ligands=["not-an-id"],
            target_ligands=["49918"],
            release=similarity_release,
        )


@pytest.mark.parametrize(
    ("query_ligands", "target_ligands"),
    [("29", ["49918"]), (["29"], "49918")],
)
def test_cross_ligand_chemical_similarity_rejects_scalar_strings(
    query_ligands,
    target_ligands,
):
    with pytest.raises(TypeError, match="must be provided as a collection"):
        scores.cross_ligand_chemical_similarity(
            query_ligands=query_ligands,
            target_ligands=target_ligands,
        )


def test_cross_ligand_chemical_similarity_maps_nodes_through_annotation(
    monkeypatch,
):
    calls = []

    def fake_query_table(table_name, *, columns, filters, release):
        calls.append((table_name, columns, filters, release))
        return pd.DataFrame(
            {
                "system_id": ["1aaa__1__1.A__1.X", "2bbb__1__1.B__1.Y"],
                "ligand_smiles_id": [0, 0],
            }
        )

    release = PlinderRelease()
    monkeypatch.setattr(ligand_module, "query_table", fake_query_table)
    result = ligand_module._map_cross_chemical_similarity(
        pd.DataFrame(
            {
                "query_ligand_id": [0],
                "target_ligand_id": [1],
                "tanimoto_similarity_ecfp4_1024": [95.0],
            }
        ),
        target_ligands={1},
        metric="tanimoto_similarity_ecfp4_1024",
        release=release,
    )

    assert set(result["system_id"]) == {
        "1aaa__1__1.A__1.X",
        "2bbb__1__1.B__1.Y",
    }
    assert calls[0] == (
        "annotation",
        ["system_id", "ligand_smiles_id"],
        [("ligand_smiles_id", "in", {0})],
        release,
    )


def test_cross_ligand_chemical_similarity_returns_empty_without_annotation(
    monkeypatch,
):
    monkeypatch.setattr(
        ligand_module,
        "query_table",
        lambda *_args, **_kwargs: pytest.fail(
            "empty similarities must not query annotation"
        ),
    )

    result = ligand_module._map_cross_chemical_similarity(
        pd.DataFrame(
            columns=[
                "query_ligand_id",
                "target_ligand_id",
                "tanimoto_similarity_ecfp4_1024",
            ]
        ),
        target_ligands={1},
        metric="tanimoto_similarity_ecfp4_1024",
        release=None,
    )

    assert result.empty
    assert result.columns.tolist() == ["system_id", "tanimoto_similarity_ecfp4_1024"]
