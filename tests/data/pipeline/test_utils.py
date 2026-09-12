# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from types import SimpleNamespace

import pandas as pd
import pytest

from plinder.core.scores.metrics import CHEMICAL_CLUSTER_SUMMARY_COLUMNS
from plinder.data.pipeline import utils


@pytest.mark.parametrize(
    "funcname, run, skip, expect",
    [
        ("a", [], [], True),
        ("a", ["a"], [], True),
        ("a", [], ["b"], True),
        ("a", [], ["a"], False),
        ("scatter_a", ["a", "b"], ["a"], False),
        ("join_a", [], [], True),
    ],
)
def test_should_run_stage(funcname, run, skip, expect):
    assert utils.should_run_stage(funcname, run, skip) == expect


def test_skipped_scatter_keeps_foreach_join_reachable() -> None:
    class Pipeline:
        cfg = SimpleNamespace(
            flow=SimpleNamespace(run_specific_stages=["other"], skip_specific_stages=[])
        )

        @utils.ingest_flow_control
        def scatter_example(self):
            raise AssertionError("skipped scatter must not execute")

    assert Pipeline().scatter_example() == [[]]


@pytest.mark.parametrize(
    "contents",
    [
        ["aaaa", "bbbb"],
        ["pdb_0000cccc", "pdb_0000dddd"],
    ],
)
def test_hash_contents(contents):
    utils.hash_contents(contents)


def test_comparability_merge_reuses_complete_collated_column(tmp_path):
    index = pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X", "1aaa__1__1.Y", None],
            "system_type": ["holo", "holo", "apo"],
            "ligand_is_proper": [True, False, False],
            "ligand_is_shape_comparable": [True, True, None],
        }
    )

    result = utils.add_ligand_3d_score_ability_column(
        index=index,
        data_dir=tmp_path,
    )

    assert result["ligand_is_shape_comparable"].tolist() == [True, False, False]
    assert str(result["ligand_is_shape_comparable"].dtype) == "boolean"


def test_finalize_index_writes_local_clusters_to_sidecar(tmp_path):
    index_dir = tmp_path / "index"
    directed_cover_file = (
        tmp_path
        / "ligand_sampling/directed_set_cover/metric=pli_qcov"
        / "threshold=100.parquet"
    )
    tanimoto_cover_file = (
        tmp_path
        / "ligand_sampling/set_cover/metric=tanimoto_similarity_ecfp4_1024"
        / "threshold=90.parquet"
    )
    index_dir.mkdir(parents=True)
    directed_cover_file.parent.mkdir(parents=True)
    tanimoto_cover_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1aaa", "1aaa"],
            "system_id": ["1aaa__1__1.A__1.X", "1aaa__2__1.A__1.X"],
            "system_id_no_biounit": ["1aaa__1.A__1.X", "1aaa__1.A__1.X"],
            "system_biounit_id": ["1", "2"],
            "system_type": ["holo", "holo"],
            "ligand_id": ["1aaa__1__1.X", "1aaa__2__1.X"],
            "ligand_is_proper": [True, False],
            "ligand_smiles": ["CCO", "CCO"],
            "ligand_tanimoto_ecfp4_1024_90_cluster": ["legacy", "legacy"],
            "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids": [1, 1],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": ["CCO"],
            "ligand_smiles_id": [0],
            "ligand_is_cofactor_like": [False],
        }
    ).to_parquet(fingerprint_dir / "ligand_similarity_annotations.parquet", index=False)
    ligand_dir = tmp_path / "ligands"
    ligand_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X"],
            "ligand_is_shape_comparable": [True],
        }
    ).to_parquet(ligand_dir / "part.parquet", index=False)
    pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X"],
            "centroid_ligand_id": ["1aaa__1__1.X"],
            "similarity_to_centroid": [100.0],
            "coverage_count": [1],
            "coverage_fraction": [1.0],
            "label": ["d0"],
            "metric": ["pli_qcov"],
            "threshold": [100],
            "directed": [True],
        }
    ).to_parquet(directed_cover_file, index=False)
    pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X"],
            "centroid_ligand_id": ["1aaa__1__1.X"],
            "label": ["t0"],
            "metric": ["tanimoto_similarity_ecfp4_1024"],
            "cluster": ["set_cover"],
            "threshold": [90],
            "directed": [False],
        }
    ).to_parquet(tanimoto_cover_file, index=False)
    utils.finalize_index(data_dir=tmp_path)
    finalized = pd.read_parquet(index_dir / "annotation_table.parquet")
    clusters = pd.read_parquet(index_dir / "ligand_clusters.parquet")
    directed_labels = clusters["pli_qcov__100__ligand__directed_set_cover"]
    assert directed_labels.iloc[0] == "d0"
    assert pd.isna(directed_labels.iloc[1])
    centroid_flags = clusters["pli_qcov__100__ligand__directed_set_cover__is_centroid"]
    assert bool(centroid_flags.iloc[0])
    assert pd.isna(centroid_flags.iloc[1])
    coverage_counts = clusters[
        "pli_qcov__100__ligand__directed_set_cover__coverage_count"
    ]
    coverage_fractions = clusters[
        "pli_qcov__100__ligand__directed_set_cover__coverage_fraction"
    ]
    assert coverage_counts.iloc[0] == 1
    assert coverage_fractions.iloc[0] == pytest.approx(1.0)
    assert pd.isna(coverage_counts.iloc[1])
    assert pd.isna(coverage_fractions.iloc[1])
    tanimoto_column = "tanimoto_similarity_ecfp4_1024__90__ligand__set_cover"
    assert clusters[tanimoto_column].tolist()[0] == "t0"
    assert pd.isna(clusters[tanimoto_column].iloc[1])
    assert bool(clusters[f"{tanimoto_column}__is_centroid"].iloc[0])
    assert clusters["ligand_tanimoto_ecfp4_1024_90_cluster"].iloc[0] == "t0"
    assert pd.isna(clusters["ligand_tanimoto_ecfp4_1024_90_cluster"].iloc[1])
    assert clusters["ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids"].iloc[0] == 1
    assert finalized.loc[0, "ligand_smiles_id"] == 0
    assert pd.isna(finalized.loc[1, "ligand_smiles_id"])
    assert finalized["ligand_is_shape_comparable"].tolist() == [True, False]
    assert "uniqueness" not in finalized.columns
    assert not any("set_cover" in column for column in finalized)
    assert "ligand_tanimoto_ecfp4_1024_90_cluster" not in finalized


def test_finalize_index_preserves_annotation_when_staging_fails(tmp_path, monkeypatch):
    index_dir = tmp_path / "index"
    index_dir.mkdir(parents=True)
    index_path = index_dir / "annotation_table.parquet"
    original = pd.DataFrame(
        {
            "entry_pdb_id": ["1aaa"],
            "system_id": ["1aaa__1__1.A__1.X"],
            "system_type": ["holo"],
            "ligand_id": ["1aaa__1__1.X"],
            "ligand_is_proper": [False],
            "ligand_smiles": [None],
            "ligand_is_shape_comparable": [False],
        }
    )
    original.to_parquet(index_path, index=False)
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": pd.Series(dtype="string"),
            "ligand_smiles_id": pd.Series(dtype="Int32"),
            "ligand_is_cofactor_like": pd.Series(dtype="boolean"),
        }
    ).to_parquet(
        fingerprint_dir / "ligand_similarity_annotations.parquet",
        index=False,
    )

    original_to_parquet = pd.DataFrame.to_parquet

    def fail_cluster_staging(frame, path, *args, **kwargs):
        if path.name == "ligand_clusters.tmp.parquet":
            raise OSError("staging failed")
        return original_to_parquet(frame, path, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_parquet", fail_cluster_staging)

    with pytest.raises(OSError, match="staging failed"):
        utils.finalize_index(data_dir=tmp_path)

    pd.testing.assert_frame_equal(pd.read_parquet(index_path), original)
    assert not (index_dir / "annotation_table.tmp.parquet").exists()


def test_cluster_index_rejects_non_tanimoto_set_cover(tmp_path):
    cover_file = (
        tmp_path / "ligand_sampling/set_cover/metric=pli_qcov" / "threshold=100.parquet"
    )
    cover_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X"],
            "centroid_ligand_id": ["1aaa__1__1.X"],
            "label": ["c0"],
        }
    ).to_parquet(cover_file, index=False)
    index = pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X"],
            "system_type": ["holo"],
            "ligand_is_proper": [True],
            "ligand_smiles_id": [0],
        }
    )

    with pytest.raises(ValueError, match="invalid ligand set-cover modes"):
        utils.build_ligand_cluster_table(index=index, data_dir=tmp_path)


@pytest.mark.parametrize(
    "metric, column", sorted(CHEMICAL_CLUSTER_SUMMARY_COLUMNS.items())
)
def test_chemical_90_set_cover_counts_distinct_pdb_ids(tmp_path, metric, column):
    cover_file = (
        tmp_path / f"ligand_sampling/set_cover/metric={metric}" / "threshold=90.parquet"
    )
    cover_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "ligand_id": ["l1", "l2", "l3"],
            "centroid_ligand_id": ["l1", "l1", "l3"],
            "label": ["c0", "c0", "c1"],
        }
    ).to_parquet(cover_file, index=False)
    index = pd.DataFrame(
        {
            "entry_pdb_id": ["1aaa", "2bbb", "1aaa"],
            "ligand_id": ["l1", "l2", "l3"],
            "system_type": ["holo", "holo", "holo"],
            "ligand_is_proper": [True, True, True],
            "ligand_smiles_id": [0, 1, 2],
        }
    )

    result = utils.build_ligand_cluster_table(index=index, data_dir=tmp_path)

    assert result[column].tolist() == [
        "c0",
        "c0",
        "c1",
    ]
    assert result[f"{column}_num_pdb_ids"].tolist() == [
        2,
        2,
        1,
    ]


def test_cluster_index_marks_only_directed_cover_centroids(tmp_path):
    cover_file = (
        tmp_path
        / "ligand_sampling/directed_set_cover/metric=pli_qcov"
        / "threshold=50.parquet"
    )
    cover_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "ligand_id": ["l1", "l2"],
            "label": ["d0", "d0"],
            "centroid_ligand_id": ["l2", "l2"],
            "coverage_count": [1, 2],
            "coverage_fraction": [0.5, 1.0],
        }
    ).to_parquet(cover_file, index=False)
    index = pd.DataFrame(
        {
            "ligand_id": ["l1", "l2", "not-eligible"],
            "system_type": ["holo", "holo", "apo"],
            "ligand_is_proper": [True, True, False],
            "ligand_smiles_id": [0, 1, pd.NA],
        }
    )

    result = utils.build_ligand_cluster_table(index=index, data_dir=tmp_path)

    label_column = "pli_qcov__50__ligand__directed_set_cover"
    centroid_column = f"{label_column}__is_centroid"
    assert result[label_column].tolist()[:2] == ["d0", "d0"]
    assert result[centroid_column].tolist()[:2] == [False, True]
    assert pd.isna(result.loc[2, centroid_column])
    assert str(result[centroid_column].dtype) == "boolean"
    count_column = f"{label_column}__coverage_count"
    fraction_column = f"{label_column}__coverage_fraction"
    assert result[count_column].tolist()[:2] == [1, 2]
    assert result[fraction_column].tolist()[:2] == [0.5, 1.0]
    assert pd.isna(result.loc[2, count_column])
    assert pd.isna(result.loc[2, fraction_column])
    assert str(result[count_column].dtype) == "Int32"
    assert str(result[fraction_column].dtype) == "Float32"


def test_cluster_index_reads_legacy_cover_during_centrality_migration(tmp_path):
    cover_file = (
        tmp_path
        / "ligand_sampling/directed_set_cover/metric=pli_qcov"
        / "threshold=50.parquet"
    )
    cover_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "ligand_id": ["l1", "l2"],
            "label": ["d0", "d0"],
            "centroid_ligand_id": ["l2", "l2"],
        }
    ).to_parquet(cover_file, index=False)
    index = pd.DataFrame(
        {
            "ligand_id": ["l1", "l2"],
            "system_type": ["holo", "holo"],
            "ligand_is_proper": [True, True],
            "ligand_smiles_id": [0, 1],
        }
    )

    result = utils.build_ligand_cluster_table(index=index, data_dir=tmp_path)

    label_column = "pli_qcov__50__ligand__directed_set_cover"
    assert result[f"{label_column}__is_centroid"].tolist() == [False, True]
    assert f"{label_column}__coverage_count" not in result
    assert f"{label_column}__coverage_fraction" not in result


def test_ligand_similarity_rejects_stale_proper_smiles_universe(tmp_path):
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": ["CC"],
            "ligand_smiles_id": [0],
        }
    ).to_parquet(
        fingerprint_dir / "ligand_similarity_annotations.parquet",
        index=False,
    )
    index = pd.DataFrame(
        {
            "system_type": ["holo", "holo"],
            "ligand_is_proper": [True, True],
            "ligand_smiles": ["CC", "CCC"],
        }
    )

    with pytest.raises(ValueError, match="proper holo SMILES universe"):
        utils.add_ligand_similarity_columns(index=index, data_dir=tmp_path)


def test_ligand_similarity_rejects_artifact_from_before_targeted_repair(
    tmp_path,
):
    fingerprint_dir = tmp_path / "fingerprints"
    index_dir = tmp_path / "index"
    fingerprint_dir.mkdir()
    index_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_smiles": ["CC"],
            "ligand_smiles_id": [0],
        }
    ).to_parquet(
        fingerprint_dir / "ligand_similarity_annotations.parquet",
        index=False,
    )
    (index_dir / "collation.json").write_text(
        '{"status": "requires_downstream_repair"}'
    )
    index = pd.DataFrame(
        {
            "system_type": ["holo"],
            "ligand_is_proper": [True],
            "ligand_rdkit_canonical_smiles": ["CC"],
        }
    )

    with pytest.raises(ValueError, match="predate"):
        utils.add_ligand_similarity_columns(index=index, data_dir=tmp_path)
