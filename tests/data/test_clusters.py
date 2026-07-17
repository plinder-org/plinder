import pandas as pd
import pytest


def test_make_components_and_communities(write_plinder_mount):
    from plinder.data.clusters import make_components_and_communities

    i = len(list((write_plinder_mount / "clusters").rglob("*")))
    make_components_and_communities(
        data_dir=write_plinder_mount,
        metric="pocket_lddt",
        threshold=50,
    )
    j = len(list((write_plinder_mount / "clusters").rglob("*")))
    assert i < j


def test_ligand_clusters_and_system_projection_are_merged(tmp_path):
    from plinder.data.clusters import make_components_and_communities
    from plinder.data.pipeline.utils import finalize_index

    index_dir = tmp_path / "index"
    score_dir = tmp_path / "scores" / "search_db=holo"
    index_dir.mkdir(parents=True)
    score_dir.mkdir(parents=True)

    system_a = "1aaa__1__1.A__1.X_1.Y"
    system_b = "2bbb__1__1.B__1.Z"
    system_c = "3ccc__1__1.C__1.W"
    system_d = "4ddd__1__1.D__1.V"
    ligand_a1 = "1aaa__1__1.X"
    ligand_a2 = "1aaa__1__1.Y"
    ligand_b = "2bbb__1__1.Z"
    ligand_c = "3ccc__1__1.W"
    ligand_d = "4ddd__1__1.V"
    annotation = pd.DataFrame(
        {
            "system_id": [system_a, system_a, system_b, system_c, system_d],
            "system_id_no_biounit": [
                "1aaa__1.A",
                "1aaa__1.A",
                "2bbb__1.B",
                "3ccc__1.C",
                "4ddd__1.D",
            ],
            "system_biounit_id": ["1"] * 5,
            "ligand_id": [ligand_a1, ligand_a2, ligand_b, ligand_c, ligand_d],
            "ligand_rdkit_canonical_smiles": [
                "CC",
                "CCC",
                "CCCC",
                "CCCCC",
                "CCCCCC",
            ],
            "system_type": ["holo"] * 5,
            "system_num_protein_chains": [1] * 5,
            "system_num_ligand_chains": [2, 2, 1, 1, 1],
            "ligand_is_proper": [True, True, True, True, False],
        }
    )
    annotation.to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": [system_a, system_b, system_a, system_a, system_d],
            "query_ligand_id": [
                ligand_a1,
                ligand_b,
                ligand_a2,
                ligand_a2,
                ligand_d,
            ],
            "target_system": [system_b, system_a, system_c, system_d, system_a],
            "target_ligand_id": [
                ligand_b,
                ligand_a1,
                ligand_c,
                ligand_d,
                ligand_a2,
            ],
            "metric": ["sucos_shape_pocket_qcov"] * 5,
            "similarity": [80, 80, 40, 90, 90],
        }
    ).to_parquet(score_dir / "scores.parquet", index=False)
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": annotation[
                "ligand_rdkit_canonical_smiles"
            ],
            "ligand_smiles_id": range(len(annotation)),
        }
    ).to_parquet(fingerprint_dir / "ligand_similarity_annotations.parquet", index=False)
    ligand_dir = tmp_path / "ligands"
    ligand_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_id": annotation["ligand_id"],
            "ligand_is_3d_score_able": [True, True, False, True, False],
        }
    ).to_parquet(ligand_dir / "part.parquet", index=False)

    make_components_and_communities(
        data_dir=tmp_path,
        metric="sucos_shape_pocket_qcov",
        threshold=50,
    )

    system_labels = pd.read_parquet(
        tmp_path
        / "clusters/cluster=components/directed=True/metric=sucos_shape_pocket_qcov/threshold=50.parquet"
    ).set_index("system_id")["label"]
    ligand_labels = pd.read_parquet(
        tmp_path
        / "ligand_clusters/cluster=components/directed=True/metric=sucos_shape_pocket_qcov/threshold=50.parquet"
    ).set_index("ligand_id")["label"]
    assert system_labels[system_a] == system_labels[system_b]
    assert system_labels[system_c] != system_labels[system_a]
    assert system_labels[system_d] != system_labels[system_a]
    assert ligand_labels[ligand_a1] == ligand_labels[ligand_b]
    assert ligand_labels[ligand_a2] != ligand_labels[ligand_a1]
    assert ligand_labels[ligand_c] != ligand_labels[ligand_a1]
    assert ligand_d not in ligand_labels

    finalized = finalize_index(data_dir=tmp_path)
    system_column = "sucos_shape_pocket_qcov__50__strong__component"
    ligand_column = "sucos_shape_pocket_qcov__50__ligand__strong__component"
    assert finalized.groupby("system_id")[system_column].nunique().max() == 1
    finalized_labels = finalized.set_index("ligand_id")[ligand_column]
    assert finalized_labels[ligand_a1] == finalized_labels[ligand_b]
    assert finalized_labels[ligand_a2] != finalized_labels[ligand_a1]
    assert not bool(
        finalized.set_index("ligand_id").loc[ligand_b, "ligand_is_3d_score_able"]
    )


def test_default_cluster_metrics_use_only_pocket_weighted_shape():
    from plinder.core.scores.metrics import LIGAND_SCORE_NAMES
    from plinder.data.pipeline.config import METRICS

    assert "shape" in LIGAND_SCORE_NAMES
    assert "color" in LIGAND_SCORE_NAMES
    assert "sucos_shape_pocket_qcov" in METRICS
    assert "tanimoto_similarity_ecfp4_1024" in METRICS
    assert "shape" not in METRICS
    assert "color" not in METRICS
    assert "shape_tanimoto" not in METRICS
    assert "sucos_shape" not in METRICS


def test_tanimoto_clusters_project_unique_smiles_without_occurrence_file(tmp_path):
    from plinder.data.clusters import make_components_and_communities

    index_dir = tmp_path / "index"
    fingerprint_dir = tmp_path / "fingerprints"
    score_dir = tmp_path / "ligand_scores"
    index_dir.mkdir()
    fingerprint_dir.mkdir()
    score_dir.mkdir()

    system_a = "1aaa__1__1.A__1.X_1.Y"
    system_b = "2bbb__1__1.B__1.Z"
    system_c = "3ccc__1__1.C__1.W"
    pd.DataFrame(
        {
            "system_id": [system_a, system_a, system_b, system_c],
            "ligand_molecular_weight": [100.0, 200.0, 150.0, 100.0],
            "ligand_rdkit_canonical_smiles": ["CC", "CCC", "CCCC", "CC"],
            "ligand_is_ion": [False] * 4,
            "ligand_is_artifact": [False] * 4,
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1, 2],
            "ligand_rdkit_canonical_smiles": ["CC", "CCC", "CCCC"],
        }
    ).to_parquet(fingerprint_dir / "ligands_per_smiles.parquet", index=False)
    pd.DataFrame(
        {
            "query_ligand_id": [0, 0, 1, 2, 2],
            "target_ligand_id": [0, 2, 1, 0, 2],
            "tanimoto_similarity_ecfp4_1024": [100.0, 95.0, 100.0, 95.0, 100.0],
        }
    ).to_parquet(score_dir / "part.parquet", index=False)

    make_components_and_communities(
        data_dir=tmp_path,
        metric="tanimoto_similarity_ecfp4_1024",
        threshold=90,
    )

    labels = pd.read_parquet(
        tmp_path
        / "clusters/cluster=components/directed=True/metric=tanimoto_similarity_ecfp4_1024/threshold=90.parquet"
    ).set_index("system_id")["label"]
    assert labels[system_b] == labels[system_c]
    assert labels[system_a] != labels[system_b]
    assert not (fingerprint_dir / "ligands_per_system.parquet").exists()


@pytest.mark.parametrize("metric", ["shape", "color", "sucos_shape"])
def test_gated_shape_diagnostics_cannot_be_clustered(tmp_path, metric):
    from plinder.data.clusters import make_components_and_communities

    with pytest.raises(ValueError, match="sucos_shape_pocket_qcov"):
        make_components_and_communities(
            data_dir=tmp_path,
            metric=metric,
            threshold=50,
        )


def test_empty_metric_graph_writes_system_and_ligand_singletons(tmp_path):
    from plinder.data.clusters import make_components_and_communities

    index_dir = tmp_path / "index"
    index_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["1aaa__1__1.A__1.X", "2bbb__1__1.B__1.Y"],
            "ligand_id": ["1aaa__1__1.X", "2bbb__1__1.Y"],
            "system_type": ["holo", "holo"],
            "system_num_protein_chains": [1, 1],
            "system_num_ligand_chains": [1, 1],
            "ligand_is_proper": [True, True],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)

    make_components_and_communities(
        data_dir=tmp_path,
        metric="sucos_shape_pocket_qcov",
        threshold=100,
    )

    system_labels = pd.read_parquet(
        tmp_path
        / "clusters/cluster=components/directed=True/metric=sucos_shape_pocket_qcov/threshold=100.parquet"
    )
    ligand_labels = pd.read_parquet(
        tmp_path
        / "ligand_clusters/cluster=components/directed=True/metric=sucos_shape_pocket_qcov/threshold=100.parquet"
    )
    assert system_labels["system_id"].nunique() == 2
    assert system_labels["label"].nunique() == 2
    assert ligand_labels["ligand_id"].nunique() == 2
    assert ligand_labels["label"].nunique() == 2
