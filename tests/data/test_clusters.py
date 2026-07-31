import json
from os import utime

import networkit as nk
import numpy as np
import pandas as pd
import pytest


def _component_partition(labels: pd.DataFrame) -> set[frozenset[str]]:
    return {
        frozenset(group["ligand_id"].astype(str))
        for _, group in labels.groupby("label")
    }


def test_exact_threshold_components_preserve_bridge_edges_across_shards():
    from plinder.data.clusters import (
        count_crossing_component_edges,
        make_exact_threshold_components,
    )

    # A is a deliberately bad representative for the 100% component {A, B}:
    # A-C is below 70, but B-C is a qualifying bridge that must not be lost.
    shards = [
        pd.DataFrame(
            {
                "query_node": ["A", "A", "D"],
                "target_node": ["B", "C", "E"],
                "similarity": [100, 60, 95],
            }
        ),
        pd.DataFrame(
            {
                "query_node": ["B", "C", "E"],
                "target_node": ["C", "D", "F"],
                "similarity": [70, 50, 30],
            }
        ),
    ]
    labels = make_exact_threshold_components(
        edge_batches=shards,
        all_nodes=["A", "B", "C", "D", "E", "F", "G"],
        thresholds=[30, 50, 70, 90, 100],
        forest_fan_in=2,
    )

    assert _component_partition(labels[100]) == {
        frozenset({"A", "B"}),
        frozenset({"C"}),
        frozenset({"D"}),
        frozenset({"E"}),
        frozenset({"F"}),
        frozenset({"G"}),
    }
    assert frozenset({"A", "B", "C"}) in _component_partition(labels[70])
    assert frozenset({"A", "B", "C", "D", "E"}) in _component_partition(labels[50])
    assert _component_partition(labels[30]) == {
        frozenset({"A", "B", "C", "D", "E", "F"}),
        frozenset({"G"}),
    }
    assert count_crossing_component_edges(
        edge_batches=shards,
        labels_by_threshold=labels,
    ) == {100: 0, 90: 0, 70: 0, 50: 0, 30: 0}


def test_exact_threshold_components_match_full_graph_and_are_order_independent():
    from plinder.data.clusters import make_exact_threshold_components

    shards = [
        pd.DataFrame(
            {
                "query_node": ["A", "B", "C", "D"],
                "target_node": ["B", "C", "A", "E"],
                "similarity": [95, 70, 50, 30],
            }
        ),
        pd.DataFrame(
            {
                "query_node": ["E", "F", "G", "A"],
                "target_node": ["F", "G", "D", "A"],
                "similarity": [100, 50, 29, 100],
            }
        ),
    ]
    nodes = list("ABCDEFGH")
    thresholds = [100, 90, 70, 50, 30]
    forward = make_exact_threshold_components(
        edge_batches=shards,
        all_nodes=nodes,
        thresholds=thresholds,
        forest_fan_in=2,
    )
    reverse = make_exact_threshold_components(
        edge_batches=reversed(shards),
        all_nodes=reversed(nodes),
        thresholds=reversed(thresholds),
        forest_fan_in=2,
    )

    for threshold in thresholds:
        qualifying = pd.concat(shards, ignore_index=True)
        qualifying = qualifying[qualifying["similarity"] >= threshold]
        expected = {frozenset({node}) for node in nodes}
        if not qualifying.empty:
            # A tiny independent union implementation is clearer than relying on
            # the production forest reducer to construct the expected partition.
            groups = {node: {node} for node in nodes}
            for edge in qualifying.itertuples(index=False):
                merged = groups[str(edge.query_node)] | groups[str(edge.target_node)]
                for node in merged:
                    groups[node] = merged
            expected = {frozenset(group) for group in groups.values()}
        assert _component_partition(forward[threshold]) == expected
        assert forward[threshold].equals(reverse[threshold])


def test_crossing_component_edge_validation_detects_split():
    from plinder.data.clusters import count_crossing_component_edges

    edges = pd.DataFrame(
        {
            "query_node": ["A", "B"],
            "target_node": ["B", "C"],
            "similarity": [80, 40],
        }
    )
    labels = {
        50: pd.DataFrame({"ligand_id": ["A", "B", "C"], "label": ["c0", "c1", "c2"]}),
        30: pd.DataFrame({"ligand_id": ["A", "B", "C"], "label": ["c0", "c0", "c1"]}),
    }

    assert count_crossing_component_edges(
        edge_batches=[edges],
        labels_by_threshold=labels,
    ) == {50: 1, 30: 1}


def test_exact_strong_components_preserve_cross_shard_directional_cycle():
    from plinder.data.clusters import make_exact_threshold_components

    shards = [
        pd.DataFrame(
            {
                "query_node": ["A", "D"],
                "target_node": ["B", "A"],
                "similarity": [90, 90],
            }
        ),
        pd.DataFrame(
            {
                "query_node": ["B"],
                "target_node": ["C"],
                "similarity": [70],
            }
        ),
        pd.DataFrame(
            {
                "query_node": ["C"],
                "target_node": ["A"],
                "similarity": [50],
            }
        ),
    ]
    labels = make_exact_threshold_components(
        edge_batches=shards,
        all_nodes=["A", "B", "C", "D"],
        thresholds=[50, 70, 90],
        forest_fan_in=2,
        directed=True,
    )

    assert _component_partition(labels[90]) == {
        frozenset({"A"}),
        frozenset({"B"}),
        frozenset({"C"}),
        frozenset({"D"}),
    }
    assert _component_partition(labels[70]) == _component_partition(labels[90])
    assert _component_partition(labels[50]) == {
        frozenset({"A", "B", "C"}),
        frozenset({"D"}),
    }


def test_weak_and_strong_components_share_one_input_scan():
    from plinder.data.clusters import (
        make_exact_weak_and_strong_threshold_components,
    )

    iterations = 0

    def batches():
        nonlocal iterations
        iterations += 1
        if iterations > 1:
            raise AssertionError("edge batches were scanned more than once")
        yield pd.DataFrame(
            {
                "query_node": ["A", "B"],
                "target_node": ["B", "A"],
                "similarity": [70, 50],
            }
        )

    labels = make_exact_weak_and_strong_threshold_components(
        edge_batches=batches(),
        all_nodes=["A", "B", "C"],
        thresholds=[50, 70],
        forest_fan_in=2,
    )

    assert iterations == 1
    assert frozenset({"A", "B"}) in _component_partition(labels[False][70])
    assert _component_partition(labels[True][70]) == {
        frozenset({"A"}),
        frozenset({"B"}),
        frozenset({"C"}),
    }
    assert frozenset({"A", "B"}) in _component_partition(labels[True][50])


def test_component_reduction_shards_are_resumable_and_merge_exactly(tmp_path):
    from plinder.data.clusters import (
        make_exact_threshold_components,
        merge_component_reduction_shards,
        write_component_reduction_shard,
    )

    nodes = list("ABCDE")
    thresholds = [50, 70, 100]
    shards = [
        pd.DataFrame(
            {
                "query_node": ["A", "B", "C"],
                "target_node": ["B", "C", "A"],
                "similarity": [100, 70, 50],
            }
        ),
        pd.DataFrame(
            {
                "query_node": ["D", "E", "C"],
                "target_node": ["E", "D", "D"],
                "similarity": [70, 50, 50],
            }
        ),
    ]
    reduction_dirs = []
    for shard_index, shard in enumerate(shards):
        source = tmp_path / f"source-{shard_index}.parquet"
        shard.to_parquet(source, index=False)
        reduction_dir = tmp_path / "reductions" / str(shard_index)
        write_component_reduction_shard(
            edge_batches=[shard],
            all_nodes=nodes,
            thresholds=thresholds,
            source_path=source,
            output_dir=reduction_dir,
            forest_fan_in=2,
        )
        reduction_dirs.append(reduction_dir)

    merged = merge_component_reduction_shards(
        reduction_dirs=reduction_dirs,
        all_nodes=nodes,
        thresholds=thresholds,
        forest_fan_in=2,
        parquet_batch_size=2,
    )
    direct = make_exact_threshold_components(
        edge_batches=shards,
        all_nodes=nodes,
        thresholds=thresholds,
        forest_fan_in=2,
    )
    for threshold in thresholds:
        assert merged[False][threshold].equals(direct[threshold])

    def should_not_be_read():
        raise AssertionError("completed reduction unexpectedly rescanned its source")
        yield pd.DataFrame()

    cached = write_component_reduction_shard(
        edge_batches=should_not_be_read(),
        all_nodes=nodes,
        thresholds=thresholds,
        source_path=tmp_path / "source-0.parquet",
        output_dir=reduction_dirs[0],
        forest_fan_in=2,
    )
    assert cached["status"] == "complete"


def test_component_score_batches_read_reciprocal_minimum_edges(tmp_path):
    from plinder.data.clusters import iter_component_score_batches

    source = tmp_path / "symmetric.parquet"
    pd.DataFrame(
        {
            "query_node": ["l1", "l2"],
            "target_node": ["l2", "l3"],
            "similarity": [70.0, 95.0],
        }
    ).to_parquet(source, index=False)

    result = pd.concat(
        iter_component_score_batches(
            source_path=source,
            metric="sucos_shape_pocket_qcov",
            batch_size=1,
        ),
        ignore_index=True,
    )

    assert result.to_dict("records") == [
        {"query_node": "l1", "target_node": "l2", "similarity": 70.0},
        {"query_node": "l2", "target_node": "l3", "similarity": 95.0},
    ]


def test_symmetric_edge_shards_take_minimum_of_directional_maxima(tmp_path):
    from plinder.data.clusters import (
        make_directed_cover_component_reduction,
        make_directed_set_cover,
        merge_directed_cover_component_reductions,
        prepare_symmetric_edge_plan,
        write_symmetric_edge_fragment_batch,
        write_symmetric_edge_shard,
    )

    metric = "sucos_shape_pocket_qcov"
    index_dir = tmp_path / "index"
    score_dir = tmp_path / "scores" / "search_db=holo"
    index_dir.mkdir(parents=True)
    score_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["s1", "s2", "s3", "s4"],
            "ligand_id": ["l1", "l2", "l3", "l4"],
            "system_type": ["holo"] * 4,
            "ligand_is_proper": [True] * 4,
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": ["s1", "s1", "s2", "s1", "s4"],
            "query_ligand_id": ["l1", "l1", "l2", "l1", "l4"],
            "target_system": ["s2", "s2", "s1", "s3", "s3"],
            "target_ligand_id": ["l2", "l2", "l1", "l3", "l3"],
            "metric": [metric] * 5,
            "similarity": [80, 75, 60, 90, 85],
        }
    ).to_parquet(score_dir / "part.parquet", index=False)
    plan = prepare_symmetric_edge_plan(
        data_dir=tmp_path,
        metrics=[metric],
        source_batch_size=1,
        bucket_count=1,
    )
    write_symmetric_edge_fragment_batch(
        data_dir=tmp_path,
        batch=plan["batches"][0],
        scratch_dir=tmp_path / "scratch-fragments",
        threads=1,
    )
    result = write_symmetric_edge_shard(
        data_dir=tmp_path,
        metric=metric,
        bucket=0,
        scratch_dir=tmp_path / "scratch-shard",
        threads=1,
    )

    assert result["rows"] == 3
    edges = pd.read_parquet(
        tmp_path
        / "ligand_clusters/symmetric_edges"
        / f"metric={metric}/bucket=000.parquet"
    ).set_index(["query_node", "target_node"])
    assert edges.loc[("l1", "l2"), "forward_similarity"] == pytest.approx(80.0)
    assert edges.loc[("l1", "l2"), "reverse_similarity"] == pytest.approx(60.0)
    assert edges.loc[("l1", "l2"), "similarity"] == pytest.approx(60.0)
    assert edges.loc[("l1", "l2"), "maximum_similarity"] == pytest.approx(80.0)
    assert edges.loc[("l1", "l3"), "forward_similarity"] == pytest.approx(90.0)
    assert pd.isna(edges.loc[("l1", "l3"), "reverse_similarity"])
    assert pd.isna(edges.loc[("l1", "l3"), "similarity"])
    assert edges.loc[("l1", "l3"), "maximum_similarity"] == pytest.approx(90.0)
    assert edges.loc[("l3", "l4"), "reverse_similarity"] == pytest.approx(85.0)
    assert pd.isna(edges.loc[("l3", "l4"), "forward_similarity"])
    assert pd.isna(edges.loc[("l3", "l4"), "similarity"])

    edge_path = (
        tmp_path
        / "ligand_clusters/symmetric_edges"
        / f"metric={metric}/bucket=000.parquet"
    )
    make_directed_cover_component_reduction(
        data_dir=tmp_path,
        metric=metric,
        thresholds=[50],
        source_path=edge_path,
    )
    merge_directed_cover_component_reductions(
        data_dir=tmp_path,
        metric=metric,
        thresholds=[50],
    )
    cover_path = make_directed_set_cover(
        data_dir=tmp_path,
        metric=metric,
        threshold=50,
        scratch_dir=tmp_path / "scratch-cover",
        threads=1,
    )
    cover = pd.read_parquet(cover_path).set_index("ligand_id")
    assert cover.loc["l1", "centroid_ligand_id"] == "l3"
    assert cover.loc["l1", "similarity_to_centroid"] == pytest.approx(90.0)
    assert cover.loc["l2", "centroid_ligand_id"] == "l2"
    assert cover.loc["l3", "centroid_ligand_id"] == "l3"
    assert cover.loc["l4", "centroid_ligand_id"] == "l3"
    assert cover.loc["l4", "similarity_to_centroid"] == pytest.approx(85.0)


def test_symmetric_edge_fragment_validation_can_be_scoped_to_one_bucket(tmp_path):
    from plinder.data.clusters import _completed_symmetric_fragment_batch

    batch = {"key": "batch-key", "kind": "score"}
    plan = {"plan_hash": "plan-hash"}
    output_dir = (
        tmp_path
        / "ligand_clusters"
        / "symmetric_edges"
        / "fragments"
        / "batch=batch-key"
    )
    selected = output_dir / "metric=pocket_qcov/bucket=0/data.parquet"
    unrelated = output_dir / "metric=pli_qcov/bucket=1/data.parquet"
    selected.parent.mkdir(parents=True)
    unrelated.parent.mkdir(parents=True)
    selected.write_bytes(b"selected")
    unrelated.write_bytes(b"unrelated")
    manifest = {
        "version": 1,
        "plan_hash": plan["plan_hash"],
        "batch": batch,
        "outputs": [
            {
                "path": selected.relative_to(output_dir).as_posix(),
                "size": selected.stat().st_size,
            },
            {
                "path": unrelated.relative_to(output_dir).as_posix(),
                "size": unrelated.stat().st_size,
            },
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest))
    unrelated.unlink()

    assert (
        _completed_symmetric_fragment_batch(
            data_dir=tmp_path,
            plan=plan,
            batch=batch,
            output_prefix="metric=pocket_qcov/bucket=0/",
        )
        == manifest
    )
    assert (
        _completed_symmetric_fragment_batch(
            data_dir=tmp_path,
            plan=plan,
            batch=batch,
        )
        is None
    )


def test_communities_use_greedy_cover_and_best_centroid_assignment(tmp_path):
    from plinder.data.clusters import (
        make_communities,
        prepare_symmetric_edge_plan,
        write_symmetric_edge_fragment_batch,
        write_symmetric_edge_shard,
    )

    metric = "sucos_shape_pocket_qcov"
    ligands = ["l1", "l2", "l3", "l4", "l5"]
    systems = [f"s{index}" for index in range(len(ligands))]
    index_dir = tmp_path / "index"
    score_dir = tmp_path / "scores" / "search_db=holo"
    component_dir = (
        tmp_path
        / "ligand_clusters"
        / "reductions"
        / f"metric={metric}"
        / "labels"
        / "directed=false"
    )
    index_dir.mkdir(parents=True)
    score_dir.mkdir(parents=True)
    component_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": systems,
            "ligand_id": ligands,
            "system_type": ["holo"] * len(ligands),
            "system_num_protein_chains": [1] * len(ligands),
            "system_num_ligand_chains": [1] * len(ligands),
            "ligand_is_proper": [True] * len(ligands),
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": [
                systems[0],
                systems[1],
                systems[1],
                systems[2],
                systems[2],
                systems[3],
            ],
            "query_ligand_id": ["l1", "l2", "l2", "l3", "l3", "l4"],
            "target_system": [
                systems[1],
                systems[0],
                systems[2],
                systems[1],
                systems[3],
                systems[2],
            ],
            "target_ligand_id": ["l2", "l1", "l3", "l2", "l4", "l3"],
            "metric": [metric] * 6,
            "similarity": [80, 70, 90, 60, 85, 75],
        }
    ).to_parquet(score_dir / "part.parquet", index=False)
    plan = prepare_symmetric_edge_plan(
        data_dir=tmp_path,
        metrics=[metric],
        source_batch_size=1,
        bucket_count=1,
    )
    write_symmetric_edge_fragment_batch(
        data_dir=tmp_path,
        batch=plan["batches"][0],
        scratch_dir=tmp_path / "scratch-fragments",
        threads=1,
    )
    write_symmetric_edge_shard(
        data_dir=tmp_path,
        metric=metric,
        bucket=0,
        scratch_dir=tmp_path / "scratch-shard",
        threads=1,
    )
    connectivity_labels = pd.DataFrame(
        {
            "ligand_id": ligands,
            "label": [
                "component0",
                "component0",
                "component0",
                "component0",
                "component1",
            ],
        }
    )
    connectivity_labels.to_parquet(component_dir / "threshold=50.parquet", index=False)

    make_communities(
        data_dir=tmp_path,
        metric=metric,
        threshold=50,
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    communities = pd.read_parquet(
        tmp_path
        / "ligand_clusters"
        / "cluster=communities"
        / "directed=False"
        / f"metric={metric}"
        / "threshold=50.parquet"
    )
    assert set(communities["ligand_id"]) == set(ligands)
    joined = communities.merge(connectivity_labels, on="ligand_id")
    assert joined.groupby("label_x")["label_y"].nunique().max() == 1
    assert _component_partition(communities) == {
        frozenset({"l1", "l2"}),
        frozenset({"l3", "l4"}),
        frozenset({"l5"}),
    }
    assert set(communities["directed"]) == {False}


def test_directed_cover_uses_query_to_centroid_scores_and_reassigns():
    from plinder.data.clusters import _greedy_directed_centroid_cover

    nodes = ["a", "b", "c", "d"]
    graph = nk.Graph(len(nodes), weighted=True, directed=True)
    graph.addEdges(
        (
            np.asarray([0.8, 0.6, 0.9]),
            (
                np.asarray([0, 2, 2], dtype=np.uint),
                np.asarray([1, 1, 3], dtype=np.uint),
            ),
        )
    )

    assignments = _greedy_directed_centroid_cover(graph, nodes)

    assert assignments == [
        ("a", "b", pytest.approx(80.0)),
        ("b", "b", 100.0),
        ("c", "d", pytest.approx(90.0)),
        ("d", "d", 100.0),
    ]


def test_community_stream_rejects_edges_crossing_components(tmp_path):
    from plinder.data.clusters import (
        make_communities,
        prepare_symmetric_edge_plan,
        write_symmetric_edge_fragment_batch,
        write_symmetric_edge_shard,
    )

    metric = "sucos_shape_pocket_qcov"
    index_dir = tmp_path / "index"
    score_dir = tmp_path / "scores" / "search_db=holo"
    component_dir = (
        tmp_path
        / "ligand_clusters"
        / "reductions"
        / f"metric={metric}"
        / "labels"
        / "directed=false"
    )
    index_dir.mkdir(parents=True)
    score_dir.mkdir(parents=True)
    component_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["s1", "s2"],
            "ligand_id": ["l1", "l2"],
            "system_type": ["holo", "holo"],
            "system_num_protein_chains": [1, 1],
            "system_num_ligand_chains": [1, 1],
            "ligand_is_proper": [True, True],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": ["s1", "s2"],
            "query_ligand_id": ["l1", "l2"],
            "target_system": ["s2", "s1"],
            "target_ligand_id": ["l2", "l1"],
            "metric": [metric, metric],
            "similarity": [80, 70],
        }
    ).to_parquet(score_dir / "part.parquet", index=False)
    plan = prepare_symmetric_edge_plan(
        data_dir=tmp_path,
        metrics=[metric],
        source_batch_size=1,
        bucket_count=1,
    )
    write_symmetric_edge_fragment_batch(
        data_dir=tmp_path,
        batch=plan["batches"][0],
        scratch_dir=tmp_path / "scratch-fragments",
        threads=1,
    )
    write_symmetric_edge_shard(
        data_dir=tmp_path,
        metric=metric,
        bucket=0,
        scratch_dir=tmp_path / "scratch-shard",
        threads=1,
    )
    pd.DataFrame(
        {"ligand_id": ["l1", "l2"], "label": ["component0", "component1"]}
    ).to_parquet(component_dir / "threshold=50.parquet", index=False)

    with pytest.raises(ValueError, match="crossing edges"):
        make_communities(
            data_dir=tmp_path,
            metric=metric,
            threshold=50,
            scratch_dir=tmp_path / "scratch",
            threads=1,
        )


def test_ligand_clusters_are_merged_without_system_projection(tmp_path):
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
    proper_annotation = annotation[annotation["ligand_is_proper"]]
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": proper_annotation[
                "ligand_rdkit_canonical_smiles"
            ],
            "ligand_smiles_id": range(len(proper_annotation)),
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

    ligand_cluster = (
        tmp_path / "ligand_clusters/cluster=components/directed=False/"
        "metric=sucos_shape_pocket_qcov/threshold=50.parquet"
    )
    ligand_cluster.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "ligand_id": [ligand_a1, ligand_a2, ligand_b, ligand_c],
            "label": ["c0", "c1", "c0", "c2"],
        }
    ).to_parquet(ligand_cluster, index=False)
    directed_cover = (
        tmp_path
        / "ligand_sampling/directed_set_cover"
        / "metric=sucos_shape_pocket_qcov/threshold=50.parquet"
    )
    directed_cover.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "ligand_id": [ligand_a1, ligand_a2, ligand_b, ligand_c],
            "label": ["d0", "d1", "d0", "d2"],
        }
    ).to_parquet(directed_cover, index=False)

    system_cluster = (
        tmp_path
        / "clusters/cluster=components/directed=True/metric=sucos_shape_pocket_qcov/threshold=50.parquet"
    )
    ligand_labels = pd.read_parquet(ligand_cluster).set_index("ligand_id")["label"]
    assert not system_cluster.exists()
    assert ligand_labels[ligand_a1] == ligand_labels[ligand_b]
    assert ligand_labels[ligand_a2] != ligand_labels[ligand_a1]
    assert ligand_labels[ligand_c] != ligand_labels[ligand_a1]
    assert ligand_d not in ligand_labels

    finalized = finalize_index(data_dir=tmp_path)
    ligand_column = "sucos_shape_pocket_qcov__50__ligand__component"
    assert "sucos_shape_pocket_qcov__50__component" not in finalized
    finalized_labels = finalized.set_index("ligand_id")[ligand_column]
    assert finalized_labels[ligand_a1] == finalized_labels[ligand_b]
    assert finalized_labels[ligand_a2] != finalized_labels[ligand_a1]
    assert not bool(
        finalized.set_index("ligand_id").loc[ligand_b, "ligand_is_3d_score_able"]
    )


def test_finalize_index_rejects_stale_ligand_cluster_universe(tmp_path):
    from plinder.data.pipeline.utils import finalize_index

    index_dir = tmp_path / "index"
    ligand_dir = tmp_path / "ligands"
    fingerprint_dir = tmp_path / "fingerprints"
    cluster_file = (
        tmp_path
        / "ligand_clusters/cluster=components/directed=False"
        / "metric=sucos_shape_pocket_qcov/threshold=30.parquet"
    )
    index_dir.mkdir(parents=True)
    ligand_dir.mkdir()
    fingerprint_dir.mkdir()
    cluster_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["s1", "s2"],
            "system_id_no_biounit": ["s1", "s2"],
            "system_type": ["holo", "holo"],
            "ligand_id": ["l1", "l2"],
            "ligand_is_proper": [True, True],
            "ligand_rdkit_canonical_smiles": ["CC", "CCC"],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "ligand_id": ["l1", "l2"],
            "ligand_is_3d_score_able": [True, True],
        }
    ).to_parquet(ligand_dir / "part.parquet", index=False)
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": ["CC", "CCC"],
            "ligand_smiles_id": [0, 1],
        }
    ).to_parquet(
        fingerprint_dir / "ligand_similarity_annotations.parquet",
        index=False,
    )
    pd.DataFrame({"ligand_id": ["l1"], "label": ["c0"]}).to_parquet(
        cluster_file,
        index=False,
    )
    directed_cover = (
        tmp_path
        / "ligand_sampling/directed_set_cover"
        / "metric=sucos_shape_pocket_qcov/threshold=30.parquet"
    )
    directed_cover.parent.mkdir(parents=True)
    pd.DataFrame({"ligand_id": ["l1"], "label": ["d0"]}).to_parquet(
        directed_cover,
        index=False,
    )

    with pytest.raises(ValueError, match="current eligible ligand universe"):
        finalize_index(data_dir=tmp_path)


def test_finalize_index_rejects_clusters_from_before_targeted_repair(tmp_path):
    from plinder.data.pipeline.utils import finalize_index

    index_dir = tmp_path / "index"
    ligand_dir = tmp_path / "ligands"
    fingerprint_dir = tmp_path / "fingerprints"
    cluster_file = (
        tmp_path
        / "ligand_clusters/cluster=components/directed=False"
        / "metric=sucos_shape_pocket_qcov/threshold=30.parquet"
    )
    index_dir.mkdir(parents=True)
    ligand_dir.mkdir()
    fingerprint_dir.mkdir()
    cluster_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["s1"],
            "system_id_no_biounit": ["s1"],
            "system_type": ["holo"],
            "ligand_id": ["l1"],
            "ligand_is_proper": [True],
            "ligand_rdkit_canonical_smiles": ["CC"],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "ligand_id": ["l1"],
            "ligand_is_3d_score_able": [True],
        }
    ).to_parquet(ligand_dir / "part.parquet", index=False)
    pd.DataFrame({"ligand_id": ["l1"], "label": ["c0"]}).to_parquet(
        cluster_file,
        index=False,
    )
    directed_cover = (
        tmp_path
        / "ligand_sampling/directed_set_cover"
        / "metric=sucos_shape_pocket_qcov/threshold=30.parquet"
    )
    directed_cover.parent.mkdir(parents=True)
    pd.DataFrame({"ligand_id": ["l1"], "label": ["d0"]}).to_parquet(
        directed_cover,
        index=False,
    )
    (index_dir / "collation.json").write_text(
        '{"status": "requires_downstream_repair"}'
    )
    fingerprint_path = fingerprint_dir / "ligand_similarity_annotations.parquet"
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": ["CC"],
            "ligand_smiles_id": [0],
        }
    ).to_parquet(fingerprint_path, index=False)
    marker_mtime_ns = (index_dir / "collation.json").stat().st_mtime_ns
    utime(
        fingerprint_path,
        ns=(marker_mtime_ns + 1_000_000_000,) * 2,
    )

    with pytest.raises(ValueError, match="cluster artifact predates"):
        finalize_index(data_dir=tmp_path)


def test_default_cluster_metrics_use_only_pocket_weighted_shape():
    from plinder.core.scores.metrics import LIGAND_SCORE_NAMES, is_ligand_level_metric
    from plinder.data.pipeline.config import METRICS

    assert "shape" in LIGAND_SCORE_NAMES
    assert "color" in LIGAND_SCORE_NAMES
    assert "sucos_shape_pocket_qcov" in METRICS
    assert "tanimoto_similarity_ecfp4_1024" in METRICS
    assert "shape" not in METRICS
    assert "color" not in METRICS
    assert "pocket_fident" not in METRICS
    assert "pocket_fident_qcov" not in METRICS
    assert "shape_tanimoto" not in METRICS
    assert "sucos_shape" not in METRICS
    assert not any(metric.startswith("protein_lddt") for metric in METRICS)
    assert not any(metric.startswith("protein_") for metric in METRICS)
    assert all(is_ligand_level_metric(metric) for metric in METRICS)


def test_component_node_universe_keeps_large_holo_targets_in_prepared_cache(
    tmp_path, monkeypatch
):
    from plinder.data import clusters

    index_dir = tmp_path / "index"
    index_dir.mkdir()
    pd.DataFrame(
        {
            "system_id": ["s1", "s2", "s3"],
            "ligand_id": ["l1", "l2", "l3"],
            "system_type": ["holo", "holo", "holo"],
            "system_num_protein_chains": [1, 6, 1],
            "system_num_ligand_chains": [1, 1, 1],
            "ligand_is_proper": [True, True, False],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)

    report = clusters.prepare_component_node_universe(tmp_path)
    assert report["ligand_count"] == 2
    assert report["system_count"] == 2
    monkeypatch.setattr(
        clusters,
        "_eligible_annotation",
        lambda data_dir: pytest.fail("prepared node universe was not reused"),
    )
    nodes, systems = clusters.component_node_universe(
        data_dir=tmp_path,
        metric="pocket_qcov",
    )
    assert nodes == ["l1", "l2"]
    assert systems == {"s1", "s2"}


def test_tanimoto_clusters_expand_unique_smiles_to_ligands(tmp_path):
    from plinder.data.clusters import expand_fingerprint_clusters_to_ligands

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
            "ligand_id": [
                "1aaa__1__1.X",
                "1aaa__1__1.Y",
                "2bbb__1__1.Z",
                "3ccc__1__1.W",
            ],
            "ligand_molecular_weight": [100.0, 200.0, 150.0, 100.0],
            "ligand_rdkit_canonical_smiles": ["CC", "CCC", "CCCC", "CC"],
            "ligand_is_proper": [True] * 4,
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

    labels = expand_fingerprint_clusters_to_ligands(
        data_dir=tmp_path,
        labeldf=pd.DataFrame(
            {
                "ligand_id": [0, 1, 2],
                "label": ["c0", "c1", "c0"],
            }
        ),
    ).set_index("ligand_id")["label"]
    assert labels["1aaa__1__1.X"] == labels["2bbb__1__1.Z"]
    assert labels["3ccc__1__1.W"] == labels["2bbb__1__1.Z"]
    assert labels["1aaa__1__1.Y"] != labels["2bbb__1__1.Z"]
    assert not (tmp_path / "clusters").exists()
    assert not (fingerprint_dir / "ligands_per_system.parquet").exists()


@pytest.mark.parametrize("metric", ["shape", "color", "sucos_shape"])
def test_gated_shape_diagnostics_cannot_be_clustered(tmp_path, metric):
    from plinder.data.pipeline.score import _cluster_parameters

    with pytest.raises(ValueError, match="sucos_shape_pocket_qcov"):
        _cluster_parameters(metrics=[metric], thresholds=[50])


def test_v3_rejects_scores_without_ligand_identifiers(tmp_path):
    from plinder.data.clusters import iter_component_score_batches

    index_dir = tmp_path / "index"
    score_dir = tmp_path / "scores" / "search_db=holo"
    index_dir.mkdir(parents=True)
    score_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["1aaa__1__1.A__1.X"],
            "ligand_id": ["1aaa__1__1.X"],
            "system_type": ["holo"],
            "system_num_protein_chains": [1],
            "system_num_ligand_chains": [1],
            "ligand_is_proper": [True],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": ["1aaa__1__1.A__1.X"],
            "target_system": ["1aaa__1__1.A__1.X"],
            "metric": ["protein_lddt_weighted_sum"],
            "similarity": [100],
        }
    ).to_parquet(score_dir / "scores.parquet", index=False)

    with pytest.raises(ValueError, match="missing"):
        list(
            iter_component_score_batches(
                source_path=score_dir / "scores.parquet",
                metric="pocket_qcov",
                eligible_systems={"1aaa__1__1.A__1.X"},
            )
        )
