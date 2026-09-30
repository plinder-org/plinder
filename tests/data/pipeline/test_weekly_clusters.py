from pathlib import Path

import pandas as pd

from plinder.data.pipeline.weekly_clusters import (
    _best_representatives,
    extend_interface_clusters,
    extend_ligand_clusters,
)


def test_weekly_assignment_keeps_qualifying_previous_representative():
    scores = pd.DataFrame(
        {
            "query_node": ["new", "new"],
            "target_node": ["first", "second"],
            "similarity": [99, 91],
        }
    )
    representatives = pd.DataFrame(
        {"centroid": ["first", "second"], "label": ["a", "b"]}
    )
    result = _best_representatives(
        scores,
        representatives,
        threshold=90,
        directed=True,
        previous_labels=pd.Series({"new": "b"}),
    )
    assert result.iloc[0]["label"] == "b"


def test_weekly_ligands_join_existing_representatives(tmp_path: Path):
    base = tmp_path / "base"
    updated = tmp_path / "updated"
    (base / "index").mkdir(parents=True)
    (updated / "index").mkdir(parents=True)
    (updated / "exports/ligand_similarity_scores").mkdir(parents=True)
    (updated / "ligand_scores").mkdir()
    (updated / "mhfp6_scores").mkdir()
    (updated / "fingerprints").mkdir()
    direct_column = "pocket_qcov__90__ligand__directed_set_cover"
    direct_component = "pocket_qcov__90__ligand__reciprocal_component"
    chemical_column = "tanimoto_similarity_ecfp4_1024__90__ligand__set_cover"
    chemical_component = (
        "tanimoto_similarity_ecfp4_1024__90__ligand__reciprocal_component"
    )
    summary = "ligand_tanimoto_ecfp4_1024_90_cluster"
    pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X", "2bbb__1__1.X"],
            direct_column: ["c0", "c1"],
            direct_component: ["r0", "r1"],
            f"{direct_column}__is_centroid": [True, True],
            f"{direct_column}__coverage_count": [2, 1],
            f"{direct_column}__coverage_fraction": [1.0, 1.0],
            chemical_column: ["t0", "t1"],
            chemical_component: ["u0", "u1"],
            f"{chemical_column}__is_centroid": [True, True],
            summary: ["t0", "t1"],
            f"{summary}_num_pdb_ids": [1, 1],
        }
    ).to_parquet(base / "index/ligand_clusters.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1aaa", "2bbb", "3ccc"],
            "ligand_id": ["1aaa__1__1.X", "2bbb__1__1.X", "3ccc__1__1.X"],
            "ligand_smiles": ["C", "CC", "CCC"],
            "ligand_is_proper": [True, True, True],
            "system_type": ["holo", "holo", "holo"],
        }
    ).to_parquet(updated / "index/annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": ["C", "CC", "CCC"],
            "ligand_smiles_id": [0, 1, 2],
        }
    ).to_parquet(
        updated / "fingerprints/ligand_similarity_annotations.parquet", index=False
    )
    pd.DataFrame(
        {
            "query_ligand_id": ["3ccc__1__1.X"],
            "target_ligand_id": ["1aaa__1__1.X"],
            "pocket_qcov": [95],
            "pli_qcov": [70],
            "sucos_shape": [80],
        }
    ).to_parquet(updated / "exports/ligand_similarity_scores/cc.parquet", index=False)
    pd.DataFrame(
        {
            "query_ligand_id": [2],
            "target_ligand_id": [0],
            "tanimoto_similarity_ecfp4_1024": [96.0],
        }
    ).to_parquet(updated / "ligand_scores/added.parquet", index=False)

    result = extend_ligand_clusters(
        base=base, data_dir=updated, affected={"3ccc"}
    ).set_index("ligand_id")

    assert result.at["3ccc__1__1.X", direct_column] == "c0"
    assert result.at["3ccc__1__1.X", chemical_column] == "t0"
    assert pd.isna(result.at["3ccc__1__1.X", direct_component])
    assert pd.isna(result.at["3ccc__1__1.X", chemical_component])
    assert result.at["1aaa__1__1.X", direct_component] == "r0"
    assert not result.at["3ccc__1__1.X", f"{direct_column}__is_centroid"]
    assert pd.isna(result.at["3ccc__1__1.X", f"{direct_column}__coverage_count"])
    assert result.at["1aaa__1__1.X", f"{summary}_num_pdb_ids"] == 2


def test_weekly_interfaces_join_existing_whole_and_side_covers(tmp_path: Path):
    base = tmp_path / "base"
    updated = tmp_path / "updated"
    (base / "index").mkdir(parents=True)
    (updated / "index").mkdir(parents=True)
    (updated / "interface_scores").mkdir()
    whole = "interface_qcov__90__directed_set_cover"
    whole_component = "interface_qcov__90__reciprocal_component"
    side = "interface_side_qcov__90__chain_1_directed_set_cover"
    side_component = "interface_side_qcov__90__chain_1_reciprocal_component"
    side_two = "interface_side_qcov__90__chain_2_directed_set_cover"
    old = "1aaa__1__1.A--1.B"
    new = "2bbb__1__1.A--1.B"
    pd.DataFrame(
        {
            "system_id": [old],
            whole: ["c0"],
            whole_component: ["r0"],
            side: ["s0"],
            side_component: ["u0"],
            side_two: ["s1"],
        }
    ).to_parquet(base / "index/interface_clusters.parquet", index=False)
    pd.DataFrame(
        {"entry_pdb_id": ["1aaa", "2bbb"], "system_id": [old, new]}
    ).to_parquet(updated / "index/interface_annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "system_id": [old, new],
            "representative_system_id": [old, new],
            "side_1_half_interface_id": [f"{old}::side=1", f"{new}::side=1"],
            "side_2_half_interface_id": [f"{old}::side=2", f"{new}::side=2"],
        }
    ).to_parquet(updated / "index/interface_membership.parquet", index=False)
    for metric, nodes, labels in (
        ("interface_qcov", [old], ["c0"]),
        (
            "interface_side_qcov",
            [f"{old}::side=1", f"{old}::side=2"],
            ["s0", "s1"],
        ),
    ):
        directory = base / f"interface_sampling/directed_set_cover/metric={metric}"
        directory.mkdir(parents=True)
        pd.DataFrame(
            {"system_id": nodes, "centroid_system_id": nodes, "label": labels}
        ).to_parquet(directory / "threshold=90.parquet", index=False)
    pd.DataFrame(
        {
            "query_system": [new, f"{new}::side=1", f"{new}::side=2"],
            "target_system": [old, f"{old}::side=1", f"{old}::side=2"],
            "metric": [
                "interface_qcov",
                "interface_side_qcov",
                "interface_side_qcov",
            ],
            "similarity": [95, 92, 93],
        }
    ).to_parquet(updated / "interface_scores/shard=bb.parquet", index=False)

    result = extend_interface_clusters(
        base=base, data_dir=updated, affected={"2bbb"}
    ).set_index("system_id")

    assert result.at[new, whole] == "c0"
    assert result.at[new, side] == "s0"
    assert result.at[new, side_two] == "s1"
    assert pd.isna(result.at[new, whole_component])
    assert pd.isna(result.at[new, side_component])
    assert result.at[old, whole_component] == "r0"
