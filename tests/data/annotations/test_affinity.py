# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path

import pandas as pd
import pytest

from plinder.data.annotations.affinity import (
    build_ligand_affinity_table,
    matched_affinity,
    publish_affinity_tables,
)


def test_matched_affinity_selects_one_target() -> None:
    candidates = [
        {"target_sequence": "ACDEFG", "endpoint": "Ki", "pchembl": 9.0, "count": 1},
        {"target_sequence": "ACDEFA", "endpoint": "Ki", "pchembl": 6.0, "count": 1},
    ]
    assert matched_affinity(candidates, {"1.A": "ACDEFG"}) == (9.0, "Ki", 1)
    assert matched_affinity(candidates, {"1.A": "ACDEFA"}) == (6.0, "Ki", 1)


def test_matched_affinity_rejects_missing_or_ambiguous_targets() -> None:
    candidates = [
        {"target_sequence": "ACDEFG", "endpoint": "Ki", "pchembl": 9.0, "count": 1},
        {"target_sequence": "ACDEFA", "endpoint": "Ki", "pchembl": 6.0, "count": 1},
    ]
    assert matched_affinity(candidates, {}) is None
    assert matched_affinity(candidates, {"1.A": "TTTTTT"}) is None
    assert matched_affinity(candidates, {"1.A": "ACDEFG", "1.B": "ACDEFA"}) is None


def test_build_ligand_affinity_table_matches_receptor_sequences(tmp_path: Path) -> None:
    index_dir = tmp_path / "index"
    affinity_dir = tmp_path / "dbs" / "affinity"
    index_dir.mkdir()
    affinity_dir.mkdir(parents=True)
    annotation_path = index_dir / "annotation_table.parquet"
    chains_path = index_dir / "entry_chains.parquet"
    candidates_path = affinity_dir / "candidates.parquet"
    pd.DataFrame(
        {
            "ligand_id": ["lig1", "lig2", "lig3"],
            "system_id": ["sys1", "sys1", "sys2"],
            "entry_pdb_id": ["1abc", "1abc", "2def"],
            "ligand_ccd_code": ["LIG", "ATP", "LIG"],
            "ligand_protein_chains_asym_id": [["1.A"], ["1.A"], ["1.B"]],
        }
    ).to_parquet(annotation_path, index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "chain_asym_id": ["A", "B"],
            "chain_sequence": ["ACDEFG", "HIJKLM"],
        }
    ).to_parquet(chains_path, index=False)
    pd.DataFrame(
        {
            "pdbid_ligid": ["1ABC_LIG", "2DEF_LIG"],
            "target_sequence": ["ACDEFG", "ACDEFG"],
            "endpoint": ["Ki", "Ki"],
            "pchembl": [8.0, 9.0],
            "count": [2, 1],
        }
    ).to_parquet(candidates_path, index=False)

    result = build_ligand_affinity_table(
        annotation_path, chains_path, candidates_path
    ).set_index("ligand_id")

    assert result.loc["lig1", "ligand_binding_affinity"] == pytest.approx(8.0)
    assert result.loc["lig1", "ligand_binding_affinity_endpoint"] == "Ki"
    assert result.loc["lig1", "ligand_binding_affinity_measurement_count"] == 2
    assert pd.isna(result.loc["lig2", "ligand_binding_affinity"])
    assert bool(result.loc["lig2", "system_has_binding_affinity"])
    assert pd.isna(result.loc["lig3", "ligand_binding_affinity"])
    assert not bool(result.loc["lig3", "system_has_binding_affinity"])

    pd.DataFrame({"source_row": [7]}).to_parquet(
        affinity_dir / "measurements.parquet", index=False
    )
    publish_affinity_tables(tmp_path)
    published = pd.read_parquet(index_dir / "ligand_affinity.parquet")
    assert published["ligand_binding_affinity"].notna().sum() == 1
    assert pd.read_parquet(index_dir / "bindingdb_measurements.parquet")[
        "source_row"
    ].tolist() == [7]


def test_build_ligand_affinity_table_keeps_empty_schema(tmp_path: Path) -> None:
    annotation_path = tmp_path / "annotation.parquet"
    chains_path = tmp_path / "chains.parquet"
    candidates_path = tmp_path / "candidates.parquet"
    pd.DataFrame(
        columns=[
            "ligand_id",
            "system_id",
            "entry_pdb_id",
            "ligand_ccd_code",
            "ligand_protein_chains_asym_id",
        ]
    ).to_parquet(annotation_path)
    pd.DataFrame(
        columns=["entry_pdb_id", "chain_asym_id", "chain_sequence"]
    ).to_parquet(chains_path)
    pd.DataFrame(
        columns=["pdbid_ligid", "target_sequence", "endpoint", "pchembl", "count"]
    ).to_parquet(candidates_path)

    result = build_ligand_affinity_table(annotation_path, chains_path, candidates_path)

    assert result.empty
    assert str(result["ligand_binding_affinity"].dtype) == "float64"
    assert str(result["ligand_binding_affinity_measurement_count"].dtype) == "Int64"
    assert str(result["system_has_binding_affinity"].dtype) == "bool"


def test_updated_affinity_matches_full_rebuild(tmp_path: Path) -> None:
    base = tmp_path / "base"
    updated = tmp_path / "updated"
    for root in (base, updated):
        (root / "index").mkdir(parents=True)
    affinity_dir = base / "dbs" / "affinity"
    affinity_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "pdbid_ligid": ["1ABC_LIG", "2DEF_LIG", "3GHI_LIG"],
            "target_sequence": ["ACDEFG", "HIJKLM", "NPQRST"],
            "endpoint": ["Ki", "Ki", "Kd"],
            "pchembl": [8.0, 9.0, 7.0],
            "count": [1, 1, 1],
        }
    ).to_parquet(affinity_dir / "candidates.parquet", index=False)
    pd.DataFrame({"source_row": [7]}).to_parquet(
        affinity_dir / "measurements.parquet", index=False
    )

    def write_index(root: Path, entries: list[tuple[str, str, str]]) -> None:
        pd.DataFrame(
            {
                "ligand_id": [ligand for _, ligand, _ in entries],
                "system_id": [f"{pdb_id}__1" for pdb_id, _, _ in entries],
                "entry_pdb_id": [pdb_id for pdb_id, _, _ in entries],
                "ligand_ccd_code": ["LIG"] * len(entries),
                "ligand_protein_chains_asym_id": [["1.A"]] * len(entries),
            }
        ).to_parquet(root / "index/annotation_table.parquet", index=False)
        pd.DataFrame(
            {
                "entry_pdb_id": [pdb_id for pdb_id, _, _ in entries],
                "chain_asym_id": ["A"] * len(entries),
                "chain_sequence": [sequence for _, _, sequence in entries],
            }
        ).to_parquet(root / "index/entry_chains.parquet", index=False)

    write_index(
        base,
        [
            ("1abc", "old1", "ACDEFG"),
            ("2def", "old2", "HIJKLM"),
            ("4jkl", "old4", "AAAAAA"),
        ],
    )
    publish_affinity_tables(base)
    write_index(
        updated,
        [
            ("1abc", "old1", "ACDEFG"),
            ("2def", "new2", "HIJKLM"),
            ("3ghi", "new3", "NPQRST"),
        ],
    )
    publish_affinity_tables(
        updated,
        affinity_dir=affinity_dir,
        base_dir=base,
        changed_pdb_ids=["2def", "3ghi"],
        removed_pdb_ids=["2def", "4jkl"],
    )

    expected = build_ligand_affinity_table(
        updated / "index/annotation_table.parquet",
        updated / "index/entry_chains.parquet",
        affinity_dir / "candidates.parquet",
    )
    pd.testing.assert_frame_equal(
        pd.read_parquet(updated / "index/ligand_affinity.parquet"), expected
    )
    assert pd.read_parquet(updated / "index/bindingdb_measurements.parquet")[
        "source_row"
    ].tolist() == [7]
