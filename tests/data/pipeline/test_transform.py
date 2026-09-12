# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path

import pandas as pd
import pytest

from plinder.data.pipeline.transform import transform_bindingdb_affinity_data


def test_bindingdb_transform_aggregates_across_chunks(tmp_path: Path) -> None:
    raw_path = tmp_path / "BindingDB_All.tsv"
    pd.DataFrame(
        {
            "Ligand HET ID in PDB": ["LIG", "LIG", "ATP", None],
            "PDB ID(s) for Ligand-Target Complex": [
                "1abc",
                "1abc",
                "2def, 3ghi",
                None,
            ],
            "Ki (nM)": [10.0, 100.0, 1.0, 5.0],
            "Kd (nM)": [None, None, None, None],
            "EC50 (nM)": [None, None, None, None],
            "BindingDB Target Chain Sequence 1": [
                "AAAA",
                "AAAA",
                "BBBB",
                "CCCC",
            ],
        }
    ).to_csv(raw_path, sep="\t", index=False)

    result = transform_bindingdb_affinity_data(
        raw_affinity_path=raw_path,
        chunksize=1,
    ).set_index("pdbid_ligid")

    assert set(result.index) == {"1ABC_LIG", "2DEF_ATP", "3GHI_ATP"}
    assert result.loc["1ABC_LIG", "pchembl"] == pytest.approx(7.5)
    assert result.loc["1ABC_LIG", "target_sequence"] == "AAAA"
    assert result.loc["2DEF_ATP", "pchembl"] == pytest.approx(9.0)


def test_bindingdb_transform_preserves_numeric_identifiers_as_strings(
    tmp_path: Path,
) -> None:
    raw_path = tmp_path / "BindingDB_numeric_ids.tsv"
    pd.DataFrame(
        {
            "Ligand HET ID in PDB": [123],
            "PDB ID(s) for Ligand-Target Complex": [1234],
            "Ki (nM)": [10.0],
            "Kd (nM)": [None],
            "EC50 (nM)": [None],
            "BindingDB Target Chain Sequence 1": ["AAAA"],
        }
    ).to_csv(raw_path, sep="\t", index=False)

    result = transform_bindingdb_affinity_data(raw_affinity_path=raw_path)

    assert result["pdbid_ligid"].tolist() == ["1234_123"]
    assert result["pchembl"].tolist() == pytest.approx([8.0])
