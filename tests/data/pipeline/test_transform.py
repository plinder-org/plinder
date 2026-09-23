# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path

import pandas as pd
import pytest
from plinder.data.pipeline.io import download_affinity_data
from plinder.data.pipeline.transform import (
    strict_bindingdb_candidates,
    transform_bindingdb_measurements,
)


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

    records = transform_bindingdb_measurements(raw_affinity_path=raw_path, chunksize=1)
    result = strict_bindingdb_candidates(records).set_index("pdbid_ligid")

    assert set(result.index) == {"1ABC_LIG", "2DEF_ATP", "3GHI_ATP"}
    assert result.loc["1ABC_LIG", "pchembl"] == pytest.approx(7.5)
    assert result.loc["1ABC_LIG", "target_sequence"] == "AAAA"
    assert result.loc["1ABC_LIG", "count"] == 2
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

    records = transform_bindingdb_measurements(raw_affinity_path=raw_path)
    result = strict_bindingdb_candidates(records)

    assert result["pdbid_ligid"].tolist() == ["1234_123"]
    assert result["pchembl"].tolist() == pytest.approx([8.0])


# Regression cases contributed by Jacob Casey in
# https://github.com/plinder-org/plinder/issues/94#issuecomment-5784851981
def test_bindingdb_measurements_keep_target_associations(tmp_path: Path) -> None:
    base = {
        "Ligand HET ID in PDB": "LIG",
        "PDB ID(s) for Ligand-Target Complex": "1abc",
        "Ki (nM)": None,
        "Kd (nM)": None,
    }
    rows = [
        {**base, "Ki (nM)": "1", "BindingDB Target Chain Sequence 1": "ACDEFG"},
        {**base, "Ki (nM)": "1000", "BindingDB Target Chain Sequence 1": "ACDEFA"},
    ]
    raw_path = tmp_path / "BindingDB_All.tsv"
    for order in (rows, list(reversed(rows))):
        pd.DataFrame(order).to_csv(raw_path, sep="\t", index=False)
        records = transform_bindingdb_measurements(
            raw_affinity_path=raw_path, chunksize=1
        )
        assert len(records) == 2
        candidates = strict_bindingdb_candidates(records).set_index("target_sequence")
        assert candidates.loc["ACDEFG", "pchembl"] == pytest.approx(9.0)
        assert candidates.loc["ACDEFA", "pchembl"] == pytest.approx(6.0)


def test_bindingdb_measurements_keep_bounds_and_endpoints(tmp_path: Path) -> None:
    raw_path = tmp_path / "BindingDB_All.tsv"
    pd.DataFrame(
        [
            {
                "Ligand HET ID in PDB": "LIG",
                "PDB ID(s) for Ligand-Target Complex": "1abc",
                "BindingDB Target Chain Sequence 1": "ACDEFG",
                "Ki (nM)": ">100",
                "Kd (nM)": None,
            },
            {
                "Ligand HET ID in PDB": "LIG",
                "PDB ID(s) for Ligand-Target Complex": "2def",
                "BindingDB Target Chain Sequence 1": "ACDEFG",
                "Ki (nM)": "1",
                "Kd (nM)": "1000",
            },
        ]
    ).to_csv(raw_path, sep="\t", index=False)
    records = transform_bindingdb_measurements(raw_affinity_path=raw_path)

    assert len(records) == 3
    bound = records.loc[records["pdbid_ligid"].eq("1ABC_LIG")].iloc[0]
    assert bound["endpoint"] == "Ki"
    assert bound["relation"] == "<"
    assert bound["pchembl"] == pytest.approx(7.0)
    assert set(records.loc[records["pdbid_ligid"].eq("2DEF_LIG"), "endpoint"]) == {
        "Ki",
        "Kd",
    }
    assert strict_bindingdb_candidates(records).empty


def test_bindingdb_multichain_target_is_not_a_scalar(tmp_path: Path) -> None:
    raw_path = tmp_path / "BindingDB_All.tsv"
    pd.DataFrame(
        [
            {
                "BindingDB Reactant_set_id": "123",
                "BindingDB MonomerID": "456",
                "Curation/DataSource": "ChEMBL",
                "Article DOI": "10.1000/example",
                "BindingDB Entry DOI": "10.7270/example",
                "Ligand HET ID in PDB": "LIG",
                "PDB ID(s) for Ligand-Target Complex": "1abc",
                "BindingDB Target Chain Sequence 1": "ACDEFG",
                "BindingDB Target Chain Sequence 2": "HIJKLM",
                "Ki (nM)": "10",
                "Kd (nM)": None,
            }
        ]
    ).to_csv(raw_path, sep="\t", index=False)
    records = transform_bindingdb_measurements(raw_affinity_path=raw_path)

    assert records.loc[0, "reactant_set_id"] == "123"
    assert records.loc[0, "monomer_id"] == "456"
    assert records.loc[0, "curation_source"] == "ChEMBL"
    assert records.loc[0, "article_doi"] == "10.1000/example"
    assert records.loc[0, "bindingdb_entry_doi"] == "10.7270/example"
    assert records.loc[0, "target_sequences"] == ["ACDEFG", "HIJKLM"]
    assert pd.isna(records.loc[0, "target_sequence"])
    assert strict_bindingdb_candidates(records).empty


def test_affinity_cache_keeps_records_and_target_candidates(tmp_path: Path) -> None:
    affinity_dir = tmp_path / "dbs" / "affinity"
    affinity_dir.mkdir(parents=True)
    raw_path = affinity_dir / "BindingDB_All.tsv"
    pd.DataFrame(
        {
            "Ligand HET ID in PDB": ["LIG", "LIG"],
            "PDB ID(s) for Ligand-Target Complex": ["1abc", "1abc"],
            "BindingDB Target Chain Sequence 1": ["ACDEFG", "HIJKLM"],
            "Ki (nM)": ["10", "100"],
            "Kd (nM)": [None, None],
        }
    ).to_csv(raw_path, sep="\t", index=False)

    candidates = download_affinity_data(data_dir=tmp_path)

    assert {row["target_sequence"] for row in candidates["1ABC_LIG"]} == {
        "ACDEFG",
        "HIJKLM",
    }
    assert len(pd.read_parquet(affinity_dir / "measurements.parquet")) == 2
    assert len(pd.read_parquet(affinity_dir / "candidates.parquet")) == 2
    assert not (affinity_dir / "affinity.json").exists()
