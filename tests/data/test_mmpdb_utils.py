from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest
from plinder.core.structure.smallmols_similarity import (
    get_mmp_similarity_dict,
    smiles2nonstereo,
)
from plinder.core.utils import schemas
from plinder.data.annotations import mmpdb_utils


def _write_ligands(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCN"],
        }
    ).to_parquet(path, index=False)


def _write_pair_file(path: Path) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write("CCO\tCCN\t0\t1\t*O>>*N\t*C\n")
        handle.write("CCO\tCCN\t0\t1\t*O>>*N\t*C\n")


def test_mmp_pair_conversion_uses_smiles_ids_and_compact_schema(tmp_path: Path):
    fingerprint_path = tmp_path / "fingerprints.parquet"
    _write_ligands(fingerprint_path)
    ligands, exclusions = mmpdb_utils._ligand_table(fingerprint_path)
    assert exclusions == []
    pair_file = tmp_path / "pairs.csv.gz"
    _write_pair_file(pair_file)
    output = tmp_path / "pairs.parquet"

    mmpdb_utils._write_pair_parquet(
        pair_files=[pair_file],
        ligands=ligands,
        output_path=output,
    )

    import pyarrow.parquet as pq

    assert pq.read_schema(output).equals(schemas.LIGAND_MMP_PAIR_SCHEMA)
    pairs = pd.read_parquet(output)
    assert len(pairs) == 1
    pair = pairs.iloc[0]
    assert pair["ligand_smiles_id_1"] == 0
    assert pair["ligand_smiles_id_2"] == 1
    assert pair["ligand_smiles_1"] == "CCO"
    assert pair["ligand_smiles_2"] == "CCN"
    assert pair["transformation"] == "*O>>*N"
    assert pair["shared_core_smiles"] == "*C"
    assert pair["num_cuts"] == 1
    assert pair["shared_core_num_heavy_atoms"] == 1
    assert pair["ligand_1_num_heavy_atoms"] == 3
    assert pair["ligand_2_num_heavy_atoms"] == 3
    assert pair["ligand_1_shared_core_fraction"] == pytest.approx(1 / 3)
    assert pair["ligand_2_shared_core_fraction"] == pytest.approx(1 / 3)

    similarities = get_mmp_similarity_dict(output, min_constant_size=1)
    ethanol = smiles2nonstereo("CCO")
    ethylamine = smiles2nonstereo("CCN")
    assert similarities[ethanol][ethylamine] == pytest.approx(100 / 3)
    assert similarities[ethylamine][ethanol] == pytest.approx(100 / 3)


def test_mmp_builder_reuses_output_for_the_same_smiles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    fingerprint_path = tmp_path / "fingerprints" / "ligands_per_smiles.parquet"
    _write_ligands(fingerprint_path)
    calls = 0

    def fake_generate_pair_files(*, work_dir: Path, **_kwargs):
        nonlocal calls
        calls += 1
        pair_file = work_dir / "partition.0000.csv.gz"
        _write_pair_file(pair_file)
        return [pair_file]

    monkeypatch.setattr(mmpdb_utils, "_generate_pair_files", fake_generate_pair_files)
    monkeypatch.setattr(mmpdb_utils.shutil, "which", lambda _name: "/mock/mmpdb")

    first = mmpdb_utils.make_ligand_mmp_pairs(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )
    second = mmpdb_utils.make_ligand_mmp_pairs(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )

    assert first == second == tmp_path / "index" / "ligand_mmp_pairs.parquet"
    assert calls == 1
    assert not list((tmp_path / "scratch").glob("plinder-mmp-*"))


def test_mmp_builder_rejects_duplicate_smiles_ids(tmp_path: Path):
    path = tmp_path / "ligands.parquet"
    pd.DataFrame(
        {
            "ligand_smiles_id": [0, 0],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCN"],
        }
    ).to_parquet(path, index=False)

    with pytest.raises(ValueError, match="duplicate ligand_smiles_id"):
        mmpdb_utils._ligand_table(path)


@pytest.mark.parametrize(
    "smiles,reason",
    [
        (
            "[B]12[B]3[B]4[B]1[B]5[B]6[B]2[B]3[B]45C6",
            "unsupported_by_mmpdb_smiles_parser",
        ),
        ("[Be](C)(C)(C)C", "unsupported_by_mmpdb_smiles_parser"),
        ("N->[Cu+2]<-N", "metal_dative_bonds"),
    ],
)
def test_mmp_input_exclusions_preserve_supported_ids(tmp_path, smiles, reason):
    path = tmp_path / "ligands.parquet"
    original = pd.DataFrame(
        {
            "ligand_smiles_id": [12, 4, 8],
            "ligand_rdkit_canonical_smiles": ["CCO", smiles, "CCN"],
        }
    )
    original.to_parquet(path, index=False)

    ligands, exclusions = mmpdb_utils._ligand_table(path)

    assert ligands["ligand_smiles_id"].tolist() == [8, 12]
    assert ligands["ligand_rdkit_canonical_smiles"].tolist() == ["CCN", "CCO"]
    assert ligands["num_heavy_atoms"].tolist() == [3, 3]
    assert exclusions == [{"ligand_smiles_id": 4, "smiles": smiles, "reason": reason}]
    pd.testing.assert_frame_equal(pd.read_parquet(path), original)


def test_mmp_invalid_smiles_still_raise(tmp_path):
    path = tmp_path / "ligands.parquet"
    pd.DataFrame(
        {
            "ligand_smiles_id": [7],
            "ligand_rdkit_canonical_smiles": ["not-a-smiles"],
        }
    ).to_parquet(path, index=False)
    with pytest.raises(ValueError, match="ligand_smiles_id 7 has invalid"):
        mmpdb_utils._ligand_table(path)


@pytest.mark.parametrize("include_supported", [False, True])
def test_mmp_build_records_exclusions_and_keeps_supported_pairs(
    tmp_path, include_supported
):
    """Exercise the real mmpdb CLI, including an entirely unsupported input."""
    path = tmp_path / "fingerprints" / "ligands_per_smiles.parquet"
    path.parent.mkdir()
    smiles = ["[Be](C)(C)(C)C"]
    if include_supported:
        smiles.extend(["CCO", "CCN"])
    pd.DataFrame(
        {
            "ligand_smiles_id": range(len(smiles)),
            "ligand_rdkit_canonical_smiles": smiles,
        }
    ).to_parquet(path, index=False)

    output = mmpdb_utils.make_ligand_mmp_pairs(
        data_dir=tmp_path, scratch_dir=tmp_path / "scratch", threads=2
    )

    pairs = pd.read_parquet(output)
    if include_supported:
        assert not pairs.empty
        assert set(pairs["ligand_smiles_id_1"]) | set(pairs["ligand_smiles_id_2"]) == {
            1,
            2,
        }
    else:
        assert pairs.empty
    manifest = json.loads(output.with_suffix(".manifest.json").read_text())
    assert manifest["input_exclusions"] == [
        {
            "ligand_smiles_id": 0,
            "smiles": smiles[0],
            "reason": "unsupported_by_mmpdb_smiles_parser",
        }
    ]
