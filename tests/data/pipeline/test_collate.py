# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from plinder.data.pipeline import collate as collate_module
from plinder.data.pipeline.collate import (
    collate_shard,
    finalize_collation,
    plan_collation,
    planned_code_batch,
    run_collation,
)


def _write_entry(
    data_dir: Path,
    pdb_id: str,
    *,
    ligand_rows: list[dict[str, object]],
    scoreability: dict[str, bool],
    ph: float | None = None,
) -> None:
    code = pdb_id[1:3]
    raw_root = data_dir / "raw_entries" / code
    raw_root.mkdir(parents=True, exist_ok=True)
    system_id = f"{pdb_id}__1__1.A__1.L"
    annotation_rows = []
    for ligand in ligand_rows:
        row: dict[str, object] = {
            "entry_pdb_id": pdb_id,
            "entry_pH": ph,
            "system_biounit_id": "1",
            "system_id": system_id,
            "system_type": "holo",
            "system_protein_chains_length": [100, 200],
            "ligand_id": ligand["ligand_id"],
            "ligand_unique_ccd_code": ligand["ccd"],
            "ligand_is_proper": ligand["proper"],
            "ligand_is_lipinski": ligand["proper"],
            "ligand_is_cofactor": ligand.get("cofactor", False),
            "ligand_is_fragment": False,
            "ligand_is_monosaccharide": False,
            "ligand_is_oligosaccharide": False,
            "ligand_is_mononucleotide": False,
            "ligand_is_oligonucleotide": False,
            "ligand_is_monopeptide": False,
            "ligand_is_oligopeptide": False,
            "ligand_is_artifact": ligand.get("artifact", False),
            "ligand_is_other": False,
            "ligand_is_covalent": False,
            "ligand_is_invalid": False,
            "ligand_is_ion": ligand.get("ion", False),
            "system_pocket_ECOD": "retired",
            "system_pocket_kinase_name": "retired",
        }
        annotation_rows.append(row)
    pd.DataFrame(annotation_rows).to_parquet(
        raw_root / f"{pdb_id}.parquet", index=False
    )

    entry_dir = raw_root / pdb_id
    entry_dir.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["A"],
            "chain_entity_id": ["1"],
            "chain_type": ["polypeptide(L)"],
            "chain_receptor_type": ["protein"],
            "chain_length": [300],
            "chain_num_unresolved_residues": [0],
            "chain_is_holo": [True],
            "chain_uniprot_ids": [["P12345"]],
        }
    ).to_parquet(entry_dir / "entry_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "biounit_id": ["1"],
            "chain_instance": ["1.A"],
            "chain_asym_id": ["A"],
            "chain_role": ["receptor"],
        }
    ).to_parquet(entry_dir / "entry_biounit_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "source_mmcif_major_revision": [1],
            "source_mmcif_minor_revision": [0],
        }
    ).to_parquet(entry_dir / "entry_source.parquet", index=False)

    proper = [row for row in ligand_rows if row["proper"]]
    ligand_table = pd.DataFrame(
        {
            "pdb_id": [pdb_id for _ in proper],
            "system_id": [system_id for _ in proper],
            "ligand_rdkit_canonical_smiles": ["C" for _ in proper],
            "ligand_ccd_code": [str(row["ccd"]) for row in proper],
            "ligand_id": [str(row["ligand_id"]) for row in proper],
            "ligand_asym_id": ["L" for _ in proper],
            "ligand_is_3d_score_able": [
                scoreability[str(row["ligand_id"])] for row in proper
            ],
        }
    )
    ligand_root = data_dir / "ligands"
    ligand_root.mkdir(exist_ok=True)
    ligand_table.to_parquet(ligand_root / f"{pdb_id}.parquet", index=False)


def _write_release(data_dir: Path) -> None:
    _write_entry(
        data_dir,
        "1abc",
        ligand_rows=[
            {
                "ligand_id": "1abc__1.L",
                "ccd": "ATP",
                "proper": True,
                "cofactor": True,
            },
            {
                "ligand_id": "1abc__1.Z",
                "ccd": "ZN",
                "proper": False,
                "ion": True,
            },
        ],
        scoreability={"1abc__1.L": True},
    )
    _write_entry(
        data_dir,
        "2def",
        ligand_rows=[{"ligand_id": "2def__1.L", "ccd": "LIG", "proper": True}],
        scoreability={"2def__1.L": False},
        ph=7.4,
    )


def test_plan_shards_and_finalize_real_v3_contract(tmp_path: Path) -> None:
    _write_release(tmp_path)

    plan = plan_collation(tmp_path)
    assert plan["entry_count"] == 2
    assert plan["codes"] == ["ab", "de"]
    assert planned_code_batch(tmp_path, batch_index=0, batch_size=1) == ["ab"]
    assert planned_code_batch(tmp_path, batch_index=1, batch_size=1) == ["de"]
    assert planned_code_batch(tmp_path, batch_index=2, batch_size=1) == []

    first = collate_shard(tmp_path, "ab", memory_limit="1GB")
    cached = collate_shard(tmp_path, "ab", memory_limit="1GB")
    collate_shard(tmp_path, "de", memory_limit="1GB")
    report = finalize_collation(tmp_path, threads=2, memory_limit="1GB")

    assert first == cached
    assert report["row_counts"] == {
        "annotation": 3,
        "entry_chains": 2,
        "entry_biounit_chains": 2,
        "entry_sources": 2,
    }
    annotation = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert not {
        column
        for column in annotation.columns
        if any(marker in column.casefold() for marker in ("ecod", "kinase"))
    }
    first_entry = annotation[annotation["entry_pdb_id"] == "1abc"]
    assert first_entry["biounit_num_ligands"].tolist() == [2, 2]
    assert first_entry["biounit_num_unique_ccd_codes"].tolist() == [2, 2]
    assert first_entry["biounit_num_proper_ligands"].tolist() == [1, 1]
    assert first_entry["system_protein_chains_total_length"].tolist() == [300, 300]
    assert first_entry["system_unique_ccd_codes"].tolist() == ["ATP-ZN", "ATP-ZN"]
    assert first_entry["system_proper_unique_ccd_codes"].tolist() == ["ATP", "ATP"]
    assert first_entry["system_ligand_has_cofactor"].all()
    assert first_entry["system_ligand_has_ion"].all()
    assert first_entry.set_index("ligand_id")["ligand_is_3d_score_able"].to_dict() == {
        "1abc__1.L": True,
        "1abc__1.Z": False,
    }
    marker = json.loads((tmp_path / "index/collation.json").read_text())
    assert marker["status"] == "complete"


def test_shard_rejects_inputs_changed_after_plan(tmp_path: Path) -> None:
    _write_release(tmp_path)
    plan_collation(tmp_path)
    annotation = tmp_path / "raw_entries/ab/1abc.parquet"
    frame = pd.read_parquet(annotation)
    frame["entry_pH"] = 6.0
    frame.to_parquet(annotation, index=False)

    with pytest.raises(RuntimeError, match="changed after planning"):
        collate_shard(tmp_path, "ab", memory_limit="1GB")


def test_finalize_rechecks_inputs_changed_after_sharding(tmp_path: Path) -> None:
    _write_release(tmp_path)
    plan_collation(tmp_path)
    collate_shard(tmp_path, "ab", memory_limit="1GB")
    collate_shard(tmp_path, "de", memory_limit="1GB")
    annotation = tmp_path / "raw_entries/ab/1abc.parquet"
    frame = pd.read_parquet(annotation)
    frame["entry_pH"] = 6.0
    frame.to_parquet(annotation, index=False)

    with pytest.raises(RuntimeError, match="changed after planning"):
        finalize_collation(tmp_path, threads=2, memory_limit="1GB")


def test_final_install_fails_closed_on_partial_replacement(
    tmp_path: Path, monkeypatch
) -> None:
    names = ("annotation", "entry_chains", "entry_biounit_chains", "entry_sources")
    temporary_paths = {name: tmp_path / f"new-{name}" for name in names}
    final_paths = {name: tmp_path / f"final-{name}" for name in names}
    marker = tmp_path / "collation.json"
    for path in temporary_paths.values():
        path.write_text("new")
    for path in final_paths.values():
        path.write_text("old")
    marker.write_text('{"status": "complete"}')
    real_replace = Path.replace

    def fail_on_second_sidecar(path: Path, target: Path) -> Path:
        if target == final_paths["entry_biounit_chains"]:
            raise OSError("simulated replacement failure")
        return real_replace(path, target)

    monkeypatch.setattr(Path, "replace", fail_on_second_sidecar)

    with pytest.raises(OSError, match="simulated replacement failure"):
        collate_module._install_final_tables_fail_closed(
            temporary_paths, final_paths, marker
        )

    assert not marker.exists()
    assert not final_paths["annotation"].exists()


def test_plan_rejects_partial_materialized_entries(tmp_path: Path) -> None:
    _write_release(tmp_path)
    (tmp_path / "ligands/1abc.parquet").unlink()

    with pytest.raises(FileNotFoundError, match="incomplete V3 entry 1abc"):
        plan_collation(tmp_path)


def test_plan_rejects_empty_release(tmp_path: Path) -> None:
    (tmp_path / "raw_entries").mkdir()

    with pytest.raises(ValueError, match="no materialized V3 entries"):
        plan_collation(tmp_path)


def test_plan_supports_bounded_parallel_inventory(tmp_path: Path) -> None:
    _write_release(tmp_path)

    plan = plan_collation(tmp_path, threads=2)

    assert plan["entry_count"] == 2
    with pytest.raises(ValueError, match="planning threads must be positive"):
        plan_collation(tmp_path, threads=0)


def test_final_validation_rejects_all_ion_or_artifact_systems(
    tmp_path: Path,
) -> None:
    _write_entry(
        tmp_path,
        "1abc",
        ligand_rows=[
            {
                "ligand_id": "1abc__1.I",
                "ccd": "ZN",
                "proper": True,
                "ion": True,
            }
        ],
        scoreability={"1abc__1.I": False},
    )

    with pytest.raises(ValueError, match="all_ion_or_artifact_systems=1"):
        run_collation(tmp_path, threads=1, memory_limit="1GB")


@pytest.mark.parametrize(
    ("length", "unresolved", "error_key"),
    [
        (0, 0, "nonpositive_lengths"),
        (300, -1, "negative_unresolved"),
        (300, 301, "unresolved_exceeds_length"),
    ],
)
def test_final_validation_rejects_invalid_chain_sequence_metadata(
    tmp_path: Path,
    length: int,
    unresolved: int,
    error_key: str,
) -> None:
    _write_release(tmp_path)
    chain_path = tmp_path / "raw_entries/ab/1abc/entry_chains.parquet"
    chains = pd.read_parquet(chain_path)
    chains.loc[0, "chain_length"] = length
    chains.loc[0, "chain_num_unresolved_residues"] = unresolved
    chains.to_parquet(chain_path, index=False)

    with pytest.raises(ValueError, match=error_key):
        run_collation(tmp_path, threads=1, memory_limit="1GB")
