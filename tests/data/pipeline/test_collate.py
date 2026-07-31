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
    finalize_repair_marker,
    plan_collation,
    planned_code_batch,
    repair_collation,
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
            "entry_pdb_id": [pdb_id, pdb_id],
            "chain_asym_id": ["A", "B"],
            "chain_auth_id": ["A", "B"],
            "chain_entity_id": ["1", "2"],
            "chain_type": ["polypeptide(L)", "polypeptide(L)"],
            "chain_receptor_type": ["protein", "protein"],
            "chain_length": [300, 200],
            "chain_num_unresolved_residues": [0, 0],
            "chain_is_holo": [True, True],
            "chain_is_ligand_like": [False, False],
            "chain_uniprot_ids": [["P12345"], ["Q12345"]],
        }
    ).to_parquet(entry_dir / "entry_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id, pdb_id],
            "biounit_id": ["1", "1"],
            "chain_instance": ["1.A", "1.B"],
            "chain_asym_id": ["A", "B"],
            "chain_role": ["receptor", "receptor"],
        }
    ).to_parquet(entry_dir / "entry_biounit_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "source_mmcif_major_revision": [1],
            "source_mmcif_minor_revision": [0],
        }
    ).to_parquet(entry_dir / "entry_source.parquet", index=False)
    pd.DataFrame({"entry_pdb_id": [pdb_id], "entry_pH": [ph]}).to_parquet(
        entry_dir / "entry_metadata.parquet", index=False
    )
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "system_id": [f"{pdb_id}__1__1.A--1.B"],
            "system_biounit_id": ["1"],
            "interface_chain_1": ["1.A"],
            "interface_chain_2": ["1.B"],
            "interface_chain_1_residue_numbers": [[1, 2, 3]],
            "interface_chain_1_residue_indices": [[0, 1, 2]],
            "interface_chain_2_residue_numbers": [[4, 5, 6]],
            "interface_chain_2_residue_indices": [[0, 1, 2]],
            "interface_num_contact_residue_pairs": [3],
        }
    ).to_parquet(entry_dir / "interfaces.parquet", index=False)

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
                "ligand_id": "1abc__1__1.L",
                "ccd": "ATP",
                "proper": True,
                "cofactor": True,
            },
            {
                "ligand_id": "1abc__1__1.Z",
                "ccd": "ZN",
                "proper": False,
                "ion": True,
            },
        ],
        scoreability={"1abc__1__1.L": True},
    )
    _write_entry(
        data_dir,
        "2def",
        ligand_rows=[{"ligand_id": "2def__1__1.L", "ccd": "LIG", "proper": True}],
        scoreability={"2def__1__1.L": False},
        ph=7.4,
    )


def _write_interface_only_entry(data_dir: Path, pdb_id: str = "3ghi") -> None:
    _write_entry(
        data_dir,
        pdb_id,
        ligand_rows=[],
        scoreability={},
    )
    (data_dir / "raw_entries" / pdb_id[1:3] / f"{pdb_id}.parquet").unlink()
    (data_dir / "ligands" / f"{pdb_id}.parquet").unlink()


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
        "entry_chains": 4,
        "entry_biounit_chains": 4,
        "entry_metadata": 2,
        "interfaces": 2,
        "entry_sources": 2,
    }
    assert report["interface_count"] == 2
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
        "1abc__1__1.L": True,
        "1abc__1__1.Z": False,
    }
    marker = json.loads((tmp_path / "index/collation.json").read_text())
    assert marker["status"] == "complete"


def test_collation_retains_entries_with_only_protein_interfaces(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    _write_interface_only_entry(tmp_path)

    report = run_collation(tmp_path, memory_limit="1GB")

    assert report["entry_count"] == 3
    assert report["interface_count"] == 3
    annotation = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert "3ghi" not in set(annotation["entry_pdb_id"])
    metadata = pd.read_parquet(tmp_path / "index/entry_metadata.parquet")
    assert set(metadata["entry_pdb_id"]) == {"1abc", "2def", "3ghi"}
    interfaces = pd.read_parquet(tmp_path / "index/interface_annotation_table.parquet")
    assert "3ghi__1__1.A--1.B" in set(interfaces["system_id"])


def test_targeted_repair_preserves_unaffected_release_only_columns(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    installed = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    installed["release_only"] = installed["entry_pdb_id"].map(
        {"1abc": "old", "2def": "keep"}
    )
    installed.to_parquet(tmp_path / "index/annotation_table.parquet", index=False)
    installed.iloc[[0]].to_parquet(
        tmp_path / "index/annotation_table_nonredundant.parquet",
        index=False,
    )
    chain_path = tmp_path / "index/entry_chains.parquet"
    before_chain_stat = chain_path.stat()

    raw = pd.read_parquet(tmp_path / "raw_entries/ab/1abc.parquet")
    raw["entry_pH"] = 6.5
    raw.to_parquet(tmp_path / "raw_entries/ab/1abc.parquet", index=False)

    report = repair_collation(tmp_path, ["1ABC", "1abc"], threads=2, memory_limit="1GB")

    assert report["mode"] == "targeted_repair"
    assert report["status"] == "requires_downstream_repair"
    assert report["repaired_entry_count"] == 1
    assert not (tmp_path / "index/annotation_table_nonredundant.parquet").exists()
    repaired = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert repaired.loc[repaired["entry_pdb_id"].eq("1abc"), "entry_pH"].eq(6.5).all()
    assert (
        repaired.loc[repaired["entry_pdb_id"].eq("1abc"), "release_only"].isna().all()
    )
    assert (
        repaired.loc[repaired["entry_pdb_id"].eq("2def"), "release_only"]
        .eq("keep")
        .all()
    )
    repaired_chains = pd.read_parquet(tmp_path / "index/entry_chains.parquet")
    lengths = (
        repaired_chains[repaired_chains["chain_asym_id"].eq("A")]
        .set_index("entry_pdb_id")["chain_length"]
        .to_dict()
    )
    assert lengths == {"1abc": 300, "2def": 300}
    after_chain_stat = chain_path.stat()
    assert (
        after_chain_stat.st_ino,
        after_chain_stat.st_size,
        after_chain_stat.st_mtime_ns,
    ) == (
        before_chain_stat.st_ino,
        before_chain_stat.st_size,
        before_chain_stat.st_mtime_ns,
    )

    with pytest.raises(FileNotFoundError, match="nonredundant"):
        finalize_repair_marker(tmp_path)
    repaired.iloc[[0]].to_parquet(
        tmp_path / "index/annotation_table_nonredundant.parquet",
        index=False,
    )
    finalized_report = finalize_repair_marker(tmp_path)
    assert finalized_report is not None
    assert finalized_report["status"] == "complete"
    assert finalized_report["downstream_repair_complete"] is True


def test_targeted_repair_rejects_chain_metadata_changes(tmp_path: Path) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    chains = pd.read_parquet(tmp_path / "raw_entries/ab/1abc/entry_chains.parquet")
    chains["chain_length"] = 301
    chains.to_parquet(
        tmp_path / "raw_entries/ab/1abc/entry_chains.parquet", index=False
    )

    with pytest.raises(ValueError, match="changed entry_chains"):
        repair_collation(tmp_path, ["1abc"], memory_limit="1GB")


def test_shard_rejects_inputs_changed_after_plan(tmp_path: Path) -> None:
    _write_release(tmp_path)
    plan_collation(tmp_path)
    annotation = tmp_path / "raw_entries/ab/1abc.parquet"
    frame = pd.read_parquet(annotation)
    frame["entry_pH"] = 6.0
    frame.to_parquet(annotation, index=False)

    with pytest.raises(RuntimeError, match="changed after planning"):
        collate_shard(tmp_path, "ab", memory_limit="1GB")


def test_shard_rejects_optional_ligand_outputs_appearing_after_plan(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    _write_interface_only_entry(tmp_path)
    plan_collation(tmp_path)
    pd.DataFrame({"entry_pdb_id": ["3ghi"]}).to_parquet(
        tmp_path / "raw_entries/gh/3ghi.parquet", index=False
    )
    pd.DataFrame({"ligand_id": ["3ghi__1__1.L"]}).to_parquet(
        tmp_path / "ligands/3ghi.parquet", index=False
    )

    with pytest.raises(RuntimeError, match="appeared after planning"):
        collate_shard(tmp_path, "gh", memory_limit="1GB")


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
    names = (
        "annotation",
        "entry_chains",
        "entry_biounit_chains",
        "entry_metadata",
        "interfaces",
        "entry_sources",
    )
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


def test_shard_cli_selects_codes_from_pdb_manifest(tmp_path: Path) -> None:
    from plinder.data.pipeline.collate import build_parser

    manifest = tmp_path / "affected.txt"
    manifest.write_text("1abc\n2abd\n3xyz\n")
    args = build_parser().parse_args(
        ["shard", str(tmp_path), "--pdb-manifest", str(manifest)]
    )

    assert args.pdb_manifest == manifest


def test_final_validation_rejects_all_ion_or_artifact_systems(
    tmp_path: Path,
) -> None:
    _write_entry(
        tmp_path,
        "1abc",
        ligand_rows=[
            {
                "ligand_id": "1abc__1__1.I",
                "ccd": "ZN",
                "proper": True,
                "ion": True,
            }
        ],
        scoreability={"1abc__1__1.I": False},
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


def test_final_validation_rejects_ligands_from_another_biounit(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    annotation_path = tmp_path / "raw_entries/ab/1abc.parquet"
    annotation = pd.read_parquet(annotation_path)
    annotation.loc[0, "ligand_id"] = "1abc__2__1.L"
    annotation.to_parquet(annotation_path, index=False)
    ligand_path = tmp_path / "ligands/1abc.parquet"
    ligands = pd.read_parquet(ligand_path)
    ligands.loc[ligands["ligand_id"].eq("1abc__1__1.L"), "ligand_id"] = "1abc__2__1.L"
    ligands.to_parquet(ligand_path, index=False)

    with pytest.raises(ValueError, match="mismatched_biounits=1"):
        run_collation(tmp_path, threads=1, memory_limit="1GB")
