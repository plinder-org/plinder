# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from __future__ import annotations

import gzip
import json
from pathlib import Path
from shutil import copy2, copytree

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from plinder.data.annotations.interface_utils import (
    INTERFACE_ANNOTATION_SCHEMA,
    MIN_INTERFACE_RESIDUES_METADATA_KEY,
)
from plinder.data.pipeline import collate as collate_module
from plinder.data.pipeline import update_entries, updates
from plinder.data.pipeline.collate import (
    COLLATION_VERSION,
    collate_shard,
    finalize_collation,
    finalize_collation_plan,
    finalize_repair_marker,
    inventory_collation_codes,
    plan_collation,
    planned_code_batch,
    planned_inventory_code_batch,
    repair_collation,
    run_collation,
    start_collation_plan,
)


def _write_entry(
    data_dir: Path,
    pdb_id: str,
    *,
    ligand_rows: list[dict[str, object]],
    comparability: dict[str, bool],
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
            "system_id_no_biounit": f"{pdb_id}__1.A__1.L",
            "system_id_legacy": f"{pdb_id}__1__1.A__1.L",
            "system_type": "holo",
            "system_protein_chains_length": [100, 200],
            "system_ligand_chains": ["1.L"],
            "system_ligand_chains_asym_id": ["1.L"],
            "system_ligand_validation_average_rsr": 0.1,
            "system_pocket_validation_average_rsr": 0.2,
            "ligand_id": ligand["ligand_id"],
            "ligand_id_legacy": ligand["ligand_id"],
            "ligand_smiles": "C",
            "ligand_rdkit_canonical_smiles": "C",
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
        for prefix in (
            "system_protein_chains_",
            "system_ligand_chains_",
            "ligand_protein_chains_",
            "ligand_neighboring_ligand_chains_",
            "ligand_interacting_ligand_chains_",
        ):
            row[f"{prefix}auth_id"] = ["A"]
            row[f"{prefix}entity_id"] = ["1"]
            row[f"{prefix}num_unresolved_residues"] = [0]
            row[f"{prefix}validation_average_rsr"] = [0.1]
            if f"{prefix}length" not in row:
                row[f"{prefix}length"] = [100]
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
            "chain_sequence": ["A" * 300, "A" * 200],
            "chain_sequence_noncanonical": ["A" * 300, "A" * 200],
            "chain_modified_residues": [[], []],
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
            "chain_num_contacting_ions": [0, 0],
            "chain_num_contacting_artifacts": [0, 0],
            "chain_num_contacting_other_ligands": [0, 0],
        }
    ).to_parquet(entry_dir / "entry_biounit_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "source_mmcif_major_revision": [1],
            "source_mmcif_minor_revision": [0],
        }
    ).to_parquet(entry_dir / "entry_source.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "entry_pH": [ph],
            "entry_pH_min": [ph],
            "entry_pH_max": [ph],
            "entry_has_ligand_of_interest": [None],
        }
    ).to_parquet(entry_dir / "entry_metadata.parquet", index=False)
    interface_row = dict.fromkeys(INTERFACE_ANNOTATION_SCHEMA.names)
    interface_row.update(
        {
            "entry_pdb_id": pdb_id,
            "system_id": f"{pdb_id}__1__1.A--1.B",
            "system_biounit_id": "1",
            "interface_chain_1": "1.A",
            "interface_chain_2": "1.B",
            "interface_chain_1_residue_numbers": [1, 2, 3, 4, 5, 6, 7],
            "interface_chain_1_residue_indices": [0, 1, 2, 3, 4, 5, 6],
            "interface_chain_2_residue_numbers": [11, 12, 13, 14, 15, 16, 17],
            "interface_chain_2_residue_indices": [0, 1, 2, 3, 4, 5, 6],
            "interface_num_contact_residue_pairs": 7,
            "prodigy_is_annotated": False,
        }
    )
    pq.write_table(
        pa.Table.from_pylist([interface_row], schema=INTERFACE_ANNOTATION_SCHEMA),
        entry_dir / "interfaces.parquet",
    )
    _set_interface_threshold(entry_dir / "interfaces.parquet", 7)

    proper = [row for row in ligand_rows if row["proper"]]
    ligand_table = pd.DataFrame(
        {
            "pdb_id": [pdb_id for _ in proper],
            "system_id": [system_id for _ in proper],
            "ligand_rdkit_canonical_smiles": ["C" for _ in proper],
            "ligand_ccd_code": [str(row["ccd"]) for row in proper],
            "ligand_id": [str(row["ligand_id"]) for row in proper],
            "ligand_asym_id": ["L" for _ in proper],
            "ligand_is_shape_comparable": [
                comparability[str(row["ligand_id"])] for row in proper
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
        comparability={"1abc__1__1.L": True},
    )
    _write_entry(
        data_dir,
        "2def",
        ligand_rows=[{"ligand_id": "2def__1__1.L", "ccd": "LIG", "proper": True}],
        comparability={"2def__1__1.L": False},
        ph=7.4,
    )


@pytest.fixture
def entry_update_case(tmp_path, monkeypatch):
    base = tmp_path / "base"
    _write_release(base)
    _write_entry(
        base,
        "3ghi",
        ligand_rows=[{"ligand_id": "3ghi__1__1.L", "ccd": "ATP", "proper": True}],
        comparability={"3ghi__1__1.L": True},
    )
    run_collation(base, memory_limit="1GB")
    for relative in update_entries.REQUIRED_REFERENCE_FILES:
        path = base / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")
    root = tmp_path / "nextgen"
    (root / "holdings").mkdir(parents=True)
    revisions = {"1abc": (2, 0), "3ghi": (1, 0), "4jkl": (1, 0)}
    for name, values in {
        "current_file_holdings": {key: {"mmcif": [key]} for key in revisions},
        "released_structures_last_modified_dates": {
            key: "2026-07-01" for key in revisions
        },
    }.items():
        with gzip.open(root / f"holdings/{name}.json.gz", "wt") as handle:
            json.dump(values, handle)
    monkeypatch.setattr(
        updates,
        "_read_revision",
        lambda path: revisions[path.parent.name.removeprefix("pdb_0000")],
    )
    obsolete = tmp_path / "obsolete.dat"
    obsolete.write_text("OBSLTE 01-JUL-26 2DEF 4JKL\n")
    plan_dir = tmp_path / "plan"
    updates.plan_update(
        base, nextgen_root=root, obsolete_path=obsolete, output_dir=plan_dir
    )
    calls = []

    def ingest(**kwargs):
        calls.append(kwargs["pdb_ids"])
        destination = kwargs["output_root"]
        for pdb_id in kwargs["pdb_ids"]:
            entry = destination / "raw_entries" / pdb_id[1:3] / pdb_id
            if entry.exists():
                continue
            _write_entry(
                destination,
                pdb_id,
                ligand_rows=[
                    {"ligand_id": f"{pdb_id}__1__1.L", "ccd": "ADP", "proper": True}
                ],
                comparability={f"{pdb_id}__1__1.L": True},
            )
            pd.DataFrame(
                {
                    "entry_pdb_id": [pdb_id],
                    "source_mmcif_major_revision": [revisions[pdb_id][0]],
                    "source_mmcif_minor_revision": [revisions[pdb_id][1]],
                }
            ).to_parquet(entry / "entry_source.parquet", index=False)
            chains = pd.read_parquet(entry / "entry_chains.parquet")
            chains["chain_sequence"] = chains.chain_length.map(
                lambda length: "C" * length
            )
            chains["chain_sequence_noncanonical"] = chains.chain_sequence
            chains.to_parquet(entry / "entry_chains.parquet", index=False)
            (entry / "ligand_files").mkdir()
            (entry / "ligand_files/L.sdf").write_text("new ligand pose")
        return destination / "metrics.json", False

    monkeypatch.setattr(update_entries, "ingest_pdb_batch", ingest)
    args = {
        "plan_dir": plan_dir,
        "output_dir": tmp_path / "updated",
        "validation_root": tmp_path,
        "annotation_cfg": {},
        "entry_cfg": {},
        "interface_cfg": {"min_interface_residues": 7, "annotate_prodigy": False},
        "threads": 1,
        "memory_limit": "1GB",
    }
    return base, args, calls


def _index_bytes(root):
    return {
        path.name: path.read_bytes()
        for path in (root / "index").glob("*")
        if path.is_file()
    }


def test_weekly_entry_update_matches_full_collation(entry_update_case, tmp_path):
    base, args, calls = entry_update_case
    before = _index_bytes(base)
    report = update_entries.apply_entry_update(**args)
    assert report["status"] == "requires_downstream_repair"
    assert report["ingested_pdb_ids"] == ["1abc", "4jkl"]
    assert report["obsolete_pdb_ids"] == ["2def"]
    assert calls == [["1abc", "4jkl"]]
    assert _index_bytes(base) == before
    output = args["output_dir"]
    assert not (output / "scores").exists()
    assert not (output / "ligand_sampling").exists()
    assert (output / ".incoming/raw_entries/ab/1abc/ligand_files/L.sdf").is_file()
    full = tmp_path / "full"
    copytree(output / ".incoming/raw_entries", full / "raw_entries")
    copytree(output / ".incoming/ligands", full / "ligands")
    copytree(base / "raw_entries/gh", full / "raw_entries/gh")
    copy2(base / "ligands/3ghi.parquet", full / "ligands/3ghi.parquet")
    run_collation(full, memory_limit="1GB")
    for filename, _ in update_entries.ENTRY_TABLES.values():
        pd.testing.assert_frame_equal(
            pd.read_parquet(output / "index" / filename),
            pd.read_parquet(full / "index" / filename),
        )
    assert update_entries.apply_entry_update(**args) == report
    assert calls == [["1abc", "4jkl"]]
    assert (
        json.loads((output / "index/collation.json").read_text())["status"]
        == "requires_downstream_repair"
    )


@pytest.mark.parametrize("with_interfaces", [True, False])
def test_weekly_update_without_incoming_ligands(
    entry_update_case, tmp_path, monkeypatch, with_interfaces
):
    base, args, _ = entry_update_case
    before = _index_bytes(base)

    def ingest(**kwargs):
        root = kwargs["output_root"]
        for pdb_id in kwargs["pdb_ids"]:
            directory = root / "raw_entries" / pdb_id[1:3] / pdb_id
            if directory.exists():
                continue
            writer = (
                _write_interface_only_entry
                if with_interfaces
                else _write_sidecar_only_entry
            )
            writer(root, pdb_id)
            source = pd.read_parquet(directory / "entry_source.parquet")
            source["source_mmcif_major_revision"] = 2 if pdb_id == "1abc" else 1
            source.to_parquet(directory / "entry_source.parquet", index=False)
            chains = pd.read_parquet(directory / "entry_chains.parquet")
            chains["chain_is_holo"] = False
            chains.to_parquet(directory / "entry_chains.parquet", index=False)
        return root / "metrics.json", False

    monkeypatch.setattr(update_entries, "ingest_pdb_batch", ingest)
    report = update_entries.apply_entry_update(**args)
    output = args["output_dir"]
    assert report["status"] == "requires_downstream_repair"
    assert _index_bytes(base) == before
    for name in ("annotation", "system_validation"):
        filename = update_entries.ENTRY_TABLES[name][0]
        incoming = output / ".incoming/index" / filename
        assert pq.ParquetFile(incoming).metadata.num_rows == 0
        assert pq.read_schema(incoming).equals(
            pq.read_schema(base / "index" / filename)
        )
    assert set(
        pd.read_parquet(output / "index/annotation_table.parquet").entry_pdb_id
    ) == {"3ghi"}

    full = tmp_path / "full"
    copytree(base / "raw_entries/gh", full / "raw_entries/gh")
    (full / "ligands").mkdir()
    copy2(base / "ligands/3ghi.parquet", full / "ligands/3ghi.parquet")
    ingest(output_root=full, pdb_ids=["1abc", "4jkl"])
    run_collation(full, memory_limit="1GB")
    for filename, _ in update_entries.ENTRY_TABLES.values():
        pd.testing.assert_frame_equal(
            pd.read_parquet(output / "index" / filename),
            pd.read_parquet(full / "index" / filename),
        )
    assert update_entries.apply_entry_update(**args) == report


def test_collation_keeps_unknown_ligand_interest_boolean(tmp_path):
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    table = pq.read_table(tmp_path / "index/entry_metadata.parquet")
    flag = table["entry_has_ligand_of_interest"]
    assert flag.type == pa.bool_()
    assert flag.null_count == table.num_rows == 2


def test_weekly_removal_only_needs_no_raw_entries(entry_update_case, tmp_path):
    base, args, calls = entry_update_case
    original_plan = json.loads((args["plan_dir"] / "plan.json").read_text())
    plan_dir = tmp_path / "removal-plan"
    updates.plan_update(
        base,
        output_dir=plan_dir,
        nextgen_root=Path(original_plan["nextgen_root"]),
        obsolete_path=Path(original_plan["inputs"]["obsolete"]["path"]),
        pdb_ids=["2def"],
    )
    args["plan_dir"] = plan_dir
    # Obsolete removals operate on the existing tables, not the old raw files.
    from shutil import rmtree

    rmtree(base / "raw_entries")
    report = update_entries.apply_entry_update(**args)
    assert report["ingested_pdb_ids"] == []
    assert calls == []
    for filename, _ in update_entries.ENTRY_TABLES.values():
        frame = pd.read_parquet(args["output_dir"] / "index" / filename)
        ids = (
            frame.entry_pdb_id
            if "entry_pdb_id" in frame
            else frame.system_id.str.split("__").str[0]
        )
        assert set(ids) == {"1abc", "3ghi"}


def test_weekly_update_failure_and_resume_preserve_base(entry_update_case, monkeypatch):
    base, args, calls = entry_update_case
    before = _index_bytes(base)
    original = collate_module._validate_final_tables

    def fail_staging(paths, **kwargs):
        if paths["annotation"].parent.name == ".entry_tables":
            raise ValueError("injected validation failure")
        return original(paths, **kwargs)

    monkeypatch.setattr(collate_module, "_validate_final_tables", fail_staging)
    with pytest.raises(ValueError, match="injected"):
        update_entries.apply_entry_update(**args)
    assert _index_bytes(base) == before
    assert not (args["output_dir"] / "index").exists()
    monkeypatch.setattr(collate_module, "_validate_final_tables", original)
    report = update_entries.apply_entry_update(**args)
    assert report["status"] == "requires_downstream_repair"
    assert _index_bytes(base) == before


@pytest.mark.parametrize(
    "change", ["base", "report", "settings", "unowned", "threshold"]
)
def test_weekly_update_rejects_changed_inputs(entry_update_case, change):
    base, args, _ = entry_update_case
    if change == "base":
        path = base / "index/entry_sources.parquet"
        frame = pd.read_parquet(path)
        frame["source_mmcif_minor_revision"] = 9
        frame.to_parquet(path, index=False)
    elif change == "report":
        path = args["plan_dir"] / "entries.parquet"
        frame = pd.read_parquet(path)
        frame["action"] = "obsolete"
        frame.to_parquet(path, index=False)
    elif change == "settings":
        update_entries.apply_entry_update(**args)
        args["annotation_cfg"] = {"min_polymer_size": 30}
    elif change == "threshold":
        args["interface_cfg"] = {"min_interface_residues": 8}
    else:
        args["output_dir"].mkdir()
        (args["output_dir"] / "user-file.txt").write_text("keep")
    before = _index_bytes(base)
    with pytest.raises(ValueError):
        update_entries.apply_entry_update(**args)
    assert _index_bytes(base) == before


def test_weekly_update_rejects_wrong_processed_revision(entry_update_case, monkeypatch):
    base, args, _ = entry_update_case
    original = update_entries.ingest_pdb_batch

    def wrong_revision(**kwargs):
        result = original(**kwargs)
        path = kwargs["output_root"] / "raw_entries/ab/1abc/entry_source.parquet"
        frame = pd.read_parquet(path)
        frame["source_mmcif_major_revision"] = 3
        frame.to_parquet(path, index=False)
        return result

    monkeypatch.setattr(update_entries, "ingest_pdb_batch", wrong_revision)
    before = _index_bytes(base)
    with pytest.raises(ValueError, match="revisions differ"):
        update_entries.apply_entry_update(**args)
    assert before == _index_bytes(base)
    assert not (args["output_dir"] / "index").exists()


def test_weekly_update_install_failure_can_resume(entry_update_case, monkeypatch):
    base, args, _ = entry_update_case
    before = _index_bytes(base)
    replace = Path.replace

    def fail_metadata(path, target):
        if (
            path.parent.name == ".entry_tables"
            and path.name == "entry_metadata.parquet"
        ):
            raise OSError("injected installation failure")
        return replace(path, target)

    monkeypatch.setattr(Path, "replace", fail_metadata)
    with pytest.raises(OSError, match="injected"):
        update_entries.apply_entry_update(**args)
    assert before == _index_bytes(base)
    assert not (args["output_dir"] / "index/annotation_table.parquet").exists()
    assert not (args["output_dir"] / "index/collation.json").exists()
    monkeypatch.setattr(Path, "replace", replace)
    assert (
        update_entries.apply_entry_update(**args)["status"]
        == "requires_downstream_repair"
    )
    assert before == _index_bytes(base)


def test_weekly_update_entry_failure_keeps_workspace_resumable(
    entry_update_case, monkeypatch
):
    base, args, _ = entry_update_case
    original = update_entries.ingest_pdb_batch
    before = _index_bytes(base)

    def failed(**kwargs):
        path, _ = original(**kwargs)
        return path, True

    monkeypatch.setattr(update_entries, "ingest_pdb_batch", failed)
    with pytest.raises(RuntimeError, match="entry update failed"):
        update_entries.apply_entry_update(**args)
    assert before == _index_bytes(base)
    assert not (args["output_dir"] / "index").exists()
    monkeypatch.setattr(update_entries, "ingest_pdb_batch", original)
    assert (
        update_entries.apply_entry_update(**args)["status"]
        == "requires_downstream_repair"
    )


def test_weekly_update_does_not_reuse_missing_marker(entry_update_case):
    _, args, _ = entry_update_case
    update_entries.apply_entry_update(**args)
    (args["output_dir"] / "index/collation.json").unlink()
    with pytest.raises(FileNotFoundError):
        update_entries.apply_entry_update(**args)


def test_weekly_update_preserves_nested_residue_mappings(
    entry_update_case, monkeypatch
):
    base, args, _ = entry_update_case
    path = base / "index/annotation_table.parquet"
    frame = pd.read_parquet(path)
    frame["ligand__members"] = [{"1.L": [3, 4]} for _ in range(len(frame))]
    frame.to_parquet(path, index=False)
    original = update_entries.ingest_pdb_batch

    def new_chains(**kwargs):
        result = original(**kwargs)
        for pdb_id in kwargs["pdb_ids"]:
            path = (
                kwargs["output_root"]
                / "raw_entries"
                / pdb_id[1:3]
                / f"{pdb_id}.parquet"
            )
            rows = pd.read_parquet(path)
            rows["ligand__members"] = [{"2.NEW": [5, 6]} for _ in range(len(rows))]
            rows.to_parquet(path, index=False)
        return result

    monkeypatch.setattr(update_entries, "ingest_pdb_batch", new_chains)
    before = _index_bytes(base)
    update_entries.apply_entry_update(**args)
    rows = pd.read_parquet(args["output_dir"] / "index/annotation_table.parquet")
    for row in rows.itertuples(index=False, name=None):
        mapping = dict(zip(rows.columns, row))
        residues = mapping["ligand__members"]
        if mapping["entry_pdb_id"] == "3ghi":
            assert list(residues["1.L"]) == [3, 4]
            assert residues["2.NEW"] is None
        else:
            assert list(residues["2.NEW"]) == [5, 6]
            assert residues["1.L"] is None
    assert _index_bytes(base) == before


def test_collation_preserves_failure_diagnostics_and_unknown_older_rows(tmp_path):
    from plinder.core import PlinderRelease, query_table

    # Mix empty lists, recorded failures, and older rows without the columns,
    # both within a shard and across shards.
    for pdb_id, failures in [("1abc", []), ("2abd", ["2"]), ("3def", None)]:
        _write_entry(
            tmp_path,
            pdb_id,
            ligand_rows=[
                {"ligand_id": f"{pdb_id}__1__1.L", "ccd": "ATP", "proper": True}
            ],
            comparability={f"{pdb_id}__1__1.L": True},
        )
        if failures is None:
            continue
        raw_root = tmp_path / "raw_entries" / pdb_id[1:3]
        metadata_path = raw_root / pdb_id / "entry_metadata.parquet"
        metadata = pd.read_parquet(metadata_path)
        metadata["entry_failed_assembly_ids"] = [failures]
        metadata.to_parquet(metadata_path, index=False)
        annotation_path = raw_root / f"{pdb_id}.parquet"
        annotation = pd.read_parquet(annotation_path)
        annotation["ligand_failed_interaction_types"] = [
            ["water_bridge"] if failures else []
        ]
        annotation.to_parquet(annotation_path, index=False)

    run_collation(tmp_path, threads=1, memory_limit="1GB")
    result = query_table(
        "annotation",
        columns=[
            "entry_pdb_id",
            "entry_failed_assembly_ids",
            "ligand_failed_interaction_types",
        ],
        release=PlinderRelease(data_dir=tmp_path),
    ).set_index("entry_pdb_id")
    assert result.loc["1abc", "entry_failed_assembly_ids"].tolist() == []
    assert result.loc["1abc", "ligand_failed_interaction_types"].tolist() == []
    assert result.loc["2abd", "entry_failed_assembly_ids"].tolist() == ["2"]
    assert result.loc["2abd", "ligand_failed_interaction_types"].tolist() == [
        "water_bridge"
    ]
    assert pd.isna(result.loc["3def", "entry_failed_assembly_ids"])
    assert pd.isna(result.loc["3def", "ligand_failed_interaction_types"])
    assert (
        "entry_failed_assembly_ids"
        not in pq.read_schema(tmp_path / "index/annotation_table.parquet").names
    )


def _write_interface_only_entry(data_dir: Path, pdb_id: str = "3ghi") -> None:
    _write_entry(
        data_dir,
        pdb_id,
        ligand_rows=[],
        comparability={},
    )
    (data_dir / "raw_entries" / pdb_id[1:3] / f"{pdb_id}.parquet").unlink()
    (data_dir / "ligands" / f"{pdb_id}.parquet").unlink()
    metrics = data_dir / "metrics" / pdb_id[1:3] / f"ingest-one-{pdb_id}.json"
    metrics.parent.mkdir(parents=True, exist_ok=True)
    metrics.write_text(
        json.dumps(
            {
                "status": "complete",
                "counts": {"annotation_rows": 0, "interface_rows": 1},
                "outputs": {
                    "entry_directory": str(
                        data_dir / "raw_entries" / pdb_id[1:3] / pdb_id
                    ),
                    "entry_parquet": None,
                    "ligand_parquet": None,
                },
            }
        )
    )


def _write_sidecar_only_entry(data_dir: Path, pdb_id: str = "4jkl") -> None:
    _write_interface_only_entry(data_dir, pdb_id)
    interface_path = (
        data_dir / "raw_entries" / pdb_id[1:3] / pdb_id / "interfaces.parquet"
    )
    pq.write_table(
        pa.Table.from_pylist([], schema=INTERFACE_ANNOTATION_SCHEMA), interface_path
    )
    _set_interface_threshold(interface_path, 7)
    metrics = data_dir / "metrics" / pdb_id[1:3] / f"ingest-one-{pdb_id}.json"
    payload = json.loads(metrics.read_text())
    payload["counts"]["interface_rows"] = 0
    metrics.write_text(json.dumps(payload))


def _set_interface_threshold(path: Path, threshold: int) -> None:
    table = pq.read_table(path)
    metadata = dict(table.schema.metadata or {})
    metadata[MIN_INTERFACE_RESIDUES_METADATA_KEY] = str(threshold).encode()
    pq.write_table(table.replace_schema_metadata(metadata), path)


def test_plan_shards_and_finalize_release_contract(tmp_path: Path) -> None:
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
        "system_validation": 2,
        "entry_chains": 4,
        "entry_biounit_chains": 4,
        "entry_metadata": 2,
        "interfaces": 2,
        "entry_sources": 2,
    }
    assert report["interface_count"] == 2
    annotation = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert [column for column in annotation.columns if column.startswith("entry_")] == [
        "entry_pdb_id"
    ]
    assert not {
        column
        for column in annotation.columns
        if any(marker in column.casefold() for marker in ("ecod", "kinase"))
    }
    retired_annotation_columns = {
        "system_id_no_biounit",
        "system_ligand_chains",
        "ligand_rdkit_canonical_smiles",
    }
    assert retired_annotation_columns.isdisjoint(annotation.columns)
    assert "ligand_smiles" in annotation.columns
    assert not any(
        column.startswith(
            (
                "system_ligand_validation_",
                "system_pocket_validation_",
            )
        )
        for column in annotation.columns
    )
    for prefix in (
        "system_protein_chains_",
        "system_ligand_chains_",
        "ligand_protein_chains_",
        "ligand_neighboring_ligand_chains_",
        "ligand_interacting_ligand_chains_",
    ):
        for suffix in (
            "auth_id",
            "entity_id",
            "length",
            "num_unresolved_residues",
        ):
            assert f"{prefix}{suffix}" not in annotation.columns
    assert not any("_chains_validation_" in column for column in annotation.columns)
    first_entry = annotation[annotation["entry_pdb_id"] == "1abc"]
    assert first_entry["biounit_num_ligands"].tolist() == [2, 2]
    assert first_entry["biounit_num_unique_ccd_codes"].tolist() == [2, 2]
    assert first_entry["biounit_num_proper_ligands"].tolist() == [1, 1]
    assert first_entry["system_protein_chains_total_length"].tolist() == [300, 300]
    assert first_entry["system_unique_ccd_codes"].tolist() == ["ATP-ZN", "ATP-ZN"]
    assert first_entry["system_proper_unique_ccd_codes"].tolist() == ["ATP", "ATP"]
    assert first_entry["system_ligand_has_cofactor"].all()
    assert first_entry["system_ligand_has_ion"].all()
    assert first_entry.set_index("ligand_id")[
        "ligand_is_shape_comparable"
    ].to_dict() == {
        "1abc__1__1.L": True,
        "1abc__1__1.Z": False,
    }
    metadata = pd.read_parquet(tmp_path / "index/entry_metadata.parquet")
    assert metadata.set_index("entry_pdb_id").loc["2def", "entry_pH"] == 7.4
    system_validation = pd.read_parquet(
        tmp_path / "index/system_validation.parquet"
    ).set_index("system_id")
    assert system_validation.loc[
        "1abc__1__1.A__1.L", "system_ligand_validation_average_rsr"
    ] == pytest.approx(0.1)
    assert system_validation.loc[
        "1abc__1__1.A__1.L", "system_pocket_validation_average_rsr"
    ] == pytest.approx(0.2)
    marker = json.loads((tmp_path / "index/collation.json").read_text())
    assert marker["status"] == "complete"
    assert marker["interface_min_residues"] == 7


def test_collate_shard_does_not_reuse_older_format(tmp_path: Path) -> None:
    _write_release(tmp_path)
    plan_collation(tmp_path)
    collate_shard(tmp_path, "ab", memory_limit="1GB")
    paths = collate_module._shard_paths(tmp_path, "ab")
    stale = pd.read_parquet(paths["entry_biounit_chains"]).drop(
        columns=[
            "chain_num_contacting_ions",
            "chain_num_contacting_artifacts",
            "chain_num_contacting_other_ligands",
        ]
    )
    stale.to_parquet(paths["entry_biounit_chains"], index=False)
    metrics = json.loads(paths["metrics"].read_text())
    metrics["version"] = COLLATION_VERSION - 1
    paths["metrics"].write_text(json.dumps(metrics))

    refreshed = collate_shard(tmp_path, "ab", memory_limit="1GB")

    assert refreshed["version"] == COLLATION_VERSION
    assert {
        "chain_num_contacting_ions",
        "chain_num_contacting_artifacts",
        "chain_num_contacting_other_ligands",
    }.issubset(pq.read_schema(paths["entry_biounit_chains"]).names)


def test_collation_rejects_unknown_biounit_contact_counts(tmp_path: Path) -> None:
    _write_release(tmp_path)
    path = tmp_path / "raw_entries/ab/1abc/entry_biounit_chains.parquet"
    frame = pd.read_parquet(path)
    frame.loc[0, "chain_num_contacting_ions"] = None
    frame.to_parquet(path, index=False)

    with pytest.raises(
        ValueError,
        match="invalid biological-assembly ligand contact counts",
    ):
        run_collation(tmp_path, memory_limit="1GB")


def test_distributed_plan_requires_and_merges_every_code_inventory(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)

    build = start_collation_plan(tmp_path)
    assert build["code_count"] == 2
    assert planned_inventory_code_batch(tmp_path, batch_index=0, batch_size=1) == ["ab"]
    assert planned_inventory_code_batch(tmp_path, batch_index=1, batch_size=1) == ["de"]

    inventory_collation_codes(tmp_path, ["ab"], threads=2)
    with pytest.raises(FileNotFoundError, match="incomplete for code de"):
        finalize_collation_plan(tmp_path)

    inventory_collation_codes(tmp_path, ["de"], threads=2)
    plan = finalize_collation_plan(tmp_path)

    assert plan["entry_count"] == 2
    assert plan["code_entry_counts"] == {"ab": 1, "de": 1}
    assert pq.read_table(plan["manifest"]).column("pdb_id").to_pylist() == [
        "1abc",
        "2def",
    ]


def test_interface_only_collation_preserves_installed_ligand_annotation(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    annotation_path = tmp_path / "index/annotation_table.parquet"
    annotation = pd.read_parquet(annotation_path)
    annotation["release_only"] = "preserve"
    annotation.to_parquet(annotation_path, index=False)
    annotation_bytes = annotation_path.read_bytes()
    system_validation_path = tmp_path / "index/system_validation.parquet"
    system_validation_bytes = system_validation_path.read_bytes()
    legacy_chain_path = tmp_path / "raw_entries/ab/1abc/entry_chains.parquet"
    legacy_chains = pd.read_parquet(legacy_chain_path).drop(
        columns="chain_is_ligand_like"
    )
    legacy_chains.to_parquet(legacy_chain_path, index=False)
    biounit_path = tmp_path / "raw_entries/ab/1abc/entry_biounit_chains.parquet"
    biounits = pd.read_parquet(biounit_path)
    biounits.loc[biounits["chain_asym_id"] == "B", "chain_role"] = "ligand"
    biounits.to_parquet(biounit_path, index=False)

    build = start_collation_plan(
        tmp_path,
        include_ligand_annotations=False,
    )
    inventory_collation_codes(tmp_path, build["codes"], threads=2)
    plan = finalize_collation_plan(tmp_path)
    for code in plan["codes"]:
        collate_shard(tmp_path, code, memory_limit="1GB")
    report = finalize_collation(tmp_path, memory_limit="1GB")

    assert report["include_ligand_annotations"] is False
    assert annotation_path.read_bytes() == annotation_bytes
    assert system_validation_path.read_bytes() == system_validation_bytes
    assert report["interface_count"] == 2
    chains = pd.read_parquet(tmp_path / "index/entry_chains.parquet")
    assert chains.set_index(["entry_pdb_id", "chain_asym_id"]).loc[
        ("1abc", "B"), "chain_is_ligand_like"
    ]


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


@pytest.mark.parametrize("mode", ["all", "ligands"])
def test_collation_retains_entries_with_only_chain_sidecars(
    tmp_path: Path, mode: str
) -> None:
    _write_release(tmp_path)
    _write_sidecar_only_entry(tmp_path)
    marker = tmp_path / "metrics/jk/ingest-one-4jkl.json"
    payload = json.loads(marker.read_text())
    payload["mode"] = mode
    if mode == "ligands":
        payload["status"] = "skipped_no_ligands"
        payload["counts"]["entry_chain_rows"] = 2
    marker.write_text(json.dumps(payload))

    report = run_collation(tmp_path, memory_limit="1GB")

    assert report["entry_count"] == 3
    assert report["interface_count"] == 2
    annotation = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert "4jkl" not in set(annotation["entry_pdb_id"])
    metadata = pd.read_parquet(tmp_path / "index/entry_metadata.parquet")
    assert set(metadata["entry_pdb_id"]) == {"1abc", "2def", "4jkl"}
    chains = pd.read_parquet(tmp_path / "index/entry_chains.parquet")
    assert "4jkl" in set(chains["entry_pdb_id"])

    # Targeted repairs must accept the same sidecar-only entries as a full build.
    repaired = repair_collation(tmp_path, ["4jkl"], memory_limit="1GB")
    assert repaired["status"] == "requires_downstream_repair"
    pd.testing.assert_frame_equal(
        pd.read_parquet(tmp_path / "index/entry_chains.parquet"), chains
    )


def test_collation_rejects_interface_sidecars_from_failed_ingest(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    _write_interface_only_entry(tmp_path)
    metrics = tmp_path / "metrics/gh/ingest-one-3ghi.json"
    metrics.write_text(json.dumps({"status": "failed"}))

    with pytest.raises(ValueError, match="no successful ingest marker"):
        plan_collation(tmp_path)


def test_collation_validates_custom_interface_residue_minimum(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    for pdb_id in ("1abc", "2def"):
        path = tmp_path / "raw_entries" / pdb_id[1:3] / pdb_id / "interfaces.parquet"
        frame = pd.read_parquet(path)
        frame["interface_chain_1_residue_numbers"] = [[1]]
        frame["interface_chain_1_residue_indices"] = [[0]]
        frame["interface_chain_2_residue_numbers"] = [[11]]
        frame["interface_chain_2_residue_indices"] = [[0]]
        frame.to_parquet(path, index=False)
        _set_interface_threshold(path, 1)

    report = run_collation(tmp_path, memory_limit="1GB")

    assert report["interface_min_residues"] == 1


def test_collation_rejects_interface_below_frozen_residue_minimum(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    path = tmp_path / "raw_entries/ab/1abc/interfaces.parquet"
    frame = pd.read_parquet(path)
    for column in (
        "interface_chain_1_residue_numbers",
        "interface_chain_1_residue_indices",
        "interface_chain_2_residue_numbers",
        "interface_chain_2_residue_indices",
    ):
        frame[column] = frame[column].map(lambda values: values[:6])
    frame.to_parquet(path, index=False)
    _set_interface_threshold(path, 7)

    with pytest.raises(ValueError, match="invalid protein interfaces"):
        run_collation(tmp_path, memory_limit="1GB")


def test_collation_rejects_mixed_interface_thresholds(tmp_path: Path) -> None:
    _write_release(tmp_path)
    _set_interface_threshold(
        tmp_path / "raw_entries/de/2def/interfaces.parquet",
        9,
    )

    with pytest.raises(ValueError, match="mixed interface.min_interface_residues"):
        plan_collation(tmp_path)


def test_collation_rejects_interface_table_without_threshold(tmp_path: Path) -> None:
    _write_release(tmp_path)
    path = tmp_path / "raw_entries/de/2def/interfaces.parquet"
    table = pq.read_table(path)
    pq.write_table(table.replace_schema_metadata(None), path)

    with pytest.raises(ValueError, match="does not record its minimum residue count"):
        plan_collation(tmp_path)


def test_targeted_repair_preserves_unaffected_release_only_columns(
    tmp_path: Path,
) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    installed = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    installed["release_only"] = installed["entry_pdb_id"].map(
        {
            "1abc": "old",
            "2def": "keep",
        }
    )
    installed.to_parquet(tmp_path / "index/annotation_table.parquet", index=False)
    chain_path = tmp_path / "index/entry_chains.parquet"
    before_chain_stat = chain_path.stat()

    metadata_path = tmp_path / "raw_entries/ab/1abc/entry_metadata.parquet"
    metadata = pd.read_parquet(metadata_path)
    metadata["entry_pH"] = 6.5
    metadata.to_parquet(metadata_path, index=False)

    report = repair_collation(tmp_path, ["1ABC", "1abc"], threads=2, memory_limit="1GB")

    assert report["mode"] == "targeted_repair"
    assert report["status"] == "requires_downstream_repair"
    assert report["repaired_entry_count"] == 1
    repaired = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert "entry_pH" not in repaired.columns
    assert (
        repaired.loc[repaired["entry_pdb_id"].eq("1abc"), "release_only"].isna().all()
    )
    assert (
        repaired.loc[repaired["entry_pdb_id"].eq("2def"), "release_only"]
        .eq("keep")
        .all()
    )
    repaired_metadata = pd.read_parquet(tmp_path / "index/entry_metadata.parquet")
    assert repaired_metadata.set_index("entry_pdb_id").loc["1abc", "entry_pH"] == 6.5
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

    finalized_report = finalize_repair_marker(tmp_path)
    assert finalized_report is not None
    assert finalized_report["status"] == "complete"
    assert finalized_report["downstream_repair_complete"] is True


@pytest.mark.parametrize("missing_installed_columns", [False, True])
def test_targeted_repair_updates_chain_annotations(
    tmp_path: Path, missing_installed_columns: bool
) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    path = tmp_path / "index/entry_chains.parquet"
    installed = pd.read_parquet(path)
    columns = ["chain_sequence_noncanonical", "chain_modified_residues"]
    if missing_installed_columns:
        installed = installed.drop(columns=columns)
    installed["release_only"] = "keep"
    installed.to_parquet(path, index=False)
    replacement_path = tmp_path / "raw_entries/ab/1abc/entry_chains.parquet"
    replacement = pd.read_parquet(replacement_path)
    replacement.loc[0, "chain_sequence_noncanonical"] = "(MSE)" + "A" * 299
    replacement.at[0, "chain_modified_residues"] = ["1:MSE:A:1"]
    replacement.to_parquet(replacement_path, index=False)
    metadata_path = tmp_path / "raw_entries/ab/1abc/entry_metadata.parquet"
    metadata = pd.read_parquet(metadata_path)
    metadata["entry_pH_min"] = 6.0
    metadata["entry_pH_max"] = 7.0
    metadata["entry_has_ligand_of_interest"] = True
    metadata.to_parquet(metadata_path, index=False)
    ligand_path = tmp_path / "raw_entries/ab/1abc.parquet"
    ligands = pd.read_parquet(ligand_path)
    ligands["ligand_contact_area"] = 42.0
    ligands["ligand_unresolved_atoms"] = [["5:LIG:L:5:C1"] for _ in range(len(ligands))]
    ligands.to_parquet(ligand_path, index=False)

    report = repair_collation(tmp_path, ["1abc"], memory_limit="1GB")

    assert report["status"] == "requires_downstream_repair"
    repaired_metadata = pd.read_parquet(tmp_path / "index/entry_metadata.parquet")
    metadata_columns = ["entry_pH_min", "entry_pH_max", "entry_has_ligand_of_interest"]
    assert repaired_metadata.loc[
        repaired_metadata.entry_pdb_id.eq("1abc"), metadata_columns
    ].to_dict("records") == metadata[metadata_columns].to_dict("records")
    repaired_ligands = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    columns_to_check = ["ligand_contact_area", "ligand_unresolved_atoms"]
    pd.testing.assert_frame_equal(
        repaired_ligands.loc[
            repaired_ligands.entry_pdb_id.eq("1abc"), columns_to_check
        ].reset_index(drop=True),
        ligands[columns_to_check].reset_index(drop=True),
    )
    repaired = pd.read_parquet(path)
    selected = repaired.loc[repaired.entry_pdb_id.eq("1abc"), replacement.columns]
    pd.testing.assert_frame_equal(selected.reset_index(drop=True), replacement)
    untouched = repaired.loc[repaired.entry_pdb_id.eq("2def"), installed.columns]
    pd.testing.assert_frame_equal(
        untouched.reset_index(drop=True),
        installed.loc[installed.entry_pdb_id.eq("2def")].reset_index(drop=True),
    )
    if missing_installed_columns:
        assert (
            repaired.loc[repaired.entry_pdb_id.eq("2def"), columns].isna().all().all()
        )


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("chain_length", 301),
        ("chain_sequence", "G" * 300),
        ("chain_asym_id", "C"),
        ("chain_auth_id", "Z"),
        ("chain_is_ligand_like", True),
    ],
)
def test_targeted_repair_rejects_chain_metadata_changes(
    tmp_path: Path, column: str, value: object
) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    chains = pd.read_parquet(tmp_path / "raw_entries/ab/1abc/entry_chains.parquet")
    chains[column] = value
    chains.to_parquet(
        tmp_path / "raw_entries/ab/1abc/entry_chains.parquet", index=False
    )

    installed_paths = list((tmp_path / "index").glob("*.parquet")) + [
        tmp_path / "index/collation.json"
    ]
    before = {path: path.read_bytes() for path in installed_paths}
    with pytest.raises(ValueError, match="changed entry_chains"):
        repair_collation(tmp_path, ["1abc"], memory_limit="1GB")
    assert {path: path.read_bytes() for path in installed_paths} == before


@pytest.mark.parametrize(
    ("filename", "column", "value", "table"),
    [
        ("entry_source.parquet", "source_mmcif_minor_revision", 99, "entry_sources"),
        (
            "entry_biounit_chains.parquet",
            "chain_instance",
            "2.A",
            "entry_biounit_chains",
        ),
    ],
)
def test_metadata_repair_preserves_release_when_scoring_sources_change(
    tmp_path: Path, filename: str, column: str, value: object, table: str
) -> None:
    _write_release(tmp_path)
    run_collation(tmp_path, memory_limit="1GB")
    entry_dir = tmp_path / "raw_entries/ab/1abc"
    chain_path = entry_dir / "entry_chains.parquet"
    chains = pd.read_parquet(chain_path)
    chains["chain_sequence_noncanonical"] = "(MSE)AA"
    chains.to_parquet(chain_path, index=False)
    source = entry_dir / filename
    frame = pd.read_parquet(source)
    frame.loc[0, column] = value
    frame.to_parquet(source, index=False)
    paths = list((tmp_path / "index").glob("*.parquet")) + [
        tmp_path / "index/collation.json"
    ]
    before = {path: path.read_bytes() for path in paths}
    with pytest.raises(ValueError, match=f"changed {table}"):
        repair_collation(tmp_path, ["1abc"], memory_limit="1GB")
    assert {path: path.read_bytes() for path in paths} == before


def test_shard_rejects_inputs_changed_after_plan(tmp_path: Path) -> None:
    _write_release(tmp_path)
    plan_collation(tmp_path)
    annotation = tmp_path / "raw_entries/ab/1abc.parquet"
    frame = pd.read_parquet(annotation)
    frame["ligand_unique_ccd_code"] = "CHANGED"
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


def test_finalize_uses_frozen_shards_after_raw_inputs_change(tmp_path: Path) -> None:
    _write_release(tmp_path)
    plan_collation(tmp_path)
    collate_shard(tmp_path, "ab", memory_limit="1GB")
    collate_shard(tmp_path, "de", memory_limit="1GB")
    annotation = tmp_path / "raw_entries/ab/1abc.parquet"
    frame = pd.read_parquet(annotation)
    frame["ligand_unique_ccd_code"] = "CHANGED"
    frame.to_parquet(annotation, index=False)

    report = finalize_collation(tmp_path, threads=2, memory_limit="1GB")

    assert report["status"] == "complete"
    installed = pd.read_parquet(tmp_path / "index/annotation_table.parquet")
    assert "CHANGED" not in set(installed["ligand_unique_ccd_code"])


def test_final_install_fails_closed_on_partial_replacement(
    tmp_path: Path, monkeypatch
) -> None:
    names = (
        "annotation",
        "system_validation",
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
    assert not final_paths["system_validation"].exists()


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
        [
            "shard",
            str(tmp_path),
            "--pdb-manifest",
            str(manifest),
        ]
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
        comparability={"1abc__1__1.I": False},
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


def test_final_validation_rejects_missing_chain_sequence(tmp_path: Path) -> None:
    _write_release(tmp_path)
    chain_path = tmp_path / "raw_entries/ab/1abc/entry_chains.parquet"
    chains = pd.read_parquet(chain_path)
    chains.loc[0, "chain_sequence"] = ""
    chains.to_parquet(chain_path, index=False)

    with pytest.raises(ValueError, match="missing_sequences"):
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
