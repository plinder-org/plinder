# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import json
from pathlib import Path

import pandas as pd
import pytest
from plinder.data.pipeline import ingest
from plinder.data.pipeline.ingest import (
    REQUIRED_REFERENCE_FILES,
    EntryInput,
    balance_entries,
    build_parser,
    check_reference_data,
    completed_entry_metrics,
    discover_entries,
    ingest_one_pdb,
    ingest_pdb_batch,
    load_manifest,
    manifest_slice,
    normalize_pdb_id,
    resolve_entry_paths,
    resolve_source_roots,
    write_manifest,
)


def _write_fake_sidecars(
    entry_dir: Path,
    pdb_id: str,
    *,
    interfaces: list[dict[str, object]] | None = None,
) -> None:
    pd.DataFrame(
        {
            "entry_pdb_id": [pdb_id],
            "chain_receptor_type": ["protein"],
            "chain_is_ligand_like": [False],
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
    pd.DataFrame({"entry_pdb_id": [pdb_id]}).to_parquet(
        entry_dir / "entry_source.parquet", index=False
    )
    pd.DataFrame({"entry_pdb_id": [pdb_id]}).to_parquet(
        entry_dir / "entry_metadata.parquet", index=False
    )
    interface_columns = [
        "entry_pdb_id",
        "system_id",
        "system_biounit_id",
        "interface_chain_1",
        "interface_chain_2",
        "interface_chain_1_residue_numbers",
        "interface_chain_1_residue_indices",
        "interface_chain_2_residue_numbers",
        "interface_chain_2_residue_indices",
        "interface_num_contact_residue_pairs",
    ]
    pd.DataFrame(interfaces or [], columns=interface_columns).to_parquet(
        entry_dir / "interfaces.parquet", index=False
    )


def test_resolve_entry_paths_uses_managed_archive_layout(tmp_path: Path) -> None:
    cif_file, validation_file = resolve_entry_paths(
        "8GRN",
        cif_root=tmp_path / "nextgen",
        validation_root=tmp_path / "validation",
    )

    assert cif_file == (
        tmp_path / "nextgen" / "gr" / "pdb_00008grn" / "pdb_00008grn_xyz-enrich.cif.gz"
    )
    assert validation_file == (
        tmp_path / "validation" / "gr" / "8grn" / "8grn_validation.xml.gz"
    )


def test_source_roots_prefer_config_then_environment_then_local_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PLINDER_PDB_NEXTGEN_ROOT", str(tmp_path / "env-cif"))
    monkeypatch.setenv("PLINDER_VALIDATION_ROOT", str(tmp_path / "env-validation"))

    assert resolve_source_roots(data_dir=tmp_path) == (
        (tmp_path / "env-cif").resolve(),
        (tmp_path / "env-validation").resolve(),
    )
    assert resolve_source_roots(
        data_dir=tmp_path,
        cif_root="configured-cif",
        validation_root="configured-validation",
    ) == (
        (tmp_path / "configured-cif").resolve(),
        (tmp_path / "configured-validation").resolve(),
    )

    monkeypatch.delenv("PLINDER_PDB_NEXTGEN_ROOT")
    monkeypatch.delenv("PLINDER_VALIDATION_ROOT")
    assert resolve_source_roots(data_dir=tmp_path) == (
        (tmp_path / "ingest").resolve(),
        (tmp_path / "reports").resolve(),
    )

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    assert resolve_source_roots(
        data_dir=tmp_path,
        cif_root="~/nextgen",
        validation_root="~/validation",
    ) == (
        (tmp_path / "home/nextgen").resolve(),
        (tmp_path / "home/validation").resolve(),
    )


@pytest.mark.parametrize("value", ["", "8gr", "abcd", "8grn-extra"])
def test_normalize_pdb_id_rejects_invalid_ids(value: str) -> None:
    with pytest.raises(ValueError, match="invalid four-character PDB ID"):
        normalize_pdb_id(value)


def test_reference_data_check_reports_every_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError) as exc_info:
        check_reference_data(tmp_path)

    message = str(exc_info.value)
    assert "components.parquet" in message
    assert "cofactors.json" in message
    assert "affinity.json" in message


def test_ingest_one_pdb_writes_entry_outputs_and_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_root = tmp_path / "output"
    cif_root = tmp_path / "nextgen"
    validation_root = tmp_path / "validation"
    cif_file, validation_file = resolve_entry_paths(
        "8grn",
        cif_root=cif_root,
        validation_root=validation_root,
    )
    for path in (cif_file, validation_file):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    for relative in REQUIRED_REFERENCE_FILES:
        path = output_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")

    class FakeAnnotation:
        def __init__(
            self,
            cif_path: Path,
            validation_path: Path,
            *,
            save_folder: Path,
        ) -> None:
            assert cif_path == cif_file
            assert validation_path == validation_file
            self.save_folder = save_folder

        def annotate(self) -> pd.DataFrame:
            ligand_dir = self.save_folder / "8grn" / "ligand_files"
            ligand_dir.mkdir(parents=True)
            (ligand_dir / "1.C.sdf").touch()
            _write_fake_sidecars(ligand_dir.parent, "8grn")
            return pd.DataFrame(
                {
                    "system_id": ["8grn__1__1.A__1.C"],
                    "ligand_id": ["8grn__1.C"],
                    "system_receptor_type": ["protein"],
                }
            )

    def fake_save_ligand_batch(
        *, data_dir: Path, annotation: pd.DataFrame, output_path: Path
    ) -> None:
        assert data_dir == output_root.resolve()
        assert len(annotation) == 1
        annotation[["ligand_id"]].assign(ligand_is_3d_score_able=True).to_parquet(
            output_path, index=False
        )

    monkeypatch.setattr(ingest, "_get_annotation_class", lambda: FakeAnnotation)
    monkeypatch.setattr(ingest, "_save_ligand_batch", fake_save_ligand_batch)

    metrics_path = ingest_one_pdb(
        pdb_id="8GRN",
        output_root=output_root,
        cif_root=cif_root,
        validation_root=validation_root,
    )

    metrics = json.loads(metrics_path.read_text())
    assert (
        metrics_path
        == output_root.resolve() / "metrics" / "gr" / "ingest-one-8grn.json"
    )
    assert metrics["status"] == "complete"
    assert metrics["counts"] == {
        "annotation_rows": 1,
        "interface_rows": 0,
        "systems": 1,
        "ligand_ids": 1,
        "canonical_ligand_sdfs": 1,
    }
    assert [stage["stage"] for stage in metrics["timings"]] == [
        "annotate_entry",
        "write_entry_parquet",
        "write_ligand_parquet",
    ]
    assert all(stage["status"] == "complete" for stage in metrics["timings"])
    assert (output_root / "raw_entries" / "gr" / "8grn.parquet").is_file()
    assert (output_root / "ligands" / "8grn.parquet").is_file()


def test_ingest_one_pdb_does_not_materialize_entries_without_systems(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_root = tmp_path / "output"
    cif_root = tmp_path / "nextgen"
    validation_root = tmp_path / "validation"
    cif_file, validation_file = resolve_entry_paths(
        "8grn",
        cif_root=cif_root,
        validation_root=validation_root,
    )
    for path in (cif_file, validation_file):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    for relative in REQUIRED_REFERENCE_FILES:
        path = output_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")

    class EmptyAnnotation:
        def __init__(self, *_args, save_folder: Path, **_kwargs) -> None:
            self.save_folder = save_folder

        def annotate(self) -> None:
            stale_dir = self.save_folder / "8grn" / "ligand_files"
            stale_dir.mkdir(parents=True)
            (stale_dir / "A.sdf").touch()
            return None

    monkeypatch.setattr(ingest, "_get_annotation_class", lambda: EmptyAnnotation)
    entry_parquet = output_root / "raw_entries" / "gr" / "8grn.parquet"
    ligand_parquet = output_root / "ligands" / "8grn.parquet"
    for path in (entry_parquet, ligand_parquet):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    metrics_path = ingest_one_pdb(
        pdb_id="8grn",
        output_root=output_root,
        cif_root=cif_root,
        validation_root=validation_root,
        force=True,
    )

    metrics = json.loads(metrics_path.read_text())
    assert metrics["status"] == "skipped_no_systems"
    assert metrics["counts"] == {
        "annotation_rows": 0,
        "interface_rows": 0,
        "systems": 0,
        "ligand_ids": 0,
        "canonical_ligand_sdfs": 0,
    }
    assert [stage["stage"] for stage in metrics["timings"]] == ["annotate_entry"]
    assert metrics["outputs"] == {
        "entry_parquet": None,
        "entry_directory": None,
        "ligand_parquet": None,
    }
    assert not entry_parquet.exists()
    assert not ligand_parquet.exists()
    assert not (output_root / "raw_entries" / "gr" / "8grn").exists()


def test_ingest_one_pdb_materializes_interface_only_entries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_root = tmp_path / "output"
    cif_root = tmp_path / "nextgen"
    validation_root = tmp_path / "validation"
    cif_file, _ = resolve_entry_paths(
        "8grn", cif_root=cif_root, validation_root=validation_root
    )
    cif_file.parent.mkdir(parents=True)
    cif_file.touch()

    class InterfaceOnlyAnnotation:
        def __init__(self, *_args, save_folder: Path, **_kwargs) -> None:
            self.save_folder = save_folder

        def annotate(self) -> pd.DataFrame:
            entry_dir = self.save_folder / "8grn"
            entry_dir.mkdir(parents=True)
            _write_fake_sidecars(
                entry_dir,
                "8grn",
                interfaces=[
                    {
                        "entry_pdb_id": "8grn",
                        "system_id": "8grn__1__1.A--1.B",
                        "system_biounit_id": "1",
                        "interface_chain_1": "1.A",
                        "interface_chain_2": "1.B",
                        "interface_chain_1_residue_numbers": [1, 2, 3],
                        "interface_chain_1_residue_indices": [0, 1, 2],
                        "interface_chain_2_residue_numbers": [4, 5, 6],
                        "interface_chain_2_residue_indices": [0, 1, 2],
                        "interface_num_contact_residue_pairs": 3,
                    }
                ],
            )
            return pd.DataFrame()

    monkeypatch.setattr(
        ingest, "_get_annotation_class", lambda: InterfaceOnlyAnnotation
    )
    metrics_path = ingest_one_pdb(
        pdb_id="8grn",
        output_root=output_root,
        cif_root=cif_root,
        validation_root=validation_root,
        check_references=False,
    )

    metrics = json.loads(metrics_path.read_text())
    assert metrics["status"] == "complete"
    assert metrics["counts"] == {
        "annotation_rows": 0,
        "interface_rows": 1,
        "systems": 0,
        "ligand_ids": 0,
        "canonical_ligand_sdfs": 0,
    }
    assert metrics["outputs"]["entry_parquet"] is None
    assert metrics["outputs"]["ligand_parquet"] is None
    assert (output_root / "raw_entries/gr/8grn/interfaces.parquet").is_file()
    assert completed_entry_metrics(output_root, "8grn") == metrics_path


def test_ingest_one_pdb_retries_partial_outputs_without_force(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "output"
    cif_root = tmp_path / "nextgen"
    validation_root = tmp_path / "validation"
    cif_file, _ = resolve_entry_paths(
        "8grn",
        cif_root=cif_root,
        validation_root=validation_root,
    )
    cif_file.parent.mkdir(parents=True)
    cif_file.touch()

    raw_entry_root = output_root / "raw_entries" / "gr"
    entry_parquet = raw_entry_root / "8grn.parquet"
    entry_directory = raw_entry_root / "8grn"
    stale_file = entry_directory / "stale.txt"
    stale_file.parent.mkdir(parents=True)
    stale_file.write_text("partial")
    entry_parquet.write_text("partial")
    metrics_path = output_root / "metrics" / "ingest-one-8grn.json"
    metrics_path.parent.mkdir(parents=True)
    metrics_path.write_text('{"status": "failed"}')
    annotation_calls = 0

    class FakeAnnotation:
        def __init__(self, *_args, save_folder: Path, **_kwargs) -> None:
            assert not stale_file.exists()
            self.save_folder = save_folder

        def annotate(self) -> pd.DataFrame:
            nonlocal annotation_calls
            annotation_calls += 1
            ligand_dir = self.save_folder / "8grn" / "ligand_files"
            ligand_dir.mkdir(parents=True)
            (ligand_dir / "C.sdf").touch()
            _write_fake_sidecars(ligand_dir.parent, "8grn")
            return pd.DataFrame(
                {
                    "system_id": ["8grn__1__1.A__1.C"],
                    "ligand_id": ["8grn__1.C"],
                    "system_receptor_type": ["protein"],
                }
            )

    def fake_save_ligand_batch(
        *, data_dir: Path, annotation: pd.DataFrame, output_path: Path
    ) -> None:
        annotation[["ligand_id"]].assign(ligand_is_3d_score_able=True).to_parquet(
            output_path, index=False
        )

    monkeypatch.setattr(ingest, "_get_annotation_class", lambda: FakeAnnotation)
    monkeypatch.setattr(ingest, "_save_ligand_batch", fake_save_ligand_batch)

    result = ingest_one_pdb(
        pdb_id="8grn",
        output_root=output_root,
        cif_root=cif_root,
        validation_root=validation_root,
        check_references=False,
    )

    assert json.loads(result.read_text())["status"] == "complete"
    assert annotation_calls == 1
    assert entry_parquet.is_file()
    assert (output_root / "ligands" / "8grn.parquet").is_file()
    with pytest.raises(FileExistsError, match="pass --force"):
        ingest_one_pdb(
            pdb_id="8grn",
            output_root=output_root,
            cif_root=cif_root,
            validation_root=validation_root,
            check_references=False,
        )
    assert annotation_calls == 1


def test_shared_ingest_cli_has_manifest_and_batch_commands() -> None:
    parser = build_parser()

    manifest_args = parser.parse_args(["manifest", "cif", "validation", "manifest.txt"])
    batch_args = parser.parse_args(
        [
            "batch",
            "manifest.txt",
            "output",
            "--batch-size",
            "10",
            "--cif-root",
            "cif",
            "--validation-root",
            "validation",
        ]
    )

    assert manifest_args.command == "manifest"
    assert manifest_args.output_path == Path("manifest.txt")
    assert batch_args.command == "batch"
    assert batch_args.batch_size == 10
    assert batch_args.interface_min_residues is None
    configured = parser.parse_args(
        [
            "batch",
            "manifest.txt",
            "output",
            "--batch-size",
            "10",
            "--cif-root",
            "cif",
            "--validation-root",
            "validation",
            "--interface-min-residues",
            "9",
        ]
    )
    assert configured.interface_min_residues == 9


def test_load_manifest_and_select_slice(tmp_path: Path) -> None:
    manifest = tmp_path / "pdb_ids.txt"
    manifest.write_text("# current release\n1abc\n\n2DEF\n3ghi\n")

    pdb_ids = load_manifest(manifest)

    assert pdb_ids == ["1abc", "2def", "3ghi"]
    assert manifest_slice(pdb_ids, batch_index=1, batch_size=2) == ["3ghi"]


def test_load_manifest_rejects_duplicates(tmp_path: Path) -> None:
    manifest = tmp_path / "pdb_ids.txt"
    manifest.write_text("1abc\n1ABC\n")

    with pytest.raises(ValueError, match="duplicate PDB IDs"):
        load_manifest(manifest)


def test_batch_continues_after_failure_and_resumes_completed_entries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_root = tmp_path / "output"
    calls: list[str] = []

    def fake_ingest_one_pdb(*, pdb_id: str, **_: object) -> Path:
        calls.append(pdb_id)
        if pdb_id == "2def":
            raise ValueError("bad entry")
        entry_path = output_root / "raw_entries" / pdb_id[1:3] / f"{pdb_id}.parquet"
        entry_directory = entry_path.with_suffix("")
        ligand_path = output_root / "ligands" / f"{pdb_id}.parquet"
        metrics_path = (
            output_root / "metrics" / pdb_id[1:3] / f"ingest-one-{pdb_id}.json"
        )
        for path in (entry_path, ligand_path, metrics_path):
            path.parent.mkdir(parents=True, exist_ok=True)
        entry_directory.mkdir()
        pd.DataFrame({"system_receptor_type": ["protein"]}).to_parquet(
            entry_path, index=False
        )
        _write_fake_sidecars(entry_directory, pdb_id)
        pd.DataFrame(
            {"ligand_id": [f"{pdb_id}__1.L"], "ligand_is_3d_score_able": [True]}
        ).to_parquet(ligand_path, index=False)
        metrics_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "counts": {"annotation_rows": 1, "interface_rows": 0},
                    "outputs": {
                        "entry_parquet": str(entry_path),
                        "entry_directory": str(entry_directory),
                        "ligand_parquet": str(ligand_path),
                    },
                }
            )
        )
        return metrics_path

    monkeypatch.setattr(ingest, "ingest_one_pdb", fake_ingest_one_pdb)

    metrics_path, had_failures = ingest_pdb_batch(
        pdb_ids=["1abc", "2def", "3ghi"],
        output_root=output_root,
        cif_root=tmp_path / "cif",
        validation_root=tmp_path / "validation",
        job_id="123",
        batch_index=4,
    )

    assert had_failures
    metrics = json.loads(metrics_path.read_text())
    assert metrics["status"] == "completed_with_failures"
    assert [entry["status"] for entry in metrics["entries"]] == [
        "complete",
        "failed",
        "complete",
    ]
    assert calls == ["1abc", "2def", "3ghi"]

    calls.clear()
    _, had_failures = ingest_pdb_batch(
        pdb_ids=["1abc", "2def", "3ghi"],
        output_root=output_root,
        cif_root=tmp_path / "cif",
        validation_root=tmp_path / "validation",
        job_id="124",
        batch_index=4,
    )

    assert had_failures
    assert calls == ["2def"]
    assert completed_entry_metrics(output_root, "1abc") is not None
    (output_root / "raw_entries/ab/1abc/entry_chains.parquet").unlink()
    assert completed_entry_metrics(output_root, "1abc") is None


def test_batch_resumes_entries_previously_skipped_without_systems(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_root = tmp_path / "output"
    entry_metrics = output_root / "metrics" / "ingest-one-1abc.json"
    entry_metrics.parent.mkdir(parents=True)
    entry_metrics.write_text(
        json.dumps(
            {
                "status": "skipped_no_systems",
                "counts": {"annotation_rows": 0, "interface_rows": 0},
            }
        )
    )

    def unexpected_ingest(**_kwargs: object) -> Path:
        pytest.fail("a completed no-system entry must not be ingested again")

    monkeypatch.setattr(ingest, "ingest_one_pdb", unexpected_ingest)

    metrics_path, had_failures = ingest_pdb_batch(
        pdb_ids=["1abc"],
        output_root=output_root,
        cif_root=tmp_path / "cif",
        validation_root=tmp_path / "validation",
    )

    assert not had_failures
    metrics = json.loads(metrics_path.read_text())
    assert metrics["entries"][0]["status"] == "skipped_complete"


def test_pre_interface_skip_is_not_considered_complete(tmp_path: Path) -> None:
    output_root = tmp_path / "output"
    entry_metrics = output_root / "metrics" / "ingest-one-1abc.json"
    entry_metrics.parent.mkdir(parents=True)
    entry_metrics.write_text(json.dumps({"status": "skipped_no_systems"}))

    assert completed_entry_metrics(output_root, "1abc") is None


def test_discover_entries_tracks_optional_validation(tmp_path: Path) -> None:
    cif_root = tmp_path / "cif"
    validation_root = tmp_path / "validation"
    for pdb_id, content in (("1abc", b"one"), ("2def", b"second")):
        cif_path = (
            cif_root
            / pdb_id[1:3]
            / f"pdb_0000{pdb_id}"
            / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
        )
        cif_path.parent.mkdir(parents=True)
        cif_path.write_bytes(content)
    validation_path = validation_root / "ab" / "1abc" / "1abc_validation.xml.gz"
    validation_path.parent.mkdir(parents=True)
    validation_path.touch()

    entries = discover_entries(cif_root, validation_root)

    assert entries == [
        EntryInput("1abc", 3, True),
        EntryInput("2def", 6, False),
    ]

    unchecked = discover_entries(
        cif_root,
        validation_root,
        check_validation=False,
    )
    assert [entry.validation_exists for entry in unchecked] == [None, None]


def test_balanced_manifest_has_fixed_slices_and_summary(tmp_path: Path) -> None:
    entries = [
        EntryInput(f"{index}abc", size, True)
        for index, size in enumerate((100, 90, 80, 20, 10), start=1)
    ]

    batches = balance_entries(entries, batch_size=2)

    assert [len(batch) for batch in batches] == [2, 2, 1]
    totals = [sum(entry.cif_size for entry in batch) for batch in batches]

    manifest = tmp_path / "pdb_ids.txt"
    summary = write_manifest(
        entries=entries,
        output_path=manifest,
        batch_size=2,
    )
    assert len(manifest.read_text().splitlines()) == len(entries)
    assert summary["batch_count"] == 3
    assert summary["maximum_batch_cif_bytes"] == max(totals)
    assert summary["cif_size_percentiles"]["maximum"] == 100
    assert summary["largest_entries"][0]["cif_size"] == 100
    assert manifest.with_suffix(".txt.json").is_file()
    inventory = manifest.with_suffix(".txt.entries.tsv")
    assert inventory.is_file()
    assert inventory.read_text().splitlines()[0] == (
        "pdb_id\tcif_size\tvalidation_exists"
    )

    manifest_ids = manifest.read_text().splitlines()
    assert [manifest_ids[index : index + 2] for index in range(0, len(entries), 2)] == [
        [entry.pdb_id for entry in batch] for batch in batches
    ]


def test_balance_entries_distributes_large_inputs_across_full_batches() -> None:
    entries = [
        EntryInput(f"{index}abc", size, True)
        for index, size in enumerate((100, 90, 80, 20, 10, 5), start=1)
    ]

    batches = balance_entries(entries, batch_size=2)
    totals = [sum(entry.cif_size for entry in batch) for batch in batches]

    assert [len(batch) for batch in batches] == [2, 2, 2]
    assert max(totals) - min(totals) <= 15
