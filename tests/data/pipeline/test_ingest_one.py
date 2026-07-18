# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import json
from pathlib import Path

import pandas as pd
import pytest
from plinder.data.pipeline import ingest_one
from plinder.data.pipeline.ingest_one import (
    REQUIRED_REFERENCE_FILES,
    check_reference_data,
    ingest_one_pdb,
    normalize_pdb_id,
    resolve_entry_paths,
    resolve_source_roots,
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
            return pd.DataFrame(
                {
                    "system_id": ["8grn__1__1.A__1.C"],
                    "ligand_id": ["8grn__1.C"],
                }
            )

    def fake_save_ligand_batch(
        *, data_dir: Path, annotation: pd.DataFrame, output_path: Path
    ) -> None:
        assert data_dir == output_root.resolve()
        assert len(annotation) == 1
        annotation[["ligand_id"]].to_parquet(output_path, index=False)

    monkeypatch.setattr(ingest_one, "GetPlinderAnnotation", FakeAnnotation)
    monkeypatch.setattr(ingest_one, "save_ligand_batch", fake_save_ligand_batch)

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

    monkeypatch.setattr(ingest_one, "GetPlinderAnnotation", EmptyAnnotation)
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
            return pd.DataFrame(
                {
                    "system_id": ["8grn__1__1.A__1.C"],
                    "ligand_id": ["8grn__1.C"],
                }
            )

    def fake_save_ligand_batch(
        *, data_dir: Path, annotation: pd.DataFrame, output_path: Path
    ) -> None:
        annotation[["ligand_id"]].to_parquet(output_path, index=False)

    monkeypatch.setattr(ingest_one, "GetPlinderAnnotation", FakeAnnotation)
    monkeypatch.setattr(ingest_one, "save_ligand_batch", fake_save_ligand_batch)

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
