# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import json
from pathlib import Path

import pytest
from plinder.data.pipeline import ingest_batch
from plinder.data.pipeline.ingest_batch import (
    ingest_pdb_batch,
    load_manifest,
    manifest_slice,
)


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
        entry_path.touch()
        ligand_path.touch()
        metrics_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "outputs": {
                        "entry_parquet": str(entry_path),
                        "entry_directory": str(entry_directory),
                        "ligand_parquet": str(ligand_path),
                    },
                }
            )
        )
        return metrics_path

    monkeypatch.setattr(ingest_batch, "ingest_one_pdb", fake_ingest_one_pdb)

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


def test_batch_resumes_entries_previously_skipped_without_systems(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_root = tmp_path / "output"
    entry_metrics = output_root / "metrics" / "ingest-one-1abc.json"
    entry_metrics.parent.mkdir(parents=True)
    entry_metrics.write_text(json.dumps({"status": "skipped_no_systems"}))

    def unexpected_ingest(**_kwargs: object) -> Path:
        pytest.fail("a completed no-system entry must not be ingested again")

    monkeypatch.setattr(ingest_batch, "ingest_one_pdb", unexpected_ingest)

    metrics_path, had_failures = ingest_pdb_batch(
        pdb_ids=["1abc"],
        output_root=output_root,
        cif_root=tmp_path / "cif",
        validation_root=tmp_path / "validation",
    )

    assert not had_failures
    metrics = json.loads(metrics_path.read_text())
    assert metrics["entries"][0]["status"] == "skipped_complete"
