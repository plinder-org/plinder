# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Run a resumable slice of per-entry V3 ingest in one Python process."""

from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from collections import Counter
from pathlib import Path
from typing import Any

from plinder.data.pipeline.ingest_one import (
    PDB_NEXTGEN_ROOT_ENV,
    VALIDATION_ROOT_ENV,
    entry_metrics_paths,
    ingest_one_pdb,
    normalize_pdb_id,
)


def load_manifest(path: Path) -> list[str]:
    """Load unique PDB IDs from a newline-delimited manifest."""
    pdb_ids = [
        normalize_pdb_id(line)
        for raw_line in path.read_text().splitlines()
        if (line := raw_line.strip()) and not line.startswith("#")
    ]
    duplicates = sorted(
        pdb_id for pdb_id, count in Counter(pdb_ids).items() if count > 1
    )
    if duplicates:
        raise ValueError(f"duplicate PDB IDs in {path}: {duplicates}")
    return pdb_ids


def manifest_slice(
    pdb_ids: list[str], *, batch_index: int, batch_size: int
) -> list[str]:
    """Select one zero-based, fixed-size manifest slice."""
    if batch_index < 0:
        raise ValueError("batch_index must be non-negative")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    start = batch_index * batch_size
    return pdb_ids[start : start + batch_size]


def _completed_metrics(output_root: Path, pdb_id: str) -> Path | None:
    for metrics_path in entry_metrics_paths(output_root, pdb_id):
        if not metrics_path.is_file():
            continue
        try:
            metrics = json.loads(metrics_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        status = metrics.get("status")
        if status == "skipped_no_systems":
            return metrics_path
        if status != "complete":
            continue
        outputs = metrics.get("outputs", {})
        required_files = (outputs.get("entry_parquet"), outputs.get("ligand_parquet"))
        entry_directory = outputs.get("entry_directory")
        if not all(path and Path(path).is_file() for path in required_files):
            continue
        if not entry_directory or not Path(entry_directory).is_dir():
            continue
        return metrics_path
    return None


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def ingest_pdb_batch(
    *,
    pdb_ids: list[str],
    output_root: Path,
    cif_root: Path,
    validation_root: Path,
    force: bool = False,
    job_id: str = "local",
    batch_index: int = 0,
) -> tuple[Path, bool]:
    """Ingest PDB IDs sequentially and continue after per-entry failures."""
    output_root = output_root.resolve()
    metrics_path = (
        output_root / "metrics" / f"ingest-batch-{job_id}-{batch_index:05d}.json"
    )
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "job_id": job_id,
        "batch_index": batch_index,
        "pdb_ids": pdb_ids,
        "status": "running",
        "entries": [],
    }
    started = time.perf_counter()
    had_failures = False
    try:
        for pdb_id in pdb_ids:
            entry_started = time.perf_counter()
            completed = None if force else _completed_metrics(output_root, pdb_id)
            if completed is not None:
                payload["entries"].append(
                    {
                        "pdb_id": pdb_id,
                        "status": "skipped_complete",
                        "metrics": str(completed),
                        "wall_seconds": time.perf_counter() - entry_started,
                    }
                )
                continue
            try:
                entry_metrics = ingest_one_pdb(
                    pdb_id=pdb_id,
                    output_root=output_root,
                    cif_root=cif_root,
                    validation_root=validation_root,
                    force=True,
                )
                entry_status = json.loads(entry_metrics.read_text()).get("status")
                payload["entries"].append(
                    {
                        "pdb_id": pdb_id,
                        "status": entry_status,
                        "metrics": str(entry_metrics),
                        "wall_seconds": time.perf_counter() - entry_started,
                    }
                )
            except Exception as exc:
                had_failures = True
                payload["entries"].append(
                    {
                        "pdb_id": pdb_id,
                        "status": "failed",
                        "error": repr(exc),
                        "traceback": traceback.format_exc(),
                        "wall_seconds": time.perf_counter() - entry_started,
                    }
                )
            finally:
                payload["total_wall_seconds"] = time.perf_counter() - started
                _write_json(metrics_path, payload)
        payload["status"] = "completed_with_failures" if had_failures else "complete"
    finally:
        payload["total_wall_seconds"] = time.perf_counter() - started
        _write_json(metrics_path, payload)
    return metrics_path, had_failures


def _default_batch_index() -> int:
    return int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))


def _default_job_id() -> str:
    return os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID", "local"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--batch-size", type=int, required=True)
    parser.add_argument("--batch-index", type=int, default=_default_batch_index())
    parser.add_argument(
        "--cif-root",
        type=Path,
        default=os.environ.get(PDB_NEXTGEN_ROOT_ENV),
        required=PDB_NEXTGEN_ROOT_ENV not in os.environ,
    )
    parser.add_argument(
        "--validation-root",
        type=Path,
        default=os.environ.get(VALIDATION_ROOT_ENV),
        required=VALIDATION_ROOT_ENV not in os.environ,
    )
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    pdb_ids = manifest_slice(
        load_manifest(args.manifest),
        batch_index=args.batch_index,
        batch_size=args.batch_size,
    )
    metrics_path, had_failures = ingest_pdb_batch(
        pdb_ids=pdb_ids,
        output_root=args.output_root,
        cif_root=args.cif_root,
        validation_root=args.validation_root,
        force=args.force,
        job_id=_default_job_id(),
        batch_index=args.batch_index,
    )
    print(metrics_path.read_text(), end="")
    if had_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
