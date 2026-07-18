# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Run the per-entry portion of batched V3 ingest."""

from __future__ import annotations

import json
import os
import re
import resource
import shutil
import time
import traceback
from pathlib import Path
from typing import Any, Callable, Mapping, TypeVar

import pandas as pd

from plinder.data.get_system_annotations import GetPlinderAnnotation
from plinder.data.pipeline.utils import save_ligand_batch

PDB_NEXTGEN_ROOT_ENV = "PLINDER_PDB_NEXTGEN_ROOT"
VALIDATION_ROOT_ENV = "PLINDER_VALIDATION_ROOT"
REQUIRED_REFERENCE_FILES = (
    Path("dbs/components/components.cif"),
    Path("dbs/components/components.parquet"),
    Path("dbs/cofactors/cofactors.json"),
    Path("dbs/affinity/affinity.json"),
)

T = TypeVar("T")


def resolve_source_roots(
    *,
    data_dir: Path,
    cif_root: str | Path | None = None,
    validation_root: str | Path | None = None,
) -> tuple[Path, Path]:
    """Resolve V3 source roots from config, environment, or local defaults.

    Explicit arguments take precedence over environment variables. Relative
    configured paths are interpreted relative to the Plinder data directory so
    the same configuration works in local and Metaflow deployments.
    """

    def resolve(
        configured: str | Path | None,
        environment_variable: str,
        fallback: str,
    ) -> Path:
        value = configured or os.environ.get(environment_variable)
        path = Path(value) if value else Path(fallback)
        if not path.is_absolute():
            path = data_dir / path
        return path.expanduser().resolve()

    return (
        resolve(cif_root, PDB_NEXTGEN_ROOT_ENV, "ingest"),
        resolve(validation_root, VALIDATION_ROOT_ENV, "reports"),
    )


def normalize_pdb_id(value: str) -> str:
    """Validate and normalize one four-character PDB ID."""
    pdb_id = value.strip().lower()
    if re.fullmatch(r"[0-9][a-z0-9]{3}", pdb_id) is None:
        raise ValueError(f"invalid four-character PDB ID: {value!r}")
    return pdb_id


def entry_metrics_paths(output_root: Path, pdb_id: str) -> tuple[Path, Path]:
    """Return the sharded metrics path followed by the legacy flat path."""
    pdb_id = normalize_pdb_id(pdb_id)
    filename = f"ingest-one-{pdb_id}.json"
    metrics_root = output_root / "metrics"
    return metrics_root / pdb_id[1:3] / filename, metrics_root / filename


def resolve_entry_paths(
    pdb_id: str,
    *,
    cif_root: Path,
    validation_root: Path,
) -> tuple[Path, Path]:
    """Resolve managed NextGen and validation inputs without staging them."""
    pdb_id = normalize_pdb_id(pdb_id)
    code = pdb_id[1:3]
    cif_file = (
        cif_root / code / f"pdb_0000{pdb_id}" / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
    )
    validation_file = validation_root / code / pdb_id / f"{pdb_id}_validation.xml.gz"
    return cif_file, validation_file


def check_reference_data(output_root: Path) -> None:
    """Require global annotation inputs prepared once for all array tasks."""
    missing = [
        output_root / relative
        for relative in REQUIRED_REFERENCE_FILES
        if not (output_root / relative).is_file()
        or (output_root / relative).stat().st_size == 0
    ]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            "missing shared ingest reference data; run the online "
            f"reference-data preparation command first:\n{formatted}"
        )


def _usage() -> dict[str, float]:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "user_cpu_seconds": float(usage.ru_utime),
        "system_cpu_seconds": float(usage.ru_stime),
        # Linux reports ru_maxrss in KiB.
        "peak_rss_mb": float(usage.ru_maxrss) / 1024.0,
    }


def _run_timed(
    stage: str,
    func: Callable[[], T],
    timings: list[dict[str, float | str]],
) -> T:
    before = _usage()
    started = time.perf_counter()
    status = "complete"
    error = None
    try:
        return func()
    except BaseException as exc:
        status = "failed"
        error = repr(exc)
        raise
    finally:
        after = _usage()
        timing: dict[str, float | str] = {
            "stage": stage,
            "status": status,
            "wall_seconds": time.perf_counter() - started,
            "user_cpu_seconds": (
                after["user_cpu_seconds"] - before["user_cpu_seconds"]
            ),
            "system_cpu_seconds": (
                after["system_cpu_seconds"] - before["system_cpu_seconds"]
            ),
            "peak_rss_mb": after["peak_rss_mb"],
        }
        if error is not None:
            timing["error"] = error
        timings.append(timing)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _entry_outputs_complete(
    *,
    metrics_path: Path,
    entry_parquet: Path,
    entry_directory: Path,
    ligand_parquet: Path,
) -> bool:
    """Return whether a prior per-entry run completed atomically."""
    if not metrics_path.is_file():
        return False
    try:
        metrics = json.loads(metrics_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return bool(
        metrics.get("status") == "complete"
        and entry_parquet.is_file()
        and entry_directory.is_dir()
        and ligand_parquet.is_file()
    )


def completed_entry_metrics(output_root: Path, pdb_id: str) -> Path | None:
    """Return the metrics file when one V3 entry has a complete output set."""
    for metrics_path in entry_metrics_paths(output_root, pdb_id):
        if not metrics_path.is_file():
            continue
        try:
            metrics = json.loads(metrics_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if metrics.get("status") == "skipped_no_systems":
            return metrics_path
        outputs = metrics.get("outputs", {})
        entry_parquet = outputs.get("entry_parquet")
        entry_directory = outputs.get("entry_directory")
        ligand_parquet = outputs.get("ligand_parquet")
        if not all((entry_parquet, entry_directory, ligand_parquet)):
            continue
        if _entry_outputs_complete(
            metrics_path=metrics_path,
            entry_parquet=Path(entry_parquet),
            entry_directory=Path(entry_directory),
            ligand_parquet=Path(ligand_parquet),
        ):
            return metrics_path
    return None


def _clear_entry_outputs(
    *,
    entry_parquet: Path,
    entry_directory: Path,
    ligand_parquet: Path,
) -> None:
    """Remove one incomplete or explicitly replaced per-entry output set."""
    entry_parquet.unlink(missing_ok=True)
    ligand_parquet.unlink(missing_ok=True)
    shutil.rmtree(entry_directory, ignore_errors=True)


def ingest_one_pdb(
    *,
    pdb_id: str,
    output_root: Path,
    cif_root: Path,
    validation_root: Path,
    force: bool = False,
    check_references: bool = True,
    annotation_cfg: Mapping[str, Any] | None = None,
    entry_cfg: Mapping[str, Any] | None = None,
) -> Path:
    """Generate all per-entry V3 Parquets and canonical ASU SDFs."""
    pdb_id = normalize_pdb_id(pdb_id)
    output_root = output_root.resolve()
    cif_file, validation_file = resolve_entry_paths(
        pdb_id,
        cif_root=cif_root,
        validation_root=validation_root,
    )
    if not cif_file.is_file():
        raise FileNotFoundError(f"missing NextGen mmCIF: {cif_file}")
    if check_references:
        check_reference_data(output_root)

    code = pdb_id[1:3]
    raw_entry_root = output_root / "raw_entries" / code
    entry_parquet = raw_entry_root / f"{pdb_id}.parquet"
    entry_directory = raw_entry_root / pdb_id
    ligand_parquet = output_root / "ligands" / f"{pdb_id}.parquet"
    metrics_path, legacy_metrics_path = entry_metrics_paths(output_root, pdb_id)
    complete = any(
        _entry_outputs_complete(
            metrics_path=candidate,
            entry_parquet=entry_parquet,
            entry_directory=entry_directory,
            ligand_parquet=ligand_parquet,
        )
        for candidate in (metrics_path, legacy_metrics_path)
    )
    if complete and not force:
        raise FileExistsError(
            f"output already exists: {entry_parquet}; pass --force to replace it"
        )
    has_partial_outputs = (
        entry_parquet.exists() or entry_directory.exists() or ligand_parquet.exists()
    )
    if force or has_partial_outputs:
        if has_partial_outputs and not force:
            print(f"clearing incomplete outputs before retrying {pdb_id}")
        _clear_entry_outputs(
            entry_parquet=entry_parquet,
            entry_directory=entry_directory,
            ligand_parquet=ligand_parquet,
        )

    timings: list[dict[str, float | str]] = []
    summary: dict[str, Any] = {
        "pdb_id": pdb_id,
        "status": "running",
        "inputs": {
            "mmcif": str(cif_file),
            "validation_xml": str(validation_file),
            "validation_xml_exists": validation_file.is_file(),
        },
        "outputs": {
            "entry_parquet": str(entry_parquet),
            "entry_directory": str(entry_directory),
            "ligand_parquet": str(ligand_parquet),
        },
        "timings": timings,
    }
    if annotation_cfg or entry_cfg:
        summary["configuration"] = {
            "annotation": dict(annotation_cfg or {}),
            "entry": {
                key: value
                for key, value in dict(entry_cfg or {}).items()
                if key != "save_folder"
            },
        }
    total_started = time.perf_counter()
    try:
        raw_entry_root.mkdir(parents=True, exist_ok=True)
        ligand_parquet.parent.mkdir(parents=True, exist_ok=True)

        def annotate() -> pd.DataFrame:
            annotation_options = dict(annotation_cfg or {})
            entry_options = dict(entry_cfg or {})
            entry_options.pop("save_folder", None)
            if entry_options:
                annotation_options["entry_cfg"] = entry_options
            annotation = GetPlinderAnnotation(
                cif_file,
                validation_file,
                save_folder=raw_entry_root,
                **annotation_options,
            ).annotate()
            return annotation if annotation is not None else pd.DataFrame()

        annotation = _run_timed("annotate_entry", annotate, timings)
        if annotation.empty:
            entry_parquet.unlink(missing_ok=True)
            ligand_parquet.unlink(missing_ok=True)
            shutil.rmtree(entry_directory, ignore_errors=True)
            summary["outputs"] = {
                "entry_parquet": None,
                "entry_directory": None,
                "ligand_parquet": None,
            }
            summary["counts"] = {
                "annotation_rows": 0,
                "systems": 0,
                "ligand_ids": 0,
                "canonical_ligand_sdfs": 0,
            }
            summary["status"] = "skipped_no_systems"
        else:

            def write_entry_parquet() -> None:
                annotation.to_parquet(entry_parquet, index=False)

            _run_timed(
                "write_entry_parquet",
                write_entry_parquet,
                timings,
            )

            def write_ligand_parquet() -> None:
                save_ligand_batch(
                    data_dir=output_root,
                    annotation=annotation,
                    output_path=ligand_parquet,
                )

            _run_timed(
                "write_ligand_parquet",
                write_ligand_parquet,
                timings,
            )

            ligand_sdfs = sorted((entry_directory / "ligand_files").glob("*.sdf"))
            summary["counts"] = {
                "annotation_rows": len(annotation),
                "systems": int(annotation["system_id"].nunique()),
                "ligand_ids": int(annotation["ligand_id"].nunique()),
                "canonical_ligand_sdfs": len(ligand_sdfs),
            }
            summary["status"] = "complete"
    except BaseException as exc:
        summary["status"] = "failed"
        summary["error"] = repr(exc)
        summary["traceback"] = traceback.format_exc()
        raise
    finally:
        summary["total_wall_seconds"] = time.perf_counter() - total_started
        summary["final_resource_usage"] = _usage()
        _write_json(metrics_path, summary)
    return metrics_path
