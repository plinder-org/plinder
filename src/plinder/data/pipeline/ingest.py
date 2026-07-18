# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Discover, ingest, and resume V3 PDB entry batches."""

from __future__ import annotations

import argparse
import heapq
import json
import math
import os
import re
import resource
import shutil
import time
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Any, Callable, Collection, Mapping, TypeVar

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


@dataclass(frozen=True)
class EntryInput:
    pdb_id: str
    cif_size: int
    validation_exists: bool | None


def discover_entries(
    cif_root: Path,
    validation_root: Path,
    *,
    check_validation: bool = True,
    threads: int = 8,
    two_char_codes: Collection[str] | None = None,
    pdb_ids: Collection[str] | None = None,
) -> list[EntryInput]:
    """Discover unique NextGen CIFs and optionally check validation reports."""
    if threads < 1:
        raise ValueError("threads must be positive")
    directory_pattern = re.compile(r"pdb_0000([0-9a-z]{4})")
    selected_codes = (
        {str(code).lower() for code in two_char_codes} if two_char_codes else None
    )
    selected_pdb_ids = {str(pdb_id).lower() for pdb_id in pdb_ids} if pdb_ids else None

    def discover_code_directory(code_directory: Path) -> list[EntryInput]:
        discovered = []
        with os.scandir(code_directory) as directory_entries:
            for directory_entry in directory_entries:
                match = directory_pattern.fullmatch(directory_entry.name)
                if match is None or not directory_entry.is_dir(follow_symlinks=False):
                    continue
                pdb_id = match.group(1)
                if selected_pdb_ids is not None and pdb_id not in selected_pdb_ids:
                    continue
                cif_path = (
                    Path(directory_entry.path) / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
                )
                try:
                    cif_size = cif_path.stat().st_size
                except FileNotFoundError:
                    continue
                discovered.append(
                    EntryInput(
                        pdb_id=pdb_id,
                        cif_size=cif_size,
                        validation_exists=(
                            (
                                validation_root
                                / pdb_id[1:3]
                                / pdb_id
                                / f"{pdb_id}_validation.xml.gz"
                            ).is_file()
                            if check_validation
                            else None
                        ),
                    )
                )
        return discovered

    code_directories = [
        path
        for path in cif_root.iterdir()
        if path.is_dir()
        and (selected_codes is None or path.name.lower() in selected_codes)
    ]
    entries: dict[str, EntryInput] = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for discovered in executor.map(discover_code_directory, code_directories):
            for entry in discovered:
                if entry.pdb_id in entries:
                    raise ValueError(f"duplicate NextGen CIF for PDB ID {entry.pdb_id}")
                entries[entry.pdb_id] = entry
    return sorted(entries.values(), key=lambda entry: entry.pdb_id)


def balance_entries(
    entries: list[EntryInput], *, batch_size: int
) -> list[list[EntryInput]]:
    """Greedily balance bytes across contiguous fixed-size manifest slices."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if not entries:
        return []
    batch_count = math.ceil(len(entries) / batch_size)
    batches: list[list[EntryInput]] = [[] for _ in range(batch_count)]
    remainder = len(entries) % batch_size
    capacities = [batch_size] * batch_count
    if remainder:
        capacities[-1] = remainder
    heap = [(0, 0, index) for index in range(batch_count)]
    heapq.heapify(heap)
    for entry in sorted(entries, key=lambda item: item.cif_size, reverse=True):
        total_size, count, index = heapq.heappop(heap)
        batches[index].append(entry)
        count += 1
        if count < capacities[index]:
            heapq.heappush(heap, (total_size + entry.cif_size, count, index))
    return [
        sorted(batch, key=lambda entry: (entry.cif_size, entry.pdb_id))
        for batch in batches
    ]


def _percentile(values: list[int], percentile: float) -> int:
    """Return a linearly interpolated percentile for integer byte counts."""
    if not values:
        return 0
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return round(ordered[lower] * (1 - weight) + ordered[upper] * weight)


def _write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value)
    temporary.replace(path)


def _validation_status(value: bool | None) -> str:
    if value is None:
        return "unknown"
    return str(value).lower()


def write_manifest(
    *,
    entries: list[EntryInput],
    output_path: Path,
    batch_size: int,
) -> dict[str, Any]:
    """Write balanced PDB IDs plus neighboring input inventory and summary."""
    batches = balance_entries(entries, batch_size=batch_size)
    ordered = [entry for batch in batches for entry in batch]
    _write_text(output_path, "".join(f"{entry.pdb_id}\n" for entry in ordered))
    inventory_path = output_path.with_suffix(output_path.suffix + ".entries.tsv")
    _write_text(
        inventory_path,
        "pdb_id\tcif_size\tvalidation_exists\n"
        + "".join(
            f"{entry.pdb_id}\t{entry.cif_size}\t"
            f"{_validation_status(entry.validation_exists)}\n"
            for entry in entries
        ),
    )
    batch_bytes = [sum(entry.cif_size for entry in batch) for batch in batches]
    cif_sizes = [entry.cif_size for entry in entries]
    thresholds = (1_000_000, 2_000_000, 4_000_000, 8_000_000, 16_000_000)
    validation_checked = all(entry.validation_exists is not None for entry in entries)
    summary: dict[str, Any] = {
        "manifest": str(output_path.resolve()),
        "input_inventory": str(inventory_path.resolve()),
        "batch_size": batch_size,
        "batch_count": len(batches),
        "entry_count": len(entries),
        "total_cif_bytes": sum(entry.cif_size for entry in entries),
        "maximum_batch_cif_bytes": max(batch_bytes, default=0),
        "median_batch_cif_bytes": median(batch_bytes) if batch_bytes else 0,
        "batch_cif_size_percentiles": {
            name: _percentile(batch_bytes, percentile)
            for name, percentile in (
                ("p50", 0.5),
                ("p90", 0.9),
                ("p95", 0.95),
                ("p99", 0.99),
                ("maximum", 1.0),
            )
        },
        "cif_size_percentiles": {
            name: _percentile(cif_sizes, percentile)
            for name, percentile in (
                ("p50", 0.5),
                ("p90", 0.9),
                ("p95", 0.95),
                ("p99", 0.99),
                ("p99_9", 0.999),
                ("maximum", 1.0),
            )
        },
        "cif_size_threshold_counts": {
            str(threshold): sum(size >= threshold for size in cif_sizes)
            for threshold in thresholds
        },
        "largest_entries": [
            asdict(entry)
            for entry in sorted(
                entries, key=lambda entry: entry.cif_size, reverse=True
            )[:100]
        ],
        "validation_checked": validation_checked,
        "missing_validation_count": (
            sum(entry.validation_exists is False for entry in entries)
            if validation_checked
            else None
        ),
        "missing_validation_entries": (
            [asdict(entry) for entry in entries if entry.validation_exists is False]
            if validation_checked
            else None
        ),
    }
    summary_path = output_path.with_suffix(output_path.suffix + ".json")
    _write_text(summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


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


def ingest_pdb_batch(
    *,
    pdb_ids: list[str],
    output_root: Path,
    cif_root: Path,
    validation_root: Path,
    force: bool = False,
    job_id: str = "local",
    batch_index: int = 0,
    annotation_cfg: Mapping[str, Any] | None = None,
    entry_cfg: Mapping[str, Any] | None = None,
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
            completed = None if force else completed_entry_metrics(output_root, pdb_id)
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
                    annotation_cfg=annotation_cfg,
                    entry_cfg=entry_cfg,
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
    """Build the shared manifest-generation and batch-ingest CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    manifest_parser = commands.add_parser(
        "manifest", help="discover source entries and write a balanced manifest"
    )
    manifest_parser.add_argument("cif_root", type=Path)
    manifest_parser.add_argument("validation_root", type=Path)
    manifest_parser.add_argument("output_path", type=Path)
    manifest_parser.add_argument("--batch-size", type=int, default=100)
    manifest_parser.add_argument("--threads", type=int, default=8)
    manifest_parser.add_argument(
        "--check-validation",
        action="store_true",
        help=(
            "stat every expected validation XML; normally omit this because "
            "each ingest task records exact availability"
        ),
    )

    batch_parser = commands.add_parser(
        "batch", help="ingest one fixed-size slice of a manifest"
    )
    batch_parser.add_argument("manifest", type=Path)
    batch_parser.add_argument("output_root", type=Path)
    batch_parser.add_argument("--batch-size", type=int, required=True)
    batch_parser.add_argument("--batch-index", type=int, default=_default_batch_index())
    batch_parser.add_argument(
        "--cif-root",
        type=Path,
        default=os.environ.get(PDB_NEXTGEN_ROOT_ENV),
        required=PDB_NEXTGEN_ROOT_ENV not in os.environ,
    )
    batch_parser.add_argument(
        "--validation-root",
        type=Path,
        default=os.environ.get(VALIDATION_ROOT_ENV),
        required=VALIDATION_ROOT_ENV not in os.environ,
    )
    batch_parser.add_argument("--force", action="store_true")
    return parser


def main() -> None:
    """Run manifest discovery or one resumable ingest batch."""
    args = build_parser().parse_args()
    if args.command == "manifest":
        entries = discover_entries(
            args.cif_root,
            args.validation_root,
            check_validation=args.check_validation,
            threads=args.threads,
        )
        if not entries:
            raise FileNotFoundError(
                f"no NextGen entry CIFs found under {args.cif_root}"
            )
        summary = write_manifest(
            entries=entries,
            output_path=args.output_path,
            batch_size=args.batch_size,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return

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
