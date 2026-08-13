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
import pyarrow.parquet as pq

PDB_NEXTGEN_ROOT_ENV = "PLINDER_PDB_NEXTGEN_ROOT"
VALIDATION_ROOT_ENV = "PLINDER_VALIDATION_ROOT"
REQUIRED_REFERENCE_FILES = (
    Path("dbs/cofactors/cofactors.json"),
    Path("dbs/affinity/affinity.json"),
)
INGEST_MODES = ("all", "ligands", "interfaces")

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
        path = Path(value).expanduser() if value else Path(fallback)
        if not path.is_absolute():
            path = data_dir / path
        return path.resolve()

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


def normalize_ingest_mode(value: str) -> str:
    """Validate one explicit entry-ingest mode."""
    mode = str(value).strip().lower()
    if mode not in INGEST_MODES:
        raise ValueError(f"ingest mode must be one of {INGEST_MODES}: {value!r}")
    return mode


def entry_metrics_paths(output_root: Path, pdb_id: str) -> tuple[Path, Path]:
    """Return the sharded metrics path followed by the legacy flat path."""
    pdb_id = normalize_pdb_id(pdb_id)
    filename = f"ingest-one-{pdb_id}.json"
    metrics_root = output_root / "metrics"
    return metrics_root / pdb_id[1:3] / filename, metrics_root / filename


def interface_metrics_path(output_root: Path, pdb_id: str) -> Path:
    """Return the independent completion marker for interface-only ingest."""
    pdb_id = normalize_pdb_id(pdb_id)
    return output_root / "metrics" / pdb_id[1:3] / f"ingest-interfaces-{pdb_id}.json"


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
    expected_interface_min_residues: int | None = None,
    expected_annotate_prodigy: bool | None = None,
    expected_ingest_mode: str | None = None,
) -> bool:
    """Return whether a prior per-entry run completed atomically."""
    if not metrics_path.is_file():
        return False
    try:
        metrics = json.loads(metrics_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if metrics.get("status") != "complete":
        return False
    if expected_ingest_mode is not None and (
        metrics.get("mode", "all") != expected_ingest_mode
    ):
        return False
    if (
        expected_annotate_prodigy is not None
        and bool(metrics.get("interface_annotate_prodigy")) != expected_annotate_prodigy
    ):
        return False
    counts = metrics.get("counts", {})
    annotation_rows = int(counts.get("annotation_rows", 0))
    sidecars = {
        "entry_chains": entry_directory / "entry_chains.parquet",
        "entry_biounit_chains": entry_directory / "entry_biounit_chains.parquet",
        "entry_metadata": entry_directory / "entry_metadata.parquet",
        "interfaces": entry_directory / "interfaces.parquet",
        "entry_source": entry_directory / "entry_source.parquet",
    }
    if not all(path.is_file() for path in sidecars.values()):
        return False
    try:
        from plinder.data.annotations.interface_utils import (
            INTERFACE_ANNOTATION_SCHEMA,
            min_interface_residues_from_schema,
        )

        interface_min_residues = min_interface_residues_from_schema(
            pq.read_schema(sidecars["interfaces"])
        )
    except (OSError, TypeError, ValueError):
        return False
    if (
        expected_interface_min_residues is not None
        and interface_min_residues != expected_interface_min_residues
    ):
        return False
    required_columns = {
        sidecars["entry_chains"]: {
            "chain_receptor_type",
            "chain_is_ligand_like",
        },
        sidecars["entry_biounit_chains"]: BIOUNIT_CONTACT_COLUMNS
        | {
            "entry_pdb_id",
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
        },
        sidecars["entry_metadata"]: {"entry_pdb_id"},
        sidecars["interfaces"]: set(INTERFACE_ANNOTATION_SCHEMA.names),
        sidecars["entry_source"]: {"entry_pdb_id"},
    }
    if annotation_rows:
        if not entry_parquet.is_file() or not ligand_parquet.is_file():
            return False
        required_columns.update(
            {
                entry_parquet: {"system_receptor_type"},
                ligand_parquet: {"ligand_id", "ligand_is_3d_score_able"},
            }
        )
    try:
        schemas_are_complete = all(
            columns.issubset(pq.read_schema(path).names)
            for path, columns in required_columns.items()
        )
        return schemas_are_complete and _biounit_contacts_are_valid(
            sidecars["entry_biounit_chains"]
        )
    except Exception:
        return False


BIOUNIT_CONTACT_COLUMNS = {
    "chain_num_contacting_ions",
    "chain_num_contacting_artifacts",
    "chain_num_contacting_other_ligands",
}


def _biounit_contacts_are_valid(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        if not BIOUNIT_CONTACT_COLUMNS.issubset(pq.read_schema(path).names):
            return False
        table = pq.read_table(path, columns=sorted(BIOUNIT_CONTACT_COLUMNS))
        return all(
            value is not None and value >= 0
            for column in table.column_names
            for value in table[column].to_pylist()
        )
    except (OSError, TypeError, ValueError):
        return False


def _get_annotation_class() -> Any:
    """Import the data-generation stack only for actual entry annotation."""
    from plinder.data.get_system_annotations import GetPlinderAnnotation

    return GetPlinderAnnotation


def _save_ligand_batch(**kwargs: Any) -> None:
    """Import ligand annotation code only for actual entry annotation."""
    from plinder.data.pipeline.utils import save_ligand_batch

    save_ligand_batch(**kwargs)


def completed_entry_metrics(
    output_root: Path,
    pdb_id: str,
    *,
    expected_interface_min_residues: int | None = None,
    expected_annotate_prodigy: bool | None = None,
    expected_ingest_mode: str | None = None,
) -> Path | None:
    """Return the metrics file when one V3 entry has a complete output set."""
    for metrics_path in entry_metrics_paths(output_root, pdb_id):
        if not metrics_path.is_file():
            continue
        try:
            metrics = json.loads(metrics_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if expected_ingest_mode is not None and (
            metrics.get("mode", "all") != expected_ingest_mode
        ):
            continue
        if metrics.get("status") == "skipped_no_ligands":
            entry_directory = metrics.get("outputs", {}).get("entry_directory")
            if (
                expected_ingest_mode == "ligands"
                and entry_directory
                and _biounit_contacts_are_valid(
                    Path(entry_directory) / "entry_biounit_chains.parquet"
                )
            ):
                return metrics_path
            continue
        if metrics.get("status") == "skipped_no_systems":
            # A pre-interface-ingest skip may actually contain a protein-only
            # interface and must be reconsidered. New skips explicitly record
            # the zero interface count.
            if "interface_rows" not in metrics.get("counts", {}):
                continue
            try:
                stored_interface_min_residues = int(metrics["interface_min_residues"])
            except (KeyError, TypeError, ValueError):
                continue
            if expected_annotate_prodigy is not None and (
                bool(metrics.get("interface_annotate_prodigy"))
                != expected_annotate_prodigy
            ):
                continue
            if (
                expected_interface_min_residues is None
                or stored_interface_min_residues == expected_interface_min_residues
            ):
                return metrics_path
            continue
        outputs = metrics.get("outputs", {})
        entry_directory = outputs.get("entry_directory")
        if not entry_directory:
            continue
        entry_parquet = outputs.get("entry_parquet") or (
            output_root / "raw_entries" / pdb_id[1:3] / f"{pdb_id}.parquet"
        )
        ligand_parquet = outputs.get("ligand_parquet") or (
            output_root / "ligands" / f"{pdb_id}.parquet"
        )
        if _entry_outputs_complete(
            metrics_path=metrics_path,
            entry_parquet=Path(entry_parquet),
            entry_directory=Path(entry_directory),
            ligand_parquet=Path(ligand_parquet),
            expected_interface_min_residues=expected_interface_min_residues,
            expected_annotate_prodigy=expected_annotate_prodigy,
            expected_ingest_mode=expected_ingest_mode,
        ):
            return metrics_path
    return None


def completed_interface_metrics(
    output_root: Path,
    pdb_id: str,
    *,
    expected_interface_min_residues: int | None = None,
    expected_annotate_prodigy: bool | None = None,
) -> Path | None:
    """Return the marker when interface-only outputs are complete and current."""
    pdb_id = normalize_pdb_id(pdb_id)
    metrics_path = interface_metrics_path(output_root, pdb_id)
    if not metrics_path.is_file():
        return None
    try:
        metrics = json.loads(metrics_path.read_text())
        stored_cutoff = int(metrics["interface_min_residues"])
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None
    if (
        expected_interface_min_residues is not None
        and stored_cutoff != expected_interface_min_residues
    ):
        return None
    if (
        expected_annotate_prodigy is not None
        and bool(metrics.get("interface_annotate_prodigy")) != expected_annotate_prodigy
    ):
        return None
    status = metrics.get("status")
    if status == "skipped_no_interfaces":
        return metrics_path
    if status != "complete":
        return None
    interface_path = (
        output_root / "raw_entries" / pdb_id[1:3] / pdb_id / "interfaces.parquet"
    )
    if not interface_path.is_file():
        return None
    try:
        from plinder.data.annotations.interface_utils import (
            INTERFACE_ANNOTATION_SCHEMA,
            min_interface_residues_from_schema,
        )

        schema = pq.read_schema(interface_path)
        if min_interface_residues_from_schema(schema) != stored_cutoff:
            return None
        if not set(INTERFACE_ANNOTATION_SCHEMA.names).issubset(schema.names):
            return None
        if pq.ParquetFile(interface_path).metadata.num_rows != int(
            metrics.get("counts", {}).get("interface_rows", -1)
        ):
            return None
    except (OSError, TypeError, ValueError):
        return None
    return metrics_path


def _ingest_interfaces(
    *,
    pdb_id: str,
    output_root: Path,
    cif_file: Path,
    validation_file: Path,
    raw_entry_root: Path,
    entry_parquet: Path,
    entry_directory: Path,
    ligand_parquet: Path,
    force: bool,
    annotation_cfg: Mapping[str, Any] | None,
    entry_cfg: Mapping[str, Any] | None,
    interface_cfg: Mapping[str, Any] | None,
) -> Path:
    """Replace only interface-derived outputs while preserving ligand assets."""
    from plinder.data.annotations.interface_utils import (
        DEFAULT_MIN_INTERFACE_RESIDUES,
    )

    interface_options = dict(interface_cfg or {})
    expected_cutoff = int(
        interface_options.get("min_interface_residues", DEFAULT_MIN_INTERFACE_RESIDUES)
    )
    annotate_prodigy = bool(interface_options.get("annotate_prodigy", True))
    completed = completed_interface_metrics(
        output_root,
        pdb_id,
        expected_interface_min_residues=expected_cutoff,
        expected_annotate_prodigy=annotate_prodigy,
    )
    if completed is not None and not force:
        raise FileExistsError(
            f"interface output already exists for {pdb_id}; pass --force to replace it"
        )
    if entry_parquet.is_file() != ligand_parquet.is_file():
        raise FileNotFoundError(
            f"cannot preserve incomplete ligand outputs for {pdb_id}: "
            f"annotation={entry_parquet.is_file()}, ligand={ligand_parquet.is_file()}"
        )
    has_ligands = entry_parquet.is_file()
    if not has_ligands and entry_directory.exists():
        # This directory is owned solely by an earlier interface-only attempt,
        # so a partial write can be rebuilt without touching ligand assets.
        shutil.rmtree(entry_directory)

    metrics_path = interface_metrics_path(output_root, pdb_id)
    timings: list[dict[str, float | str]] = []
    summary: dict[str, Any] = {
        "pdb_id": pdb_id,
        "status": "running",
        "mode": "interfaces",
        "interface_min_residues": expected_cutoff,
        "interface_annotate_prodigy": annotate_prodigy,
        "inputs": {
            "mmcif": str(cif_file),
            "validation_xml": str(validation_file),
            "validation_xml_exists": validation_file.is_file(),
        },
        "outputs": {
            "entry_parquet": str(entry_parquet) if entry_parquet.is_file() else None,
            "entry_directory": str(entry_directory),
            "interface_parquet": str(entry_directory / "interfaces.parquet"),
            "ligand_parquet": str(ligand_parquet) if ligand_parquet.is_file() else None,
        },
        "timings": timings,
    }
    total_started = time.perf_counter()
    try:
        raw_entry_root.mkdir(parents=True, exist_ok=True)

        def annotate_interfaces() -> Any:
            annotation_options = dict(annotation_cfg or {})
            entry_options = dict(entry_cfg or {})
            entry_options.pop("save_folder", None)
            if entry_options:
                annotation_options["entry_cfg"] = entry_options
            annotation_options["interface_cfg"] = interface_options
            return _get_annotation_class()(
                cif_file,
                validation_file,
                save_folder=raw_entry_root,
                **annotation_options,
            ).annotate_interfaces()

        interface_table = _run_timed(
            "annotate_interfaces", annotate_interfaces, timings
        )
        interface_rows = int(interface_table.num_rows)
        if not interface_rows and not has_ligands:
            shutil.rmtree(entry_directory, ignore_errors=True)
            summary["status"] = "skipped_no_interfaces"
            summary["outputs"]["entry_directory"] = None
            summary["outputs"]["interface_parquet"] = None
        else:
            summary["status"] = "complete"
        summary["counts"] = {
            "interface_rows": interface_rows,
            "prodigy_annotated_rows": int(
                interface_table.column("prodigy_is_annotated").to_numpy().sum()
            ),
            "preserved_ligand_annotation_rows": (
                pq.ParquetFile(entry_parquet).metadata.num_rows if has_ligands else 0
            ),
        }
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


def _clear_ligand_outputs(
    *,
    entry_parquet: Path,
    entry_directory: Path,
    ligand_parquet: Path,
) -> None:
    """Remove ligand-derived outputs without touching interface sidecars."""
    entry_parquet.unlink(missing_ok=True)
    ligand_parquet.unlink(missing_ok=True)
    shutil.rmtree(entry_directory / "ligand_files", ignore_errors=True)


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
    interface_cfg: Mapping[str, Any] | None = None,
    mode: str = "all",
) -> Path:
    """Ingest ligands, protein interfaces, or both for one PDB entry."""
    pdb_id = normalize_pdb_id(pdb_id)
    mode = normalize_ingest_mode(mode)
    include_interfaces = mode == "all"
    output_root = output_root.resolve()
    cif_file, validation_file = resolve_entry_paths(
        pdb_id,
        cif_root=cif_root,
        validation_root=validation_root,
    )
    if not cif_file.is_file():
        raise FileNotFoundError(f"missing NextGen mmCIF: {cif_file}")
    if check_references and mode != "interfaces":
        check_reference_data(output_root)

    code = pdb_id[1:3]
    raw_entry_root = output_root / "raw_entries" / code
    entry_parquet = raw_entry_root / f"{pdb_id}.parquet"
    entry_directory = raw_entry_root / pdb_id
    ligand_parquet = output_root / "ligands" / f"{pdb_id}.parquet"
    metrics_path, legacy_metrics_path = entry_metrics_paths(output_root, pdb_id)
    if mode == "interfaces":
        return _ingest_interfaces(
            pdb_id=pdb_id,
            output_root=output_root,
            cif_file=cif_file,
            validation_file=validation_file,
            raw_entry_root=raw_entry_root,
            entry_parquet=entry_parquet,
            entry_directory=entry_directory,
            ligand_parquet=ligand_parquet,
            force=force,
            annotation_cfg=annotation_cfg,
            entry_cfg=entry_cfg,
            interface_cfg=interface_cfg,
        )
    from plinder.data.annotations.interface_utils import (
        DEFAULT_MIN_INTERFACE_RESIDUES,
    )

    expected_interface_min_residues = int(
        (interface_cfg or {}).get(
            "min_interface_residues", DEFAULT_MIN_INTERFACE_RESIDUES
        )
    )
    expected_annotate_prodigy = include_interfaces and bool(
        (interface_cfg or {}).get("annotate_prodigy", True)
    )
    complete = any(
        _entry_outputs_complete(
            metrics_path=candidate,
            entry_parquet=entry_parquet,
            entry_directory=entry_directory,
            ligand_parquet=ligand_parquet,
            expected_interface_min_residues=(
                expected_interface_min_residues if include_interfaces else None
            ),
            expected_annotate_prodigy=(
                expected_annotate_prodigy if include_interfaces else None
            ),
            expected_ingest_mode=mode,
        )
        for candidate in (metrics_path, legacy_metrics_path)
    )
    if complete and not force:
        raise FileExistsError(
            f"output already exists: {entry_parquet}; pass --force to replace it"
        )
    has_partial_outputs = entry_parquet.exists() or ligand_parquet.exists()
    if mode == "all":
        has_partial_outputs = has_partial_outputs or entry_directory.exists()
    else:
        has_partial_outputs = (
            has_partial_outputs or (entry_directory / "ligand_files").exists()
        )
    if force or has_partial_outputs:
        if has_partial_outputs and not force:
            print(f"clearing incomplete outputs before retrying {pdb_id}")
        clear = _clear_entry_outputs if mode == "all" else _clear_ligand_outputs
        clear(
            entry_parquet=entry_parquet,
            entry_directory=entry_directory,
            ligand_parquet=ligand_parquet,
        )

    timings: list[dict[str, float | str]] = []
    summary: dict[str, Any] = {
        "pdb_id": pdb_id,
        "status": "running",
        "mode": mode,
        "interface_min_residues": expected_interface_min_residues,
        "interface_annotate_prodigy": expected_annotate_prodigy,
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
    if annotation_cfg or entry_cfg or interface_cfg:
        summary["configuration"] = {
            "annotation": dict(annotation_cfg or {}),
            "entry": {
                key: value
                for key, value in dict(entry_cfg or {}).items()
                if key != "save_folder"
            },
            "interface": dict(interface_cfg or {}),
        }
    total_started = time.perf_counter()
    try:
        raw_entry_root.mkdir(parents=True, exist_ok=True)
        ligand_parquet.parent.mkdir(parents=True, exist_ok=True)

        def annotate() -> pd.DataFrame:
            annotation_options = dict(annotation_cfg or {})
            entry_options = dict(entry_cfg or {})
            entry_options.pop("save_folder", None)
            entry_options.pop("data_dir", None)
            if entry_options:
                annotation_options["entry_cfg"] = entry_options
            if interface_cfg:
                annotation_options["interface_cfg"] = dict(interface_cfg)
            annotator = _get_annotation_class()(
                cif_file,
                validation_file,
                save_folder=raw_entry_root,
                data_dir=output_root,
                **annotation_options,
            )
            annotation = (
                annotator.annotate()
                if include_interfaces
                else annotator.annotate(include_interfaces=False)
            )
            return annotation if annotation is not None else pd.DataFrame()

        annotation = _run_timed("annotate_entry", annotate, timings)
        interface_path = entry_directory / "interfaces.parquet"
        interface_rows = (
            pq.ParquetFile(interface_path).metadata.num_rows
            if interface_path.is_file()
            else 0
        )
        if annotation.empty and mode == "ligands":
            entry_parquet.unlink(missing_ok=True)
            ligand_parquet.unlink(missing_ok=True)
            shutil.rmtree(entry_directory / "ligand_files", ignore_errors=True)
            summary["outputs"]["entry_parquet"] = None
            summary["outputs"]["ligand_parquet"] = None
            if not entry_directory.exists():
                summary["outputs"]["entry_directory"] = None
            summary["counts"] = {
                "annotation_rows": 0,
                "interface_rows": interface_rows,
                "interface_rows_generated": 0,
                "systems": 0,
                "ligand_ids": 0,
                "canonical_ligand_sdfs": 0,
            }
            summary["status"] = "skipped_no_ligands"
        elif annotation.empty and interface_rows == 0:
            entry_parquet.unlink(missing_ok=True)
            ligand_parquet.unlink(missing_ok=True)
            shutil.rmtree(entry_directory / "ligand_files", ignore_errors=True)
            summary["outputs"]["entry_parquet"] = None
            summary["outputs"]["ligand_parquet"] = None
            summary["counts"] = {
                "annotation_rows": 0,
                "interface_rows": 0,
                "systems": 0,
                "ligand_ids": 0,
                "canonical_ligand_sdfs": 0,
            }
            summary["status"] = "complete"
        else:
            if not annotation.empty:

                def write_entry_parquet() -> None:
                    annotation.to_parquet(entry_parquet, index=False)

                _run_timed(
                    "write_entry_parquet",
                    write_entry_parquet,
                    timings,
                )

                def write_ligand_parquet() -> None:
                    _save_ligand_batch(
                        data_dir=output_root,
                        annotation=annotation,
                        output_path=ligand_parquet,
                    )

                _run_timed(
                    "write_ligand_parquet",
                    write_ligand_parquet,
                    timings,
                )
            else:
                entry_parquet.unlink(missing_ok=True)
                ligand_parquet.unlink(missing_ok=True)
                summary["outputs"]["entry_parquet"] = None
                summary["outputs"]["ligand_parquet"] = None

            ligand_sdfs = sorted((entry_directory / "ligand_files").glob("*.sdf"))
            summary["counts"] = {
                "annotation_rows": len(annotation),
                "interface_rows": interface_rows,
                "systems": (
                    int(annotation["system_id"].nunique())
                    if not annotation.empty
                    else 0
                ),
                "ligand_ids": (
                    int(annotation["ligand_id"].nunique())
                    if not annotation.empty
                    else 0
                ),
                "canonical_ligand_sdfs": len(ligand_sdfs),
            }
            if mode == "ligands":
                summary["counts"]["interface_rows_generated"] = 0
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
    interface_cfg: Mapping[str, Any] | None = None,
    mode: str = "all",
) -> tuple[Path, bool]:
    """Ingest PDB IDs sequentially and continue after per-entry failures."""
    mode = normalize_ingest_mode(mode)
    output_root = output_root.resolve()
    batch_prefix = {
        "all": "ingest-batch",
        "ligands": "ligand-ingest-batch",
        "interfaces": "interface-ingest-batch",
    }[mode]
    metrics_path = (
        output_root / "metrics" / (f"{batch_prefix}-{job_id}-{batch_index:05d}.json")
    )
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "job_id": job_id,
        "batch_index": batch_index,
        "pdb_ids": pdb_ids,
        "status": "running",
        "mode": mode,
        "entries": [],
    }
    started = time.perf_counter()
    had_failures = False
    from plinder.data.annotations.interface_utils import (
        DEFAULT_MIN_INTERFACE_RESIDUES,
    )

    expected_interface_min_residues = int(
        (interface_cfg or {}).get(
            "min_interface_residues", DEFAULT_MIN_INTERFACE_RESIDUES
        )
    )
    expected_annotate_prodigy = bool(
        (interface_cfg or {}).get("annotate_prodigy", True)
    )
    try:
        for pdb_id in pdb_ids:
            entry_started = time.perf_counter()
            completed = (
                None
                if force
                else (
                    completed_interface_metrics(
                        output_root,
                        pdb_id,
                        expected_interface_min_residues=(
                            expected_interface_min_residues
                        ),
                        expected_annotate_prodigy=expected_annotate_prodigy,
                    )
                    if mode == "interfaces"
                    else completed_entry_metrics(
                        output_root,
                        pdb_id,
                        expected_interface_min_residues=(
                            expected_interface_min_residues if mode == "all" else None
                        ),
                        expected_annotate_prodigy=(
                            expected_annotate_prodigy if mode == "all" else None
                        ),
                        expected_ingest_mode=mode,
                    )
                )
            )
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
                    interface_cfg=interface_cfg,
                    mode=mode,
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
    batch_parser.add_argument(
        "--interface-min-residues",
        type=int,
        default=os.environ.get("PLINDER_INTERFACE_MIN_RESIDUES"),
        help=(
            "minimum resolved contact residues required on each interface side; "
            "defaults to interface.min_interface_residues"
        ),
    )
    batch_parser.add_argument("--force", action="store_true")
    batch_parser.add_argument(
        "--mode",
        choices=INGEST_MODES,
        default="all",
        help=(
            "ingest ligands and interfaces (all), only ligands, or only "
            "interfaces while preserving ligand assets"
        ),
    )
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
    if args.interface_min_residues is not None and args.interface_min_residues < 1:
        raise ValueError("interface minimum residues must be positive")
    metrics_path, had_failures = ingest_pdb_batch(
        pdb_ids=pdb_ids,
        output_root=args.output_root,
        cif_root=args.cif_root,
        validation_root=args.validation_root,
        force=args.force,
        job_id=_default_job_id(),
        batch_index=args.batch_index,
        interface_cfg=(
            {"min_interface_residues": args.interface_min_residues}
            if args.interface_min_residues is not None
            else None
        ),
        mode=args.mode,
    )
    print(metrics_path.read_text(), end="")
    if had_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
