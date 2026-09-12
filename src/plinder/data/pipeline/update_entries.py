# Copyright (c) 2026, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Prepare updated entry tables in a separate, resumable workspace."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from shutil import copy2
from typing import Any, Mapping

import duckdb
import pandas as pd
import pyarrow as pa
from omegaconf import OmegaConf

from plinder.data.pipeline import collate, config
from plinder.data.pipeline.ingest import REQUIRED_REFERENCE_FILES, ingest_pdb_batch
from plinder.data.pipeline.updates import REVISION_COLUMNS, _signature

ENTRY_TABLES = {
    "annotation": ("annotation_table.parquet", "entry_pdb_id, system_id, ligand_id"),
    "system_validation": ("system_validation.parquet", "system_id"),
    "entry_chains": ("entry_chains.parquet", "entry_pdb_id, chain_asym_id"),
    "entry_biounit_chains": (
        "entry_biounit_chains.parquet",
        "entry_pdb_id, biounit_id, chain_instance",
    ),
    "entry_metadata": ("entry_metadata.parquet", "entry_pdb_id"),
    "interfaces": ("interface_annotation_table.parquet", "entry_pdb_id, system_id"),
    "entry_sources": ("entry_sources.parquet", "entry_pdb_id"),
}


def _tables(root: Path) -> dict[str, Path]:
    return {
        name: root / "index" / filename for name, (filename, _) in ENTRY_TABLES.items()
    }


def _file_stats(paths: Mapping[str, Path]) -> dict[str, dict[str, int]]:
    return {
        name: {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
        for name, path in paths.items()
    }


def _check_plan(plan_dir: Path) -> tuple[dict[str, Any], pd.DataFrame]:
    plan = json.loads((plan_dir / "plan.json").read_text())
    if plan["status"] != "ready":
        raise ValueError("entry updates require a ready plan with changes")
    for source in plan["inputs"].values():
        if _signature(Path(source["path"])) != source:
            raise ValueError(f"update source changed after planning: {source['path']}")
    for name in ("entries.parquet", "snapshot.parquet"):
        if plan.get("reports", {}).get(name) != _signature(plan_dir / name)["sha256"]:
            raise ValueError(
                f"update report changed or lacks a signature; regenerate the plan: {name}"
            )
    return plan, pd.read_parquet(plan_dir / "entries.parquet")


def _combine_tables(
    base: Path,
    incoming: Path | None,
    stage: Path,
    removed: list[str],
    *,
    threads: int,
    memory_limit: str,
    scratch_dir: Path | None,
) -> dict[str, int]:
    """Replace complete entries, including their system-validation rows."""
    stage.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect()
    try:
        collate._configure_duckdb(
            connection,
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
        )
        connection.register(
            "removed_entries",
            pa.table({"entry_pdb_id": pa.array(removed, type=pa.string())}),
        )
        counts = {}
        for name, (filename, order_by) in ENTRY_TABLES.items():
            connection.read_parquet(str(base / "index" / filename)).create_view(
                "old_rows", replace=True
            )
            condition = (
                "split_part(old_rows.system_id, '__', 1) = removed_entries.entry_pdb_id"
                if name == "system_validation"
                else "old_rows.entry_pdb_id = removed_entries.entry_pdb_id"
            )
            retained = f"SELECT old_rows.* FROM old_rows ANTI JOIN removed_entries ON {condition}"
            query = retained
            if incoming is not None:
                connection.read_parquet(str(incoming / "index" / filename)).create_view(
                    "new_rows", replace=True
                )
                query += " UNION ALL BY NAME SELECT * FROM new_rows"
            counts[name] = int(
                collate._fetch_scalar(connection, f"SELECT count(*) FROM ({query})")
            )
            collate._copy_query(
                connection,
                f"SELECT * FROM ({query}) ORDER BY {order_by}",
                stage / filename,
                row_group_size=100_000,
            )
        return counts
    finally:
        connection.close()


def apply_entry_update(
    plan_dir: Path,
    output_dir: Path,
    *,
    validation_root: Path,
    annotation_cfg: Mapping[str, Any],
    entry_cfg: Mapping[str, Any],
    interface_cfg: Mapping[str, Any],
    threads: int = 4,
    memory_limit: str = "8GB",
    scratch_dir: Path | None = None,
) -> dict[str, Any]:
    """Prepare entry tables; leave scores, archives, and clusters pending.

    Use the release's original annotation settings. ``threads`` limits table
    processing; entry annotation uses the existing sequential batch runner.
    The existing release is read-only. Failed entry batches can be resumed in
    the same workspace with the same plan and inputs.
    """
    if threads < 1:
        raise ValueError("threads must be positive")
    plan_dir = plan_dir.resolve(strict=True)
    plan, entries = _check_plan(plan_dir)
    base = Path(plan["data_dir"]).resolve(strict=True)
    output_dir = output_dir.resolve()
    validation_root = validation_root.resolve(strict=True)
    if output_dir == base or base in output_dir.parents:
        raise ValueError("write entry updates outside the existing release")
    base_tables = _tables(base)
    base_marker_path = base / "index" / collate.FINAL_MARKER_NAME
    base_marker = json.loads(base_marker_path.read_text())
    if base_marker["status"] != "complete":
        raise ValueError("finish the existing release before preparing an update")
    min_residues = int(base_marker["interface_min_residues"])
    if int(interface_cfg["min_interface_residues"]) != min_residues:
        raise ValueError("update interface threshold differs from the existing release")
    changed = entries.loc[entries.action.isin(["added", "revised"])].copy()
    remove = sorted(entries.loc[entries.action.isin(["revised", "obsolete"]), "pdb_id"])
    reference_files = (
        {str(path): base / path for path in REQUIRED_REFERENCE_FILES}
        if len(changed)
        else {}
    )
    binding = {
        "plan": _signature(plan_dir / "plan.json"),
        "base_tables": _file_stats(base_tables),
        "base_marker": _signature(base_marker_path),
        "validation_root": str(validation_root),
        "annotation": dict(annotation_cfg),
        "entry": dict(entry_cfg),
        "interface": dict(interface_cfg),
        "references": {
            name: _signature(path) for name, path in reference_files.items()
        },
    }
    marker_path = output_dir / "entry_update.json"
    prepared_paths = {
        **_tables(output_dir),
        "collation": output_dir / "index" / collate.FINAL_MARKER_NAME,
    }
    if output_dir.exists():
        if not marker_path.is_file():
            raise ValueError("output directory is not an entry-update workspace")
        previous: dict[str, Any] = json.loads(marker_path.read_text())
        if previous["inputs"] != binding:
            raise ValueError("update inputs changed; use a new workspace")
        if previous["status"] == collate.REPAIR_REQUIRED_STATUS:
            if previous["outputs"] != _file_stats(prepared_paths):
                raise ValueError("prepared entry tables changed; use a new workspace")
            return previous
    output_dir.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {
        "status": "preparing_entries",
        "inputs": binding,
        "base_release": str(base),
    }
    collate._write_json_atomic(marker_path, report)
    incoming = output_dir / ".incoming"
    if len(changed):
        for relative, source in reference_files.items():
            destination = incoming / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            copy2(source, destination)
        _, failed = ingest_pdb_batch(
            pdb_ids=changed.pdb_id.tolist(),
            output_root=incoming,
            cif_root=Path(plan["nextgen_root"]) / "data/entries/divided",
            validation_root=validation_root,
            annotation_cfg=annotation_cfg,
            entry_cfg=entry_cfg,
            interface_cfg=interface_cfg,
            mode="all",
        )
        if failed:
            raise RuntimeError(
                f"entry update failed; inspect {incoming / 'metrics'} and rerun"
            )
        collate.run_collation(
            incoming,
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir,
        )
        observed = pd.read_parquet(incoming / "index/entry_sources.parquet").set_index(
            "entry_pdb_id"
        )[REVISION_COLUMNS]
        expected = changed.set_index("pdb_id")[
            ["current_major_revision", "current_minor_revision"]
        ]
        expected.columns = REVISION_COLUMNS
        if (
            not observed.astype("Int64")
            .sort_index()
            .equals(expected.astype("Int64").sort_index())
        ):
            raise ValueError("processed entry revisions differ from the update plan")
    stage = output_dir / ".entry_tables"
    counts = _combine_tables(
        base,
        incoming if len(changed) else None,
        stage,
        remove,
        threads=threads,
        memory_limit=memory_limit,
        scratch_dir=scratch_dir,
    )
    staged_tables = {
        name: stage / filename for name, (filename, _) in ENTRY_TABLES.items()
    }
    validation = collate._validate_final_tables(
        staged_tables,
        expected_counts=counts,
        min_interface_residues=min_residues,
        threads=threads,
        memory_limit=memory_limit,
        scratch_dir=scratch_dir,
    )
    _check_plan(plan_dir)
    if binding["base_tables"] != _file_stats(base_tables) or binding[
        "base_marker"
    ] != _signature(base_marker_path):
        raise RuntimeError("existing release changed during entry processing")
    (output_dir / "index").mkdir(exist_ok=True)
    collate._install_final_tables_fail_closed(
        staged_tables,
        _tables(output_dir),
        output_dir / "index" / collate.FINAL_MARKER_NAME,
    )
    report.update(
        status=collate.REPAIR_REQUIRED_STATUS,
        incoming_entries=str(incoming) if len(changed) else None,
        ingested_pdb_ids=changed.pdb_id.tolist(),
        obsolete_pdb_ids=entries.loc[entries.action == "obsolete", "pdb_id"].tolist(),
        **validation,
    )
    collate._write_json_atomic(
        output_dir / "index" / collate.FINAL_MARKER_NAME,
        {
            "status": collate.REPAIR_REQUIRED_STATUS,
            "interface_min_residues": min_residues,
            "entry_update": str(marker_path),
            **validation,
        },
    )
    report["outputs"] = _file_stats(prepared_paths)
    collate._write_json_atomic(marker_path, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--validation-root", type=Path, required=True)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="original release's pipeline configuration",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-limit", default="8GB")
    parser.add_argument("--scratch-dir", type=Path)
    args = parser.parse_args()
    cfg = config.get_config(
        config=OmegaConf.load(args.config), config_args=[], cached=False
    )
    report = apply_entry_update(
        args.plan_dir,
        args.output_dir,
        validation_root=args.validation_root,
        annotation_cfg=dict(cfg.annotation),
        entry_cfg=dict(cfg.entry),
        interface_cfg=dict(cfg.interface),
        threads=args.threads,
        memory_limit=args.memory_limit,
        scratch_dir=args.scratch_dir,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
