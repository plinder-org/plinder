# Copyright (c) 2026, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Plan entry updates from a local NextGen snapshot and PDB obsolete.dat."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Collection

import pandas as pd

from plinder.data.pipeline.ingest import load_manifest, normalize_pdb_id

REVISION_COLUMNS = ["source_mmcif_major_revision", "source_mmcif_minor_revision"]
SNAPSHOT_COLUMNS = ["entry_pdb_id", "last_modified", *REVISION_COLUMNS, "error"]


def read_obsolete(path: Path) -> dict[str, list[str]]:
    """Read obsolete IDs and all their successors; successors are not aliases."""
    replacements: dict[str, set[str]] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if not fields or fields[0] != "OBSLTE":
            continue
        if len(fields) < 3:
            raise ValueError(f"incomplete OBSLTE record at {path}:{number}")
        pdb_id = normalize_pdb_id(fields[2])
        replacements.setdefault(pdb_id, set()).update(
            normalize_pdb_id(value) for value in fields[3:]
        )
    return {key: sorted(values) for key, values in sorted(replacements.items())}


def _holdings(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt") as handle:
        contents = json.load(handle)
    return {
        normalize_pdb_id(key.removeprefix("pdb_0000")): value
        for key, value in contents.items()
    }


def _read_revision(path: Path) -> tuple[int, int]:
    from plinder.data.annotations.cif_utils import (
        get_mmcif_revision,
        read_mmcif_container,
    )

    return get_mmcif_revision(read_mmcif_container(path))


def read_snapshot(
    nextgen_root: Path,
    *,
    threads: int = 8,
    previous_snapshot: Path | None = None,
    pdb_ids: Collection[str] | None = None,
) -> pd.DataFrame:
    """Read revisions, reusing a prior inventory for unchanged holdings timestamps.

    ``nextgen_root`` contains ``holdings/`` and ``data/entries/divided/``.
    The first comparison reads the CIFs to establish their revisions. A saved
    inventory is only a read cache: it never proves an update was applied.
    """
    if threads < 1:
        raise ValueError("threads must be positive")
    root = nextgen_root.resolve(strict=True)
    files = _holdings(root / "holdings/current_file_holdings.json.gz")
    modified = _holdings(
        root / "holdings/released_structures_last_modified_dates.json.gz"
    )
    selected = (
        {normalize_pdb_id(value) for value in pdb_ids} if pdb_ids is not None else None
    )
    current = sorted(
        key
        for key, value in files.items()
        if value.get("mmcif") and (selected is None or key in selected)
    )
    cached = {}
    if previous_snapshot is not None:
        cache = pd.read_parquet(previous_snapshot, columns=SNAPSHOT_COLUMNS)
        cached = cache.set_index("entry_pdb_id", verify_integrity=True).to_dict(
            orient="index"
        )

    def read(pdb_id: str) -> dict[str, Any]:
        timestamp = modified.get(pdb_id)
        previous = cached.get(pdb_id)
        if (
            previous is not None
            and timestamp is not None
            and previous["last_modified"] == timestamp
            and pd.isna(previous["error"])
            and all(pd.notna(previous[column]) for column in REVISION_COLUMNS)
        ):
            return {"entry_pdb_id": pdb_id, **previous}
        row: dict[str, Any] = dict.fromkeys(SNAPSHOT_COLUMNS)
        row.update(entry_pdb_id=pdb_id, last_modified=timestamp)
        path = (
            root
            / "data/entries/divided"
            / pdb_id[1:3]
            / f"pdb_0000{pdb_id}"
            / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
        )
        try:
            major, minor = _read_revision(path)
            row.update(zip(REVISION_COLUMNS, (major, minor)))
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
        return row

    with ThreadPoolExecutor(max_workers=threads) as executor:
        rows = list(executor.map(read, current))
    return pd.DataFrame(rows, columns=SNAPSHOT_COLUMNS).astype(
        {column: "Int64" for column in REVISION_COLUMNS}
    )


def compare_entries(
    sources: pd.DataFrame,
    snapshot: pd.DataFrame,
    obsolete: dict[str, list[str]],
) -> pd.DataFrame:
    """Classify each entry without treating unexplained absences as removals."""
    old = sources.set_index("entry_pdb_id", verify_integrity=True)
    new = snapshot.set_index("entry_pdb_id", verify_integrity=True)
    rows = []
    for pdb_id in sorted(set(old.index) | set(new.index)):
        before = (
            tuple(old.loc[pdb_id, REVISION_COLUMNS]) if pdb_id in old.index else None
        )
        after = (
            tuple(new.loc[pdb_id, REVISION_COLUMNS]) if pdb_id in new.index else None
        )
        error = new.at[pdb_id, "error"] if after is not None else None
        reason = ""
        if pdb_id in obsolete and after is not None:
            action, reason = "blocked", "entry is both current and obsolete"
        elif pdb_id in obsolete:
            action = "obsolete"
        elif after is None:
            action, reason = "blocked", "absent from holdings but not in obsolete.dat"
        elif pd.notna(error):
            action, reason = "blocked", str(error)
        elif any(pd.isna(value) for value in (*(before or ()), *after)):
            action, reason = "blocked", "missing source revision"
        elif before is None:
            action = "added"
        elif after < before:
            action, reason = "blocked", "snapshot revision is older than the release"
        elif after > before:
            action = "revised"
        else:
            action = "unchanged"
        rows.append(
            {
                "pdb_id": pdb_id,
                "action": action,
                "reason": reason,
                "previous_major_revision": before[0] if before else None,
                "previous_minor_revision": before[1] if before else None,
                "current_major_revision": after[0] if after else None,
                "current_minor_revision": after[1] if after else None,
                "replacement_pdb_ids": obsolete.get(pdb_id, []),
            }
        )
    columns = [
        "pdb_id",
        "action",
        "reason",
        "previous_major_revision",
        "previous_minor_revision",
        "current_major_revision",
        "current_minor_revision",
        "replacement_pdb_ids",
    ]
    return pd.DataFrame(rows, columns=columns).astype(
        {column: "Int64" for column in columns if column.endswith("_revision")}
    )


def _signature(path: Path) -> dict[str, str]:
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def plan_update(
    data_dir: Path,
    *,
    nextgen_root: Path,
    obsolete_path: Path,
    output_dir: Path,
    threads: int = 8,
    previous_snapshot: Path | None = None,
    pdb_ids: Collection[str] | None = None,
) -> dict[str, Any]:
    """Write an update report and work lists into a new, separate directory.

    Nothing is applied to the release. Blocked plans are review reports, not
    executable work lists. A partial ``pdb_ids`` scope is useful for testing;
    omit it for a complete release comparison. Coordinate revisions determine
    changes here; validation-only or NextGen enrichment changes without a PDB
    revision need a separate annotation refresh.
    """
    data_dir = data_dir.resolve(strict=True)
    nextgen_root = nextgen_root.resolve(strict=True)
    obsolete_path = obsolete_path.resolve(strict=True)
    output_dir = output_dir.resolve()
    if output_dir == data_dir or data_dir in output_dir.parents:
        raise ValueError("write the update plan outside the existing release")
    if output_dir.exists():
        raise FileExistsError(output_dir)
    sources_path = data_dir / "index/entry_sources.parquet"
    inputs = {
        "entry_sources": sources_path,
        "obsolete": obsolete_path,
        "holdings": nextgen_root / "holdings/current_file_holdings.json.gz",
        "modified": nextgen_root
        / "holdings/released_structures_last_modified_dates.json.gz",
    }
    signatures = {name: _signature(path) for name, path in inputs.items()}
    sources = pd.read_parquet(sources_path, columns=["entry_pdb_id", *REVISION_COLUMNS])
    scope = (
        sorted({normalize_pdb_id(value) for value in pdb_ids})
        if pdb_ids is not None
        else None
    )
    if scope is not None:
        sources = sources.loc[sources.entry_pdb_id.isin(scope)]
    snapshot = read_snapshot(
        nextgen_root,
        threads=threads,
        previous_snapshot=previous_snapshot,
        pdb_ids=scope,
    )
    entries = compare_entries(sources, snapshot, read_obsolete(obsolete_path))
    counts = dict(sorted(Counter(entries.action).items()))
    changed = bool(entries.action.isin(["added", "revised", "obsolete"]).any())
    status = (
        "blocked" if counts.get("blocked") else "ready" if changed else "no_changes"
    )
    plan = {
        "status": status,
        "data_dir": str(data_dir),
        "nextgen_root": str(nextgen_root),
        "scope": scope,
        "inputs": signatures,
        "counts": counts,
        "work": {
            "ingest": "ingest.txt",
            "remove": "remove.txt",
            "invalidate_queries_and_targets": "invalidate.txt",
            "searches": [
                {"queries": "ingest.txt", "targets": "active.txt"},
                {"queries": "unchanged.txt", "targets": "ingest.txt"},
            ],
            "refresh": [
                "entry tables and canonical ligand archives",
                "protein clusters, search databases, and representative membership",
                "ligand fingerprints, chemical similarities, CCD matches, and MMP pairs",
                "ligand-pocket, PLI, SuCOS, interface, and linked-apo scores",
                "set covers, cluster summaries, apo links, and final release tables",
            ]
            if changed
            else [],
        },
        "notes": [
            "Do not apply a blocked plan; resolve its source errors first.",
            "Remove stale rows on both query and target sides before replacing scores.",
            "Search directions are candidate work, not a complete search-cache repair plan: "
            "changed representatives and capped hit lists can require full-query searches.",
            "Recompute covers globally; their representatives and labels can change.",
            "This plan compares PDB coordinate revisions, not independent changes to "
            "validation reports, NextGen enrichment, CCD, or annotation settings.",
        ],
    }
    # Check small catalogues again in case a source changed during revision reads.
    if signatures != {name: _signature(path) for name, path in inputs.items()}:
        raise RuntimeError(
            "source catalogues changed while planning; use a fixed snapshot"
        )
    output_dir.mkdir(parents=True)
    snapshot.to_parquet(output_dir / "snapshot.parquet", index=False)
    entries.to_parquet(output_dir / "entries.parquet", index=False)
    entries.to_csv(output_dir / "entries.tsv", sep="\t", index=False)
    if status != "blocked":
        for name, actions in {
            "ingest": ["added", "revised"],
            "remove": ["obsolete"],
            "invalidate": ["revised", "obsolete"],
            "unchanged": ["unchanged"],
            "active": ["added", "revised", "unchanged"],
        }.items():
            ids = entries.loc[entries.action.isin(actions), "pdb_id"]
            (output_dir / f"{name}.txt").write_text(
                "".join(f"{value}\n" for value in ids)
            )
    # Written last so an interrupted report is not mistaken for a finished plan.
    (output_dir / "plan.json").write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n"
    )
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--nextgen-root", type=Path, required=True)
    parser.add_argument("--obsolete", type=Path, required=True)
    parser.add_argument("--previous-snapshot", type=Path)
    parser.add_argument("--pdb-manifest", type=Path)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    plan = plan_update(
        args.data_dir,
        nextgen_root=args.nextgen_root,
        obsolete_path=args.obsolete,
        output_dir=args.output_dir,
        threads=args.threads,
        previous_snapshot=args.previous_snapshot,
        pdb_ids=load_manifest(args.pdb_manifest) if args.pdb_manifest else None,
    )
    print(json.dumps(plan, indent=2, sort_keys=True))
    if plan["status"] == "blocked":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
