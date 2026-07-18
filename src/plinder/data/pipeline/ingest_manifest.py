# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Build a size-balanced manifest for per-entry V3 ingest."""

from __future__ import annotations

import argparse
import heapq
import json
import math
import os
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Any


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
) -> list[EntryInput]:
    """Discover unique NextGen CIFs and optionally check validation reports."""
    if threads < 1:
        raise ValueError("threads must be positive")
    directory_pattern = re.compile(r"pdb_0000([0-9a-z]{4})")

    def discover_code_directory(code_directory: Path) -> list[EntryInput]:
        discovered = []
        with os.scandir(code_directory) as directory_entries:
            for directory_entry in directory_entries:
                match = directory_pattern.fullmatch(directory_entry.name)
                if match is None or not directory_entry.is_dir(follow_symlinks=False):
                    continue
                pdb_id = match.group(1)
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

    code_directories = [path for path in cif_root.iterdir() if path.is_dir()]
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif_root", type=Path)
    parser.add_argument("validation_root", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--check-validation",
        action="store_true",
        help=(
            "stat every expected validation XML; normally omit this because "
            "each ingest task records exact availability"
        ),
    )
    args = parser.parse_args()

    entries = discover_entries(
        args.cif_root,
        args.validation_root,
        check_validation=args.check_validation,
        threads=args.threads,
    )
    if not entries:
        raise FileNotFoundError(f"no NextGen entry CIFs found under {args.cif_root}")
    summary = write_manifest(
        entries=entries,
        output_path=args.output_path,
        batch_size=args.batch_size,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
