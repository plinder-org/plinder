# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Build a clean public release tree from a working ingest."""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path

from plinder.core.release import RELEASE_PATHS, RELEASE_TABLES

_FILE_ARTIFACTS = (
    *(table["artifact"] for table in RELEASE_TABLES.values()),
    "ligand_similarity_scores",
    "interface_similarity_scores",
    "interface_half_similarity_scores",
)
_PARQUET_DIRECTORIES = (
    "protein_similarity_scores",
    "monomer_similarity_scores",
    "ligand_archives",
    "ligand_scores",
    "alignment_cigars",
)
_SAMPLING_DIRECTORIES = ("ligand_sampling", "interface_sampling")
_SEARCH_DATABASES = (
    "holo_mmseqs",
    "holo_foldseek",
    "monomer_mmseqs",
    "monomer_foldseek",
)


def _public_files(source: Path) -> list[Path]:
    files = {Path(RELEASE_PATHS[name]) for name in _FILE_ARTIFACTS}
    for name in _PARQUET_DIRECTORIES:
        directory = source / RELEASE_PATHS.get(name, name)
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        files.update(path.relative_to(source) for path in directory.rglob("*.parquet"))
    for name in _SAMPLING_DIRECTORIES:
        directory = source / RELEASE_PATHS.get(name, name)
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        for cover in ("set_cover", "directed_set_cover"):
            files.update(
                path.relative_to(source)
                for path in (directory / cover).glob("metric=*/threshold=*.parquet")
            )

    database_root = source / RELEASE_PATHS["search_databases"]
    for name in _SEARCH_DATABASES:
        directory = database_root / name
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        files.update(
            path.relative_to(source)
            for path in directory.rglob("*")
            if path.is_file() and path.suffix != ".json"
        )
    shadowed = database_root / "shadowed_entries.parquet"
    if shadowed.is_file():
        files.add(shadowed.relative_to(source))
    overlay = database_root / "weekly_delta"
    if overlay.is_dir():
        files.update(
            path.relative_to(source)
            for path in overlay.rglob("*")
            if path.is_file() and path.suffix != ".json"
        )

    missing = [str(source / path) for path in files if not (source / path).is_file()]
    if missing:
        raise FileNotFoundError(f"missing public release files: {missing[:10]}")
    included = {source / item for item in files}
    for relative in files:
        path = source / relative
        if path.is_symlink():
            if Path(os.readlink(path)).is_absolute():
                raise ValueError(f"absolute link in public release: {path}")
            if path.resolve() not in included:
                raise ValueError(
                    f"search database link leaves the public release: {path}"
                )
    return sorted(files)


def prepare_public_release(
    source: Path, destination: Path, *, replace: bool = False, copy: bool = False
) -> int:
    """Stage only public files beside an ingest; return the number of files."""
    source = source.resolve(strict=True)
    if destination.is_symlink():
        raise ValueError("public release destination cannot be a symlink")
    destination = destination.parent.resolve() / destination.name
    if (
        source == destination
        or source in destination.parents
        or destination in source.parents
    ):
        raise ValueError("working and public release directories must be separate")
    if destination.exists() and not replace:
        raise FileExistsError(destination)
    if destination.exists() and not destination.is_dir():
        raise NotADirectoryError(destination)
    files = _public_files(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}-", dir=destination.parent)
    )
    previous = staging.with_name(staging.name + "-previous")
    try:
        for relative in files:
            original = source / relative
            target = staging / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            if original.is_symlink():
                target.symlink_to(os.readlink(original))
            elif copy:
                shutil.copy2(original, target)
            else:
                os.link(original, target)
        if destination.exists():
            destination.rename(previous)
        try:
            staging.rename(destination)
        except OSError:
            if previous.exists():
                previous.rename(destination)
            raise
        if previous.exists():
            shutil.rmtree(previous)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return len(files)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("working", type=Path, help="completed working ingest")
    parser.add_argument("public", type=Path, help="separate public release directory")
    parser.add_argument(
        "--replace", action="store_true", help="replace an existing public tree"
    )
    parser.add_argument(
        "--copy", action="store_true", help="copy bytes instead of hard-linking"
    )
    args = parser.parse_args()
    count = prepare_public_release(
        args.working, args.public, replace=args.replace, copy=args.copy
    )
    print(f"Prepared {count} public files in {args.public}")


if __name__ == "__main__":
    main()
