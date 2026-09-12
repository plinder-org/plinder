"""Small file helpers shared by release builders and readers."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from errno import EACCES, EPERM, EXDEV
from pathlib import Path
from typing import Any
from uuid import uuid4


def link_or_copy_file(source: Path, destination: Path) -> None:
    """Reuse immutable files on one filesystem, copying across filesystems."""
    if source.is_symlink():
        destination.symlink_to(source.readlink())
        return
    try:
        os.link(source, destination)
    except OSError as exc:
        if exc.errno not in {EACCES, EPERM, EXDEV}:
            raise
        shutil.copy2(source, destination)


def file_sha256(path: Path) -> str:
    """Hash a file without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_json_cache(path: Path) -> dict[str, Any] | None:
    """Read an optional cache marker; missing or malformed markers are cache misses."""
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    """Replace a JSON file only after writing its complete contents."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
