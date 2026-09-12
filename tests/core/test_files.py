import hashlib
import json
import os
from errno import EXDEV
from pathlib import Path

import pytest

from plinder.core.utils.files import (
    file_sha256,
    link_or_copy_file,
    read_json_cache,
    write_json_atomic,
)


def test_file_sha256_streams_multiple_chunks(tmp_path):
    content = b"plinder" * 400_000
    path = tmp_path / "source"
    path.write_bytes(content)
    assert file_sha256(path) == hashlib.sha256(content).hexdigest()


@pytest.mark.parametrize("content", [None, "broken", "[]", "null", '"text"'])
def test_unusable_json_cache_is_a_miss(tmp_path, content):
    path = tmp_path / "marker.json"
    if content is not None:
        path.write_text(content)
    assert read_json_cache(path) is None


def test_write_json_atomic_replaces_marker(tmp_path):
    path = tmp_path / "nested" / "marker.json"
    write_json_atomic(path, {"state": "pending"})
    write_json_atomic(path, {"state": "complete"})
    assert read_json_cache(path) == {"state": "complete"}
    assert list(path.parent.iterdir()) == [path]


def test_failed_json_install_preserves_marker(tmp_path, monkeypatch):
    path = tmp_path / "marker.json"
    write_json_atomic(path, {"state": "complete"})

    def fail_replace(self, target):
        raise OSError("install failed")

    monkeypatch.setattr(Path, "replace", fail_replace)
    with pytest.raises(OSError, match="install failed"):
        write_json_atomic(path, {"state": "pending"})
    assert json.loads(path.read_text()) == {"state": "complete"}
    assert list(tmp_path.iterdir()) == [path]


def test_link_or_copy_file_copies_across_filesystems(tmp_path, monkeypatch):
    source, target = tmp_path / "source", tmp_path / "target"
    source.write_bytes(b"immutable archive")

    def cross_device(*args):
        raise OSError(EXDEV, "cross-device link")

    monkeypatch.setattr(os, "link", cross_device)
    link_or_copy_file(source, target)
    assert target.read_bytes() == source.read_bytes()
    assert target.stat().st_ino != source.stat().st_ino


def test_link_or_copy_file_does_not_overwrite_existing_file(tmp_path):
    source, target = tmp_path / "source", tmp_path / "target"
    source.write_bytes(b"new")
    target.write_bytes(b"old")
    with pytest.raises(FileExistsError):
        link_or_copy_file(source, target)
    assert target.read_bytes() == b"old"
