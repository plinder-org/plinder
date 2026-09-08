# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from plinder.core import PlinderSystem
from plinder.core.utils import config
from plinder.core.utils import io as core_io
from plinder.core.utils.cpl import is_offline


class _Response:
    def __init__(self, content: bytes, status_code: int = 200) -> None:
        self.content = content
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int):
        yield from (
            self.content[offset : offset + chunk_size]
            for offset in range(0, len(self.content), chunk_size)
        )


@pytest.fixture
def release_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PLINDER_MOUNT", tmp_path.as_posix())
    monkeypatch.setenv("PLINDER_BUCKET", "plinder")
    monkeypatch.setenv("PLINDER_RELEASE", "test")
    monkeypatch.setenv("PLINDER_RELEASE_NUMBER", "")
    monkeypatch.delenv("PLINDER_OFFLINE", raising=False)
    monkeypatch.delenv("PLINDER_OFFLINE_MODE", raising=False)
    config._config._clear()
    yield Path(config.get_config().data.plinder_dir)
    config._config._clear()


@pytest.fixture
def source_manifest(release_cache):
    path = release_cache / "index" / "entry_sources.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "entry_pdb_id": ["2y4i", "9zzz"],
            "source_mmcif_major_revision": [1, 1],
            "source_mmcif_minor_revision": [5, 0],
        }
    ).to_parquet(path, index=False)
    return path


def test_source_mmcif_download_cache_and_system_resolution(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    calls = []

    def get(url, *, stream, timeout):
        calls.append((url, stream, timeout))
        return _Response(cif_2y4i.read_bytes())

    monkeypatch.setattr(core_io.requests, "get", get)

    path = core_io.get_pdb_mmcif("2Y4I", manifest_path=source_manifest)
    expected = release_cache / "source_mmcifs" / "y4" / "2y4i_v1-5.cif.gz"
    assert path == expected
    assert path.is_file()
    assert calls == [
        (
            "https://files-versioned.wwpdb.org/pdb_versioned/data/entries/"
            "y4/pdb_00002y4i/pdb_00002y4i_xyz_v1-5.cif.gz",
            True,
            (10, 120),
        )
    ]

    # A system ID resolves to the same valid cache entry without another request.
    assert (
        core_io.get_pdb_mmcif(
            "2y4i__1__1.B__1.E_1.F",
            manifest_path=source_manifest,
        )
        == expected
    )
    assert len(calls) == 1

    manifest_requests = []

    def fetch_release_artifact(_release, name, **parameters):
        manifest_requests.append((name, parameters))
        return source_manifest

    monkeypatch.setattr(
        core_io.PlinderRelease,
        "fetch",
        fetch_release_artifact,
    )
    assert (
        PlinderSystem(system_id="2y4i__1__1.B__1.E_1.F").source_mmcif_path == expected
    )
    assert manifest_requests == [("entry_sources", {})]
    assert len(calls) == 1


def test_source_mmcif_offline_mode_requires_valid_cache(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    monkeypatch.setenv("PLINDER_OFFLINE_MODE", "true")
    assert is_offline()
    expected = core_io.pdb_mmcif_cache_path(
        "2y4i",
        manifest_path=source_manifest,
    )

    with pytest.raises(FileNotFoundError, match="online node first") as exc:
        core_io.get_pdb_mmcif("2y4i", manifest_path=source_manifest)
    assert str(expected) in str(exc.value)
    assert not expected.parent.exists()

    expected.parent.mkdir(parents=True)
    expected.write_bytes(cif_2y4i.read_bytes())
    monkeypatch.setattr(
        core_io.requests,
        "get",
        lambda *args, **kwargs: pytest.fail("offline mode attempted a request"),
    )
    assert (
        core_io.get_pdb_mmcif(
            "2y4i",
            force_update=True,
            manifest_path=source_manifest,
        )
        == expected
    )


def test_source_mmcif_replaces_corrupt_online_cache(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    expected = core_io.pdb_mmcif_cache_path(
        "2y4i",
        manifest_path=source_manifest,
    )
    expected.parent.mkdir(parents=True)
    expected.write_bytes(b"not a gzip file")
    monkeypatch.setattr(
        core_io.requests,
        "get",
        lambda *args, **kwargs: _Response(cif_2y4i.read_bytes()),
    )

    assert core_io.get_pdb_mmcif("2y4i", manifest_path=source_manifest) == expected
    core_io._validate_gzipped_mmcif(expected, expected_revision=(1, 5))
    assert not list(expected.parent.glob("*.tmp"))


def test_source_mmcif_reports_missing_versioned_entry(
    release_cache, source_manifest, monkeypatch
):
    calls = []

    def get(*args, **kwargs):
        calls.append(args[0])
        return _Response(b"", status_code=404)

    monkeypatch.setattr(
        core_io.requests,
        "get",
        get,
    )
    monkeypatch.setattr(
        core_io,
        "sleep",
        lambda seconds: pytest.fail("404 response was retried"),
    )

    with pytest.raises(FileNotFoundError, match="revision 1.0"):
        core_io.get_pdb_mmcif("9zzz", manifest_path=source_manifest)
    assert len(calls) == 1
    assert not core_io.pdb_mmcif_cache_path(
        "9zzz",
        manifest_path=source_manifest,
    ).exists()


def test_download_pdb_mmcifs_deduplicates_system_ids(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    calls = []

    def get(url, *, stream, timeout):
        calls.append(url)
        return _Response(cif_2y4i.read_bytes())

    monkeypatch.setattr(core_io.requests, "get", get)
    paths = core_io.download_pdb_mmcifs(
        ["2y4i", "2y4i__1__1.B__1.E_1.F"],
        manifest_path=source_manifest,
    )

    assert set(paths) == {"2y4i"}
    assert paths["2y4i"].is_file()
    assert calls == [
        "https://files-versioned.wwpdb.org/pdb_versioned/data/entries/"
        "y4/pdb_00002y4i/pdb_00002y4i_xyz_v1-5.cif.gz"
    ]

    with pytest.raises(ValueError, match="max_workers"):
        core_io.download_pdb_mmcifs(
            ["2y4i"],
            max_workers=0,
            manifest_path=source_manifest,
        )


def test_source_mmcif_retries_transient_http_status(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    responses = [
        _Response(b"temporarily unavailable", status_code=503),
        _Response(cif_2y4i.read_bytes()),
    ]
    waits = []
    monkeypatch.setattr(
        core_io.requests,
        "get",
        lambda *args, **kwargs: responses.pop(0),
    )
    monkeypatch.setattr(core_io, "sleep", waits.append)

    path = core_io.get_pdb_mmcif("2y4i", manifest_path=source_manifest)

    assert path.is_file()
    assert responses == []
    assert waits == [1]


def test_batch_source_download_retries_timeout(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    calls = 0
    waits = []

    def get(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise core_io.requests.exceptions.Timeout("temporary timeout")
        return _Response(cif_2y4i.read_bytes())

    monkeypatch.setattr(core_io.requests, "get", get)
    monkeypatch.setattr(core_io, "sleep", waits.append)

    paths = core_io.download_pdb_mmcifs(
        ["2y4i"],
        manifest_path=source_manifest,
    )

    assert paths["2y4i"].is_file()
    assert calls == 2
    assert waits == [1]


def test_source_mmcif_transient_retries_are_bounded(
    release_cache, source_manifest, monkeypatch
):
    calls = 0
    waits = []

    def get(*args, **kwargs):
        nonlocal calls
        calls += 1
        return _Response(b"temporarily unavailable", status_code=503)

    monkeypatch.setattr(core_io.requests, "get", get)
    monkeypatch.setattr(core_io, "sleep", waits.append)

    with pytest.raises(RuntimeError, match="HTTP 503"):
        core_io.get_pdb_mmcif("2y4i", manifest_path=source_manifest)

    assert calls == core_io.SOURCE_MMCIF_DOWNLOAD_ATTEMPTS == 3
    assert waits == [1, 2]


def test_parser_failure_replaces_cache_online(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    expected = core_io.pdb_mmcif_cache_path(
        "2y4i",
        manifest_path=source_manifest,
    )
    expected.parent.mkdir(parents=True)
    expected.write_bytes(cif_2y4i.read_bytes())
    original_read = core_io.CIFFile.read
    read_calls = 0

    def read(handle):
        nonlocal read_calls
        read_calls += 1
        if read_calls == 1:
            raise core_io.DeserializationError("malformed cached CIF")
        return original_read(handle)

    monkeypatch.setattr(core_io.CIFFile, "read", read)
    monkeypatch.setattr(
        core_io.requests,
        "get",
        lambda *args, **kwargs: _Response(cif_2y4i.read_bytes()),
    )

    assert core_io.get_pdb_mmcif("2y4i", manifest_path=source_manifest) == expected
    assert read_calls == 2
    assert not list(expected.parent.glob("*.tmp"))


def test_parser_failure_is_missing_cache_offline(
    release_cache, source_manifest, cif_2y4i, monkeypatch
):
    expected = core_io.pdb_mmcif_cache_path(
        "2y4i",
        manifest_path=source_manifest,
    )
    expected.parent.mkdir(parents=True)
    expected.write_bytes(cif_2y4i.read_bytes())
    monkeypatch.setenv("PLINDER_OFFLINE", "true")

    def read(handle):
        raise core_io.DeserializationError("malformed cached CIF")

    monkeypatch.setattr(core_io.CIFFile, "read", read)
    monkeypatch.setattr(
        core_io.requests,
        "get",
        lambda *args, **kwargs: pytest.fail("offline mode attempted a request"),
    )

    with pytest.raises(FileNotFoundError, match="online node first"):
        core_io.get_pdb_mmcif("2y4i", manifest_path=source_manifest)


def test_source_manifest_pins_revision(source_manifest):
    assert core_io.get_pdb_revision(
        "PDB_00002Y4I__1__1.B__1.E_1.F",
        manifest_path=source_manifest,
    ) == (1, 5)
    with pytest.raises(KeyError, match="missing"):
        core_io.get_pdb_revision("1abc", manifest_path=source_manifest)
    assert core_io._get_pdb_revisions(
        ["2y4i", "9zzz"],
        manifest_path=source_manifest,
    ) == {"2y4i": (1, 5), "9zzz": (1, 0)}


@pytest.mark.parametrize(
    "value, expected",
    [("true", True), ("1", True), ("false", False), ("0", False), ("", False)],
)
def test_offline_mode_alias_parsing(monkeypatch, value, expected):
    monkeypatch.delenv("PLINDER_OFFLINE", raising=False)
    monkeypatch.setenv("PLINDER_OFFLINE_MODE", value)
    assert is_offline() is expected


def test_source_mmcif_rejects_wrong_version(release_cache, cif_2y4i, monkeypatch):
    monkeypatch.setattr(
        core_io.requests,
        "get",
        lambda *args, **kwargs: _Response(cif_2y4i.read_bytes()),
    )

    with pytest.raises(ValueError, match="revision mismatch"):
        core_io.get_pdb_mmcif("2y4i", revision=(2, 0))
    assert not core_io.pdb_mmcif_cache_path(
        "2y4i",
        revision=(2, 0),
    ).exists()


@pytest.mark.parametrize("value", ["", "../../etc/passwd", "not-a-pdb"])
def test_source_mmcif_rejects_invalid_ids(value):
    with pytest.raises(ValueError, match="invalid PDB ID"):
        core_io.pdb_mmcif_cache_path(value, revision=(1, 0))
