import base64
import gzip
import hashlib
import json
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread
from types import SimpleNamespace

import pytest
from plinder.core.utils import cpl
from plinder.core.utils import r2 as dataset


@pytest.fixture
def mirror(tmp_path, monkeypatch):
    origin = tmp_path / "origin"
    release = origin / "2024-06/v2"
    release.mkdir(parents=True)
    records = []
    for key, data in [
        ("systems/ab.zip", b"archive"),
        ("scores/a.parquet", b"one"),
        ("scores/b.parquet", b"two"),
    ]:
        path = release / key
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(data)
        records.append(
            {
                "key": key,
                "size": len(data),
                "md5": base64.b64encode(hashlib.md5(data).digest()).decode(),
            }
        )
    (release / "manifest.jsonl.gz").write_bytes(
        gzip.compress(("\n".join(map(json.dumps, records))).encode())
    )
    server = ThreadingHTTPServer(
        ("127.0.0.1", 0), partial(SimpleHTTPRequestHandler, directory=str(origin))
    )
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    cfg = SimpleNamespace(
        data=SimpleNamespace(
            plinder_bucket="plinder",
            plinder_release="2024-06",
            plinder_iteration="v2",
            plinder_dir=str(tmp_path / "cache"),
            plinder_remote=f"http://127.0.0.1:{server.server_port}/2024-06/v2",
            force_update=False,
        )
    )
    monkeypatch.setattr(cpl, "get_config", lambda: cfg)
    monkeypatch.delenv("PLINDER_OFFLINE", raising=False)
    cpl.manifest.cache_clear()
    yield cfg, release
    server.shutdown()
    server.server_close()
    thread.join()
    cpl.manifest.cache_clear()


def test_directory_download_and_truncated_cache_repair(mirror):
    cfg, _ = mirror
    path = cpl.get_plinder_path(rel="scores")
    assert (path / "a.parquet").read_bytes() == b"one"
    (path / "a.parquet").write_bytes(b"short")
    cpl.get_plinder_path(rel="scores")
    assert (path / "a.parquet").read_bytes() == b"one"
    assert not (Path(cfg.data.plinder_dir) / "systems/ab.zip").exists()


def test_archive_listing_without_download(mirror):
    cfg, _ = mirror
    root = cpl.get_plinder_path(rel="systems", download=False)
    paths = cpl.list_zip_paths("systems")
    assert paths == [root / "ab.zip"]
    assert not paths[0].exists()
    assert cpl.get_plinder_paths(paths=paths) == paths
    assert paths[0].read_bytes() == b"archive"


def test_small_batch_failure_propagates(mirror):
    cfg, origin = mirror
    (origin / "systems/ab.zip").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        cpl.get_plinder_path(rel="systems/ab.zip")
    assert not (Path(cfg.data.plinder_dir) / "systems/ab.zip").exists()


def test_offline_does_not_access_network(mirror, monkeypatch):
    monkeypatch.setenv("PLINDER_OFFLINE", "true")
    monkeypatch.setattr(cpl, "_get_client", lambda: pytest.fail("network"))
    assert isinstance(cpl.get_plinder_path(rel="missing"), Path)
    assert cpl.list_zip_paths("systems") == []


@pytest.mark.parametrize(
    "rel", ["../escape", "/absolute", "systems/../../escape", "a\\b"]
)
def test_path_traversal_rejected(mirror, rel):
    with pytest.raises(ValueError):
        cpl.get_plinder_path(rel=rel)


def test_unsupported_release_rejected(mirror):
    cfg, _ = mirror
    cfg.data.plinder_iteration = "v1"
    with pytest.raises(ValueError, match="2024-06/v2"):
        cpl.get_plinder_path(rel="systems")


def test_missing_file_fails(mirror):
    with pytest.raises(FileNotFoundError):
        cpl.get_plinder_path(rel="missing")


def test_archive_extraction_and_repair(mirror):
    from zipfile import ZipFile

    from plinder.core.utils import unpack

    cfg, origin = mirror
    cfg.data.systems = "systems"
    cfg.context = SimpleNamespace(system_ids=[], pdb_ids=[], two_char_codes=["ab"])
    archive = origin / "systems/ab.zip"
    with ZipFile(archive, "w") as stream:
        stream.writestr("example/protein.cif", "structure")
    payload = archive.read_bytes()
    record = {
        "key": "systems/ab.zip",
        "size": len(payload),
        "md5": base64.b64encode(hashlib.md5(payload).digest()).decode(),
    }
    (origin / "manifest.jsonl.gz").write_bytes(
        gzip.compress(json.dumps(record).encode())
    )
    unpack.get_zips_to_unpack(kind="systems", cfg=cfg)
    root = Path(cfg.data.plinder_dir) / "systems"
    assert (root / "example/protein.cif").read_text() == "structure"
    assert (root / "ab_done").exists()
    (root / "ab.zip").write_bytes(b"broken")
    cpl.download_paths(paths=[root / "ab.zip"])
    assert not (root / "ab_done").exists()


def test_cached_same_size_edit_requires_force_update(mirror, monkeypatch):
    cfg, _ = mirror
    path = cpl.get_plinder_path(rel="scores/a.parquet")
    path.write_bytes(b"bad")
    assert cpl.get_plinder_path(rel="scores/a.parquet").read_bytes() == b"bad"
    cfg.data.force_update = True
    assert cpl.get_plinder_path(rel="scores/a.parquet").read_bytes() == b"one"


def test_failed_checksum_preserves_existing_cache(mirror):
    cfg, origin = mirror
    path = cpl.get_plinder_path(rel="scores/a.parquet")
    (origin / "scores/a.parquet").write_bytes(b"bad")
    cfg.data.force_update = True
    with pytest.raises(ValueError, match="checksum"):
        cpl.get_plinder_path(rel="scores/a.parquet")
    assert path.read_bytes() == b"one"
    assert list(path.parent.iterdir()) == [path]


def test_native_cloudpath_listing_and_download(mirror, tmp_path):
    from cloudpathlib import CloudPath

    client = cpl._get_client()
    root = client.path()
    assert isinstance(root, CloudPath)
    assert root.is_dir() and not root.is_file()
    assert {p.name for p in root.iterdir()} == {"scores", "systems"}
    assert {p.name for p in root.rglob("*.parquet")} == {"a.parquet", "b.parquet"}
    path = root / "scores" / "a.parquet"
    assert path.client is client
    assert not (root / "absent").exists()
    assert path.download_to(tmp_path / "explicit").read_bytes() == b"one"


def test_cache_symlink_escape_is_rejected(mirror, tmp_path):
    cfg, _ = mirror
    root = Path(cfg.data.plinder_dir)
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (root / "scores").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError):
        cpl.get_plinder_path(rel="scores/a.parquet")
    assert not (outside / "a.parquet").exists()


def test_force_update_fetches_manifest_once(mirror, monkeypatch):
    cfg, _ = mirror
    actual = dataset.manifest.__wrapped__
    calls = []

    def tracked(remote):
        calls.append(remote)
        return actual(remote)

    monkeypatch.setattr(dataset, "manifest", tracked)
    cfg.data.force_update = True
    cpl.get_plinder_path(rel="scores")
    assert len(calls) == 1


def test_cached_file_is_not_read_or_transferred(mirror, monkeypatch):
    path = cpl.get_plinder_path(rel="scores/a.parquet")
    actual_open = Path.open

    def guarded_open(self, *args, **kwargs):
        if self == path:
            pytest.fail("cached file must not be opened")
        return actual_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded_open)
    monkeypatch.setattr(
        dataset.ReleaseClient, "_download_file", lambda *a: pytest.fail("download")
    )
    assert cpl.get_plinder_path(rel="scores/a.parquet") == path


def test_interrupted_native_transfer_restarts(tmp_path, monkeypatch):
    import socket
    from http.server import BaseHTTPRequestHandler

    payload = b"abcdef" * 500000
    row = {
        "key": "archive.zip",
        "size": len(payload),
        "md5": base64.b64encode(hashlib.md5(payload).digest()).decode(),
    }
    listing = gzip.compress(json.dumps(row).encode())
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_HEAD(self):
            self.send_response(200)
            self.end_headers()

        def do_GET(self):
            if self.path.endswith("manifest.jsonl.gz"):
                data = listing
            else:
                requests.append(self.headers.get("Range"))
                data = payload
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            if data is payload and len(requests) == 1:
                self.wfile.write(data[:1500000])
                self.wfile.flush()
                self.connection.shutdown(socket.SHUT_WR)
            else:
                self.wfile.write(data)

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    monkeypatch.setattr(dataset, "sleep", lambda _: None)
    try:
        client = dataset.ReleaseClient(f"http://127.0.0.1:{server.server_port}")
        target = client.path("archive.zip").download_to(tmp_path / "archive.zip")
        assert target.read_bytes() == payload
        assert requests == [None, None]
    finally:
        server.shutdown()
        thread.join()
        server.server_close()


def test_full_download_includes_database_metadata(tmp_path, monkeypatch):
    from dataclasses import asdict

    from omegaconf import DictConfig
    from plinder.core.index import utils
    from plinder.core.utils.config import DataConfig

    cfg = SimpleNamespace(
        data=DictConfig(asdict(DataConfig(plinder_mount=str(tmp_path))))
    )
    requested = []

    def download(*, rel, **kwargs):
        requested.append(rel)
        return tmp_path / rel

    monkeypatch.setattr(utils, "get_config", lambda **kwargs: cfg)
    monkeypatch.setattr(cpl, "get_plinder_path", download)
    monkeypatch.setattr(utils, "get_zips_to_unpack", lambda **kwargs: {})
    utils.download_plinder_cmd(["--yes"])
    assert "dbs" in requested
