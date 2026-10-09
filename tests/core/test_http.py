import os
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread
from types import SimpleNamespace

import pytest

from plinder.core.utils import cpl
from plinder.core.utils import http as dataset


@pytest.fixture
def mirror(tmp_path, monkeypatch):
    origin = tmp_path / "origin"
    release = origin / "PLINDER-2026-09"
    release.mkdir(parents=True)
    for key, data in [
        ("systems/ab.zip", b"archive"),
        ("scores/a.parquet", b"one"),
        ("scores/b.parquet", b"two"),
    ]:
        path = release / key
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(data)
        os.utime(path, (1700000000, 1700000000))
    requests = []

    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_HEAD(self):
            requests.append(("HEAD", self.path))
            super().do_HEAD()

        def do_GET(self):
            requests.append(("GET", self.path))
            super().do_GET()

    server = ThreadingHTTPServer(
        ("127.0.0.1", 0), partial(Handler, directory=str(origin))
    )
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    cfg = SimpleNamespace(
        data=SimpleNamespace(
            plinder_bucket="plinder",
            plinder_release="2026-09",
            plinder_release_number="",
            plinder_dir=str(tmp_path / "cache"),
            plinder_remote=f"http://127.0.0.1:{server.server_port}/PLINDER-2026-09",
            force_update=False,
        )
    )
    monkeypatch.setattr(cpl, "get_config", lambda: cfg)
    monkeypatch.delenv("PLINDER_OFFLINE", raising=False)
    cfg.requests = requests
    monkeypatch.delenv("PLINDER_OFFLINE_MODE", raising=False)
    yield cfg, release
    server.shutdown()
    server.server_close()
    thread.join()


def test_directory_download_and_truncated_cache_repair(mirror):
    cfg, _ = mirror
    path = cpl.get_plinder_path(rel="scores")
    assert (path / "a.parquet").read_bytes() == b"one"
    (path / "a.parquet").write_bytes(b"short")
    cpl.get_plinder_path(rel="scores")
    assert (path / "a.parquet").read_bytes() == b"one"
    assert not (Path(cfg.data.plinder_dir) / "systems/ab.zip").exists()


def test_directory_download_prunes_parquet_files_the_release_dropped(mirror):
    cfg, _ = mirror
    cache = Path(cfg.data.plinder_dir)
    stale = cache / "scores/old.parquet"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"stale")
    (cache / "scores/notes.txt").write_text("kept")
    elsewhere = cache / "other/x.parquet"
    elsewhere.parent.mkdir()
    elsewhere.write_bytes(b"kept")

    path = cpl.get_plinder_path(rel="scores")

    assert sorted(p.name for p in path.glob("*.parquet")) == ["a.parquet", "b.parquet"]
    assert (path / "notes.txt").read_text() == "kept"
    assert elsewhere.read_bytes() == b"kept"


def test_small_batch_failure_propagates(mirror):
    cfg, _ = mirror
    with pytest.raises(FileNotFoundError):
        cpl.download_paths(paths=[Path(cfg.data.plinder_dir) / "systems/missing.zip"])
    assert not (Path(cfg.data.plinder_dir) / "systems/ab.zip").exists()


def test_offline_does_not_access_network(mirror, monkeypatch):
    monkeypatch.setenv("PLINDER_OFFLINE", "true")
    monkeypatch.setattr(cpl, "_get_client", lambda: pytest.fail("network"))
    assert isinstance(cpl.get_plinder_path(rel="missing"), Path)


@pytest.mark.parametrize(
    "rel", ["../escape", "/absolute", "systems/../../escape", "a\\b"]
)
def test_path_traversal_rejected(mirror, rel):
    with pytest.raises(ValueError):
        cpl.get_plinder_path(rel=rel)


@pytest.mark.parametrize(
    "release", [("2024-06", "v2"), ("2026-09", "1"), ("2026-9", ""), ("2026-13", "")]
)
def test_unsupported_release_rejected(mirror, release):
    cfg, _ = mirror
    cfg.data.plinder_release, cfg.data.plinder_release_number = release
    with pytest.raises(ValueError, match="YYYY-MM"):
        cpl.get_plinder_path(rel="systems")


def test_missing_file_fails(mirror):
    with pytest.raises(FileNotFoundError):
        cpl.get_plinder_path(rel="missing")


def test_same_size_local_edit_preserving_mtime_requires_force_update(mirror):
    cfg, _ = mirror
    path = cpl.get_plinder_path(rel="scores/a.parquet")
    path.write_bytes(b"bad")
    os.utime(path, (1700000000, 1700000000))
    assert cpl.get_plinder_path(rel="scores/a.parquet").read_bytes() == b"bad"
    cfg.data.force_update = True
    assert cpl.get_plinder_path(rel="scores/a.parquet").read_bytes() == b"one"


def test_same_size_server_hotfix_is_fetched_without_force_update(mirror):
    cfg, origin = mirror
    path = cpl.get_plinder_path(rel="scores/a.parquet")
    source = origin / "scores/a.parquet"
    source.write_bytes(b"new")
    os.utime(source, (1700000010, 1700000010))

    assert not cfg.data.force_update
    assert cpl.get_plinder_path(rel="scores/a.parquet").read_bytes() == b"new"
    assert path.stat().st_mtime == source.stat().st_mtime


def test_directory_listing_refreshes_added_and_removed_files(mirror):
    _, origin = mirror
    path = cpl.get_plinder_path(rel="scores")
    (origin / "scores/a.parquet").unlink()
    (origin / "scores/c.parquet").write_bytes(b"new")

    assert cpl.get_plinder_path(rel="scores") == path
    assert sorted(p.name for p in path.glob("*.parquet")) == ["b.parquet", "c.parquet"]


def test_encoded_nested_paths_and_empty_directory(mirror):
    _, origin = mirror
    shard = origin / "alignments" / "alignment_type=mmseqs" / "shard=a b%.parquet"
    shard.parent.mkdir(parents=True)
    shard.write_bytes(b"data")
    (origin / "empty").mkdir()

    assert (
        cpl.get_plinder_path(rel="alignments")
        / shard.relative_to(origin / "alignments")
    ).read_bytes() == b"data"
    assert cpl.get_plinder_path(rel="empty").is_dir()


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


def test_force_update_fetches_metadata_once_per_file(mirror):
    cfg, _ = mirror
    cfg.data.force_update = True
    cpl.get_plinder_path(rel="scores")
    for name in ["a", "b"]:
        assert (
            cfg.requests.count(("HEAD", f"/PLINDER-2026-09/scores/{name}.parquet")) == 1
        )
    assert all("manifest" not in url for _, url in cfg.requests)


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
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_HEAD(self):
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()

        def do_GET(self):
            requests.append(self.headers.get("Range"))
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            if len(requests) == 1:
                self.wfile.write(payload[:1500000])
                self.wfile.flush()
                self.connection.shutdown(socket.SHUT_WR)
            else:
                self.wfile.write(payload)

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


@pytest.mark.parametrize(
    "resource,failure",
    [
        (resource, failure)
        for resource in ["metadata", "listing", "object"]
        for failure in [
            "503",
            "stall_headers",
            "stall_body",
            "exhausted",
            "404",
            "truncated",
        ]
        if resource != "metadata" or failure not in {"stall_body", "truncated"}
    ],
)
def test_download_recovery(tmp_path, monkeypatch, resource, failure):
    from http.server import BaseHTTPRequestHandler
    from threading import Event

    payload = b"dataset"
    listing = b'<title>Index of /</title><a href="data">data</a>'
    release_stall = Event()
    requests = []
    original_timeout = dataset._TimeoutHandler.http_request

    def short_timeout(self, request):
        request = original_timeout(self, request)
        assert request.timeout == 60
        request.timeout = 0.1
        return request

    monkeypatch.setattr(dataset._TimeoutHandler, "http_request", short_timeout)
    monkeypatch.setattr(dataset, "sleep", lambda _: None)

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def respond(self, head=False):
            is_listing = self.path == "/"
            selected = (
                head
                if resource == "metadata"
                else not head and is_listing == (resource == "listing")
            )
            if selected:
                requests.append(self.path)
            fail = selected and (
                len(requests) == 1 or failure in {"404", "exhausted", "truncated"}
            )
            if fail and failure in {"503", "404", "exhausted"}:
                self.send_error(503 if failure == "exhausted" else int(failure))
                return
            if fail and failure == "stall_headers":
                release_stall.wait(5)
                return
            data = listing if is_listing else payload
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            if head:
                return
            if fail and failure == "stall_body":
                release_stall.wait(5)
                return
            self.wfile.write(data[:-1] if fail and failure == "truncated" else data)

        def do_HEAD(self):
            self.respond(head=True)

        def do_GET(self):
            self.respond()

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    destination = tmp_path / "data"
    destination.write_bytes(b"existing")

    def download():
        client = dataset.ReleaseClient(f"http://127.0.0.1:{server.server_port}")
        if resource == "listing":
            list(client.path().iterdir())
        client.path("data").download_to(destination)

    try:
        if failure in {"404", "exhausted", "truncated"}:
            from cloudpathlib.exceptions import CloudPathNotExistsError

            error = (OSError, dataset.HTTPException, CloudPathNotExistsError)
            with pytest.raises(error):
                download()
            assert len(requests) == (1 if failure == "404" else 3)
            assert destination.read_bytes() == b"existing"
            assert list(tmp_path.iterdir()) == [destination]
        else:
            download()
            assert len(requests) == 2
            assert destination.read_bytes() == payload
    finally:
        release_stall.set()
        server.shutdown()
        server.server_close()
        thread.join()


def test_listing_ignores_navigation_and_escaping_links(mirror, monkeypatch):
    cfg, _ = mirror
    client = cpl._get_client()
    original_open = client.opener.open

    def open_index(request):
        response = original_open(request)
        if request.full_url.endswith("/scores/"):
            from io import BytesIO

            response.read = BytesIO(
                b"<title>Index of /scores/</title>"
                b'<a href="../">parent</a><a href="?C=N;O=D">sort</a>'
                b'<a href="https://example.com/file">external</a>'
                b'<a href="/file">absolute</a><a href="%2e%2e/escape">escape</a>'
                b'<a href="nested%2fescape">encoded escape</a>'
                b'<a href="a.parquet">a</a><a href="b.parquet">b</a>'
            ).read
        return response

    monkeypatch.setattr(client.opener, "open", open_index)
    assert {p.name for p in client.path("scores").iterdir()} == {
        "a.parquet",
        "b.parquet",
    }
    assert all(
        "escape" not in url and "example.com" not in url for _, url in cfg.requests
    )


def test_invalid_listing_preserves_cached_parquet(mirror, monkeypatch):
    cfg, _ = mirror
    path = cpl.get_plinder_path(rel="scores")
    monkeypatch.setattr(dataset._DirectoryIndex, "handle_data", lambda *args: None)
    with pytest.raises(ValueError, match="directory index"):
        cpl.get_plinder_path(rel="scores")
    assert (path / "a.parquet").read_bytes() == b"one"
    assert (path / "b.parquet").read_bytes() == b"two"


def test_empty_directory_cannot_prune_through_a_cache_symlink(mirror, tmp_path):
    cfg, origin = mirror
    (origin / "empty").mkdir()
    root = Path(cfg.data.plinder_dir)
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "kept.parquet").write_bytes(b"kept")
    (root / "empty").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError):
        cpl.get_plinder_path(rel="empty")
    assert (outside / "kept.parquet").read_bytes() == b"kept"
