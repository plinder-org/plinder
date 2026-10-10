# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Storage routing must never change directed score results."""

import json
from contextlib import contextmanager
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from threading import Thread
from types import SimpleNamespace
from urllib.parse import unquote

import pandas as pd
import pytest

from plinder.core.release import PlinderRelease
from plinder.core.scores.query import query_score_table
from plinder.core.scores.store import score_bucket
from plinder.data.pipeline.score_store import export_score_store


@pytest.fixture
def stored_scores(tmp_path):
    source = tmp_path / "source"
    path = source / "exports/interface_similarity_scores.parquet"
    path.parent.mkdir(parents=True)
    frame = pd.DataFrame(
        {
            "query_system": ["1abc__1", "2def__1", "3ghi__1", "1abc__1"],
            "target_system": ["2def__1", "1abc__1", "2def__1", "2def__1"],
            "similarity": pd.Series([91, 32, 73, 91], dtype="int8"),
        }
    )
    frame.to_parquet(path, index=False)
    destination = tmp_path / "optimized"
    manifest = export_score_store(source, destination, "interface_similarity_scores")
    return source, destination, frame, manifest


def _canonical(frame):
    return frame.sort_values(list(frame.columns)).reset_index(drop=True)


@pytest.mark.parametrize(
    "filters",
    [
        [("query_system", "==", "1abc__1")],
        [("target_system", "==", "2def__1")],
        [("target_system", "in", {"2def__1", "1abc__1"})],
        [("query_system", "==", "not-present")],
        [("target_system", "==", "not-present")],
        [("query_system", "!=", "1abc__1")],
        [("similarity", ">", 50)],
        [
            [("query_system", "==", "1abc__1"), ("similarity", ">", 50)],
            [("target_system", "==", "1abc__1")],
        ],
        [
            [("query_system", "==", "1abc__1")],
            [("similarity", "<", 50)],
        ],
        [("query_system", "==", "1abc__1"), ("query_system", "==", "2def__1")],
    ],
)
@pytest.mark.parametrize("columns", [None, ["similarity"]])
def test_sharded_queries_match_legacy(stored_scores, filters, columns):
    source, destination, _, _ = stored_scores
    expected = query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(source),
        filters=filters,
        columns=columns,
    )
    actual = query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(destination),
        filters=filters,
        columns=columns,
    )
    pd.testing.assert_frame_equal(_canonical(actual), _canonical(expected))


def test_only_routed_files_are_opened(stored_scores, monkeypatch):
    import plinder.core.scores.store as store

    _, destination, _, manifest = stored_scores
    opened = []
    original = store._file

    def record(root, relative, release):
        opened.append(relative)
        return original(root, relative, release)

    monkeypatch.setattr(store, "_file", record)
    query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(destination),
        filters=[("query_system", "==", "1abc__1")],
    )
    assert opened == manifest["data"][score_bucket("1abc__1")]
    opened.clear()
    result = query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(destination),
        filters=[("target_system", "==", "1abc__1")],
    )
    assert result["similarity"].tolist() == [32]  # Not the reverse score, 91.
    assert set(opened) == set(
        manifest["targets"][score_bucket("1abc__1")]
        + manifest["data"][score_bucket("2def__1")]
    )


def test_reject_unbounded_query_before_fetch(monkeypatch):
    monkeypatch.setattr(PlinderRelease, "fetch", lambda *a, **kw: pytest.fail("fetch"))
    for filters in [None, [], [[]]]:
        with pytest.raises(ValueError, match="filter"):
            query_score_table(
                "interface_similarity_scores", release=PlinderRelease(), filters=filters
            )


def test_missing_shard_does_not_fall_back(stored_scores):
    _, destination, _, manifest = stored_scores
    root = destination / "score_store/interface_similarity_scores"
    (root / manifest["data"][score_bucket("1abc__1")][0]).unlink()
    with pytest.raises(FileNotFoundError, match="score-store"):
        query_score_table(
            "interface_similarity_scores",
            release=PlinderRelease(destination),
            filters=[("query_system", "==", "1abc__1")],
        )


def test_bad_manifest_is_not_legacy(stored_scores):
    _, destination, _, _ = stored_scores
    path = destination / "score_store/interface_similarity_scores/manifest.json"
    path.write_text(json.dumps({"version": 999}))
    with pytest.raises(ValueError, match="manifest"):
        query_score_table(
            "interface_similarity_scores",
            release=PlinderRelease(destination),
            filters=[("similarity", ">", 50)],
        )


def test_export_never_overwrites_existing_store(stored_scores):
    source, destination, _, _ = stored_scores
    with pytest.raises(FileExistsError):
        export_score_store(source, destination, "interface_similarity_scores")


def test_export_isolates_engines_and_builds_exact_target_index(tmp_path, monkeypatch):
    import plinder.data.pipeline.score_store as exporter

    source = tmp_path / "source"
    path = source / "exports/interface_similarity_scores.parquet"
    path.parent.mkdir(parents=True)
    frame = pd.DataFrame(
        [(f"1a{q}b__1", f"2c{t}d__1", 80) for q in "0123" for t in "0123456789ab"],
        columns=["query_system", "target_system", "similarity"],
    )
    frame = pd.concat([frame, frame.iloc[:1]], ignore_index=True)
    frame.to_parquet(path, index=False)
    original = exporter._connection
    connections = []
    active = 0

    @contextmanager
    def isolated(*args):
        nonlocal active
        assert active == 0
        active += 1
        try:
            with original(*args) as con:
                connections.append(con)
                yield con
        finally:
            active -= 1

    monkeypatch.setattr(exporter, "_connection", isolated)
    destination = tmp_path / "out"
    manifest = export_score_store(source, destination, "interface_similarity_scores")
    root = destination / "score_store/interface_similarity_scores"
    index = pd.concat(
        [
            pd.read_parquet(root / p)
            for paths in manifest["targets"].values()
            for p in paths
        ]
    )
    expected = frame.assign(
        column="target_system",
        value=frame["target_system"],
        shard=frame["query_system"].map(score_bucket),
    )[["column", "value", "shard"]].drop_duplicates()
    pd.testing.assert_frame_equal(_canonical(index), _canonical(expected))
    assert manifest["rows"] == len(frame)  # Index deduplication must not drop scores.
    assert not list(root.glob("index-parts-*"))
    # One initial partitioning engine, one sort and one index engine per query
    # bucket, then between one and eight coarse target-index engines.
    assert (
        1 + 2 * len(manifest["data"])
        < len(connections)
        <= 9 + 2 * len(manifest["data"])
    )
    for con in connections:
        with pytest.raises(Exception, match="closed"):
            con.execute("SELECT 1")


def test_export_failure_removes_staging_and_preserves_source(tmp_path, monkeypatch):
    import plinder.data.pipeline.score_store as exporter

    source = tmp_path / "source"
    path = source / "exports/interface_similarity_scores.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame(
        {"query_system": ["1abc__1"], "target_system": ["2def__1"], "similarity": [80]}
    ).to_parquet(path, index=False)
    before = path.read_bytes()

    def fail(*args, **kwargs):
        raise RuntimeError("injected sort failure")

    monkeypatch.setattr(exporter, "_sort_bucket", fail)
    destination = tmp_path / "out"
    with pytest.raises(RuntimeError, match="injected sort failure"):
        export_score_store(source, destination, "interface_similarity_scores")
    assert path.read_bytes() == before
    assert not list((destination / "score_store").iterdir())


def test_numeric_buckets_and_duplicate_rows(tmp_path):
    source = tmp_path / "source"
    (source / "ligand_scores").mkdir(parents=True)
    frame = pd.DataFrame(
        {
            "query_ligand_id": [29, 1025, 29],
            "target_ligand_id": [1025, 29, 1025],
            "score": [91, 32, 91],
        }
    )
    frame.to_parquet(source / "ligand_scores/part.parquet", index=False)
    export_score_store(source, tmp_path / "out", "ligand_scores")
    result = query_score_table(
        "ligand_scores",
        release=PlinderRelease(tmp_path / "out"),
        filters=[("query_ligand_id", "==", 29)],
    )
    assert result["score"].tolist() == [91, 91]
    reverse = query_score_table(
        "ligand_scores",
        release=PlinderRelease(tmp_path / "out"),
        filters=[("target_ligand_id", "==", 29.0)],
    )
    assert reverse["score"].tolist() == [32]


def test_query_aliases_route_to_the_same_entry_bucket(tmp_path):
    source = tmp_path / "source"
    path = source / "exports/ligand_similarity_scores.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_system": ["query__1__A__B", "1abc__1__A__B"],
            "query_ligand_id": ["query__1__B", "1abc__1__B"],
            "target_system": ["2def__1__C__D"] * 2,
            "target_ligand_id": ["2def__1__D"] * 2,
            "score": [91, 32],
        }
    ).to_parquet(path, index=False)
    output = tmp_path / "output"
    export_score_store(source, output, "ligand_similarity_scores")
    for column, value in [
        ("query_ligand_id", "query__1__B"),
        ("target_ligand_id", "2def__1__D"),
    ]:
        filters = [(column, "==", value)]
        result = query_score_table(
            "ligand_similarity_scores", release=PlinderRelease(output), filters=filters
        )
        pd.testing.assert_frame_equal(
            _canonical(result), _canonical(pd.read_parquet(path, filters=filters))
        )


@pytest.fixture
def score_mirror(stored_scores, tmp_path, monkeypatch):
    from plinder.core import release as release_module
    from plinder.core.utils import cpl

    _, origin, _, manifest = stored_scores
    requests = []

    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            requests.append(unquote(self.path))
            super().do_GET()

    server = ThreadingHTTPServer(
        ("127.0.0.1", 0), partial(Handler, directory=str(origin))
    )
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    cfg = SimpleNamespace(
        data=SimpleNamespace(
            plinder_release="2026-09",
            plinder_release_number="",
            plinder_dir=str(tmp_path / "cache"),
            plinder_remote=f"http://127.0.0.1:{server.server_port}",
            force_update=False,
        )
    )
    monkeypatch.setattr(cpl, "get_config", lambda: cfg)
    monkeypatch.setattr(release_module, "get_config", lambda: cfg)
    monkeypatch.delenv("PLINDER_OFFLINE", raising=False)
    monkeypatch.delenv("PLINDER_OFFLINE_MODE", raising=False)
    yield cfg, requests, manifest
    server.shutdown()
    server.server_close()
    thread.join()


def test_http_downloads_only_manifest_and_needed_shard(score_mirror, monkeypatch):
    cfg, requests, manifest = score_mirror
    filters = [("query_system", "==", "1abc__1")]
    expected_requests = ["/score_store/interface_similarity_scores/manifest.json"] + [
        "/score_store/interface_similarity_scores/" + p
        for p in manifest["data"][score_bucket("1abc__1")]
    ]
    result = query_score_table(
        "interface_similarity_scores", release=PlinderRelease(), filters=filters
    )
    assert result["similarity"].tolist() == [91, 91]
    assert requests == expected_requests
    requests.clear()
    query_score_table(
        "interface_similarity_scores", release=PlinderRelease(), filters=filters
    )
    assert requests == []  # Metadata checks don't re-download cached files.
    monkeypatch.setenv("PLINDER_OFFLINE", "true")
    result = query_score_table(
        "interface_similarity_scores", release=PlinderRelease(), filters=filters
    )
    assert len(result) == 2
    assert requests == []
    # A different bucket is absent, not a valid empty result in a partial cache.
    with pytest.raises(FileNotFoundError, match="score-store"):
        query_score_table(
            "interface_similarity_scores",
            release=PlinderRelease(),
            filters=[("query_system", "==", "2def__1")],
        )
    assert requests == []


def test_http_target_only_lookup_preserves_direction(score_mirror):
    _, requests, manifest = score_mirror
    result = query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(),
        filters=[("target_system", "==", "1abc__1")],
    )
    assert result["similarity"].tolist() == [32]
    expected = (
        ["manifest.json"]
        + manifest["targets"][score_bucket("1abc__1")]
        + manifest["data"][score_bucket("2def__1")]
    )
    assert requests == [
        "/score_store/interface_similarity_scores/" + p for p in expected
    ]


def test_http_unknown_query_downloads_no_score_file(score_mirror):
    _, requests, _ = score_mirror
    result = query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(),
        filters=[("query_system", "==", "missing")],
    )
    assert result.empty
    assert requests == ["/score_store/interface_similarity_scores/manifest.json"]


def test_export_preserves_hive_columns_and_uint8(tmp_path):
    source = tmp_path / "source"
    path = (
        source / "exports/protein_similarity_scores/alignment_type=mmseqs/part.parquet"
    )
    path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "target_entry": ["2def"],
            "qcov": pd.Series([80], dtype="uint8"),
        }
    ).to_parquet(path, index=False)
    output = tmp_path / "output"
    export_score_store(source, output, "protein_similarity_scores")
    filters = [("query_entry", "==", "1abc")]
    expected = pd.read_parquet(path.parents[1], filters=filters)
    result = query_score_table(
        "protein_similarity_scores", release=PlinderRelease(output), filters=filters
    )
    pd.testing.assert_frame_equal(expected, result)


def test_empty_export_retains_schema(tmp_path):
    source = tmp_path / "source"
    path = source / "exports/interface_similarity_scores.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_system": pd.Series(dtype="str"),
            "target_system": pd.Series(dtype="str"),
            "similarity": pd.Series(dtype="int8"),
        }
    ).to_parquet(path, index=False)
    export_score_store(source, tmp_path / "output", "interface_similarity_scores")
    result = query_score_table(
        "interface_similarity_scores",
        release=PlinderRelease(tmp_path / "output"),
        filters=[("similarity", ">", 0)],
    )
    assert result.empty
    assert str(result["similarity"].dtype) == "int8"


def test_unsafe_manifest_path_is_rejected_before_download(stored_scores):
    _, destination, _, manifest = stored_scores
    path = destination / "score_store/interface_similarity_scores/manifest.json"
    manifest["data"][score_bucket("1abc__1")] = ["../../../outside.parquet"]
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="unsafe"):
        query_score_table(
            "interface_similarity_scores",
            release=PlinderRelease(destination),
            filters=[("query_system", "==", "1abc__1")],
        )
