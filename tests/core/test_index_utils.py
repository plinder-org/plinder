# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os
from pathlib import Path

import pytest
from plinder.core.index import utils


def mock_path(*, rel: str = "", download: bool = False, force_progress: bool = False):
    obj = Path(
        "/".join(
            [
                str(os.getenv("PLINDER_MOUNT")),
                str(os.getenv("PLINDER_BUCKET")),
                str(os.getenv("PLINDER_RELEASE")),
            ]
        )
    )
    return obj / rel if rel else obj


@pytest.fixture
def mock_cpl(read_plinder_mount, monkeypatch):
    # patch cpl at core.utils not core.index.utils because of unpack
    monkeypatch.setattr(
        "plinder.core.utils.cpl.get_plinder_path",
        mock_path,
    )
    monkeypatch.setattr(
        "plinder.core.utils.cpl.download_paths",
        lambda **kws: None,
    )
    yield

    from plinder.core.utils import config

    config._config._clear()


def test_get_plindex(mock_cpl):
    df = utils.get_plindex()
    assert len(df.index) == 57
    assert "pli_unique_qcov__50__strong__component" in df.columns


def test_get_manifest(mock_cpl):
    df = utils.get_manifest()
    assert len(df.index) == 57


@pytest.mark.parametrize(
    "args",
    [
        [],
        ["--release", "2024-04"],
        ["--release-number", "2"],
        ["--release", "2026-07", "--release-number", "2"],
    ],
)
def test_download_cmd(args, mock_cpl):
    utils.download_plinder_cmd(args=args + ["-y"])


def test_download_cmd_forwards_release_identity(mock_cpl, monkeypatch):
    calls = []
    original_get_config = utils.get_config

    def track_config(**kwargs):
        calls.append(kwargs.get("config"))
        return original_get_config(**kwargs)

    monkeypatch.setattr(utils, "get_config", track_config)

    utils.download_plinder_cmd(
        args=["--release", "2026-07", "--release-number", "2", "-y"]
    )

    assert calls == [
        {
            "data": {
                "plinder_release": "2026-07",
                "plinder_release_number": "2",
            }
        }
    ]


def test_download_cmd_does_not_fetch_source_mmcif_cache(mock_cpl, monkeypatch):
    from plinder.core.utils import cpl

    requested = []

    def track_path(**kwargs):
        requested.append(kwargs.get("rel", ""))
        return mock_path(**kwargs)

    monkeypatch.setattr(cpl, "get_plinder_path", track_path)

    utils.download_plinder_cmd(args=["-y"])

    assert "index/annotation_table.parquet" in requested
    assert "index/linked_apo_structures.parquet" in requested
    assert not any(path.startswith("source_mmcifs") for path in requested)


def test_download_uses_only_release_artifacts(mock_cpl, monkeypatch):
    from plinder.core.utils import cpl

    requested = []

    def track_path(**kwargs):
        requested.append(kwargs.get("rel", ""))
        return mock_path(**kwargs)

    monkeypatch.setattr(cpl, "get_plinder_path", track_path)

    utils.download_plinder_cmd(args=["-y"])

    assert "alignments" in requested
    assert "search_databases" in requested
    assert not any(path == "scores" or path.startswith("scores/") for path in requested)
    assert "entries" not in requested
    assert "systems" not in requested
    assert "links" not in requested
    assert "linked_structures" not in requested
    assert set(requested) <= {
        path for path in utils.RELEASE_PATHS.values() if "{" not in path
    }


def test_download_can_skip_large_artifact_groups(mock_cpl, monkeypatch):
    from plinder.core.utils import cpl

    requested = []

    def track_path(**kwargs):
        requested.append(kwargs.get("rel", ""))
        return mock_path(**kwargs)

    monkeypatch.setattr(cpl, "get_plinder_path", track_path)
    monkeypatch.setattr("builtins.input", lambda prompt: "n")

    utils.download_plinder_cmd(args=[])

    large_paths = {
        utils.RELEASE_PATHS[artifact_name]
        for _, artifact_names, large_download in utils._DOWNLOAD_GROUPS
        if large_download
        for artifact_name in artifact_names
    }
    assert requested
    assert large_paths.isdisjoint(requested)
