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
        ["--iteration", "v1"],
        ["--release", "2024-06", "--iteration", "v2"],
    ],
)
def test_download_cmd(args, mock_cpl):
    utils.download_plinder_cmd(args=args + ["-y"])


def test_download_cmd_does_not_fetch_source_mmcif_cache(mock_cpl, monkeypatch):
    from plinder.core.utils import cpl

    requested = []

    def track_path(**kwargs):
        requested.append(kwargs.get("rel", ""))
        return mock_path(**kwargs)

    monkeypatch.setattr(cpl, "get_plinder_path", track_path)

    utils.download_plinder_cmd(args=["-y"])

    assert "index" in requested
    assert "source_mmcifs" not in requested


def test_v3_download_uses_alignments_instead_of_legacy_scores(mock_cpl, monkeypatch):
    from plinder.core.utils import cpl

    requested = []

    def track_path(**kwargs):
        requested.append(kwargs.get("rel", ""))
        return mock_path(**kwargs)

    monkeypatch.setattr(cpl, "get_plinder_path", track_path)

    utils.download_plinder_cmd(args=["--iteration", "v3", "-y"])

    assert "alignments" in requested
    assert "search_databases" in requested
    assert not any(path == "scores" or path.startswith("scores/") for path in requested)
    assert "entries" not in requested
    assert "systems" not in requested


def test_v2_download_does_not_request_search_databases(mock_cpl, monkeypatch):
    from plinder.core.utils import cpl

    requested = []

    def track_path(**kwargs):
        requested.append(kwargs.get("rel", ""))
        return mock_path(**kwargs)

    monkeypatch.setattr(cpl, "get_plinder_path", track_path)

    utils.download_plinder_cmd(args=["--iteration", "v2", "-y"])

    assert "search_databases" not in requested
