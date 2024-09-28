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


def test_get_plindex(mock_cpl):
    df = utils.get_plindex()
    assert len(df.index) == 57
    assert "pli_unique_qcov__50__strong__component" in df.columns


def test_get_manifest(mock_cpl):
    df = utils.get_manifest()
    assert len(df.index) == 57


def test_load_entries(mock_cpl):
    blob = utils.load_entries(two_char_codes=["9h"])
    assert len(blob) == 1


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
