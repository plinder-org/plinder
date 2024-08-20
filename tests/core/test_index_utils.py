# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os
from pathlib import Path

import pytest
from plinder.core.index import utils
from plinder.core.split import utils as split_utils


def mock_path(*, rel: str = "", download: bool = False):
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
    assert len(df.index) == 10


def test_get_manifest(mock_cpl):
    df = utils.get_manifest()
    assert len(df.index) == 10


def test_load_entries(mock_cpl):
    blob = utils.load_entries(two_char_codes=["9h"])
    assert len(blob) == 1


def test_download_cmd(mock_cpl):
    utils.download_plinder_cmd()


def test_get_extended_plindex(mock_cpl):
    df = split_utils.get_extended_plindex()
    assert len(df.index) == 10
    assert "pli_unique_qcov__50__components" in df.columns
