import os
from pathlib import Path

import pytest

from plinder.core.index import utils

def mock_path(rel):
    class obj:
        def __init__(self, rel):
            self.fspath = Path("/".join([
                str(os.getenv("PLINDER_MOUNT")),
                str(os.getenv("PLINDER_BUCKET")),
                str(os.getenv("PLINDER_RELEASE")),
                str(rel),
            ]))
        def __truediv__(self, other):
            return self.fspath / other
    return obj(rel)

@pytest.fixture
def mock_cpl(read_plinder_mount, monkeypatch):

    # patch cpl at core.utils not core.index.utils because of unpack
    monkeypatch.setattr(
        "plinder.core.utils.cpl.get_plinder_path",
        mock_path,
    )
    # monkeypatch.setattr(
    #     "plinder.core.index.utils.cpl.get_plinder_path",
    #     mock_path,
    # )
    monkeypatch.setattr(
        "plinder.core.utils.cpl.download_many",
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
