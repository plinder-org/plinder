# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import os

import pytest
from omegaconf import DictConfig
from plinder.core.utils import gcs


class _Blob:
    def __init__(self, name):
        self.name = name

    def download_as_string(self):
        return b"test"

    def download_to_filename(self, local_path):
        pass


class _Bucket:
    def __init__(self, name):
        self.name = name

    def blob(self, name):
        return _Blob(name)

    def list_blobs(self, prefix):
        return []


CONF = DictConfig(
    {
        "data": {
            "plinder_bucket": "plinder",
        }
    }
)


@pytest.fixture
def mock_buckets(monkeypatch):
    monkeypatch.setattr(
        "plinder.core.utils.gcs.BUCKETS",
        {
            "plinder": _Bucket("plinder"),
            "plinder-test": _Bucket("plinder-test"),
        },
    )


def test_download_as_str(mock_buckets):
    assert gcs.download_as_str(gcs_path="gs://plinder/test", cfg=CONF) == "test"


def test_download_to_file(mock_buckets, tmp_path):
    gcs.download_to_file(
        gcs_path="gs://plinder/test",
        local_path=(tmp_path / "afile.txt").as_posix(),
        cfg=CONF,
    )


def test_download_many(mock_buckets, tmp_path):
    gcs_paths = [
        "gs://plinder/afile.txt",
        "gs://plinder/bfile.txt",
    ]
    local_paths = [
        (tmp_path / "afile.txt").as_posix(),
        (tmp_path / "bfile.txt").as_posix(),
    ]
    gcs.download_many(gcs_paths=gcs_paths, local_paths=local_paths, cfg=CONF)


def test_list_dir(mock_buckets):
    assert isinstance(gcs.list_dir(gcs_path="gs://plinder/test", cfg=CONF), list)


def test_real_download():
    gcs.download_as_str(
        gcs_path="gs://plinder/2024-04/v1/README.md", cfg=CONF
    ).startswith("plinder-data")


def test_real_list_dir():
    assert len(gcs.list_dir(gcs_path="gs://plinder/2024-04/v1/", cfg=CONF))


def test_real_download_many(tmp_path):
    gcs_paths = [
        "gs://plinder/2024-04/v1/README.md",
        "gs://plinder/2024-04/v1/README.md",
    ]
    local_paths = [
        (tmp_path / "afile.txt").as_posix(),
        (tmp_path / "bfile.txt").as_posix(),
    ]
    for path in local_paths:
        assert not os.path.exists(path)
    gcs.download_many(gcs_paths=gcs_paths, local_paths=local_paths, cfg=CONF)
    for path in local_paths:
        assert os.path.exists(path)


def test_default_behavior():
    assert gcs.download_as_str(gcs_path="gs://plinder/2024-04/v1/README.md").startswith(
        "plinder-data"
    )
