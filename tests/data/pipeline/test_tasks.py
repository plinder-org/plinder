# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest
from plinder.data.pipeline import io, tasks


@pytest.mark.parametrize(
    "inputs, expected",
    [
        ({"batch_size": 4, "two_char_codes": []}, [4, 4]),
        ({"batch_size": 4, "two_char_codes": []}, [4, 3]),
        ({"two_char_codes": ["xx"], "batch_size": 4}, [1]),
    ],
)
def test_scatter_download_rcsb_files(inputs, expected, tmp_path):
    codes = ["aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh"]
    _orig_rsync_rcsb = io.rsync_rcsb
    _orig_list_rcsb = io.list_rcsb
    io.rsync_rcsb = lambda **kws: codes
    io.list_rcsb = lambda **kws: codes if expected[0] == expected[1] else codes[:-1]
    chunks = tasks.scatter_download_rcsb_files(data_dir=tmp_path, **inputs)
    for chunk, expect in zip(chunks, expected):
        assert len(chunk) == expect
    io.rsync_rcsb = _orig_rsync_rcsb
    io.list_rcsb = _orig_list_rcsb


def test_download_rcsb_files(tmp_path):
    _orig_rsync_rcsb = io.rsync_rcsb
    io.rsync_rcsb = lambda *args, **kws: None
    tasks.download_rcsb_files(data_dir=tmp_path, two_char_codes=["aa"])
    io.rsync_rcsb = _orig_rsync_rcsb
