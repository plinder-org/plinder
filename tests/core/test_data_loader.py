# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest

try:
    import atom3d
except ImportError:
    pytest.skip(allow_module_level=True)


def test_data_loader(read_plinder_mount):
    from plinder.core.loader import PlinderDataset

    ds = PlinderDataset(file_with_system_ids=read_plinder_mount / "index/annotation_table.parquet")
    assert len(ds[0])
