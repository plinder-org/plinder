# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
def test_data_loader(read_plinder_mount):
    from plinder.core.loader import PlinderDataset

    ds = PlinderDataset(split="removed", use_alternate_structures=False)
    assert len(ds[0])
