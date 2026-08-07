# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
def test_data_loader(cached_plinder_system):
    from plinder.core.loader import PlinderDataset

    ds = PlinderDataset(
        split="removed",
        use_alternate_structures=False,
        system_factory=cached_plinder_system,
    )
    assert len(ds[0])
