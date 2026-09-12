# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
def test_data_loader(cached_plinder_system):
    from plinder.core.loader import PlinderDataset
    from plinder.core.loader.dataset import get_torch_loader

    system = cached_plinder_system("19hc__1__1.B__1.T")
    ds = PlinderDataset(
        filters=[("system_id", "==", "19hc__1__1.B__1.T")],
        system_factory=lambda _system_id: system,
    )
    assert len(ds) == 1
    item = ds[0]
    assert set(item) == {
        "system_id",
        "holo_structure",
        "features_and_coords",
        "path",
    }

    batch = next(iter(get_torch_loader(ds, batch_size=1, num_workers=0)))
    assert set(batch) == {
        "system_ids",
        "holo_structures",
        "features_and_coords",
        "paths",
    }
    assert batch["system_ids"] == ["19hc__1__1.B__1.T"]
