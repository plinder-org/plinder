# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
def test_data_loader(read_plinder_mount):
    from plinder.core.loader import PlinderDataset

    FILTERS = [
        (
            "system_id",
            "in",
            ["19hc__1__1.B__1.T", "19hc__1__1.A__1.I", "1avd__1__1.A__1.C"],
        )
    ]
    ds = PlinderDataset(
        split="removed", filters=FILTERS, use_alternate_structures=False
    )
    assert len(ds[0])
