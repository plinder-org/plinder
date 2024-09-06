def test_make_components_and_communities(write_plinder_mount):
    from plinder.data.clusters import make_components_and_communities

    i = len(list((write_plinder_mount / "clusters").rglob("*")))
    make_components_and_communities(
        data_dir=write_plinder_mount,
        metric="pocket_lddt",
        threshold=50,
    )
    j = len(list((write_plinder_mount / "clusters").rglob("*")))
    assert i < j
