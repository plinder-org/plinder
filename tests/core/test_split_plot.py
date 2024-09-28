def test_split_plot(write_plinder_mount, split_plot_split_file, tmp_path):
    from plinder.core.split.plot import SplitPropertiesPlotter

    output_dir = tmp_path / "split_plots"
    output_dir.mkdir(exist_ok=True)

    plotter = SplitPropertiesPlotter.from_files(
        data_dir=write_plinder_mount,
        split_file=split_plot_split_file,
        output_dir=output_dir,
        # stratified_train_test_file=write_plinder_mount / "strat" / "train_vs_test_data" / "test_set.parquet",
        # stratified_train_val_file=write_plinder_mount / "strat" / "train_vs_val_data" / "val_set.parquet",
        stratified_val_test_file=write_plinder_mount
        / "strat"
        / "val_vs_test_data"
        / "test_set.parquet",
        make_plots=False,
    )
    try:
        plotter.plot_all()
    except Exception as e:
        print(e)
        pass

    assert len(list(output_dir.rglob("*"))) > 0
