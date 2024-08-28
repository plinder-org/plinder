
def test_split_plot(monkeypatch, write_plinder_mount, tmp_path):

    def update_plindex():
        from plinder.core.index import utils

        df = utils.get_plindex()
        utils._PLINDEX = None
        df["pli_unique_qcov__50__communities"] = 0
        df["pocket_qcov__50__communities"] = 0
        df["tanimoto_similarity_max__50__communities"] = 0
        df["system_proper_ligand_unique_ccd_codes"] = "IDK"
        df["ligand_is_proper"] = True
        df.to_parquet(write_plinder_mount / "index" / "annotation_table.parquet", index=False)
        return df

    monkeypatch.setattr("plinder.core.split.plot.get_extended_plindex", update_plindex)

    from plinder.core.split.plot import SplitPropertiesPlotter

    output_dir = tmp_path / "split_plots"
    output_dir.mkdir(exist_ok=True)

    plotter = SplitPropertiesPlotter.from_files(
        data_dir=write_plinder_mount,
        split_file=write_plinder_mount / "splits" / "split.parquet",
        output_dir=output_dir,
        # stratified_train_test_file=write_plinder_mount / "strat" / "train_vs_test_data" / "test_set.parquet",
        # stratified_train_val_file=write_plinder_mount / "strat" / "train_vs_val_data" / "val_set.parquet",
        stratified_val_test_file=write_plinder_mount / "strat" / "val_vs_test_data" / "test_set.parquet",
        make_plots=False,
    )
    try:
        plotter.plot_all()
    except Exception:
        pass

    assert len(list(output_dir.rglob("*"))) > 0
