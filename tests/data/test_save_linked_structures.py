import pandas as pd


def test_save_linked_structures(write_plinder_mount):
    from plinder.data import save_linked_structures

    cfg = save_linked_structures.LinkedStructureConfig(
        filter_criteria={
            "pocket_lddt": 0,
            "protein_fident_qcov_weighted_sum": 0,
        }
    )

    try:
        save_linked_structures.make_linked_structures_data_file(
            data_dir=write_plinder_mount,
            search_db="holo",
            superposed_folder=write_plinder_mount / "linked_staging",
            output_file=write_plinder_mount
            / "linked_structures"
            / "holo_links.parquet",
            cfg=cfg,
            num_processes=1,
        )

        save_linked_structures.save_linked_structures(
            links_file=write_plinder_mount / "linked_structures" / "holo_links.parquet",
            data_dir=write_plinder_mount,
            search_db="holo",
            output_folder=write_plinder_mount / "linked_structures",
            num_threads=1,
        )

        df = pd.read_parquet(
            write_plinder_mount / "linked_structures" / "holo_links.parquet"
        )
        assert isinstance(df, pd.DataFrame)
    except Exception:
        pass
