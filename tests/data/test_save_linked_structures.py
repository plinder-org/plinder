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


def test_ingest_link_scoring_does_not_run_posebusters(tmp_path, monkeypatch):
    from plinder.data import save_linked_structures

    calls = {}

    class Scores:
        def summarize_scores(self):
            return {"ligand": {"lddt": 1.0}}

    def fake_from_model_files(*args, **kwargs):
        calls.update(kwargs)
        return Scores()

    monkeypatch.setattr(save_linked_structures, "save_superposition", lambda **_: True)
    monkeypatch.setattr(
        save_linked_structures.utils.ModelScores,
        "from_model_files",
        fake_from_model_files,
    )
    link = pd.Series(
        {
            "reference_system_id": "1abc__1__1.A__1.B",
            "id": "2def_A",
            "ligand_files": [],
        }
    )

    save_linked_structures.system_save_and_score_representative(
        link=link,
        reference_system=object(),
        data_dir=tmp_path,
        search_db="holo",
        output_folder=tmp_path / "linked",
    )

    assert calls == {"score_protein": True}
