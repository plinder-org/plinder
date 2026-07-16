import pandas as pd


def test_make_holo_links_uses_full_system_ids_and_no_ligand_paths(
    tmp_path, monkeypatch
):
    from plinder.data import save_linked_structures

    query_system = "1abc__1__1.A__1.B"
    target_system = "2def__1__1.C__1.D"
    cfg = save_linked_structures.LinkedStructureConfig(
        filter_criteria={
            "pocket_lddt": 0,
            "protein_fident_qcov_weighted_sum": 0,
        }
    )
    similarities = pd.DataFrame(
        {
            "query_system": [query_system, query_system],
            "target_system": [target_system, target_system],
            "metric": ["pocket_lddt", "protein_fident_qcov_weighted_sum"],
            "similarity": [100, 100],
        }
    )
    monkeypatch.setattr(
        save_linked_structures.scores,
        "query_protein_similarity",
        lambda **_: similarities,
    )
    monkeypatch.setattr(save_linked_structures, "get_resolution", lambda _: 2.0)
    output_file = tmp_path / "linked_structures" / "holo_links.parquet"

    save_linked_structures.make_linked_structures_data_file(
        data_dir=tmp_path,
        search_db="holo",
        superposed_folder=tmp_path / "linked_staging",
        output_file=output_file,
        cfg=cfg,
        num_processes=1,
    )

    links = pd.read_parquet(output_file)
    assert links["reference_system_id"].tolist() == [query_system]
    assert links["id"].tolist() == [target_system]
    assert "ligand_files" not in links
    assert links["receptor_file"].tolist() == [
        (
            tmp_path
            / "linked_staging"
            / "holo"
            / query_system
            / target_system
            / "superposed.cif"
        ).as_posix()
    ]


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

    class ReferenceSystem:
        ligand_sdfs = {"1.B": (tmp_path / "canonical.sdf").as_posix()}

    save_linked_structures.system_save_and_score_representative(
        link=link,
        reference_system=ReferenceSystem(),
        data_dir=tmp_path,
        search_db="holo",
        output_folder=tmp_path / "linked",
    )

    assert calls == {"score_protein": True}


def test_holo_cif_is_reconstructed_from_original_pdb_mmcif(tmp_path, monkeypatch):
    from plinder.data import save_linked_structures

    captured = {}

    class FakeSystem:
        def __init__(self, **kwargs):
            captured.update(kwargs)
            self.receptor_cif = tmp_path / "reconstructed" / "receptor.cif"

    monkeypatch.setattr(save_linked_structures, "PlinderSystem", FakeSystem)
    system_id = "1abc__1__1.A__1.B"
    receptor = save_linked_structures.get_cif_file(tmp_path, "holo", system_id)

    assert receptor == tmp_path / "reconstructed" / "receptor.cif"
    assert captured == {
        "system_id": system_id,
        "source_mmcif": (
            tmp_path
            / "ingest"
            / "ab"
            / "pdb_00001abc"
            / "pdb_00001abc_xyz-enrich.cif.gz"
        ),
        "reconstruction_dir": tmp_path / "reconstructed_systems" / system_id,
        "canonical_ligand_dir": (
            tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files"
        ),
    }
