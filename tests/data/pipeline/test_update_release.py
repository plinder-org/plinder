import json
import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from omegaconf import OmegaConf
from plinder.core.structure.smallmols_similarity import mol2morgan_fp
from plinder.core.utils import schemas
from plinder.core.utils.files import file_sha256, write_json_atomic
from plinder.data.annotations import get_similarity_scores
from plinder.data.pipeline import update_release
from rdkit import DataStructs


def _fingerprint(smiles: str) -> bytes:
    return DataStructs.BitVectToBinaryText(
        mol2morgan_fp(
            smiles,
            radius=get_similarity_scores.ECFP4_RADIUS,
            nbits=get_similarity_scores.ECFP4_NBITS,
        )
    )


def test_chemical_update_keeps_old_pairs_and_adds_both_directions(tmp_path):
    fingerprints = tmp_path / "fingerprints"
    fingerprints.mkdir()
    fingerprint_path = fingerprints / "ligands_per_smiles.parquet"
    get_similarity_scores.write_ecfp4_fingerprint_table(
        pd.DataFrame(
            {
                "ligand_smiles_id": [0, 1],
                "ligand_rdkit_canonical_smiles": ["CCO", "CCC"],
                "fingerprint": [_fingerprint("CCO"), _fingerprint("CCC")],
            }
        ),
        fingerprint_path,
    )
    source = tmp_path / "base-score.parquet"
    get_similarity_scores.ligand_scores(
        ligand_ids=[0, 1],
        data_dir=tmp_path,
        output_path=source,
        minimum_similarity=0,
    )
    scores = tmp_path / "ligand_scores"
    scores.mkdir()
    os.link(source, scores / "base.parquet")
    prior_manifest = get_similarity_scores.ligand_score_manifest_payload(
        fingerprint_path=fingerprint_path,
        fingerprint_col="fingerprint",
        metric="tanimoto_similarity_ecfp4_1024",
        minimum_similarity=0,
        number_id_col="ligand_smiles_id",
    )
    write_json_atomic(
        scores / get_similarity_scores.LIGAND_SCORE_MANIFEST,
        prior_manifest,
    )
    source_hash = file_sha256(source)
    get_similarity_scores.write_ecfp4_fingerprint_table(
        pd.DataFrame(
            {
                "ligand_smiles_id": [0, 2],
                "ligand_rdkit_canonical_smiles": ["CCO", "c1ccccc1"],
                "fingerprint": [_fingerprint("CCO"), _fingerprint("c1ccccc1")],
            }
        ),
        fingerprint_path,
    )
    current_manifest = get_similarity_scores.ligand_score_manifest_payload(
        fingerprint_path=fingerprint_path,
        fingerprint_col="fingerprint",
        metric="tanimoto_similarity_ecfp4_1024",
        minimum_similarity=0,
        number_id_col="ligand_smiles_id",
    )

    report = update_release._refresh_chemical_score_shards(
        data_dir=tmp_path,
        directory="ligand_scores",
        metric="tanimoto_similarity_ecfp4_1024",
        batch_size=10,
        number_id_col="ligand_smiles_id",
        minimum_similarity=0,
        scratch_dir=tmp_path / "scratch",
        prior_ligands=pd.DataFrame(
            {
                "ligand_smiles_id": [0, 1],
                "ligand_rdkit_canonical_smiles": ["CCO", "CCC"],
                "fingerprint": [_fingerprint("CCO"), _fingerprint("CCC")],
            }
        ),
        prior_manifest=prior_manifest,
        current_manifest=current_manifest,
    )

    result = pd.read_parquet(scores).sort_values(
        ["query_ligand_id", "target_ligand_id"], ignore_index=True
    )
    assert result[["query_ligand_id", "target_ligand_id"]].to_records(
        index=False
    ).tolist() == [(0, 0), (0, 2), (2, 0), (2, 2)]
    assert report["new_query_scores"] == 1
    assert report["removed_ligands"] == 1
    assert file_sha256(source) == source_hash
    assert json.loads((scores / "_manifest.json").read_text())["metric"] == (
        "tanimoto_similarity_ecfp4_1024"
    )

    resumed = update_release._refresh_chemical_score_shards(
        data_dir=tmp_path,
        directory="ligand_scores",
        metric="tanimoto_similarity_ecfp4_1024",
        batch_size=10,
        number_id_col="ligand_smiles_id",
        minimum_similarity=0,
        scratch_dir=tmp_path / "scratch",
        prior_ligands=pd.DataFrame(
            {
                "ligand_smiles_id": [0, 2],
                "ligand_rdkit_canonical_smiles": ["CCO", "c1ccccc1"],
                "fingerprint": [_fingerprint("CCO"), _fingerprint("c1ccccc1")],
            }
        ),
        prior_manifest=current_manifest,
        current_manifest=current_manifest,
    )
    assert resumed["new_query_scores"] == 0
    assert (
        pd.read_parquet(scores)
        .sort_values(["query_ligand_id", "target_ligand_id"], ignore_index=True)
        .equals(result)
    )

    (scores / get_similarity_scores.LIGAND_SCORE_MANIFEST).unlink()
    rebuilt = update_release._refresh_chemical_score_shards(
        data_dir=tmp_path,
        directory="ligand_scores",
        metric="tanimoto_similarity_ecfp4_1024",
        batch_size=10,
        number_id_col="ligand_smiles_id",
        minimum_similarity=0,
        scratch_dir=tmp_path / "scratch",
        prior_ligands=pd.DataFrame(
            {
                "ligand_smiles_id": [0, 2],
                "ligand_rdkit_canonical_smiles": ["CCO", "c1ccccc1"],
                "fingerprint": [_fingerprint("CCO"), _fingerprint("c1ccccc1")],
            }
        ),
        prior_manifest=current_manifest,
        current_manifest=current_manifest,
    )
    assert rebuilt["new_query_scores"] == 2
    assert (
        pd.read_parquet(scores)
        .sort_values(["query_ligand_id", "target_ligand_id"], ignore_index=True)
        .equals(result)
    )


def test_chemical_update_rescores_changed_fingerprint(tmp_path):
    fingerprints = tmp_path / "fingerprints"
    fingerprints.mkdir()
    fingerprint_path = fingerprints / "ligands_per_smiles.parquet"
    prior_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": [0],
            "ligand_rdkit_canonical_smiles": ["CCO"],
            "fingerprint": [_fingerprint("CCC")],
        }
    )
    get_similarity_scores.write_ecfp4_fingerprint_table(prior_ligands, fingerprint_path)
    scores = tmp_path / "ligand_scores"
    scores.mkdir()
    get_similarity_scores.ligand_scores(
        ligand_ids=[0],
        data_dir=tmp_path,
        output_path=scores / "base.parquet",
        minimum_similarity=0,
    )
    prior_manifest = get_similarity_scores.ligand_score_manifest_payload(
        fingerprint_path=fingerprint_path,
        fingerprint_col="fingerprint",
        metric="tanimoto_similarity_ecfp4_1024",
        minimum_similarity=0,
        number_id_col="ligand_smiles_id",
    )
    write_json_atomic(
        scores / get_similarity_scores.LIGAND_SCORE_MANIFEST,
        prior_manifest,
    )
    get_similarity_scores.write_ecfp4_fingerprint_table(
        pd.DataFrame(
            {
                "ligand_smiles_id": [0],
                "ligand_rdkit_canonical_smiles": ["CCO"],
                "fingerprint": [_fingerprint("CCO")],
            }
        ),
        fingerprint_path,
    )
    current_manifest = get_similarity_scores.ligand_score_manifest_payload(
        fingerprint_path=fingerprint_path,
        fingerprint_col="fingerprint",
        metric="tanimoto_similarity_ecfp4_1024",
        minimum_similarity=0,
        number_id_col="ligand_smiles_id",
    )

    report = update_release._refresh_chemical_score_shards(
        data_dir=tmp_path,
        directory="ligand_scores",
        metric="tanimoto_similarity_ecfp4_1024",
        batch_size=10,
        number_id_col="ligand_smiles_id",
        minimum_similarity=0,
        scratch_dir=tmp_path / "scratch",
        prior_ligands=prior_ligands,
        prior_manifest=prior_manifest,
        current_manifest=current_manifest,
    )

    assert report["new_query_scores"] == 1
    assert pd.read_parquet(scores)[["query_ligand_id", "target_ligand_id"]].to_records(
        index=False
    ).tolist() == [(0, 0)]


def test_remove_affected_shape_scores(tmp_path):
    cache = tmp_path / "scores/ligand_3d_by_query/ab.parquet"
    cache.parent.mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "query_entry": "1abc",
                    "query_ligand_asym_id": "L",
                    "target_entry": "2def",
                    "target_ligand_asym_id": "M",
                    "shape": 0.5,
                    "color": 0.4,
                    "sucos_shape": 0.45,
                },
                {
                    "query_entry": "2def",
                    "query_ligand_asym_id": "M",
                    "target_entry": "1abc",
                    "target_ligand_asym_id": "L",
                    "shape": 0.5,
                    "color": 0.4,
                    "sucos_shape": 0.45,
                },
                {
                    "query_entry": "2def",
                    "query_ligand_asym_id": "M",
                    "target_entry": "3ghi",
                    "target_ligand_asym_id": "N",
                    "shape": 0.6,
                    "color": 0.5,
                    "sucos_shape": 0.55,
                },
            ],
            schema=schemas.LIGAND_3D_SCORE_SCHEMA,
        ),
        cache,
    )

    removed = update_release._remove_affected_shape_scores(
        tmp_path,
        shards=["ab"],
        affected={"1abc"},
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    assert removed == 2
    assert pd.read_parquet(cache)[["query_entry", "target_entry"]].to_dict(
        "records"
    ) == [{"query_entry": "2def", "target_entry": "3ghi"}]


def test_holo_repair_discards_interrupted_shape_batches(tmp_path, monkeypatch):
    partial = tmp_path / "scores/ligand_3d_pair_repairs/0.parquet"
    partial.parent.mkdir(parents=True)
    partial.write_text("partial")
    monkeypatch.setattr(
        update_release.score,
        "plan_score_batches",
        lambda *_args, **_kwargs: None,
    )

    def no_active_queries(*args, **kwargs):
        raise ValueError("no active score queries are affected by the repair")

    monkeypatch.setattr(update_release.score, "plan_score_repair", no_active_queries)

    report = update_release.repair_holo_scores(
        tmp_path,
        base_data_dir=tmp_path / "base",
        affected={"1abc"},
        full_alignment_queries=set(),
        scorer_cfg=OmegaConf.create(
            {
                "max_query_protein_chains": 40,
                "max_query_proper_ligand_chains": 20,
            }
        ),
        scratch_dir=tmp_path / "scratch",
        threads=1,
        memory_limit="1GB",
        score_batch_size=10,
        ligand_batch_size=10,
    )

    assert report == {"status": "complete", "query_count": 0, "pair_count": 0}
    assert not partial.parent.exists()


def test_ligand_score_export_removes_stale_shards(tmp_path, monkeypatch):
    shard_dir = tmp_path / "exports" / "ligand_similarity_score_shards"
    shard_dir.mkdir(parents=True)
    for suffix in ("parquet", "json"):
        (shard_dir / f"stale.{suffix}").write_text("stale")

    monkeypatch.setattr(
        update_release.score,
        "published_scoring_query_ids",
        lambda _data_dir: ["1abc"],
    )

    def export(*args, **kwargs):
        del args
        assert kwargs["shards"] == ["ab"]
        assert not (shard_dir / "stale.parquet").exists()
        assert not (shard_dir / "stale.json").exists()

    monkeypatch.setattr(
        update_release.score,
        "export_ligand_similarity_scores_batch",
        export,
    )
    monkeypatch.setattr(
        update_release.score,
        "finalize_ligand_similarity_scores",
        lambda *args, **kwargs: {"status": "complete"},
    )

    assert update_release.refresh_ligand_score_export(
        tmp_path,
        scratch_dir=tmp_path / "scratch",
        threads=1,
        memory_limit="1GB",
    ) == {"status": "complete"}


def test_weekly_stage_report_is_strictly_ordered(tmp_path):
    path = tmp_path / "weekly_update.json"
    report = {"status": "running", "inputs": {"plan": "x"}, "completed_stages": []}
    update_release._finish_stage(
        path, report, "entries_and_archives", {"status": "complete"}
    )
    loaded = update_release._read_weekly_report(path, inputs={"plan": "x"})
    assert loaded is not None
    assert loaded["completed_stages"] == ["entries_and_archives"]


def test_base_reuse_skips_transient_staging_trees(tmp_path):
    base = tmp_path / "base"
    workspace = tmp_path / "workspace"
    durable = base / "index" / "annotation_table.parquet"
    stale = base / "index" / ".staging" / "old" / "entries.parquet"
    temporary = base / "index" / "annotation_table.parquet.tmp"
    installing = base / "dbs" / ".foldseek.installing" / "foldseek"
    partial_shape = base / "scores/ligand_3d_pair_repairs/0.parquet"
    for path in (durable, stale, temporary, installing, partial_shape):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(path.name)

    update_release._reuse_base_artifacts(base, workspace)

    assert (workspace / "index" / durable.name).read_text() == durable.name
    assert not (workspace / "index" / ".staging").exists()
    assert not (workspace / "index" / temporary.name).exists()
    assert not (workspace / "dbs" / ".foldseek.installing").exists()
    assert not (workspace / "scores/ligand_3d_pair_repairs").exists()


def test_alignment_repair_queries_are_saved_before_restart(tmp_path, monkeypatch):
    state_path = tmp_path / "weekly_update.json"
    report = {
        "status": "running",
        "inputs": {"plan": "x"},
        "completed_stages": ["entries_and_archives", "search_inputs"],
    }
    planned = {"holo": {"1abc", "2def"}, "pred": {"1abc"}}
    calls: list[str] = []

    def plan(*_args, search_db, **_kwargs):
        calls.append(search_db)
        return planned[search_db]

    monkeypatch.setattr(update_release, "_plan_alignment_repairs", plan)
    result = update_release._load_or_plan_alignment_repairs(
        tmp_path,
        search_databases=["holo", "pred"],
        affected={"1abc"},
        scorer_cfg=OmegaConf.create({}),
        foldseek_cfg=OmegaConf.create({}),
        mmseqs_cfg=OmegaConf.create({}),
        scratch_dir=tmp_path / "scratch",
        threads=1,
        batch_size=10,
        state_path=state_path,
        report=report,
    )

    assert result == planned
    assert calls == ["holo", "pred"]
    saved = json.loads(state_path.read_text())
    assert saved["alignment_repair_queries"] == {
        "holo": ["1abc", "2def"],
        "pred": ["1abc"],
    }

    monkeypatch.setattr(
        update_release,
        "_plan_alignment_repairs",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("restart must use the saved repair queries")
        ),
    )
    assert (
        update_release._load_or_plan_alignment_repairs(
            tmp_path,
            search_databases=["holo", "pred"],
            affected=set(),
            scorer_cfg=OmegaConf.create({}),
            foldseek_cfg=OmegaConf.create({}),
            mmseqs_cfg=OmegaConf.create({}),
            scratch_dir=tmp_path / "scratch",
            threads=1,
            batch_size=10,
            state_path=state_path,
            report=saved,
        )
        == planned
    )


def test_pred_alignment_repairs_only_changed_queries(tmp_path, monkeypatch):
    manifest = tmp_path / update_release.score.MANIFEST_RELATIVE
    manifest.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": ["1abc", "2def"]}).to_parquet(manifest, index=False)

    def unexpected(*_args, **_kwargs):
        raise AssertionError("pred targets do not change with a PDB update")

    monkeypatch.setattr(update_release, "_queries_with_old_target_hits", unexpected)
    monkeypatch.setattr(update_release, "_changed_target_database", unexpected)
    monkeypatch.setattr(update_release, "_search_changed_targets", unexpected)

    assert update_release._plan_alignment_repairs(
        tmp_path,
        search_db="pred",
        affected={"1abc", "9xyz"},
        scorer_cfg=OmegaConf.create({}),
        foldseek_cfg=OmegaConf.create({}),
        mmseqs_cfg=OmegaConf.create({}),
        scratch_dir=tmp_path / "scratch",
        threads=1,
        batch_size=10,
    ) == {"1abc"}


def test_release_update_runs_only_configured_search_databases(tmp_path, monkeypatch):
    plan_dir = tmp_path / "plan"
    plan_dir.mkdir()
    (plan_dir / "plan.json").write_text("{}")
    entries = pd.DataFrame({"pdb_id": ["1abc"], "action": ["revised"]})
    entries.to_parquet(plan_dir / "entries.parquet", index=False)
    config_path = tmp_path / "config.yaml"
    config_path.write_text("{}")
    validation_root = tmp_path / "validation"
    validation_root.mkdir()
    base = tmp_path / "base"
    base.mkdir()
    workspace = tmp_path / "workspace"
    (workspace / "index").mkdir(parents=True)
    state_path = workspace / "weekly_update.json"
    write_json_atomic(
        state_path,
        {
            "status": "running",
            "inputs": update_release._weekly_inputs(
                plan_dir, config_path, validation_root
            ),
            "base_release": str(base),
            "affected_pdb_ids": ["1abc"],
            "completed_stages": ["entries_and_archives"],
        },
    )
    cfg = OmegaConf.create(
        {
            "scorer": {
                "sub_databases": ["holo", "pred"],
                "max_query_protein_chains": 30,
                "max_query_proper_ligand_chains": 30,
            },
            "foldseek": {"max_seqs": 100},
            "mmseqs": {},
            "flow": {
                "protein_sequence_cluster_identity": 0.4,
                "protein_structure_cluster_lddt": 0.7,
                "protein_cluster_coverage": 0.8,
                "make_batch_scores_batch_size": 10,
                "make_ligand_3d_scores_batch_size": 10,
            },
        }
    )
    calls: dict[str, list] = {
        "make_dbs": [],
        "make_sub_dbs": [],
        "plan": [],
        "alignments": [],
        "scores": [],
    }

    monkeypatch.setattr(update_release, "_load_configuration", lambda _path: cfg)
    monkeypatch.setattr(
        update_release,
        "load_update_plan",
        lambda _path: (
            {"data_dir": str(base), "nextgen_root": str(tmp_path / "nextgen")},
            entries,
        ),
    )
    monkeypatch.setattr(
        update_release.tasks,
        "make_alignment_chain_lookup",
        lambda **_kwargs: tmp_path / "lookup.parquet",
    )
    monkeypatch.setattr(
        update_release.score, "plan_protein_scoring", lambda *_args, **_kwargs: {}
    )
    monkeypatch.setattr(
        update_release.score,
        "make_foldseek_input_manifest",
        lambda *_args, **_kwargs: tmp_path / "foldseek.tsv",
    )
    monkeypatch.setattr(
        update_release.score,
        "make_mmseqs_input_fasta",
        lambda *_args, **_kwargs: tmp_path / "mmseqs.fasta",
    )
    monkeypatch.setattr(
        update_release.protein_clusters,
        "make_protein_sequence_clusters",
        lambda **_kwargs: tmp_path / "sequence.parquet",
    )
    monkeypatch.setattr(
        update_release.protein_clusters,
        "make_protein_structure_clusters",
        lambda **_kwargs: tmp_path / "structure.parquet",
    )
    monkeypatch.setattr(
        update_release.tasks,
        "make_dbs",
        lambda **kwargs: calls["make_dbs"].append(kwargs["sub_databases"]),
    )
    monkeypatch.setattr(
        update_release.tasks,
        "make_sub_dbs",
        lambda **kwargs: calls["make_sub_dbs"].append(kwargs["sub_databases"]),
    )
    monkeypatch.setattr(
        update_release.score,
        "plan_linked_apo_scoring",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("apo is not configured")
        ),
    )

    def plan_alignments(*_args, search_db, **_kwargs):
        calls["plan"].append(search_db)
        return {"1abc"}

    def repair_alignments(*_args, search_db, full_queries, **_kwargs):
        calls["alignments"].append(search_db)
        return full_queries

    monkeypatch.setattr(update_release, "_plan_alignment_repairs", plan_alignments)
    monkeypatch.setattr(update_release, "repair_alignments", repair_alignments)
    monkeypatch.setattr(
        update_release.score,
        "finalize_alignment_artifacts",
        lambda *_args, **_kwargs: {"status": "complete"},
    )
    monkeypatch.setattr(
        update_release,
        "repair_holo_scores",
        lambda *_args, **_kwargs: calls["scores"].append("holo") or {},
    )
    monkeypatch.setattr(
        update_release,
        "repair_interface_scores",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        update_release,
        "refresh_ligand_score_export",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        update_release,
        "repair_apo_scores",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("apo is not configured")
        ),
    )

    def repair_non_holo(*_args, search_db, **_kwargs):
        calls["scores"].append(search_db)
        return {}

    monkeypatch.setattr(update_release, "repair_non_holo_scores", repair_non_holo)
    monkeypatch.setattr(
        update_release,
        "refresh_ligand_chemistry",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        update_release, "rebuild_similarity_covers", lambda *_args, **_kwargs: {}
    )

    def finalize_index(*, data_dir):
        write_json_atomic(data_dir / "index" / "collation.json", {"status": "complete"})

    monkeypatch.setattr(update_release.tasks, "finalize_index", finalize_index)

    result = update_release.apply_release_update(
        plan_dir,
        workspace,
        validation_root=validation_root,
        config_path=config_path,
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    assert calls == {
        "make_dbs": [["holo"]],
        "make_sub_dbs": [["holo"]],
        "plan": ["holo", "pred"],
        "alignments": ["holo", "pred"],
        "scores": ["holo", "pred"],
    }
    assert result["alignments"]["full_queries"] == {
        "holo": ["1abc"],
        "pred": ["1abc"],
    }
    assert set(result["scores"]) == {"ligand", "interface", "ligand_export", "pred"}
