import gzip
import json
import os
import threading
from types import SimpleNamespace

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from omegaconf import OmegaConf
from rdkit import DataStructs

from plinder.core.structure.smallmols_similarity import mol2morgan_fp
from plinder.core.utils import schemas
from plinder.core.utils.files import file_sha256, write_json_atomic
from plinder.data.annotations import get_similarity_scores
from plinder.data.pipeline import update_entries, update_release
from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig


def _fingerprint(smiles: str) -> bytes:
    return DataStructs.BitVectToBinaryText(
        mol2morgan_fp(
            smiles,
            radius=get_similarity_scores.ECFP4_RADIUS,
            nbits=get_similarity_scores.ECFP4_NBITS,
        )
    )


def test_weekly_table_merge_counts_written_rows(tmp_path, monkeypatch):
    monkeypatch.setattr(
        update_entries,
        "ENTRY_TABLES",
        {
            "entry_metadata": ("entry_metadata.parquet", "entry_pdb_id"),
            "entry_sources": ("entry_sources.parquet", "entry_pdb_id"),
        },
    )
    base = tmp_path / "base/index"
    incoming = tmp_path / "incoming/index"
    base.mkdir(parents=True)
    incoming.mkdir(parents=True)
    for name in ("entry_metadata", "entry_sources"):
        pd.DataFrame({"entry_pdb_id": ["1abc", "2def"], "value": [1, 2]}).to_parquet(
            base / f"{name}.parquet"
        )
        pd.DataFrame({"entry_pdb_id": ["1abc"], "value": [3]}).to_parquet(
            incoming / f"{name}.parquet"
        )
    stage = tmp_path / "stage"

    counts = update_entries._combine_tables(
        base.parent,
        incoming.parent,
        stage,
        ["1abc"],
        threads=2,
        memory_limit="1GB",
        scratch_dir=tmp_path / "scratch",
    )

    assert counts == {"entry_metadata": 2, "entry_sources": 2}
    for name in counts:
        assert pd.read_parquet(stage / f"{name}.parquet").to_dict("records") == [
            {"entry_pdb_id": "1abc", "value": 3},
            {"entry_pdb_id": "2def", "value": 2},
        ]


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
    other_cache = cache.with_name("cd.parquet")
    pq.write_table(pq.read_table(cache), other_cache)

    removed = update_release._remove_affected_shape_scores(
        tmp_path,
        shards=["ab", "cd"],
        affected={"1abc"},
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )

    assert removed == 4
    for path in (cache, other_cache):
        assert pd.read_parquet(path)[["query_entry", "target_entry"]].to_dict(
            "records"
        ) == [{"query_entry": "2def", "target_entry": "3ghi"}]


def test_shape_repair_uses_changed_ligand_archives(tmp_path):
    base = tmp_path / "base"
    updated = tmp_path / "updated"
    for root, changed_sdf in ((base, b"old"), (updated, b"new")):
        archive = root / "ligand_archives/qz.parquet"
        archive.parent.mkdir(parents=True)
        pd.DataFrame(
            {
                "pdb_id": ["1qz5", "2qz6"],
                "ligand_asym_id": ["B", "B"],
                "sdf": [b"same", changed_sdf],
                "pharmacophore_features": [None, None],
            }
        ).to_parquet(archive, index=False)

    assert update_release._changed_ligand_archive_entries(
        base, updated, {"1qz5", "2qz6", "3qz7"}
    ) == {"2qz6"}


def test_weekly_metadata_only_revision_reuses_scoring_inputs(tmp_path):
    source = tmp_path / "pdb_00001abc/pdb_00001abc_xyz-enrich.cif.gz"
    source.parent.mkdir()
    source.write_bytes(gzip.compress(b"coordinates", mtime=1))
    roots = [tmp_path / "base", tmp_path / "updated"]
    for root in roots:
        (root / "index").mkdir(parents=True)
        (root / "manifests").mkdir()
        (root / "manifests/foldseek_createdb_inputs.tsv").write_text(f"{source}\n")
        pd.DataFrame(
            {
                "entry_pdb_id": ["1abc"],
                "chain_asym_id": ["A"],
                "chain_auth_id": ["A"],
                "chain_type": ["polypeptide(L)"],
                "chain_receptor_type": ["protein"],
                "chain_sequence": ["ACD"],
                "chain_sequence_noncanonical": ["ACD"],
                "chain_is_holo": [True],
                "chain_is_ligand_like": [False],
            }
        ).to_parquet(root / "index/entry_chains.parquet")
        pd.DataFrame(
            {
                "entry_pdb_id": ["1abc"],
                "biounit_id": ["1"],
                "chain_instance": ["1.A"],
                "chain_asym_id": ["A"],
                "chain_role": ["receptor"],
                "chain_num_contacting_ions": [0],
                "chain_num_contacting_artifacts": [0],
                "chain_num_contacting_other_ligands": [1],
                "chain_num_contacting_proteins": [0],
            }
        ).to_parquet(root / "index/entry_biounit_chains.parquet")
        pd.DataFrame(
            {
                "entry_pdb_id": ["1abc"],
                "ligand_id": ["1abc__1__L"],
                "system_id": ["1abc__1__A__L"],
                "ligand_smiles": ["CCO"],
                "ligand_is_proper": [True],
                "ligand_neighboring_residues": [["A_1"]],
                "ligand_interactions": [[]],
            }
        ).to_parquet(root / "index/annotation_table.parquet")
        pd.DataFrame(
            {
                "entry_pdb_id": ["1abc"],
                "system_id": ["1abc__1__A__B"],
                "interface_chain_1_residue_numbers": [[1]],
                "interface_chain_2_residue_numbers": [[2]],
            }
        ).to_parquet(root / "index/interface_annotation_table.parquet")

    base, updated = roots
    assert update_release._entries_with_interfaces(base, updated, {"1abc"}) == {"1abc"}
    assert not update_release._entries_with_interfaces(base, updated, {"2def"})
    assert update_release._unchanged_scoring_entries(base, updated, {"1abc"}) == {
        "1abc"
    }
    recompressed = tmp_path / "other/pdb_00001abc/pdb_00001abc_xyz-enrich.cif.gz"
    recompressed.parent.mkdir(parents=True)
    recompressed.write_bytes(gzip.compress(b"coordinates", mtime=2))
    (updated / "manifests/foldseek_createdb_inputs.tsv").write_text(f"{recompressed}\n")
    assert source.read_bytes() != recompressed.read_bytes()
    assert update_release._unchanged_scoring_entries(base, updated, {"1abc"}) == {
        "1abc"
    }
    recompressed.write_bytes(gzip.compress(b"different coordinates", mtime=2))
    assert not update_release._unchanged_scoring_entries(base, updated, {"1abc"})
    recompressed.write_bytes(gzip.compress(b"coordinates", mtime=2))
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["A"],
            "chain_type": ["polypeptide(L)"],
            "chain_receptor_type": ["protein"],
            "chain_sequence": ["ACE"],
            "chain_sequence_noncanonical": ["ACE"],
            "chain_is_holo": [True],
            "chain_is_ligand_like": [False],
        }
    ).to_parquet(updated / "index/entry_chains.parquet")
    assert not update_release._unchanged_scoring_entries(base, updated, {"1abc"})


def test_weekly_reused_foldseek_manifest_uses_current_snapshot(tmp_path):
    old = (
        tmp_path
        / "old/data/entries/divided/ab/pdb_00001abc/pdb_00001abc_xyz-enrich.cif.gz"
    )
    current_root = tmp_path / "current"
    current = current_root / old.relative_to(tmp_path / "old")
    old.parent.mkdir(parents=True)
    current.parent.mkdir(parents=True)
    old.write_bytes(gzip.compress(b"same cif", mtime=1))
    current.write_bytes(gzip.compress(b"same cif", mtime=2))
    data_dir = tmp_path / "release"
    (data_dir / "manifests").mkdir(parents=True)
    manifest = data_dir / "manifests/foldseek_createdb_inputs.tsv"
    manifest.write_text(f"{old}\n")

    update_release._refresh_foldseek_source_manifest(
        data_dir, nextgen_root=current_root, affected={"1abc"}
    )

    assert manifest.read_text() == f"{current}\n"


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
        targeted_alignment_queries=set(),
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


def test_weekly_score_repair_runs_independent_batches(tmp_path):
    for pdb_id in ("1abc", "2def"):
        score_file = tmp_path / "dbs/subdbs/search_db=holo" / f"{pdb_id}.parquet"
        score_file.parent.mkdir(parents=True, exist_ok=True)
        score_file.touch()

    update_release._repair_score_batches(
        tmp_path,
        pd.DataFrame(
            {
                "repair_batch_index": [0, 1],
                "pdb_id": ["1abc", "2def"],
                "repair_mode": ["drop", "drop"],
            }
        ),
        OmegaConf.create({"sub_databases": ["holo"]}),
        tmp_path / "scratch",
        threads=8,
    )

    assert not list((tmp_path / "dbs/subdbs/search_db=holo").glob("*.parquet"))


def test_completed_weekly_query_repair_survives_a_resume(tmp_path):
    manifests = tmp_path / "manifests"
    affected = update_release._write_pdb_manifest(
        manifests / "weekly_affected_entries.parquet", {"1abc"}
    )
    full = update_release._write_pdb_manifest(
        manifests / "weekly_full_score_queries.parquet", {"1abc"}
    )
    target = update_release._write_pdb_manifest(
        manifests / "weekly_target_score_queries.parquet", {"2def"}
    )
    repair = manifests / "score_repair.parquet"
    pd.DataFrame({"pdb_id": ["1abc"]}).to_parquet(repair, index=False)
    plan = {
        "affected_manifest": update_release.score._source_signature(affected),
        "additional_full_query_manifest": update_release.score._source_signature(full),
        "additional_target_query_manifest": update_release.score._source_signature(
            target
        ),
        "output": update_release.score._source_signature(repair),
    }
    write_json_atomic(repair.with_suffix(".json"), plan)
    write_json_atomic(
        manifests / "score_repair_incomplete_queries.json", {"status": "complete"}
    )

    old_mtime = target.stat().st_mtime_ns
    update_release._write_pdb_manifest(target, {"2def"})
    assert target.stat().st_mtime_ns == old_mtime
    assert (
        update_release._completed_holo_query_repair(
            tmp_path,
            affected_manifest=affected,
            full_manifest=full,
            target_manifest=target,
        )
        == plan
    )

    update_release._write_pdb_manifest(target, {"3ghi"})
    assert (
        update_release._completed_holo_query_repair(
            tmp_path,
            affected_manifest=affected,
            full_manifest=full,
            target_manifest=target,
        )
        is None
    )


def test_weekly_alignment_shards_can_map_concurrently(tmp_path, monkeypatch):
    barrier = threading.Barrier(2)
    observed = set()

    def map_shard(*, shards, **_kwargs):
        observed.update(shards)
        barrier.wait(timeout=10)

    monkeypatch.setattr(update_release.tasks, "map_batch_alignments", map_shard)
    update_release._refresh_alignment_shards(
        tmp_path,
        search_db="holo",
        shards={"ab", "cd"},
        replacement_query_ids=set(),
        scorer_cfg=OmegaConf.create({"sub_databases": ["holo"]}),
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )

    assert observed == {"ab", "cd"}


def test_weekly_candidate_shards_can_pack_concurrently(tmp_path, monkeypatch):
    barrier = threading.Barrier(2)
    observed = {}

    def collate(*, shards, scratch_dir, threads, **_kwargs):
        observed[shards[0]] = (scratch_dir, threads)
        barrier.wait(timeout=10)

    monkeypatch.setattr(update_release.tasks, "collate_ligand_3d_candidates", collate)
    update_release._collate_repaired_candidate_shards(
        tmp_path,
        ["ab", "cd"],
        tmp_path / "scratch",
        threads=8,
    )

    assert observed == {
        "ab": (tmp_path / "scratch/ab", 4),
        "cd": (tmp_path / "scratch/cd", 4),
    }


def test_weekly_target_repair_restores_only_needed_query_caches(tmp_path):
    base = tmp_path / "base"
    workspace = tmp_path / "workspace"
    old_score = base / "dbs/subdbs/search_db=holo/1abc.parquet"
    old_score.parent.mkdir(parents=True)
    pd.DataFrame({"query_system": ["1abc__1"]}).to_parquet(old_score, index=False)
    for folder, schema in (
        ("ligand_3d_candidate_shards", schemas.LIGAND_3D_CANDIDATE_SCHEMA),
        ("ligand_pair_score_shards", schemas.LIGAND_PAIR_SCORE_SCHEMA),
    ):
        packed = base / "scores" / folder / "shard=ab.parquet"
        packed.parent.mkdir(parents=True)
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {"query_entry": "1abc", "target_entry": "2def"},
                    {"query_entry": "2abc", "target_entry": "3ghi"},
                ],
                schema=schema,
            ),
            packed,
        )

    update_release._restore_score_query_caches(
        base_data_dir=base,
        data_dir=workspace,
        query_ids={"1abc", "2abc"},
    )

    assert pd.read_parquet(workspace / "dbs/subdbs/search_db=holo/1abc.parquet")[
        "query_system"
    ].tolist() == ["1abc__1"]
    for folder in ("ligand_3d_candidates", "ligand_pair_scores"):
        query_file = (
            workspace / "scores" / folder / "search_db=holo/shard=ab/1abc.parquet"
        )
        assert pd.read_parquet(query_file)["query_entry"].tolist() == ["1abc"]
        assert not (
            workspace / "scores" / folder / "search_db=holo/shard=ab/2abc.parquet"
        ).exists()


def test_ligand_score_export_removes_stale_shards(tmp_path, monkeypatch):
    shard_dir = tmp_path / "exports" / "ligand_similarity_scores"
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
    old_shape_batch = base / "scores/ligand_3d_pairs/0.parquet"
    old_raw_alignment = base / "dbs/subdbs/holo_mmseqs/aln/1abc.parquet"
    old_holo_score = base / "dbs/subdbs/search_db=holo/1abc.parquet"
    old_apo_score = base / "dbs/subdbs/search_db=apo/1abc.parquet"
    for path in (
        durable,
        stale,
        temporary,
        installing,
        partial_shape,
        old_shape_batch,
        old_raw_alignment,
        old_holo_score,
        old_apo_score,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(path.name)

    update_release._reuse_base_artifacts(base, workspace)

    assert (workspace / "index" / durable.name).read_text() == durable.name
    assert not (workspace / "index" / ".staging").exists()
    assert not (workspace / "index" / temporary.name).exists()
    assert not (workspace / "dbs" / ".foldseek.installing").exists()
    assert not (workspace / "scores/ligand_3d_pair_repairs").exists()
    assert not (workspace / "scores/ligand_3d_pairs").exists()
    assert not (workspace / "dbs/subdbs/holo_mmseqs/aln").exists()
    assert not (workspace / "dbs/subdbs/search_db=holo").exists()
    assert not (workspace / "dbs/subdbs/search_db=apo").exists()


def test_weekly_apo_scores_replace_only_changed_queries(tmp_path):
    scores = tmp_path / "scores/search_db=apo/apo.parquet"
    scores.parent.mkdir(parents=True)
    source_dir = tmp_path / "dbs/subdbs/search_db=apo"
    source_dir.mkdir(parents=True)

    def row(query: str, similarity: int) -> dict[str, object]:
        return {
            "query_system": f"{query}__1__1.A__1.B",
            "query_ligand_id": f"{query}__1__1.B",
            "target_system": "9xyz__1__1.A__1.B",
            "target_ligand_id": "9xyz__1__1.B",
            "protein_mapping": "1.A:1.A",
            "mapping": "1.B:1.B",
            "protein_mapper": "mmseqs",
            "source": "mmseqs",
            "metric": "pocket_fident",
            "similarity": similarity,
        }

    pd.DataFrame([row("1abc", 40), row("2def", 80)]).assign(search_db="apo").to_parquet(
        scores, index=False
    )
    pd.DataFrame([row("1abc", 90)]).to_parquet(source_dir / "1abc.parquet", index=False)
    update_release._merge_apo_scores(
        tmp_path,
        output=scores,
        new_queries=["1abc"],
        replaced_queries={"1abc", "3ghi"},
        scratch_dir=tmp_path / "scratch",
        threads=1,
        memory_limit="1GB",
    )
    result = pq.ParquetFile(scores).read().to_pandas().sort_values("query_system")
    assert result["query_system"].tolist() == [
        "1abc__1__1.A__1.B",
        "2def__1__1.A__1.B",
    ]
    assert result["similarity"].tolist() == [90, 80]
    assert result["search_db"].tolist() == ["apo", "apo"]


def test_weekly_chemistry_reuses_unchanged_ligands(tmp_path):
    base = tmp_path / "base"
    updated = tmp_path / "updated"
    rows = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "ligand_id": ["1abc__1__1.B", "2def__1__1.B"],
            "ligand_smiles": ["CCO", "CCN"],
            "ligand_is_proper": [True, True],
            "system_type": ["holo", "holo"],
        }
    )
    for root in (base, updated):
        (root / "index").mkdir(parents=True)
        rows.to_parquet(root / "index/annotation_table.parquet", index=False)

    assert update_release._ligand_chemistry_unchanged(base, updated, {"1abc"})
    changed = rows.copy()
    changed.loc[0, "ligand_smiles"] = "CCC"
    changed.to_parquet(updated / "index/annotation_table.parquet", index=False)
    assert not update_release._ligand_chemistry_unchanged(base, updated, {"1abc"})


def test_old_target_discovery_reads_retained_cigars(tmp_path):
    cigar = (
        tmp_path
        / "alignment_cigars/search_db=holo/alignment_type=mmseqs/shard=ab.parquet"
    )
    cigar.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc", "2abc", "3abc"],
            "target_entry": ["9xyz", "8xyz", "9xyz"],
        }
    ).to_parquet(cigar, index=False)

    assert update_release._queries_with_old_target_hits(
        tmp_path, search_db="holo", affected={"9xyz"}
    ) == {"1abc", "3abc"}


def test_changed_target_probe_searches_in_reverse_once_per_backend(
    tmp_path, monkeypatch
):
    for backend in ("mmseqs", "foldseek"):
        marker = tmp_path / "changed" / f"holo_{backend}" / "exact_cluster.json"
        marker.parent.mkdir(parents=True)
        marker.write_text("{}")
    scorer = SimpleNamespace(
        source_to_full_db_file={
            "holo_mmseqs": tmp_path / "full_mmseqs",
            "holo_foldseek": tmp_path / "full_foldseek",
        },
        get_config=lambda _db, backend: (
            MMSeqsConfig() if backend == "mmseqs" else FoldseekConfig()
        ),
    )
    monkeypatch.setattr(
        update_release.utils,
        "get_scorer",
        lambda **_kwargs: (scorer, [], tmp_path / "probe"),
    )
    monkeypatch.setattr(
        update_release.databases,
        "exact_search_database_paths",
        lambda *_args: (tmp_path / "changed_target", tmp_path / "unused", None),
    )
    calls = []

    def run_alignment(**kwargs):
        calls.append(kwargs["aln_type"])
        assert kwargs["query_ids_only"]
        assert kwargs["id_column"] == "target"
        assert kwargs["query_db"] == tmp_path / "changed_target"
        assert kwargs["target_db"] == tmp_path / f"full_{kwargs['aln_type']}"
        assert kwargs["search_target_db"] == kwargs["target_db"]
        assert kwargs["alignment_config"].evalue == 2.0
        assert kwargs["alignment_config"].max_seqs == 2_000_000
        return (
            {"1abc_A", "9xyz_X"}
            if kwargs["aln_type"] == "mmseqs"
            else {"pdb_00002def_xyz-enrich_B"}
        )

    monkeypatch.setattr(
        update_release.get_similarity_scores, "run_alignment", run_alignment
    )
    hits = update_release._search_changed_targets(
        tmp_path,
        search_db="holo",
        query_ids=["1abc", "2def"],
        target_database_dir=tmp_path / "changed",
        scorer_cfg=OmegaConf.create({}),
        foldseek_cfg=OmegaConf.create({}),
        mmseqs_cfg=OmegaConf.create({}),
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )

    assert calls == ["mmseqs", "foldseek"]
    assert hits == {"1abc", "2def"}


def test_alignment_query_only_mode_skips_residue_output(tmp_path, monkeypatch):
    monkeypatch.setattr(
        get_similarity_scores,
        "_alignment_search_command",
        lambda **_kwargs: ["mmseqs", "search"],
    )

    def check_call(command, **_kwargs):
        if command[1] == "convertalis":
            assert command[command.index("--format-output") + 1] == "query"
            (tmp_path / "hits.tsv").write_text("query\n1abc_A\n1abc_A\n2def_B\n")

    monkeypatch.setattr(get_similarity_scores.subprocess, "check_call", check_call)
    identifiers = get_similarity_scores.run_alignment(
        aln_type="mmseqs",
        query_db=tmp_path / "query",
        target_db=tmp_path / "target",
        search_target_db=tmp_path / "target",
        search_db=tmp_path / "search",
        aln_file=tmp_path / "hits.tsv",
        alignment_config=MMSeqsConfig(),
        tmp_dir=tmp_path / "tmp",
        remove_tmp=False,
        query_ids_only=True,
    )

    assert identifiers == {"1abc_A", "2def_B"}
    assert not (tmp_path / "hits.parquet").exists()


def test_alignment_target_only_mode_returns_target_ids(tmp_path, monkeypatch):
    def search_command(**kwargs):
        assert kwargs["expand_exact_clusters"] is False
        return ["foldseek", "search"]

    monkeypatch.setattr(
        get_similarity_scores,
        "_alignment_search_command",
        search_command,
    )

    def check_call(command, **_kwargs):
        if command[1] == "convertalis":
            assert command[command.index("--format-output") + 1] == "target"
            (tmp_path / "hits.tsv").write_text("target\n1abc_A\n1abc_A\n2def_B\n")

    monkeypatch.setattr(get_similarity_scores.subprocess, "check_call", check_call)
    identifiers = get_similarity_scores.run_alignment(
        aln_type="foldseek",
        query_db=tmp_path / "query",
        target_db=tmp_path / "target",
        search_target_db=tmp_path / "target",
        search_db=tmp_path / "search",
        aln_file=tmp_path / "hits.tsv",
        alignment_config=FoldseekConfig(),
        tmp_dir=tmp_path / "tmp",
        query_ids_only=True,
        id_column="target",
    )

    assert identifiers == {"1abc_A", "2def_B"}


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


def test_unchanged_alignment_rebase_accepts_public_mapping_schema(
    tmp_path,
    monkeypatch,
):
    inputs = {"foldseek": [{"name": "1abc.parquet"}], "mmseqs": []}
    output = update_release.tasks._alignment_release_path(
        data_dir=tmp_path,
        search_db="holo",
        alignment_type="foldseek",
        shard="ab",
    )
    output.parent.mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist(
            [],
            schema=schemas.release_alignment_mapping_schema(alignment_type="foldseek"),
        ),
        output,
    )
    stat = output.stat()
    manifest = tmp_path / "alignments/manifests/shard=ab.json"
    manifest.parent.mkdir(parents=True)
    write_json_atomic(
        manifest,
        {
            "inputs": inputs,
            "outputs": {
                "foldseek": {
                    "name": output.name,
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                },
                "mmseqs": None,
            },
        },
    )
    lookup = {"name": "alignment_chain_lookup.parquet"}
    monkeypatch.setattr(
        update_release.tasks,
        "_completed_alignment_chain_lookup",
        lambda _data_dir: lookup,
    )
    monkeypatch.setattr(
        update_release.tasks,
        "_alignment_input_signatures",
        lambda **_kwargs: inputs,
    )

    update_release._rebase_unchanged_alignment_manifests(
        tmp_path,
        search_db="holo",
        repaired_shards=set(),
    )

    assert json.loads(manifest.read_text())["alignment_chain_lookup"] == lookup


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
        "overlay": [],
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
        update_release,
        "extend_protein_clusters",
        lambda **kwargs: tmp_path / f"{kwargs['backend']}.parquet",
    )
    monkeypatch.setattr(
        update_release,
        "_prepare_weekly_search_overlay",
        lambda **kwargs: calls["overlay"].append(kwargs["search_databases"])
        or {"1abc"},
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
        "_isolated_ligand_score_export",
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
        update_release, "_ligand_chemistry_unchanged", lambda *_args: False
    )
    monkeypatch.setattr(
        update_release,
        "refresh_ligand_chemistry",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        update_release,
        "extend_ligand_clusters",
        lambda **_kwargs: pd.DataFrame({"ligand_id": ["1abc__1__1.X"]}),
    )
    monkeypatch.setattr(
        update_release,
        "extend_interface_clusters",
        lambda **_kwargs: pd.DataFrame({"system_id": []}),
    )

    def finalize_index(*, data_dir, **_kwargs):
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
        "overlay": [["holo"]],
        "plan": ["holo", "pred"],
        "alignments": ["holo", "pred"],
        "scores": ["holo", "pred"],
    }
    assert result["alignments"]["full_queries"] == {
        "holo": ["1abc"],
        "pred": ["1abc"],
    }
    assert set(result["scores"]) == {"ligand", "interface", "ligand_export", "pred"}
