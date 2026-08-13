# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import ast
import json
import os
from pathlib import Path
from shutil import copyfile
from types import SimpleNamespace

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from plinder.core.utils import schemas
from plinder.data.pipeline import io, tasks
from plinder.data.pipeline.config import LigandConfig


def _write_alignment_chain_lookup(data_dir: Path) -> None:
    index = data_dir / "index"
    index.mkdir(exist_ok=True, parents=True)
    for name in ["annotation_table.parquet", "entry_chains.parquet"]:
        path = index / name
        if not path.is_file():
            pd.DataFrame({"entry_pdb_id": ["1abc"]}).to_parquet(path, index=False)
    lookup = data_dir / tasks.ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["A"],
            "pocket_residue_numbers": [[10]],
            "pocket_residue_indices": [[9]],
        }
    ).to_parquet(lookup, index=False)
    tasks._write_alignment_chain_lookup_manifest(data_dir)


def _write_alignment_mapping_manifest(
    data_dir: Path,
    shard: str,
    *,
    skipped_queries: dict[str, dict[str, object]] | None = None,
) -> None:
    if tasks._completed_alignment_chain_lookup(data_dir) is None:
        _write_alignment_chain_lookup(data_dir)
    inputs = tasks._alignment_input_signatures(data_dir=data_dir, shard=shard)
    outputs = {}
    for alignment_type, signatures in inputs.items():
        output = tasks._alignment_release_path(
            data_dir=data_dir,
            search_db="holo",
            alignment_type=alignment_type,
            shard=shard,
        )
        if not signatures:
            outputs[alignment_type] = None
            continue
        stat = output.stat()
        outputs[alignment_type] = {
            "name": output.name,
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    manifest = tasks._alignment_mapping_manifest_path(data_dir=data_dir, shard=shard)
    manifest.parent.mkdir(exist_ok=True, parents=True)
    manifest.write_text(
        json.dumps(
            {
                "shard": shard,
                "alignment_chain_lookup": tasks._completed_alignment_chain_lookup(
                    data_dir
                ),
                "inputs": inputs,
                "outputs": outputs,
                "skipped_queries": skipped_queries or {},
            }
        )
    )


def test_make_batch_scores_threads_node_scratch_to_score_writer(tmp_path, monkeypatch):
    calls = []

    class FakeScorer:
        shape_score_threads = 1

        def get_score_df(self, *args, **kwargs):
            calls.append((args, kwargs))

    fake_scorer = FakeScorer()
    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: (fake_scorer, ["1abc"], tmp_path / "batch"),
    )
    scratch = tmp_path / "node-scratch"

    tasks.make_batch_scores(
        data_dir=tmp_path,
        pdb_ids=["1abc"],
        scorer_cfg=SimpleNamespace(sub_databases=["holo"]),
        force_update=False,
        scratch_dir=scratch,
        threads=3,
    )

    assert fake_scorer.shape_score_threads == 3
    assert calls == [
        (
            (tmp_path, "1abc"),
            {
                "search_db": "holo",
                "overwrite": False,
                "map_alignments": False,
                "scratch_dir": scratch,
                "source_to_aln_file": {
                    "holo_foldseek": tmp_path
                    / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet",
                    "holo_mmseqs": tmp_path
                    / "alignments/search_db=holo/alignment_type=mmseqs/shard=ab.parquet",
                },
                "defer_ligand_3d": True,
            },
        )
    ]


def test_make_entries_uses_shared_v3_batch(tmp_path, monkeypatch):
    calls = []

    def capture(**kwargs):
        calls.append(kwargs)
        metrics_path = tmp_path / "metrics" / "batch.json"
        metrics_path.parent.mkdir(parents=True)
        metrics_path.write_text(
            json.dumps(
                {
                    "entries": [
                        {"pdb_id": "1abc", "status": "complete"},
                        {"pdb_id": "2def", "status": "failed"},
                    ]
                }
            )
        )
        return metrics_path, True

    monkeypatch.setattr(tasks, "ingest_pdb_batch", capture)
    failed = tasks.make_entries(
        data_dir=tmp_path,
        pdb_ids=["1ABC", "2def"],
        cif_root=tmp_path / "nextgen",
        validation_root=tmp_path / "validation",
        force_update=False,
        annotation_cfg={"min_polymer_size": 10},
        entry_cfg={"plip_complex_threshold": 8},
        cpu=1,
    )

    assert failed == ["2def"]
    assert calls[0]["pdb_ids"] == ["1abc", "2def"]
    assert calls[0]["cif_root"] == tmp_path / "nextgen"
    assert calls[0]["validation_root"] == tmp_path / "validation"
    assert calls[0]["annotation_cfg"] == {"min_polymer_size": 10}
    assert calls[0]["entry_cfg"] == {"plip_complex_threshold": 8}


def test_scatter_make_entries_discovers_configured_source_root(tmp_path):
    cif_root = tmp_path / "external-nextgen"
    cif_path = cif_root / "ab" / "pdb_00001abc" / "pdb_00001abc_xyz-enrich.cif.gz"
    cif_path.parent.mkdir(parents=True)
    cif_path.write_bytes(b"1234")

    chunks = tasks.scatter_make_entries(
        data_dir=tmp_path / "release",
        cif_root=cif_root,
        validation_root=tmp_path / "external-validation",
        batch_size=10,
        two_char_codes=[],
        pdb_ids=["1ABC"],
        force_update=False,
        discovery_threads=1,
    )

    assert chunks == [["1abc"]]


def test_scatter_make_entries_pdb_ids_override_two_char_codes(tmp_path):
    cif_root = tmp_path / "external-nextgen"
    for pdb_id in ("1abc", "2def"):
        cif_path = (
            cif_root
            / pdb_id[1:3]
            / f"pdb_0000{pdb_id}"
            / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
        )
        cif_path.parent.mkdir(parents=True)
        cif_path.touch()

    chunks = tasks.scatter_make_entries(
        data_dir=tmp_path / "release",
        cif_root=cif_root,
        validation_root=tmp_path / "external-validation",
        batch_size=10,
        two_char_codes=["ab"],
        pdb_ids=["2def"],
        force_update=False,
        discovery_threads=1,
    )

    assert chunks == [["2def"]]


def test_entry_collation_tasks_use_shared_core(tmp_path, monkeypatch):
    shard_calls = []
    finalize_calls = []
    monkeypatch.setattr(
        tasks.collate,
        "plan_collation",
        lambda data_dir: {"codes": ["aa", "ab", "ac"]},
    )
    monkeypatch.setattr(
        tasks.collate,
        "collate_shard",
        lambda *args, **kwargs: shard_calls.append((args, kwargs)),
    )
    monkeypatch.setattr(
        tasks.collate,
        "finalize_collation",
        lambda *args, **kwargs: finalize_calls.append((args, kwargs)) or {"ok": True},
    )

    chunks = tasks.scatter_collate_entries(data_dir=tmp_path, batch_size=2)
    tasks.collate_entries(
        data_dir=tmp_path,
        two_char_codes=chunks[0],
        cpu=2,
        memory_limit="7GB",
    )
    result = tasks.finalize_entry_collation(
        data_dir=tmp_path,
        cpu=4,
        memory_limit="32GB",
    )

    assert chunks == [["aa", "ab"], ["ac"]]
    assert [args[1] for args, _ in shard_calls] == ["aa", "ab"]
    assert all(kwargs["threads"] == 2 for _, kwargs in shard_calls)
    assert result == {"ok": True}
    assert finalize_calls[0][1]["threads"] == 4


def test_archive_scatter_pdb_ids_override_two_char_codes(tmp_path):
    for code in ("ab", "de"):
        (tmp_path / "raw_entries" / code).mkdir(parents=True)

    chunks = tasks.scatter_make_canonical_ligand_archives(
        data_dir=tmp_path,
        batch_size=1,
        two_char_codes=["ab"],
        pdb_ids=["2DEF"],
    )

    assert chunks == [["de"]]


@pytest.mark.parametrize(
    "inputs, expected",
    [
        ({"batch_size": 4, "two_char_codes": []}, [4, 4]),
        ({"batch_size": 4, "two_char_codes": []}, [4, 3]),
        ({"two_char_codes": ["xx"], "batch_size": 4}, [1]),
    ],
)
def test_scatter_download_rcsb_files(inputs, expected, tmp_path):
    codes = ["aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh"]
    _orig_rsync_rcsb = io.rsync_rcsb
    _orig_list_rcsb = io.list_rcsb
    io.rsync_rcsb = lambda **kws: codes
    io.list_rcsb = lambda **kws: codes if expected[0] == expected[1] else codes[:-1]
    chunks = tasks.scatter_download_rcsb_files(data_dir=tmp_path, **inputs)
    for chunk, expect in zip(chunks, expected):
        assert len(chunk) == expect
    io.rsync_rcsb = _orig_rsync_rcsb
    io.list_rcsb = _orig_list_rcsb


def test_download_rcsb_files(tmp_path):
    _orig_rsync_rcsb = io.rsync_rcsb
    io.rsync_rcsb = lambda *args, **kws: None
    tasks.download_rcsb_files(data_dir=tmp_path, two_char_codes=["aa"])
    io.rsync_rcsb = _orig_rsync_rcsb


def test_alternative_downloads_use_only_approved_sources(tmp_path, monkeypatch):
    calls = []
    approved = {
        "download_cofactors",
        "download_seqres_data",
        "refresh_bundled_ccd",
        "download_affinity_data",
    }

    for name in approved:
        monkeypatch.setattr(io, name, lambda *, _name=name, **_: calls.append(_name))

    tasks.download_alternative_datasets(
        data_dir=tmp_path,
        threads=1,
        force_update=False,
    )

    assert set(calls) == approved
    assert not hasattr(io, "download_ecod_data")
    assert not hasattr(io, "download_panther_data")
    assert not hasattr(io, "download_kinase_data")


def test_scoring_finalization_stage_order_and_partitions():
    assert tasks.STAGES.index("collate_entries") < tasks.STAGES.index("make_dbs")
    assert tasks.STAGES.index("make_canonical_ligand_archives") < tasks.STAGES.index(
        "finalize_ligand_archives"
    )
    assert tasks.STAGES.index("finalize_ligand_archives") < tasks.STAGES.index(
        "compute_ligand_fingerprints"
    )
    assert tasks.STAGES.index("make_ligand_scores") < tasks.STAGES.index(
        "annotate_ligand_similarity"
    )
    assert tasks.STAGES.index("annotate_ligand_similarity") < tasks.STAGES.index(
        "make_sub_dbs"
    )
    assert tasks.STAGES.index("map_batch_alignments") < tasks.STAGES.index(
        "collate_alignments"
    )
    assert tasks.STAGES.index("collate_alignments") < tasks.STAGES.index(
        "finalize_alignments"
    )
    assert tasks.STAGES.index("finalize_alignments") < tasks.STAGES.index(
        "make_batch_scores"
    )
    assert tasks.STAGES.index("make_batch_scores") < tasks.STAGES.index(
        "collate_ligand_3d_candidates"
    )
    assert tasks.STAGES.index("collate_ligand_3d_candidates") < tasks.STAGES.index(
        "plan_ligand_3d_scores"
    )
    assert tasks.STAGES.index("plan_ligand_3d_scores") < tasks.STAGES.index(
        "make_ligand_3d_scores"
    )
    assert tasks.STAGES.index("make_ligand_3d_scores") < tasks.STAGES.index(
        "collate_ligand_3d_scores"
    )
    assert tasks.STAGES.index("collate_ligand_3d_scores") < tasks.STAGES.index(
        "merge_ligand_3d_scores"
    )
    assert tasks.STAGES.index("merge_ligand_3d_scores") < tasks.STAGES.index(
        "finalize_scores"
    )
    assert tasks.STAGES.index("finalize_scores") < tasks.STAGES.index(
        "export_sucos_shape_pocket_qcov"
    )
    assert tasks.STAGES.index("export_sucos_shape_pocket_qcov") < tasks.STAGES.index(
        "finalize_sucos_export"
    )
    assert tasks.STAGES.index("finalize_sucos_export") < tasks.STAGES.index(
        "collate_partitions"
    )
    assert tasks.STAGES.index("collate_partitions") < tasks.STAGES.index(
        "make_component_reductions"
    )
    assert tasks.STAGES.index("make_component_reductions") < tasks.STAGES.index(
        "merge_component_reductions"
    )
    assert tasks.STAGES.index("merge_component_reductions") < tasks.STAGES.index(
        "make_communities"
    )
    assert tasks.STAGES.index("make_communities") < tasks.STAGES.index(
        "make_directed_set_covers"
    )
    assert tasks.STAGES.index("make_directed_set_covers") < tasks.STAGES.index(
        "summarize_clusters"
    )
    assert tasks.STAGES.index("summarize_clusters") < tasks.STAGES.index(
        "finalize_index"
    )
    assert tasks.STAGES.index("finalize_index") < tasks.STAGES.index("make_mmp_index")
    partitions = tasks.scatter_collate_partitions()
    assert len(partitions) == 38
    assert ["0"] in partitions
    assert ["z"] in partitions
    assert ["apo"] in partitions
    assert ["pred"] in partitions


def test_directed_set_cover_scatter_skips_only_complete_outputs(tmp_path):
    work = tasks.scatter_make_directed_set_covers(
        data_dir=tmp_path,
        metrics=["pocket_qcov"],
        thresholds=[100, 50],
        stop_on_cluster=0,
        skip_existing=True,
    )
    assert work == [[("pocket_qcov", 100)], [("pocket_qcov", 50)]]

    output = (
        tmp_path
        / "ligand_sampling/directed_set_cover/metric=pocket_qcov"
        / "threshold=100.parquet"
    )
    output.parent.mkdir(parents=True)
    output.touch()
    work = tasks.scatter_make_directed_set_covers(
        data_dir=tmp_path,
        metrics=["pocket_qcov"],
        thresholds=[100, 50],
        stop_on_cluster=0,
        skip_existing=True,
    )
    assert work == [[("pocket_qcov", 50)]]

    output.with_name("threshold=50.parquet").touch()
    assert tasks.scatter_make_directed_set_covers(
        data_dir=tmp_path,
        metrics=["pocket_qcov"],
        thresholds=[100, 50],
        stop_on_cluster=0,
        skip_existing=True,
    ) == [[]]


def test_ligand_score_threshold_must_cover_frequency_clustering() -> None:
    with pytest.raises(ValueError, match="must not exceed"):
        LigandConfig(minimum_similarity=95)


def test_scatter_protein_scoring_uses_v3_chain_index(tmp_path) -> None:
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi", "4jkl"],
            "chain_receptor_type": [
                "protein",
                "dna",
                "protein",
                "protein",
            ],
            "chain_is_holo": [True, True, True, False],
        }
    ).to_parquet(index_dir / "entry_chains.parquet", index=False)

    assert tasks.scatter_protein_scoring(
        data_dir=tmp_path,
        batch_size=1,
        two_char_codes=[],
        pdb_ids=[],
    ) == [["1abc"], ["3ghi"]]
    assert tasks.scatter_protein_scoring(
        data_dir=tmp_path,
        batch_size=10,
        two_char_codes=["ab"],
        pdb_ids=["3GHI"],
    ) == [["3ghi"]]


def test_protein_scoring_plan_and_alignment_finalization(tmp_path, monkeypatch) -> None:
    from plinder.data.pipeline.score import (
        finalize_alignment_artifacts,
        make_foldseek_input_manifest,
        plan_protein_scoring,
    )

    index_dir = tmp_path / "index"
    index_dir.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "2def"],
            "chain_receptor_type": [
                "protein",
                "protein",
                "dna",
            ],
            "chain_is_holo": [True, False, True],
        }
    ).to_parquet(index_dir / "entry_chains.parquet", index=False)
    plan = plan_protein_scoring(tmp_path)

    assert plan["query_count"] == 1
    assert plan["protein_chain_count"] == 1
    assert plan["max_seqs"] == 10_000
    foldseek_inputs = make_foldseek_input_manifest(tmp_path, Path("/nextgen"))
    assert foldseek_inputs.read_text().splitlines() == [
        "/nextgen/ab/pdb_00001abc/pdb_00001abc_xyz-enrich.cif.gz"
    ]
    for alignment_type in ["foldseek", "mmseqs"]:
        backend = tmp_path / "dbs/subdbs" / f"holo_{alignment_type}"
        backend.mkdir(parents=True)
        source_index = backend / f"holo_{alignment_type}.index"
        source_index.write_text("0\t0\t1\n")
        source_identifier = (
            "pdb_00001abc_xyz-enrich_A" if alignment_type == "foldseek" else "1abc_A"
        )
        (backend / f"holo_{alignment_type}.lookup").write_text(
            f"0\t{source_identifier}\t0\n"
        )
        search_target = (
            "clustered" if alignment_type == "foldseek" else "representatives"
        )
        conversion_target = (
            "clustered" if alignment_type == "foldseek" else f"holo_{alignment_type}"
        )
        for database in {search_target, conversion_target}:
            (backend / f"{database}.dbtype").touch()
        (backend / "exact_cluster.json").write_text(
            json.dumps(
                {
                    "alignment_type": alignment_type,
                    "identity": 1.0,
                    "coverage": 1.0,
                    "coverage_mode": 0,
                    "source_index": {
                        "name": source_index.name,
                        "size": source_index.stat().st_size,
                        "mtime_ns": source_index.stat().st_mtime_ns,
                    },
                    "search_target": search_target,
                    "conversion_target": conversion_target,
                }
            )
        )
        output = backend / "aln"
        output.mkdir()
        (output / "1abc.parquet").touch()
        (output / "4qp3.tmp.parquet").touch()
        release = (
            tmp_path / "alignments/search_db=holo" / f"alignment_type={alignment_type}"
        )
        release.mkdir(parents=True)
        columns = {
            "query_entry": ["1abc"],
            "target_entry": ["1abc"],
            "query_chain_mapped": ["A"],
            "target_chain_mapped": ["A"],
            "source": [alignment_type],
            "query_pocket_residue_numbers": [[1]],
            "target_pocket_residue_numbers": [[1]],
            "pocket_residue_identity": [bytes([1])],
            "qcov": [1.0],
            "fident": [1.0],
            "seqsim": [1.0],
        }
        if alignment_type == "foldseek":
            columns["lddt"] = [1.0]
        pd.DataFrame(columns).to_parquet(release / "shard=ab.parquet", index=False)
    _write_alignment_mapping_manifest(tmp_path, "ab")
    with monkeypatch.context() as patch:
        patch.setattr(
            tasks,
            "_alignment_input_signatures",
            lambda **_kwargs: pytest.fail("finalization rescanned raw alignments"),
        )
        patch.setattr(
            tasks,
            "alignment_mapping_shard_is_current",
            lambda **_kwargs: pytest.fail("finalization used the per-shard validator"),
        )
        report = finalize_alignment_artifacts(tmp_path)
    assert report["status"] == "complete"
    assert report["artifact_counts"]["foldseek_release_shards"] == 1
    assert report["target_clustering"]["expand_to_chain_level"] is True
    assert report["skipped_queries"] == {}
    assert (tmp_path / "alignments/manifest.json").is_file()

    skipped = {
        "1abc": {
            "reason": "raw_alignment_row_budget_exceeded",
            "maximum_rows": 1,
            "total_rows": 2,
            "rows_by_alignment_type": {"foldseek": 1, "mmseqs": 1},
        }
    }
    for alignment_type in ["foldseek", "mmseqs"]:
        release = tasks._alignment_release_path(
            data_dir=tmp_path,
            search_db="holo",
            alignment_type=alignment_type,
            shard="ab",
        )
        pq.write_table(
            pa.Table.from_pylist(
                [],
                schema=schemas.mapped_alignment_schema(alignment_type=alignment_type),
            ),
            release,
        )
    _write_alignment_mapping_manifest(tmp_path, "ab", skipped_queries=skipped)

    with monkeypatch.context() as patch:
        patch.setattr(
            tasks,
            "_alignment_input_signatures",
            lambda **_kwargs: pytest.fail("finalization rescanned raw alignments"),
        )
        skipped_report = finalize_alignment_artifacts(tmp_path)
    assert skipped_report["skipped_queries"] == skipped
    assert skipped_report["artifact_counts"]["foldseek_skipped_queries"] == 1
    assert skipped_report["artifact_counts"]["mmseqs_skipped_queries"] == 1


def test_protein_scoring_plan_groups_queries_by_two_character_shard(tmp_path):
    from plinder.data.pipeline.score import (
        MANIFEST_RELATIVE,
        plan_protein_scoring,
    )

    index_dir = tmp_path / "index"
    index_dir.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1zzz", "2aaa", "3aab"],
            "chain_receptor_type": ["protein"] * 3,
            "chain_is_holo": [True] * 3,
        }
    ).to_parquet(index_dir / "entry_chains.parquet", index=False)

    plan_protein_scoring(tmp_path)

    manifest = pd.read_parquet(tmp_path / MANIFEST_RELATIVE)
    assert manifest["pdb_id"].tolist() == ["2aaa", "3aab", "1zzz"]
    assert manifest["shard"].tolist() == ["aa", "aa", "zz"]


def test_score_work_plan_balances_expensive_queries_into_fixed_batches(tmp_path):
    from plinder.data.pipeline.score import (
        SCORE_WORK_RELATIVE,
        _score_batch,
        plan_protein_scoring,
        plan_score_batches,
    )

    pdb_ids = ["1aaa", "2aab", "3aac", "4aad", "5aae"]
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": pdb_ids,
            "chain_asym_id": ["A"] * len(pdb_ids),
            "chain_receptor_type": ["protein"] * len(pdb_ids),
            "chain_is_holo": [True] * len(pdb_ids),
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    annotation_rows = [
        {
            "entry_pdb_id": pdb_id,
            "system_id": f"{pdb_id}__1",
            "ligand_id": f"{pdb_id}__1__1.{ligand_asym_id}",
            "ligand_asym_id": ligand_asym_id,
            "ligand_is_proper": True,
            "ligand_is_3d_score_able": True,
            "system_type": "holo",
            "system_protein_chains_asym_id": ["1.A"],
            "system_num_protein_chains": 1,
            "system_num_ligand_chains": 6 if pdb_id == "5aae" else 1,
        }
        for pdb_id, ligand_asym_id in zip(pdb_ids, ["L", "M", "N", "O", "P"])
    ]
    annotation_rows.extend(
        {
            **annotation_rows[-1],
            "ligand_id": f"5aae__1__1.{ligand_asym_id}",
            "ligand_asym_id": ligand_asym_id,
        }
        for ligand_asym_id in ["Q", "R", "S", "T", "U"]
    )
    pd.DataFrame(annotation_rows).to_parquet(
        index / "annotation_table.parquet", index=False
    )
    plan_protein_scoring(tmp_path)
    release = (
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=aa.parquet"
    )
    release.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1aaa"] * 10 + ["2aab"] * 8 + ["3aac"] + ["4aad"],
            "target_entry": ["4aad"] * 20,
        }
    ).to_parquet(release, index=False)

    report = plan_score_batches(
        tmp_path,
        batch_size=2,
        threads=1,
        scratch_dir=tmp_path / "scratch",
        max_query_protein_chains=5,
        max_query_proper_ligand_chains=5,
    )

    assert report["batch_count"] == 2
    work = pd.read_parquet(tmp_path / SCORE_WORK_RELATIVE)
    assert set(work["pdb_id"]) == set(pdb_ids[:4])
    assert work.groupby("score_batch_index").size().tolist() == [2, 2]
    assert work.groupby("score_batch_index")["shard"].nunique().max() == 1
    batches = work.set_index("pdb_id")["score_batch_index"]
    assert batches["1aaa"] != batches["2aab"]
    assert sorted(
        _score_batch(tmp_path, 0, 2) + _score_batch(tmp_path, 1, 2)
    ) == sorted(pdb_ids[:4])


def test_score_work_plan_ignores_nonproper_ligands_in_query_cap(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        SCORE_WORK_RELATIVE,
        plan_protein_scoring,
        plan_score_batches,
    )

    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"] * 6,
            "chain_asym_id": list("ABCDEF"),
            "chain_receptor_type": ["protein", "dna", "dna", "rna", "rna", "rna"],
            "chain_is_holo": [True] * 6,
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"] * 6,
            "system_id": ["1abc__1"] * 6,
            "ligand_id": [f"1abc__1__1.{value}" for value in "LMNOPQ"],
            "ligand_asym_id": list("LMNOPQ"),
            "ligand_is_proper": [True, False, False, False, False, False],
            "ligand_is_3d_score_able": [True] * 6,
            "system_type": ["holo"] * 6,
            "system_protein_chains_asym_id": [[f"1.{chain}" for chain in "ABCDEF"]] * 6,
            "system_num_protein_chains": [6] * 6,
            "system_num_ligand_chains": [6] * 6,
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)
    plan_protein_scoring(tmp_path)
    release = (
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet"
    )
    release.parent.mkdir(parents=True)
    pd.DataFrame({"query_entry": ["1abc"], "target_entry": ["1abc"]}).to_parquet(
        release, index=False
    )

    plan_score_batches(
        tmp_path,
        batch_size=1,
        threads=1,
        scratch_dir=tmp_path / "scratch",
        max_query_protein_chains=5,
        max_query_proper_ligand_chains=5,
    )

    assert pd.read_parquet(tmp_path / SCORE_WORK_RELATIVE)["pdb_id"].tolist() == [
        "1abc"
    ]


def test_score_work_plan_keeps_over_cap_holo_systems_as_targets(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        SCORE_WORK_RELATIVE,
        plan_protein_scoring,
        plan_score_batches,
    )

    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "chain_asym_id": ["A", "A"],
            "chain_receptor_type": ["protein", "protein"],
            "chain_is_holo": [True, True],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    rows = [
        {
            "entry_pdb_id": "1abc",
            "system_id": "1abc__1",
            "ligand_id": "1abc__1__1.L",
            "ligand_asym_id": "L",
            "ligand_is_proper": True,
            "ligand_is_3d_score_able": True,
            "system_type": "holo",
            "system_protein_chains_asym_id": ["1.A"],
            "system_num_protein_chains": 1,
            "system_num_ligand_chains": 1,
        }
    ]
    rows.extend(
        {
            **rows[0],
            "entry_pdb_id": "2def",
            "system_id": "2def__1",
            "ligand_id": f"2def__1__1.{asym_id}",
            "ligand_asym_id": asym_id,
            "system_num_ligand_chains": 6,
        }
        for asym_id in "LMNOPQ"
    )
    pd.DataFrame(rows).to_parquet(index / "annotation_table.parquet", index=False)
    plan_protein_scoring(tmp_path)
    release = (
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet"
    )
    release.parent.mkdir(parents=True)
    pd.DataFrame({"query_entry": ["1abc"], "target_entry": ["2def"]}).to_parquet(
        release, index=False
    )

    plan_score_batches(
        tmp_path,
        batch_size=1,
        threads=1,
        scratch_dir=tmp_path / "scratch",
        max_query_protein_chains=5,
        max_query_proper_ligand_chains=5,
    )

    work = pd.read_parquet(tmp_path / SCORE_WORK_RELATIVE).set_index("pdb_id")
    assert work.index.tolist() == ["1abc"]
    assert work.loc["1abc", "ligand_pair_work"] == 6
    assert work.loc["1abc", "canonical_3d_work"] == 6


def test_score_manifest_batch_supports_single_pdb_retries(tmp_path) -> None:
    from plinder.data.pipeline.score import _score_manifest_batch

    manifest = tmp_path / "retry.parquet"
    pd.DataFrame({"pdb_id": ["1abc", "2def", "3ghi"]}).to_parquet(manifest, index=False)

    assert _score_manifest_batch(tmp_path, manifest, 0, 1) == ["1abc"]
    assert _score_manifest_batch(tmp_path, manifest, 1, 1) == ["2def"]
    assert _score_manifest_batch(tmp_path, manifest, 1, 2) == ["3ghi"]

    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2def", "3ghi"],
            "retry_batch_index": [1, 0, 1],
        }
    ).to_parquet(manifest, index=False)
    assert _score_manifest_batch(tmp_path, manifest, 0, 2) == ["2def"]
    assert _score_manifest_batch(tmp_path, manifest, 1, 2) == ["1abc", "3ghi"]


def test_repair_batch_assignment_balances_a_remainder() -> None:
    from plinder.data.pipeline.score import _assign_repair_batches

    records = [
        {"pdb_id": f"query-{index:03}", "estimated_work": 1} for index in range(101)
    ]

    batch_count = _assign_repair_batches(
        records,
        batch_size=100,
        index_offset=7,
    )

    counts = pd.Series(
        [record["repair_batch_index"] for record in records]
    ).value_counts()
    assert batch_count == 2
    assert set(counts.index) == {7, 8}
    assert sorted(counts.tolist()) == [50, 51]


def test_score_repair_plans_full_and_target_only_queries(
    tmp_path: Path,
) -> None:
    from plinder.data.pipeline.score import (
        MANIFEST_RELATIVE,
        PLAN_RELATIVE,
        SCORE_WORK_RELATIVE,
        _score_repair_batch,
        _source_signature,
        plan_score_repair,
    )

    query_manifest = tmp_path / MANIFEST_RELATIVE
    query_manifest.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": ["1abc", "2def", "3ghi", "4jkl"]}).to_parquet(
        query_manifest, index=False
    )
    (tmp_path / PLAN_RELATIVE).write_text(
        json.dumps(
            {
                "manifest": _source_signature(query_manifest),
                "score_max_query_protein_chains": 5,
                "score_max_query_proper_ligand_chains": 5,
            }
        )
    )
    work = tmp_path / SCORE_WORK_RELATIVE
    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2def", "3ghi"],
            "estimated_work": [10, 20, 30],
        }
    ).to_parquet(work, index=False)
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "system_id": ["1abc__1", "2def__1", "3ghi__1"],
            "ligand_id": ["1abc__1__1.L", "2def__1__1.L", "3ghi__1__1.L"],
            "system_type": ["holo", "holo", "holo"],
            "system_protein_chains_asym_id": [["1.A"], ["1.A"], ["1.A"]],
            "ligand_is_proper": [True, True, True],
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "chain_asym_id": ["A", "A", "A"],
            "chain_receptor_type": ["protein", "protein", "protein"],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    alignment = (
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet"
    )
    alignment.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc", "1abc", "1abc", "3ghi", "9zzz"],
            "target_entry": ["2def", "2def", "4jkl", "8nop", "2def"],
        }
    ).to_parquet(alignment, index=False)
    inactive_score = tmp_path / "dbs/subdbs/search_db=holo/4jkl.parquet"
    inactive_score.parent.mkdir(parents=True)
    inactive_score.touch()
    affected = tmp_path / "affected.txt"
    affected.write_text("2def\n4jkl\n")
    additional_full = tmp_path / "additional_full.txt"
    additional_full.write_text("3ghi\n9zzz\n")

    report = plan_score_repair(
        tmp_path,
        affected_manifest=affected,
        additional_full_query_manifest=additional_full,
        batch_size=2,
    )

    assert report["full_query_count"] == 2
    assert report["additional_full_query_count"] == 1
    assert report["ignored_additional_full_query_count"] == 1
    assert report["dropped_query_count"] == 1
    assert report["target_only_query_count"] == 1
    assert report["batch_count"] == 2
    assert report["full_query_batch_count"] == 1
    assert report["target_query_batch_count"] == 1
    repairs = [
        *_score_repair_batch(
            tmp_path / "manifests/score_repair.parquet",
            batch_index=0,
            batch_size=2,
        ),
        *_score_repair_batch(
            tmp_path / "manifests/score_repair.parquet",
            batch_index=1,
            batch_size=2,
        ),
    ]
    repair_by_query = {str(repair["pdb_id"]): repair for repair in repairs}
    assert repair_by_query["2def"]["repair_mode"] == "full"
    assert repair_by_query["3ghi"]["repair_mode"] == "full"
    assert repair_by_query["4jkl"]["repair_mode"] == "drop"
    assert list(repair_by_query["1abc"]["target_pdb_ids"]) == ["2def", "4jkl"]


def test_repair_batch_scores_uses_full_and_target_only_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[str, str, object]] = []

    class FakeScorer:
        shape_score_threads = 0

        def get_score_df(self, _data_dir, pdb_id, **kwargs):
            calls.append(("full", pdb_id, kwargs["defer_ligand_3d"]))

        def repair_score_df_targets(
            self, _data_dir, pdb_id, *, affected_target_entries, **_kwargs
        ):
            calls.append(("targets", pdb_id, affected_target_entries))

    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: (FakeScorer(), ["1abc", "2def"], tmp_path / "db"),
    )

    tasks.repair_batch_scores(
        data_dir=tmp_path,
        repairs=[
            {
                "pdb_id": "1abc",
                "repair_mode": "targets",
                "target_pdb_ids": ["2def"],
            },
            {"pdb_id": "2def", "repair_mode": "full", "target_pdb_ids": []},
        ],
        scorer_cfg=SimpleNamespace(sub_databases=["holo"]),
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )

    assert calls == [
        ("targets", "1abc", {"2def"}),
        ("full", "2def", True),
    ]


def test_repair_batch_scores_drops_queries_that_are_no_longer_eligible(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    score = tmp_path / "dbs/subdbs/search_db=holo/1abc.parquet"
    candidate = (
        tmp_path / "scores/ligand_3d_candidates/search_db=holo/shard=ab/1abc.parquet"
    )
    score.parent.mkdir(parents=True)
    candidate.parent.mkdir(parents=True)
    score.touch()
    candidate.touch()

    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: pytest.fail("drop-only repair must not create a scorer"),
    )

    tasks.repair_batch_scores(
        data_dir=tmp_path,
        repairs=[{"pdb_id": "1abc", "repair_mode": "drop", "target_pdb_ids": []}],
        scorer_cfg=SimpleNamespace(sub_databases=["holo"]),
        scratch_dir=tmp_path / "scratch",
    )

    assert not score.exists()
    assert not candidate.exists()


def test_repair_target_remainders_drops_stalled_query_and_resumes_followers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from plinder.data.pipeline import score as score_pipeline

    repair_manifest = tmp_path / "repair.parquet"
    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2def", "3ghi"],
            "repair_mode": ["targets", "targets", "targets"],
            "target_pdb_ids": [["9xyz"], ["9xyz"], ["9xyz"]],
            "repair_batch_index": [0, 0, 0],
        }
    ).to_parquet(repair_manifest, index=False)
    score, candidate = score_pipeline._score_repair_query_paths(tmp_path, "1abc")
    score.parent.mkdir(parents=True)
    candidate.parent.mkdir(parents=True)
    score.touch()
    candidate.touch()
    completed_ns = repair_manifest.stat().st_mtime_ns + 1_000_000_000
    os.utime(score, ns=(completed_ns, completed_ns))
    os.utime(candidate, ns=(completed_ns, completed_ns))
    captured: list[dict[str, object]] = []

    def fake_repair_batch_scores(**kwargs) -> None:
        captured.extend(kwargs["repairs"])

    monkeypatch.setattr(tasks, "repair_batch_scores", fake_repair_batch_scores)

    report = score_pipeline.repair_target_score_remainders(
        tmp_path,
        repair_manifest=repair_manifest,
        batch_index=0,
        batch_size=3,
        scorer_cfg=SimpleNamespace(),
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    assert report["dropped_query"] == "2def"
    assert report["resumed_query_count"] == 1
    assert [repair["pdb_id"] for repair in captured] == ["3ghi"]
    marker_dir = score_pipeline._score_repair_marker_dir(tmp_path, repair_manifest)
    marker = json.loads((marker_dir / "2def.json").read_text())
    assert marker["pdb_id"] == "2def"
    (marker_dir / "concurrent.tmp.json").write_text("{")

    second_report = score_pipeline.repair_target_score_remainders(
        tmp_path,
        repair_manifest=repair_manifest,
        batch_index=0,
        batch_size=3,
        scorer_cfg=SimpleNamespace(),
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    assert second_report["dropped_query"] == "3ghi"
    assert second_report["resumed_query_count"] == 0
    assert score_pipeline._score_repair_marked_targets(tmp_path, repair_manifest) == {
        "2def",
        "3ghi",
    }


def test_repair_target_remainders_removes_planned_drops_without_marking_them(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from plinder.data.pipeline import score as score_pipeline

    repair_manifest = tmp_path / "repair.parquet"
    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2def"],
            "repair_mode": ["drop", "targets"],
            "target_pdb_ids": [[], ["9xyz"]],
            "repair_batch_index": [0, 0],
        }
    ).to_parquet(repair_manifest, index=False)
    for pdb_id in ["1abc", "2def"]:
        score, candidate = score_pipeline._score_repair_query_paths(tmp_path, pdb_id)
        score.parent.mkdir(parents=True, exist_ok=True)
        candidate.parent.mkdir(parents=True, exist_ok=True)
        score.touch()
        candidate.touch()
    completed_ns = repair_manifest.stat().st_mtime_ns + 1_000_000_000
    for path in score_pipeline._score_repair_query_paths(tmp_path, "2def"):
        os.utime(path, ns=(completed_ns, completed_ns))

    monkeypatch.setattr(
        tasks,
        "repair_batch_scores",
        lambda **_kwargs: pytest.fail("no target repair should be needed"),
    )

    report = score_pipeline.repair_target_score_remainders(
        tmp_path,
        repair_manifest=repair_manifest,
        batch_index=0,
        batch_size=2,
        scorer_cfg=SimpleNamespace(),
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    assert report["dropped_query"] is None
    assert report["resumed_query_count"] == 0
    assert not any(
        path.exists()
        for path in score_pipeline._score_repair_query_paths(tmp_path, "1abc")
    )
    assert not score_pipeline._score_repair_marker_dir(
        tmp_path, repair_manifest
    ).exists()


def test_finalize_score_repair_records_marked_targets_and_incomplete_full_queries(
    tmp_path: Path,
) -> None:
    from plinder.data.pipeline import score as score_pipeline

    repair_manifest = tmp_path / "repair.parquet"
    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2def", "3ghi"],
            "repair_mode": ["targets", "targets", "full"],
        }
    ).to_parquet(repair_manifest, index=False)
    current_score, current_candidate = score_pipeline._score_repair_query_paths(
        tmp_path, "1abc"
    )
    current_score.parent.mkdir(parents=True)
    current_candidate.parent.mkdir(parents=True)
    current_score.touch()
    current_candidate.touch()
    completed_ns = repair_manifest.stat().st_mtime_ns + 1_000_000_000
    os.utime(current_score, ns=(completed_ns, completed_ns))
    os.utime(current_candidate, ns=(completed_ns, completed_ns))
    marker_dir = score_pipeline._score_repair_marker_dir(tmp_path, repair_manifest)
    marker = marker_dir / "2def.json"
    marker.parent.mkdir(parents=True)
    marker.write_text(
        json.dumps(
            {
                "pdb_id": "2def",
                "repair_run_id": score_pipeline._score_repair_run_id(repair_manifest),
                "repair_batch_index": 4,
                "repair_mode": "targets",
                "reason": "resource_limit",
            }
        )
    )
    existing = tmp_path / score_pipeline.DROPPED_QUERY_RELATIVE
    existing.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "pdb_id": ["4jkl"],
            "stage": ["derived_scoring"],
            "reason": ["resource_limit"],
            "details": ["{}"],
        }
    ).to_parquet(existing, index=False)
    plan = tmp_path / score_pipeline.PLAN_RELATIVE
    plan.write_text("{}")

    report = score_pipeline.finalize_score_repair_queries(
        tmp_path,
        repair_manifest=repair_manifest,
        max_new_drops=2,
    )

    assert report["new_target_query_drops"] == 1
    assert report["new_full_query_drops"] == 1
    dropped = pd.read_parquet(tmp_path / score_pipeline.DROPPED_QUERY_RELATIVE)
    assert set(dropped["pdb_id"].astype(str)) == {"2def", "3ghi", "4jkl"}


def test_score_repair_scores_only_new_canonical_ligand_pairs(
    tmp_path: Path,
) -> None:
    from plinder.data.pipeline.score import (
        SCORE_REPAIR_LIGAND_3D_RELATIVE,
        merge_score_repair_ligand_3d,
        plan_score_repair_ligand_3d,
    )

    repair_manifest = tmp_path / "repair.parquet"
    pd.DataFrame({"pdb_id": ["1abc"], "shard": ["ab"]}).to_parquet(
        repair_manifest, index=False
    )
    candidate = tmp_path / "scores/ligand_3d_pair_candidate_shards/shard=ab.parquet"
    candidate.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc", "1abc"],
            "query_ligand_asym_id": ["B", "B"],
            "target_entry": ["2def", "3ghi"],
            "target_ligand_asym_id": ["Y", "Z"],
            "pocket_qcov": [0.5, 0.7],
        }
    ).to_parquet(candidate, index=False)
    cached = tmp_path / "scores/ligand_3d_by_query/ab.parquet"
    cached.parent.mkdir(parents=True)
    cached_row = {
        "query_entry": "1abc",
        "query_ligand_asym_id": "B",
        "target_entry": "2def",
        "target_ligand_asym_id": "Y",
        "shape": 0.5,
        "color": 0.4,
        "sucos_shape": 0.3,
    }
    pd.DataFrame([cached_row]).to_parquet(
        cached, index=False, schema=schemas.LIGAND_3D_SCORE_SCHEMA
    )
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "ligand_asym_id": ["B", "Y", "Z"],
            "ligand_num_heavy_atoms": [10, 20, 30],
            "ligand_is_proper": [True, True, True],
            "ligand_is_3d_score_able": [True, True, True],
            "system_type": ["holo", "holo", "holo"],
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)

    report = plan_score_repair_ligand_3d(
        tmp_path,
        repair_manifest=repair_manifest,
        batch_size=10,
        threads=1,
        scratch_dir=tmp_path / "plan-scratch",
        memory_limit="1GB",
    )

    assert report["pair_count"] == 1
    work = pd.read_parquet(tmp_path / SCORE_REPAIR_LIGAND_3D_RELATIVE)
    assert work[["query_entry", "target_entry"]].to_dict("records") == [
        {"query_entry": "1abc", "target_entry": "3ghi"}
    ]
    repaired = work[
        [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    ].copy()
    repaired["shape"] = 0.8
    repaired["color"] = 0.7
    repaired["sucos_shape"] = 0.6
    repair_score = tmp_path / "scores/ligand_3d_pair_repairs/0.parquet"
    repair_score.parent.mkdir()
    repaired.to_parquet(
        repair_score, index=False, schema=schemas.LIGAND_3D_SCORE_SCHEMA
    )

    merge_score_repair_ligand_3d(
        tmp_path,
        shards=["ab"],
        scratch_dir=tmp_path / "merge-scratch",
        threads=1,
    )

    merged = pd.read_parquet(cached)
    assert set(merged["target_entry"]) == {"2def", "3ghi"}


def test_dropped_queries_combine_mapping_and_scoring_stages(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        DROPPED_QUERY_RELATIVE,
        MANIFEST_RELATIVE,
        PLAN_RELATIVE,
        SCORE_WORK_RELATIVE,
        _score_batch,
        _source_signature,
        active_scoring_query_ids,
        record_dropped_queries,
    )

    query_ids = ["1abc", "2def", "3ghi", "4jkl"]
    query_manifest = tmp_path / MANIFEST_RELATIVE
    query_manifest.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": query_ids}).to_parquet(query_manifest, index=False)
    score_work = tmp_path / SCORE_WORK_RELATIVE
    pd.DataFrame(
        {
            "pdb_id": query_ids,
            "score_batch_index": [0, 0, 0, 0],
            "estimated_work": [4.0, 3.0, 2.0, 1.0],
        }
    ).to_parquet(score_work, index=False)
    (tmp_path / PLAN_RELATIVE).write_text(
        json.dumps(
            {
                "manifest": _source_signature(query_manifest),
                "score_work": _source_signature(score_work),
                "score_batch_size": 4,
            }
        )
    )
    alignment_manifest = tmp_path / "alignments" / "manifest.json"
    alignment_manifest.parent.mkdir()
    alignment_manifest.write_text(
        json.dumps(
            {
                "skipped_queries": {
                    "1abc": {
                        "reason": "raw_alignment_row_budget_exceeded",
                        "total_rows": 6_000_000,
                    }
                }
            }
        )
    )
    score_drops = tmp_path / "score-drops.parquet"
    pd.DataFrame({"pdb_id": ["2def", "3ghi"]}).to_parquet(score_drops, index=False)

    report = record_dropped_queries(
        tmp_path,
        score_query_manifest=score_drops,
    )

    assert report["alignment_mapping"] == 1
    assert report["derived_scoring"] == 2
    assert report["total"] == 3
    dropped = pd.read_parquet(tmp_path / DROPPED_QUERY_RELATIVE)
    assert dropped[["pdb_id", "stage"]].to_dict("records") == [
        {"pdb_id": "1abc", "stage": "alignment_mapping"},
        {"pdb_id": "2def", "stage": "derived_scoring"},
        {"pdb_id": "3ghi", "stage": "derived_scoring"},
    ]
    assert active_scoring_query_ids(tmp_path) == {"4jkl"}
    assert _score_batch(tmp_path, 0, 4) == ["4jkl"]


def test_ligand_3d_retries_require_and_reassemble_every_retry(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        LIGAND_3D_RETRY_WORK_RELATIVE,
        LIGAND_3D_WORK_RELATIVE,
        MANIFEST_RELATIVE,
        PLAN_RELATIVE,
        _source_signature,
        finalize_ligand_3d_retries,
        plan_ligand_3d_retries,
    )

    query_manifest = tmp_path / MANIFEST_RELATIVE
    query_manifest.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": ["1abc"]}).to_parquet(query_manifest, index=False)
    pair_rows = pd.DataFrame(
        {
            "query_entry": ["1abc", "1abc", "1abc"],
            "query_ligand_asym_id": ["A", "A", "A"],
            "target_entry": ["2def", "3ghi", "4jkl"],
            "target_ligand_asym_id": ["B", "C", "D"],
            "estimated_work": [3, 2, 1],
            "ligand_3d_batch_index": [0, 1, 1],
        }
    )
    work_path = tmp_path / LIGAND_3D_WORK_RELATIVE
    pair_rows.to_parquet(work_path, index=False)
    (tmp_path / PLAN_RELATIVE).write_text(
        json.dumps(
            {
                "manifest": _source_signature(query_manifest),
                "ligand_3d_plan_complete": True,
                "ligand_3d_batch_count": 2,
                "ligand_3d_work": _source_signature(work_path),
            }
        )
    )
    canonical = tmp_path / "scores" / "ligand_3d_pairs"
    canonical.mkdir(parents=True)
    completed = pair_rows.iloc[[0]][
        [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    ].copy()
    completed["shape"] = 0.5
    completed["color"] = 0.4
    completed["sucos_shape"] = 0.3
    completed.to_parquet(
        canonical / "0.parquet",
        schema=schemas.LIGAND_3D_SCORE_SCHEMA,
        index=False,
    )

    report = plan_ligand_3d_retries(tmp_path, batch_size=1)
    assert report["missing_original_batch_count"] == 1
    assert report["retry_batch_count"] == 2
    retry_work = pd.read_parquet(tmp_path / LIGAND_3D_RETRY_WORK_RELATIVE)
    first_retry = retry_work[retry_work["retry_batch_index"].eq(0)]
    scored = first_retry[
        [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    ].copy()
    scored["shape"] = 0.8
    scored["color"] = 0.7
    scored["sucos_shape"] = 0.6
    retry_dir = tmp_path / "scores" / "ligand_3d_pair_retries"
    retry_dir.mkdir()
    scored.to_parquet(
        retry_dir / "0.parquet",
        schema=schemas.LIGAND_3D_SCORE_SCHEMA,
        index=False,
    )

    with pytest.raises(ValueError, match="retries remain incomplete"):
        finalize_ligand_3d_retries(tmp_path)

    second_retry = retry_work[retry_work["retry_batch_index"].eq(1)]
    second_scored = second_retry[
        [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    ].copy()
    second_scored["shape"] = 0.2
    second_scored["color"] = 0.1
    second_scored["sucos_shape"] = 0.1
    second_scored.to_parquet(
        retry_dir / "1.parquet",
        schema=schemas.LIGAND_3D_SCORE_SCHEMA,
        index=False,
    )

    finalized = finalize_ligand_3d_retries(tmp_path)
    assert finalized["completed_retry_batch_count"] == 2
    combined = pd.read_parquet(canonical / "1.parquet")
    assert len(combined) == 2
    assert combined["shape"].notna().sum() == 2


def test_ligand_3d_plan_deduplicates_positive_pocket_candidates(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        LIGAND_3D_WORK_RELATIVE,
        SCORE_WORK_RELATIVE,
        _ligand_3d_batch,
        plan_ligand_3d_batches,
        plan_protein_scoring,
    )

    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "chain_receptor_type": ["protein"] * 3,
            "chain_is_holo": [True, True, False],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi"],
            "ligand_asym_id": ["B", "Y", "Z"],
            "ligand_num_heavy_atoms": [10, 20, 30],
            "ligand_is_proper": [True] * 3,
            "ligand_is_3d_score_able": [True] * 3,
            "system_type": ["holo"] * 3,
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)
    plan_protein_scoring(tmp_path)
    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2def"],
            "score_batch_index": [0, 0],
            "estimated_work": [1.0, 1.0],
        }
    ).to_parquet(tmp_path / SCORE_WORK_RELATIVE, index=False)

    candidate_dir = tmp_path / "scores/ligand_3d_candidates/search_db=holo/shard=ab"
    candidate_dir.mkdir(parents=True)
    repeated = {
        "query_system": "1abc_system",
        "query_ligand_id": "1abc__1__1.B",
        "query_entry": "1abc",
        "query_ligand_asym_id": "B",
        "target_system": "2def_system",
        "target_ligand_id": "2def__1__1.Y",
        "target_entry": "2def",
        "target_ligand_asym_id": "Y",
        "protein_mapping": "1.A:1.X",
        "protein_mapper": "foldseek",
        "pocket_qcov": 0.5,
    }
    pq.write_table(
        pa.Table.from_pylist(
            [
                repeated,
                {**repeated, "query_system": "1abc_other", "pocket_qcov": 0.75},
                {
                    **repeated,
                    "target_system": "3ghi_system",
                    "target_ligand_id": "3ghi__1__1.Z",
                    "target_entry": "3ghi",
                    "target_ligand_asym_id": "Z",
                    "pocket_qcov": 0.29,
                },
            ],
            schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA,
        ),
        candidate_dir / "1abc.parquet",
    )
    target_candidate_dir = candidate_dir.parent / "shard=de"
    target_candidate_dir.mkdir()
    pq.write_table(
        pa.Table.from_pylist([], schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA),
        target_candidate_dir / "2def.parquet",
    )
    tasks.collate_ligand_3d_candidates(
        data_dir=tmp_path,
        shards=["ab", "de"],
        scratch_dir=tmp_path / "candidate-scratch",
        threads=1,
    )
    compact = pd.read_parquet(
        tmp_path / "scores/ligand_3d_pair_candidate_shards/shard=ab.parquet"
    )
    assert len(compact) == 2
    assert compact.loc[
        compact["target_entry"] == "2def", "pocket_qcov"
    ].item() == pytest.approx(0.75)
    # Packed candidate shards are self-contained planning inputs. The
    # per-query score files may be archived or unavailable after collation.
    (candidate_dir / "1abc.parquet").unlink()
    (target_candidate_dir / "2def.parquet").unlink()

    report = plan_ligand_3d_batches(
        tmp_path,
        batch_size=2,
        threads=1,
        scratch_dir=tmp_path / "scratch",
    )

    assert report["pair_count"] == 2
    assert report["batch_count"] == 1
    work = pd.read_parquet(tmp_path / LIGAND_3D_WORK_RELATIVE)
    assert not work.duplicated(
        [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    ).any()
    assert set(work["estimated_work"]) == {200, 300}
    assert set(work["target_entry"]) == {"2def", "3ghi"}
    assert len(_ligand_3d_batch(tmp_path, 0, 2)) == 2


def test_make_ligand_3d_scores_writes_and_reuses_complete_batch(
    tmp_path, monkeypatch
) -> None:
    pairs = pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "query_ligand_asym_id": ["B"],
            "target_entry": ["2def"],
            "target_ligand_asym_id": ["Y"],
        }
    )
    calls = 0

    class FakeScorer:
        shape_score_threads = 1

        def score_canonical_ligand_pairs(self, _data_dir, observed_pairs):
            nonlocal calls
            calls += 1
            result = observed_pairs.copy()
            result["shape"] = 0.8
            result["color"] = 0.6
            result["sucos_shape"] = 0.7
            return result

    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: (FakeScorer(), [], tmp_path / "unused"),
    )
    kwargs = {
        "data_dir": tmp_path,
        "pairs": pairs,
        "batch_index": 3,
        "scorer_cfg": SimpleNamespace(),
        "force_update": False,
        "scratch_dir": tmp_path / "scratch",
        "threads": 2,
    }

    output = tasks.make_ligand_3d_scores(**kwargs)
    assert tasks.make_ligand_3d_scores(**kwargs) == output
    assert calls == 1
    assert pd.read_parquet(output)["sucos_shape"].item() == pytest.approx(0.7)


def test_collate_ligand_3d_scores_repartitions_balanced_batches(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        LIGAND_3D_WORK_RELATIVE,
        SCORE_WORK_RELATIVE,
    )

    work = pd.DataFrame(
        {
            "query_entry": ["1abc", "2abc"],
            "query_ligand_asym_id": ["B", "C"],
            "target_entry": ["3def", "4ghi"],
            "target_ligand_asym_id": ["Y", "Z"],
            "estimated_work": [10, 20],
            "ligand_3d_batch_index": [0, 1],
        }
    )
    work_path = tmp_path / LIGAND_3D_WORK_RELATIVE
    work_path.parent.mkdir(parents=True)
    work.to_parquet(work_path, index=False)
    pd.DataFrame({"pdb_id": ["1abc", "2abc"]}).to_parquet(
        tmp_path / SCORE_WORK_RELATIVE, index=False
    )
    pair_dir = tmp_path / "scores/ligand_3d_pairs"
    pair_dir.mkdir(parents=True)
    for row in work.itertuples(index=False):
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {
                        "query_entry": row.query_entry,
                        "query_ligand_asym_id": row.query_ligand_asym_id,
                        "target_entry": row.target_entry,
                        "target_ligand_asym_id": row.target_ligand_asym_id,
                        "shape": 0.8,
                        "color": 0.6,
                        "sucos_shape": 0.7,
                    }
                ],
                schema=schemas.LIGAND_3D_SCORE_SCHEMA,
            ),
            pair_dir / f"{row.ligand_3d_batch_index}.parquet",
        )

    outputs = tasks.collate_ligand_3d_scores(
        data_dir=tmp_path,
        shards=["ab"],
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    assert outputs == [tmp_path / "scores/ligand_3d_by_query/ab.parquet"]
    observed = pd.read_parquet(outputs[0])
    assert observed["query_entry"].tolist() == ["1abc", "2abc"]
    assert observed["shape"].tolist() == [0.8, 0.8]


def test_finalize_ligand_3d_scores_validates_pair_and_packed_shards(
    tmp_path, monkeypatch
) -> None:
    from plinder.data.pipeline.score import (
        SCORE_WORK_RELATIVE,
        finalize_ligand_3d_artifacts,
        plan_ligand_3d_batches,
        plan_protein_scoring,
    )

    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "chain_receptor_type": ["protein", "protein"],
            "chain_is_holo": [True, False],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "ligand_asym_id": ["B", "Y"],
            "ligand_num_heavy_atoms": [10, 20],
            "ligand_is_proper": [True, True],
            "ligand_is_3d_score_able": [True, True],
            "system_type": ["holo", "holo"],
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)
    plan_protein_scoring(tmp_path)
    pd.DataFrame(
        {
            "pdb_id": ["1abc"],
            "score_batch_index": [0],
            "estimated_work": [1.0],
        }
    ).to_parquet(tmp_path / SCORE_WORK_RELATIVE, index=False)
    candidate_dir = tmp_path / "scores/ligand_3d_candidates/search_db=holo/shard=ab"
    candidate_dir.mkdir(parents=True)
    candidate = {
        "query_system": "1abc_system",
        "query_ligand_id": "1abc__1__1.B",
        "query_entry": "1abc",
        "query_ligand_asym_id": "B",
        "target_system": "2def_system",
        "target_ligand_id": "2def__1__1.Y",
        "target_entry": "2def",
        "target_ligand_asym_id": "Y",
        "protein_mapping": "1.A:1.X",
        "protein_mapper": "foldseek",
        "pocket_qcov": 0.5,
    }
    pq.write_table(
        pa.Table.from_pylist([candidate], schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA),
        candidate_dir / "1abc.parquet",
    )
    tasks.collate_ligand_3d_candidates(
        data_dir=tmp_path,
        shards=["ab"],
        scratch_dir=tmp_path / "candidate-scratch",
        threads=1,
    )
    plan_ligand_3d_batches(
        tmp_path,
        batch_size=10,
        threads=1,
        scratch_dir=tmp_path / "plan-scratch",
    )
    pair_dir = tmp_path / "scores/ligand_3d_pairs"
    pair_dir.mkdir(parents=True)
    pair = {
        key: candidate[key]
        for key in [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    }
    pair.update({"shape": 0.8, "color": 0.6, "sucos_shape": 0.7})
    pq.write_table(
        pa.Table.from_pylist([pair], schema=schemas.LIGAND_3D_SCORE_SCHEMA),
        pair_dir / "0.parquet",
    )
    final_score = tmp_path / "scores/search_db=holo/ab.parquet"
    final_score.parent.mkdir(parents=True)
    pd.DataFrame(columns=schemas.PROTEIN_SIMILARITY_SCHEMA.names).to_parquet(
        final_score,
        index=False,
        schema=schemas.PROTEIN_SIMILARITY_SCHEMA,
    )

    report = finalize_ligand_3d_artifacts(tmp_path)

    assert report == {
        "status": "complete",
        "canonical_pair_count": 1,
        "pair_batch_count": 1,
        "query_shard_count": 1,
        "query_count": 1,
    }
    assert (tmp_path / "scores/ligand_3d_pair_validation.json").is_file()

    import duckdb

    monkeypatch.setattr(
        duckdb,
        "connect",
        lambda: pytest.fail("signature-matched pair validation was not reused"),
    )
    assert finalize_ligand_3d_artifacts(tmp_path) == report


def test_finalize_score_repair_accepts_refreshed_candidate_shards(tmp_path) -> None:
    from plinder.data.pipeline.score import finalize_score_repair_artifacts

    repair_manifest = tmp_path / "manifests/repair.parquet"
    repair_manifest.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": ["1abc"], "shard": ["ab"]}).to_parquet(
        repair_manifest, index=False
    )
    candidate = {
        "query_system": "1abc_system",
        "query_ligand_id": "1abc__1__1.B",
        "query_entry": "1abc",
        "query_ligand_asym_id": "B",
        "target_system": "2def_system",
        "target_ligand_id": "2def__1__1.Y",
        "target_entry": "2def",
        "target_ligand_asym_id": "Y",
        "protein_mapping": "1.A:1.X",
        "protein_mapper": "foldseek",
        "pocket_qcov": 0.5,
    }
    candidate_path = tmp_path / "scores/ligand_3d_candidate_shards/shard=ab.parquet"
    candidate_path.parent.mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist([candidate], schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA),
        candidate_path,
    )
    pair_candidate_path = (
        tmp_path / "scores/ligand_3d_pair_candidate_shards/shard=ab.parquet"
    )
    pair_candidate_path.parent.mkdir(parents=True)
    pair_candidate = {
        key: candidate[key]
        for key in [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
            "pocket_qcov",
        ]
    }
    pq.write_table(
        pa.Table.from_pylist(
            [pair_candidate], schema=schemas.LIGAND_3D_PAIR_CANDIDATE_SCHEMA
        ),
        pair_candidate_path,
    )
    candidate_stat = candidate_path.stat()
    pair_candidate_stat = pair_candidate_path.stat()
    candidate_path.with_suffix(".json").write_text(
        json.dumps(
            {
                "shard": "ab",
                "inputs": [],
                "output": {
                    "path": str(candidate_path.resolve()),
                    "size": candidate_stat.st_size,
                    "mtime_ns": candidate_stat.st_mtime_ns,
                    "rows": 1,
                },
                "pair_output": {
                    "path": str(pair_candidate_path.resolve()),
                    "size": pair_candidate_stat.st_size,
                    "mtime_ns": pair_candidate_stat.st_mtime_ns,
                    "rows": 1,
                },
            }
        )
    )

    pair_path = tmp_path / "scores/ligand_3d_by_query/ab.parquet"
    pair_path.parent.mkdir(parents=True)
    pair = {
        key: candidate[key]
        for key in [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    }
    pair.update({"shape": 0.8, "color": 0.6, "sucos_shape": 0.7})
    pq.write_table(
        pa.Table.from_pylist([pair], schema=schemas.LIGAND_3D_SCORE_SCHEMA),
        pair_path,
    )
    score_path = tmp_path / "scores/search_db=holo/ab.parquet"
    score_path.parent.mkdir(parents=True)
    pd.DataFrame(columns=schemas.PROTEIN_SIMILARITY_SCHEMA.names).to_parquet(
        score_path,
        index=False,
        schema=schemas.PROTEIN_SIMILARITY_SCHEMA,
    )

    report = finalize_score_repair_artifacts(
        tmp_path,
        repair_manifest=repair_manifest,
    )

    assert report["status"] == "complete"
    assert report["shard_count"] == 1
    assert report["candidate_rows"] == 1
    assert report["pair_rows"] == 1
    assert (tmp_path / "scores/score_repair_manifest.json").is_file()

    score_mtime_ns = score_path.stat().st_mtime_ns
    os.utime(
        pair_path,
        ns=(score_mtime_ns + 1_000_000_000, score_mtime_ns + 1_000_000_000),
    )
    with pytest.raises(ValueError, match="predates its current inputs"):
        finalize_score_repair_artifacts(
            tmp_path,
            repair_manifest=repair_manifest,
        )


def test_merge_ligand_3d_scores_uses_full_precision_pocket_coverage(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        LIGAND_3D_WORK_RELATIVE,
        SCORE_WORK_RELATIVE,
    )

    score_path = tmp_path / "dbs/subdbs/search_db=holo/1abc.parquet"
    score_path.parent.mkdir(parents=True)
    base = pd.DataFrame(
        [
            {
                "query_system": "1abc_system",
                "query_ligand_id": "1abc__1__1.B",
                "target_system": "2def_system",
                "target_ligand_id": "2def__1__1.Y",
                "protein_mapping": "1.A:1.X",
                "mapping": "1.A:1.X",
                "protein_mapper": "foldseek",
                "source": "foldseek",
                "metric": "pocket_qcov",
                "similarity": 67,
            }
        ]
    )
    base.to_parquet(
        score_path,
        index=False,
        schema=schemas.PROTEIN_SIMILARITY_SCHEMA.with_metadata(
            {b"plinder.ligand_3d": b"deferred"}
        ),
    )
    candidate_path = tmp_path / "scores/ligand_3d_candidate_shards/shard=ab.parquet"
    candidate_path.parent.mkdir(parents=True)
    candidate = {
        "query_system": "1abc_system",
        "query_ligand_id": "1abc__1__1.B",
        "query_entry": "1abc",
        "query_ligand_asym_id": "B",
        "target_system": "2def_system",
        "target_ligand_id": "2def__1__1.Y",
        "target_entry": "2def",
        "target_ligand_asym_id": "Y",
        "protein_mapping": "1.A:1.X",
        "protein_mapper": "foldseek",
        "pocket_qcov": 2 / 3,
    }
    pq.write_table(
        pa.Table.from_pylist([candidate], schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA),
        candidate_path,
    )
    pair_path = tmp_path / "scores/ligand_3d_by_query/ab.parquet"
    pair_path.parent.mkdir(parents=True)
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "query_entry": "1abc",
                    "query_ligand_asym_id": "B",
                    "target_entry": "2def",
                    "target_ligand_asym_id": "Y",
                    "shape": 0.8,
                    "color": 0.6,
                    "sucos_shape": 0.46,
                }
            ],
            schema=schemas.LIGAND_3D_SCORE_SCHEMA,
        ),
        pair_path,
    )
    work = pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "query_ligand_asym_id": ["B"],
            "target_entry": ["2def"],
            "target_ligand_asym_id": ["Y"],
            "estimated_work": [200],
            "ligand_3d_batch_index": [0],
        }
    )
    work_path = tmp_path / LIGAND_3D_WORK_RELATIVE
    work_path.parent.mkdir(exist_ok=True, parents=True)
    work.to_parquet(work_path, index=False)
    pd.DataFrame({"pdb_id": ["1abc"]}).to_parquet(
        tmp_path / SCORE_WORK_RELATIVE, index=False
    )
    pair_batch_dir = tmp_path / "scores/ligand_3d_pairs"
    pair_batch_dir.mkdir(parents=True)
    copyfile(pair_path, pair_batch_dir / "0.parquet")

    tasks.merge_ligand_3d_scores(
        data_dir=tmp_path,
        shards=["ab"],
        scorer_cfg=SimpleNamespace(minimum_threshold=0.3, minimum_thresholds={}),
        force_update=False,
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    merged_path = tmp_path / "scores/search_db=holo/ab.parquet"
    merged = pd.read_parquet(merged_path).set_index("metric")
    assert merged["similarity"].to_dict() == {
        "shape": 80,
        "pocket_qcov": 67,
        "color": 60,
        "sucos_shape": 46,
        "sucos_shape_pocket_qcov": 31,
    }
    assert pd.read_parquet(score_path)["metric"].tolist() == ["pocket_qcov"]


def test_merge_ligand_3d_scores_fails_when_query_shard_is_not_ready(
    tmp_path, monkeypatch
) -> None:
    score_work = tmp_path / "manifests" / "protein_scoring_work.parquet"
    score_work.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": ["1abc"]}).to_parquet(score_work, index=False)
    monkeypatch.setattr(
        tasks,
        "_ligand_3d_query_shard_is_ready",
        lambda **_kwargs: False,
    )

    with pytest.raises(RuntimeError, match="retry this merge task"):
        tasks.merge_ligand_3d_scores(
            data_dir=tmp_path,
            shards=["ab"],
            scorer_cfg=SimpleNamespace(
                minimum_threshold=0.3,
                minimum_thresholds={},
            ),
            force_update=False,
            scratch_dir=tmp_path / "scratch",
            threads=1,
        )


def test_protein_scoring_finalizer_rejects_changed_chain_index(tmp_path) -> None:
    from plinder.data.pipeline.score import (
        finalize_alignment_artifacts,
        plan_protein_scoring,
    )

    index_dir = tmp_path / "index"
    index_dir.mkdir()
    chain_path = index_dir / "entry_chains.parquet"
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_receptor_type": ["protein"],
            "chain_is_holo": [True],
        }
    ).to_parquet(chain_path, index=False)
    plan_protein_scoring(tmp_path)
    chain_path.write_bytes(chain_path.read_bytes() + b"changed")

    with pytest.raises(ValueError, match="entry chain index changed"):
        finalize_alignment_artifacts(tmp_path)


def test_v3_source_roots_are_configurable() -> None:
    from plinder.data.pipeline.config import get_config

    cfg = get_config(
        config={
            "source": {
                "pdb_nextgen_root": "/archives/nextgen",
                "validation_root": "/archives/validation",
                "discovery_threads": 3,
            }
        },
        cached=False,
    )

    assert cfg.source.pdb_nextgen_root == "/archives/nextgen"
    assert cfg.source.validation_root == "/archives/validation"
    assert cfg.source.discovery_threads == 3


def test_make_dbs_uses_configured_source_files(tmp_path, monkeypatch) -> None:
    cif_root = tmp_path / "nextgen"
    cif_root.mkdir()
    seqres_path = tmp_path / "pdb_seqres.txt.gz"
    seqres_path.touch()
    create_calls = []
    index_calls = []
    monkeypatch.setattr(
        tasks.databases,
        "create_db",
        lambda source, output, kind, threads: create_calls.append(
            (source, output, kind, threads)
        ),
    )
    monkeypatch.setattr(
        tasks.databases,
        "create_db_index",
        lambda output, kind, tmp_dir, threads: index_calls.append(
            (output, kind, tmp_dir, threads)
        ),
    )

    tasks.make_dbs(
        data_dir=tmp_path,
        sub_databases=["holo"],
        cpu=3,
        cif_root=cif_root,
        seqres_path=seqres_path,
    )

    assert create_calls == [
        (cif_root, tmp_path / "dbs/foldseek", "foldseek", 3),
        (seqres_path, tmp_path / "dbs/mmseqs", "mmseqs", 3),
    ]
    assert index_calls == [
        (
            tmp_path / "dbs/foldseek",
            "foldseek",
            tmp_path / "scratch/databases/foldseek",
            3,
        ),
        (
            tmp_path / "dbs/mmseqs",
            "mmseqs",
            tmp_path / "scratch/databases/mmseqs",
            3,
        ),
    ]


def test_make_dbs_reuses_completed_createdb_output(tmp_path, monkeypatch) -> None:
    cif_root = tmp_path / "nextgen"
    cif_root.mkdir()
    seqres_path = tmp_path / "pdb_seqres.txt.gz"
    seqres_path.touch()
    monkeypatch.setattr(
        tasks.databases,
        "created_database_is_complete",
        lambda *_args: True,
    )
    monkeypatch.setattr(
        tasks.databases,
        "create_db",
        lambda *_args, **_kwargs: pytest.fail("completed database was rebuilt"),
    )

    tasks.make_dbs(
        data_dir=tmp_path,
        sub_databases=["holo"],
        cpu=2,
        cif_root=cif_root,
        seqres_path=seqres_path,
        index=False,
    )


def test_get_scorer_uses_configured_search_limits(tmp_path) -> None:
    from plinder.data.pipeline import utils
    from plinder.data.pipeline.config import get_config

    cfg = get_config(
        config={
            "foldseek": {"max_seqs": 12_345},
            "mmseqs": {"max_seqs": 6_789},
            "scorer": {"sub_databases": "holo"},
        },
        cached=False,
    )
    scorer, _, _ = utils.get_scorer(
        data_dir=tmp_path,
        pdb_ids=[],
        scorer_cfg=cfg.scorer,
        load_entries=False,
        foldseek_cfg=cfg.foldseek,
        mmseqs_cfg=cfg.mmseqs,
        scratch_dir=tmp_path / "scratch",
    )

    assert scorer.get_config("holo", "foldseek").max_seqs == 12_345
    assert scorer.get_config("holo", "mmseqs").max_seqs == 6_789
    assert scorer.get_config("holo", "foldseek").min_seq_id == 0.2
    assert scorer.get_config("holo", "mmseqs").min_seq_id == 0.2
    assert scorer.get_config("apo", "foldseek").min_seq_id == 0.9
    assert scorer.get_config("pred", "mmseqs").min_seq_id == 0.9
    assert scorer.get_config("apo", "foldseek").coverage == 0.9
    assert scorer.get_config("pred", "mmseqs").coverage == 0.9

    from plinder.data.pipeline.score import _scoring_config

    v3_scoring = _scoring_config(tmp_path, 10_000)
    assert v3_scoring.foldseek.min_seq_id == 0.0
    assert v3_scoring.mmseqs.min_seq_id == 0.0


def test_run_batch_searches_skips_completed_backend_queries(
    tmp_path, monkeypatch
) -> None:
    completed = tmp_path / "dbs/subdbs/holo_foldseek/aln/1abc.parquet"
    completed.parent.mkdir(parents=True)
    completed.touch()
    calls = []

    class FakeScorer:
        def run_alignments(self, **kwargs):
            calls.append(kwargs)
            alignment_type = kwargs["alignment_types"][0]
            output = tmp_path / f"dbs/subdbs/holo_{alignment_type}/aln"
            output.mkdir(parents=True, exist_ok=True)
            for pdb_id in kwargs["entry_ids"]:
                (output / f"{pdb_id}.parquet").touch()

    def fake_get_scorer(*, pdb_ids, **_kwargs):
        scratch = tmp_path / "scratch" / "-".join(pdb_ids)
        scratch.mkdir(parents=True)
        return FakeScorer(), pdb_ids, scratch

    monkeypatch.setattr(tasks.utils, "get_scorer", fake_get_scorer)
    monkeypatch.setattr(
        tasks.databases,
        "database_identifiers",
        lambda path: (
            {"pdb_00002def_A"}
            if path.name.endswith("foldseek")
            else {"1abc_A", "2def_A"}
        ),
    )
    cfg = SimpleNamespace(sub_databases=["holo"])

    tasks.run_batch_searches(
        data_dir=tmp_path,
        pdb_ids=["1abc", "2def"],
        scorer_cfg=cfg,
        foldseek_cfg=SimpleNamespace(),
        mmseqs_cfg=SimpleNamespace(),
        cpu=4,
        scratch_dir=tmp_path / "node-scratch",
        alignment_types=["foldseek", "mmseqs"],
    )

    assert [call["entry_ids"] for call in calls] == [
        ["2def"],
        ["1abc", "2def"],
    ]
    assert [call["alignment_types"] for call in calls] == [
        ["foldseek"],
        ["mmseqs"],
    ]


def test_run_batch_searches_rejects_missing_eligible_output(tmp_path, monkeypatch):
    class FakeScorer:
        def run_alignments(self, **_kwargs):
            return None

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: (FakeScorer(), ["1abc"], scratch),
    )
    monkeypatch.setattr(
        tasks.databases,
        "database_identifiers",
        lambda _path: {"1abc_A"},
    )

    with pytest.raises(RuntimeError, match="no output.*1 eligible"):
        tasks.run_batch_searches(
            data_dir=tmp_path,
            pdb_ids=["1abc"],
            scorer_cfg=SimpleNamespace(sub_databases=["holo"]),
            foldseek_cfg=SimpleNamespace(),
            mmseqs_cfg=SimpleNamespace(),
            cpu=1,
            alignment_types=["mmseqs"],
        )


@pytest.mark.parametrize("alignment_type", ["foldseek", "mmseqs"])
def test_alignment_release_shard_preserves_typed_schema_when_all_hits_are_empty(
    tmp_path, alignment_type
) -> None:
    schema = schemas.mapped_alignment_schema(alignment_type=alignment_type)
    source = tmp_path / "mapped" / "1abc.parquet"
    source.parent.mkdir()
    pq.write_table(pa.Table.from_pylist([], schema=schema), source)
    target = tmp_path / f"shard-{alignment_type}.parquet"

    tasks._write_alignment_release_shard(
        sources=[source],
        target=target,
        alignment_type=alignment_type,
        temp_dir=tmp_path / "scratch" / alignment_type,
        threads=1,
        memory_limit="1GB",
    )

    assert pq.ParquetFile(target).metadata.num_rows == 0
    assert pq.read_schema(target).equals(schema)


def test_alignment_release_shard_unifies_empty_and_populated_list_types(
    tmp_path,
) -> None:
    columns = {
        "query_entry": ["1abc"],
        "target_entry": ["2def"],
        "query_chain_mapped": ["A"],
        "target_chain_mapped": ["B"],
        "source": ["foldseek"],
        "qcov": [1.0],
        "fident": [1.0],
        "seqsim": [1.0],
        "query_pocket_residue_numbers": [[]],
        "target_pocket_residue_numbers": [[]],
        "pocket_residue_identity": [b""],
        "lddt": [1.0],
    }
    empty_lists = tmp_path / "mapped" / "1abc.parquet"
    populated_lists = tmp_path / "mapped" / "2abc.parquet"
    empty_lists.parent.mkdir()
    pd.DataFrame(columns).to_parquet(empty_lists, index=False)
    populated = dict(columns)
    populated["query_entry"] = ["2abc"]
    populated["query_pocket_residue_numbers"] = [[1]]
    populated["target_pocket_residue_numbers"] = [[2]]
    populated["pocket_residue_identity"] = [bytes([1])]
    pd.DataFrame(populated).to_parquet(populated_lists, index=False)
    target = tmp_path / "foldseek.parquet"

    tasks._write_alignment_release_shard(
        sources=[empty_lists, populated_lists],
        target=target,
        alignment_type="foldseek",
        temp_dir=tmp_path / "scratch",
        threads=1,
        memory_limit="1GB",
    )

    result = pd.read_parquet(target)
    assert result["query_entry"].tolist() == ["1abc", "2abc"]
    assert result.loc[0, "query_pocket_residue_numbers"].tolist() == []
    assert result.loc[1, "query_pocket_residue_numbers"].tolist() == [1]


def test_ligand_score_threshold_must_cover_requested_chemical_clusters() -> None:
    from plinder.data.pipeline.config import get_config

    with pytest.raises(ValueError, match=r"clustering threshold \(30\)"):
        get_config(
            config={"ligand": {"minimum_similarity": 60}},
            cached=False,
        )


def test_derived_score_threshold_must_cover_requested_clusters() -> None:
    from plinder.data.pipeline.config import get_config

    with pytest.raises(ValueError, match=r"clustering threshold \(30\)"):
        get_config(
            config={"scorer": {"minimum_threshold": 0.31}},
            cached=False,
        )


def test_collate_alignments_writes_query_addressable_shards(tmp_path):
    columns = {
        "query_entry": ["1abc", "2abd", "3xyz"],
        "target_entry": ["9zzz", "8yyy", "7xxx"],
        "query_chain_mapped": ["A", "B", "C"],
        "target_chain_mapped": ["D", "E", "F"],
        "source": ["foldseek"] * 3,
        "qcov": [1.0, 0.9, 0.8],
        "fident": [1.0, 0.9, 0.8],
        "seqsim": [1.0, 0.9, 0.8],
        "query_pocket_residue_numbers": [[1], [2], [3]],
        "target_pocket_residue_numbers": [[11], [12], [13]],
        "pocket_residue_identity": [bytes([1])] * 3,
        "lddt": [1.0, 0.9, 0.8],
    }
    mapped_dir = tmp_path / "dbs/subdbs/holo_foldseek/mapped_aln"
    mapped_dir.mkdir(parents=True)
    for index, pdb_id in enumerate(columns["query_entry"]):
        pd.DataFrame(
            {column: [values[index]] for column, values in columns.items()}
        ).to_parquet(mapped_dir / f"{pdb_id}.parquet", index=False)
    # An interrupted atomic install must not become a query or shard.
    pd.DataFrame(
        {column: [values[0]] for column, values in columns.items()}
    ).to_parquet(mapped_dir / "4qp3.tmp.parquet", index=False)

    assert tasks.scatter_collate_alignments(data_dir=tmp_path) == [["ab"], ["xy"]]
    tasks.collate_alignments(data_dir=tmp_path, partition=["ab"])

    shard = pd.read_parquet(
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet"
    )
    assert shard["query_entry"].tolist() == ["1abc", "2abd"]
    assert "3xyz" not in set(shard["query_entry"])
    assert set(shard.columns) == set(columns)


def test_empty_alignment_scatter_has_noop_branch(tmp_path):
    assert tasks.scatter_collate_alignments(data_dir=tmp_path) == [[]]
    tasks.collate_alignments(data_dir=tmp_path, partition=[])


def test_mapping_scatter_requires_current_shard_manifest(tmp_path):
    raw_dir = tmp_path / "dbs/subdbs/holo_foldseek/aln"
    raw_dir.mkdir(parents=True)
    pd.DataFrame({"query": ["1abc_A"]}).to_parquet(
        raw_dir / "1abc.parquet", index=False
    )

    assert tasks.scatter_missing_alignment_mappings(
        data_dir=tmp_path, batch_size=10
    ) == [["ab"]]

    current_columns = {
        "query_entry": ["1abc"],
        "target_entry": ["1abc"],
        "query_chain_mapped": ["A"],
        "target_chain_mapped": ["A"],
        "source": ["foldseek"],
        "qcov": [1.0],
        "fident": [1.0],
        "seqsim": [1.0],
        "query_pocket_residue_numbers": [[1]],
        "target_pocket_residue_numbers": [[1]],
        "pocket_residue_identity": [bytes([1])],
        "lddt": [1.0],
    }
    release = tasks._alignment_release_path(
        data_dir=tmp_path,
        search_db="holo",
        alignment_type="foldseek",
        shard="ab",
    )
    release.parent.mkdir(parents=True)
    pd.DataFrame(current_columns).to_parquet(release, index=False)
    _write_alignment_mapping_manifest(tmp_path, "ab")
    assert tasks.scatter_missing_alignment_mappings(
        data_dir=tmp_path, batch_size=10
    ) == [[]]


def test_map_batch_alignments_publishes_atomic_shard(tmp_path, monkeypatch):
    raw_dir = tmp_path / "dbs/subdbs/holo_foldseek/aln"
    raw_dir.mkdir(parents=True)
    pd.DataFrame({"query": ["1abc_A"], "target_pdb_id": ["1abc"]}).to_parquet(
        raw_dir / "1abc.parquet", index=False
    )
    _write_alignment_chain_lookup(tmp_path)
    calls = []

    class FakeScorer:
        entries = {}

        def map_alignment_files(self, *_args, **kwargs):
            calls.append(kwargs)
            output = kwargs["mapped_db_dir"] / "holo_foldseek/mapped_aln/1abc.parquet"
            output.parent.mkdir(parents=True)
            pd.DataFrame(
                {
                    "query_entry": ["1abc"],
                    "target_entry": ["2def"],
                    "query_chain_mapped": ["A"],
                    "target_chain_mapped": ["B"],
                    "source": ["foldseek"],
                    "qcov": [1.0],
                    "fident": [1.0],
                    "seqsim": [1.0],
                    "query_pocket_residue_numbers": [[1]],
                    "target_pocket_residue_numbers": [[2]],
                    "pocket_residue_identity": [bytes([1])],
                    "lddt": [1.0],
                }
            ).to_parquet(output, index=False)
            return [output]

    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: (FakeScorer(), ["1abc"], tmp_path / "unused"),
    )
    scratch = tmp_path / "scratch"
    tasks.map_batch_alignments(
        data_dir=tmp_path,
        shards=["ab"],
        scorer_cfg=SimpleNamespace(sub_databases=["holo"]),
        force_update=False,
        scratch_dir=scratch,
    )

    release = tasks._alignment_release_path(
        data_dir=tmp_path,
        search_db="holo",
        alignment_type="foldseek",
        shard="ab",
    )
    assert pd.read_parquet(release)["query_entry"].tolist() == ["1abc"]
    assert tasks.alignment_mapping_shard_is_current(data_dir=tmp_path, shard="ab")
    assert not (tmp_path / "dbs/subdbs/holo_foldseek/mapped_aln").exists()
    assert not any(scratch.iterdir())

    tasks.map_batch_alignments(
        data_dir=tmp_path,
        shards=["ab"],
        scorer_cfg=SimpleNamespace(sub_databases=["holo"]),
        force_update=False,
        scratch_dir=scratch,
    )
    assert len(calls) == 1


def test_map_batch_alignments_records_and_skips_oversized_queries(
    tmp_path, monkeypatch
) -> None:
    raw_dir = tmp_path / "dbs/subdbs/holo_foldseek/aln"
    raw_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "query": ["1abc_A", "1abc_A"],
            "target_pdb_id": ["1abc", "1abc"],
        }
    ).to_parquet(raw_dir / "1abc.parquet", index=False)
    _write_alignment_chain_lookup(tmp_path)

    class FakeScorer:
        entries = {}

        def map_alignment_files(self, *_args, **_kwargs):
            raise AssertionError("oversized query must not be mapped")

    monkeypatch.setattr(
        tasks.utils,
        "get_scorer",
        lambda **_kwargs: (FakeScorer(), [], tmp_path / "unused"),
    )

    tasks.map_batch_alignments(
        data_dir=tmp_path,
        shards=["ab"],
        scorer_cfg=SimpleNamespace(
            sub_databases=["holo"], max_alignment_rows_per_query=1
        ),
        force_update=True,
        scratch_dir=tmp_path / "scratch",
    )

    manifest = json.loads(
        tasks._alignment_mapping_manifest_path(
            data_dir=tmp_path, shard="ab"
        ).read_text()
    )
    assert manifest["skipped_queries"] == {
        "1abc": {
            "reason": "raw_alignment_row_budget_exceeded",
            "maximum_rows": 1,
            "total_rows": 2,
            "rows_by_alignment_type": {"foldseek": 2},
        }
    }
    release = tasks._alignment_release_path(
        data_dir=tmp_path,
        search_db="holo",
        alignment_type="foldseek",
        shard="ab",
    )
    assert pd.read_parquet(release).empty
    assert tasks.alignment_mapping_shard_is_current(data_dir=tmp_path, shard="ab")


def test_alignment_chain_lookup_compacts_mapping_inputs(tmp_path) -> None:
    from plinder.core.scores.entries import load_alignment_entry_views

    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc"],
            "ligand_is_proper": [True, False],
            "ligand_neighboring_residues": [
                ["1.A_10_9_10", "1.A_11_10_11"],
                ["1.A_99_98_99"],
            ],
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc"],
            "chain_asym_id": ["A", "B"],
            "chain_auth_id": ["X", "Y"],
            "chain_receptor_type": ["protein", "dna"],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)

    lookup = tasks.make_alignment_chain_lookup(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    frame = pd.read_parquet(lookup)
    assert frame["chain_asym_id"].tolist() == ["A"]
    assert frame.loc[0, "pocket_residue_numbers"].tolist() == [10, 11]
    assert frame.loc[0, "pocket_residue_indices"].tolist() == [9, 10]
    assert tasks._completed_alignment_chain_lookup(tmp_path) is not None
    manifest = tmp_path / tasks.ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE
    assert manifest.is_file()
    views = load_alignment_entry_views(lookup_path=lookup, pdb_ids=["1abc"])
    assert views["1abc"].author_to_asym == {"X": "A"}
    assert views["1abc"].pocket_index_to_number_per_chain == {"A": {9: 10, 10: 11}}

    annotation = index / "annotation_table.parquet"
    annotation_frame = pd.read_parquet(annotation)
    annotation_frame.loc[0, "ligand_is_proper"] = False
    annotation_frame.to_parquet(annotation, index=False)
    assert tasks._completed_alignment_chain_lookup(tmp_path) is None


def test_alignment_chain_lookup_keeps_identity_when_only_system_rows_change(
    tmp_path: Path,
) -> None:
    index = tmp_path / "index"
    index.mkdir()
    annotation = index / "annotation_table.parquet"
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "system_id": ["1abc__1__1.A__1.B"],
            "ligand_is_proper": [True],
            "ligand_neighboring_residues": [["1.A_10_9_10"]],
        }
    ).to_parquet(annotation, index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["X"],
            "chain_receptor_type": ["protein"],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    lookup = tasks.make_alignment_chain_lookup(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch-1",
        threads=1,
    )
    original_stat = lookup.stat()
    frame = pd.read_parquet(annotation)
    frame["system_id"] = "1abc__2__1.A__1.B"
    frame.to_parquet(annotation, index=False)

    tasks.make_alignment_chain_lookup(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch-2",
        threads=1,
        force_update=True,
    )

    refreshed_stat = lookup.stat()
    assert (refreshed_stat.st_size, refreshed_stat.st_mtime_ns) == (
        original_stat.st_size,
        original_stat.st_mtime_ns,
    )
    assert tasks._completed_alignment_chain_lookup(tmp_path) is not None


def test_alignment_chain_lookup_replaces_a_legacy_schema(tmp_path: Path) -> None:
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "ligand_is_proper": [True],
            "ligand_neighboring_residues": [["1.A_10_9_10"]],
        }
    ).to_parquet(index / "annotation_table.parquet", index=False)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["X"],
            "chain_receptor_type": ["protein"],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    lookup = tasks.make_alignment_chain_lookup(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch-1",
        threads=1,
    )
    legacy = pd.read_parquet(lookup)
    legacy["legacy_extra"] = "obsolete"
    legacy.to_parquet(lookup, index=False)

    tasks.make_alignment_chain_lookup(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch-2",
        threads=1,
        force_update=True,
    )

    assert pq.read_schema(lookup).names == [
        "entry_pdb_id",
        "chain_asym_id",
        "chain_auth_id",
        "pocket_residue_numbers",
        "pocket_residue_indices",
    ]
    assert tasks._completed_alignment_chain_lookup(tmp_path) is not None


def test_finalize_index_preserves_current_alignment_lookup(
    tmp_path, monkeypatch
) -> None:
    _write_alignment_chain_lookup(tmp_path)
    annotation = tmp_path / "index" / "annotation_table.parquet"

    def enrich_index(*, data_dir: Path) -> None:
        frame = pd.read_parquet(data_dir / "index" / "annotation_table.parquet")
        frame["example_cluster"] = "c0"
        frame.to_parquet(data_dir / "index" / "annotation_table.parquet", index=False)

    monkeypatch.setattr(tasks.utils, "finalize_index", enrich_index)
    monkeypatch.setattr(
        tasks.utils, "create_nonredundant_dataset", lambda *, data_dir: None
    )
    tasks.finalize_index(data_dir=tmp_path)

    assert "example_cluster" in pd.read_parquet(annotation).columns
    assert tasks._completed_alignment_chain_lookup(tmp_path) is not None


def test_component_reduction_scatter_reads_each_physical_source_once(
    tmp_path, monkeypatch
):
    score_a = tmp_path / "scores-a.parquet"
    score_b = tmp_path / "scores-b.parquet"
    ligand = tmp_path / "ligands.parquet"

    def sources(*, data_dir, metric):
        del data_dir
        if metric == "tanimoto_similarity_ecfp4_1024":
            return [ligand]
        return [score_a, score_b]

    monkeypatch.setattr(tasks.clusters, "component_score_sources", sources)

    assert tasks.scatter_component_reduction_sources(
        data_dir=tmp_path,
        metrics=[
            "protein_lddt_max",
            "sucos_shape_pocket_qcov",
            "tanimoto_similarity_ecfp4_1024",
        ],
        batch_size=2,
    ) == [[str(ligand), str(score_a)], [str(score_b)]]


def test_score_partition_collation_stages_atomically_from_scratch(tmp_path):
    source = tmp_path / "dbs" / "subdbs" / "search_db=holo"
    source.mkdir(parents=True)
    pd.DataFrame({"query_system": ["a"], "similarity": [30]}).to_parquet(
        source / "abc1.parquet", index=False
    )
    pd.DataFrame({"query_system": ["b"], "similarity": [50]}).to_parquet(
        source / "def1.parquet", index=False
    )

    tasks.collate_partitions(
        data_dir=tmp_path,
        partition=["1"],
        scratch_dir=tmp_path / "scratch",
        threads=1,
    )

    output = pd.read_parquet(tmp_path / "scores" / "search_db=holo" / "1.parquet")
    assert output.sort_values("query_system").to_dict("records") == [
        {"query_system": "a", "similarity": 30, "search_db": "holo"},
        {"query_system": "b", "similarity": 50, "search_db": "holo"},
    ]
    assert not (tmp_path / "scores" / "search_db=holo" / "1.parquet.tmp").exists()


def test_component_reduction_task_copies_generic_source_once(tmp_path, monkeypatch):
    source = (
        tmp_path
        / "ligand_clusters/symmetric_edges"
        / "metric=sucos_shape_pocket_qcov"
        / "bucket=000.parquet"
    )
    source.parent.mkdir(parents=True)
    source.write_bytes(b"scores")
    scratch = tmp_path / "scratch"
    calls = []

    monkeypatch.setattr(
        tasks.clusters,
        "component_node_universe",
        lambda **kwargs: (["l1", "l2"], {"s1", "s2"}),
    )
    monkeypatch.setattr(
        tasks.clusters,
        "score_component_reduction_is_complete",
        lambda **kwargs: False,
    )
    monkeypatch.setattr(
        tasks.clusters,
        "directed_cover_component_reduction_is_complete",
        lambda **kwargs: False,
    )

    def record_reciprocal_reduction(**kwargs):
        assert kwargs["read_path"].read_bytes() == b"scores"
        calls.append(("reciprocal", kwargs["metric"]))
        return {"outputs": []}

    monkeypatch.setattr(
        tasks.clusters,
        "make_score_component_reduction",
        record_reciprocal_reduction,
    )

    def record_cover_reduction(**kwargs):
        assert kwargs["read_path"].read_bytes() == b"scores"
        calls.append(("directed_cover", kwargs["metric"]))
        return {"outputs": []}

    monkeypatch.setattr(
        tasks.clusters,
        "make_directed_cover_component_reduction",
        record_cover_reduction,
    )

    tasks.make_component_reductions(
        data_dir=tmp_path,
        source_paths=[str(source)],
        metrics=["protein_lddt_max", "sucos_shape_pocket_qcov"],
        thresholds=[30, 50, 70, 90, 100],
        scratch_dir=scratch,
        force_update=False,
    )

    assert calls == [
        ("reciprocal", "sucos_shape_pocket_qcov"),
        ("directed_cover", "sucos_shape_pocket_qcov"),
    ]
    assert list(scratch.iterdir()) == []


def test_component_reduction_metric_workers_must_be_positive(tmp_path):
    with pytest.raises(ValueError, match="metric workers must be positive"):
        tasks.make_component_reductions(
            data_dir=tmp_path,
            source_paths=["unused.parquet"],
            metrics=["pocket_qcov"],
            thresholds=[30],
            scratch_dir=tmp_path / "scratch",
            force_update=False,
            metric_workers=0,
        )


def test_clustering_plan_matches_slurm_array_bounds(tmp_path, monkeypatch):
    from plinder.data.pipeline.score import _community_batch, plan_clustering

    monkeypatch.setattr(
        "plinder.data.pipeline.score.clusters.prepare_component_node_universe",
        lambda data_dir: {"status": "complete"},
    )
    monkeypatch.setattr(
        "plinder.data.pipeline.score.clusters.prepare_symmetric_edge_plan",
        lambda **kwargs: {"batches": [{}, {}, {}], "plan_hash": "test-plan"},
    )

    plan = plan_clustering(
        tmp_path,
        metrics=["pocket_qcov", "tanimoto_similarity_ecfp4_1024"],
        thresholds=[30, 100],
        source_batch_size=2,
        community_batch_size=3,
        symmetric_bucket_count=4,
    )

    assert plan["symmetric_fragment_batch_count"] == 3
    assert plan["symmetric_edge_shard_count"] == 8
    assert plan["component_reduction_batch_count"] == 8
    assert plan["community_task_count"] == 4
    assert plan["community_batch_count"] == 2
    assert plan["directed_cover_batch_count"] == 2
    assert _community_batch(
        metrics=plan["metrics"],
        thresholds=plan["thresholds"],
        batch_index=1,
        batch_size=3,
    ) == [("tanimoto_similarity_ecfp4_1024", 30)]


def test_sucos_release_export_retains_scores_below_cluster_cutoff(tmp_path):
    from plinder.data.pipeline.score import (
        MANIFEST_RELATIVE,
        PLAN_RELATIVE,
        SCORE_WORK_RELATIVE,
        _source_signature,
        export_sucos_shape_pocket_qcov_batch,
        finalize_sucos_shape_pocket_qcov_export,
    )

    manifest = tmp_path / MANIFEST_RELATIVE
    manifest.parent.mkdir(parents=True)
    pd.DataFrame({"pdb_id": ["1abc", "2def"], "shard": ["ab", "de"]}).to_parquet(
        manifest, index=False
    )
    score_work = tmp_path / SCORE_WORK_RELATIVE
    pd.DataFrame({"pdb_id": ["1abc", "2def"]}).to_parquet(score_work, index=False)
    (tmp_path / PLAN_RELATIVE).write_text(
        json.dumps({"manifest": _source_signature(manifest)})
    )

    rows = [
        ("ab", "1abc", "B", "2def", "Y", 0.2, 0.5),
        ("de", "2def", "Y", "1abc", "B", 0.8, 0.5),
    ]
    for shard, query, query_asym, target, target_asym, qcov, sucos in rows:
        candidate = {
            "query_system": f"{query}_system",
            "query_ligand_id": f"{query}__1__1.{query_asym}",
            "query_entry": query,
            "query_ligand_asym_id": query_asym,
            "target_system": f"{target}_system",
            "target_ligand_id": f"{target}__1__1.{target_asym}",
            "target_entry": target,
            "target_ligand_asym_id": target_asym,
            "protein_mapping": "1.A:1.X",
            "protein_mapper": "foldseek",
            "pocket_qcov": qcov,
        }
        candidate_path = (
            tmp_path
            / "scores"
            / "ligand_3d_candidate_shards"
            / f"shard={shard}.parquet"
        )
        candidate_path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(
            pa.Table.from_pylist(
                [candidate], schema=schemas.LIGAND_3D_CANDIDATE_SCHEMA
            ),
            candidate_path,
        )
        pair_path = tmp_path / "scores" / "ligand_3d_by_query" / f"{shard}.parquet"
        pair_path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {
                        "query_entry": query,
                        "query_ligand_asym_id": query_asym,
                        "target_entry": target,
                        "target_ligand_asym_id": target_asym,
                        "shape": 0.5,
                        "color": 0.5,
                        "sucos_shape": sucos,
                    }
                ],
                schema=schemas.LIGAND_3D_SCORE_SCHEMA,
            ),
            pair_path,
        )

    shard_dir = tmp_path / "release-shards"
    for batch_index in range(2):
        export_sucos_shape_pocket_qcov_batch(
            tmp_path,
            output_dir=shard_dir,
            batch_index=batch_index,
            batch_size=1,
            scratch_dir=tmp_path / "scratch" / str(batch_index),
            threads=1,
            memory_limit="1GB",
        )
    output = tmp_path / "exports" / "all_sucos_shape_pocket_qcov.parquet"
    report = finalize_sucos_shape_pocket_qcov_export(
        tmp_path,
        source_dir=shard_dir,
        output=output,
        scratch_dir=tmp_path / "finalize-scratch",
        threads=1,
        memory_limit="1GB",
    )

    exported = pd.read_parquet(output)
    assert sorted(exported["similarity"].tolist()) == [10, 40]
    assert report["row_count"] == 2
    assert (
        finalize_sucos_shape_pocket_qcov_export(
            tmp_path,
            source_dir=shard_dir,
            output=output,
            scratch_dir=tmp_path / "finalize-scratch",
            threads=1,
            memory_limit="1GB",
        )
        == report
    )


def test_clustering_statistics_validate_complete_monotonic_artifacts(tmp_path):
    from plinder.data.pipeline.score import summarize_clustering_artifacts
    from plinder.data.pipeline.utils import _read_local_cluster_rows

    metric = "pocket_qcov"
    for threshold, labels in [(100, ["c0", "c1"]), (50, ["c0", "c0"])]:
        for cluster, directed in [
            ("components", False),
            ("communities", False),
        ]:
            path = (
                tmp_path
                / "ligand_clusters"
                / f"cluster={cluster}"
                / f"directed={directed}"
                / f"metric={metric}"
                / f"threshold={threshold}.parquet"
            )
            path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                {
                    "ligand_id": ["l1", "l2"],
                    "label": labels,
                    "metric": [metric, metric],
                    "threshold": [threshold, threshold],
                    "cluster": [cluster, cluster],
                    "directed": [directed, directed],
                }
            ).to_parquet(path, index=False)
        directed_cover = (
            tmp_path
            / "ligand_sampling"
            / "directed_set_cover"
            / f"metric={metric}"
            / f"threshold={threshold}.parquet"
        )
        directed_cover.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "ligand_id": ["l1", "l2"],
                "centroid_ligand_id": ["l1", "l2"],
                "similarity_to_centroid": [100.0, 100.0],
                "label": labels,
                "metric": [metric, metric],
                "threshold": [threshold, threshold],
                "directed": [True, True],
            }
        ).to_parquet(directed_cover, index=False)

    report = summarize_clustering_artifacts(
        tmp_path,
        metrics=[metric],
        thresholds=[100, 50],
    )

    assert report["status"] == "complete"
    assert report["artifact_count"] == 6
    assert report["issues"] == []
    assert (tmp_path / "ligand_clusters" / "stats.json").is_file()
    stats = pd.read_parquet(tmp_path / "ligand_clusters" / "stats.parquet")
    components = stats[stats["cluster"].eq("components")].set_index("threshold")
    assert components.loc[100, "cluster_count"] == 2
    assert components.loc[50, "cluster_count"] == 1
    # The diagnostics live under the same root but are not cluster-label rows.
    cluster_rows = _read_local_cluster_rows(
        root=tmp_path / "ligand_clusters",
        node_column="ligand_id",
    )
    assert len(cluster_rows) == 8

    component_50 = (
        tmp_path
        / "ligand_clusters"
        / "cluster=components"
        / "directed=False"
        / f"metric={metric}"
        / "threshold=50.parquet"
    )
    pd.DataFrame(
        {
            "ligand_id": ["l1", "l2"],
            "label": ["c0", "c1"],
            "metric": [metric, metric],
            "threshold": [50, 50],
            "cluster": ["components", "components"],
            "directed": [False, False],
        }
    ).to_parquet(component_50, index=False)
    component_100 = component_50.with_name("threshold=100.parquet")
    pd.DataFrame(
        {
            "ligand_id": ["l1", "l2"],
            "label": ["c0", "c0"],
            "metric": [metric, metric],
            "threshold": [100, 100],
            "cluster": ["components", "components"],
            "directed": [False, False],
        }
    ).to_parquet(component_100, index=False)
    with pytest.raises(ValueError, match="invalid clustering artifacts"):
        summarize_clustering_artifacts(
            tmp_path,
            metrics=[metric],
            thresholds=[100, 50],
        )
    invalid_report = json.loads(
        (tmp_path / "ligand_clusters" / "stats.json").read_text()
    )
    assert invalid_report["status"] == "invalid"
    assert any("gains clusters" in issue for issue in invalid_report["issues"])


def test_v3_score_slurm_exposes_exact_clustering_stages():
    repository = Path(__file__).resolve().parents[3]
    script_path = repository / "scripts" / "slurm" / "score_v3.sbatch"
    if not script_path.is_file():
        pytest.skip("Slurm assets are not installed in wheel-only test layouts")
    script = script_path.read_text()

    assert "plan-clusters)" in script
    assert "symmetric-edge-fragments|symmetric-edge-shards" in script
    assert "component-reductions|communities|directed-covers)" in script
    assert (
        "finalize-alignments|finalize-ligands|finalize-scores|"
        "finalize-ligand-3d-retries|finalize-index|merge-components|cluster-stats)"
        in script
    )
    assert "PLINDER_CLUSTER_METRICS" in script
    assert "PLINDER_CLUSTER_THRESHOLDS" in script


def test_metaflow_graph_uses_canonical_ligand_archive_stage():
    repository = Path(__file__).resolve().parents[3]
    flow_path = repository / "flows" / "data_ingest.py"
    if not flow_path.is_file():
        pytest.skip("Metaflow sources are not installed in wheel-only test layouts")
    flow = flow_path.read_text()

    assert "def scatter_structure_qc" not in flow
    assert "self.pipeline.structure_qc" not in flow
    assert "self.next(self.scatter_make_entries)" in flow
    assert "def join_make_entries" in flow
    assert "self.next(self.scatter_collate_entries)" in flow
    assert "def collate_entries" in flow
    assert "self.pipeline.collate_entries(self.input)" in flow
    assert "self.pipeline.join_collate_entries" in flow
    assert "self.next(self.scatter_make_canonical_ligand_archives)" in flow
    assert "def make_canonical_ligand_archives" in flow
    assert "self.next(self.finalize_ligand_archives)" in flow
    assert "self.pipeline.finalize_ligand_archives()" in flow
    assert "self.next(self.annotate_ligand_similarity)" in flow
    assert "self.pipeline.annotate_ligand_similarity()" in flow
    assert "self.next(self.scatter_collate_partitions)" in flow
    assert "self.next(self.scatter_collate_alignments)" in flow
    assert "self.pipeline.collate_alignments(self.input)" in flow
    assert "self.next(self.finalize_alignments)" in flow
    assert "self.pipeline.finalize_alignments()" in flow
    assert "self.next(self.scatter_map_batch_alignments)" in flow
    assert "self.pipeline.map_batch_alignments(self.input)" in flow
    assert "self.pipeline.collate_ligand_3d_candidates(self.input)" in flow
    assert "self.next(self.plan_ligand_3d_scores)" in flow
    assert "self.pipeline.plan_ligand_3d_scores()" in flow
    assert "self.pipeline.make_ligand_3d_scores(self.input)" in flow
    assert "self.pipeline.collate_ligand_3d_scores(self.input)" in flow
    assert "self.pipeline.merge_ligand_3d_scores(self.input)" in flow
    assert "self.pipeline.finalize_scores()" in flow
    assert "self.next(self.scatter_export_sucos_shape_pocket_qcov)" in flow
    assert "self.pipeline.export_sucos_shape_pocket_qcov(self.input)" in flow
    assert "self.pipeline.finalize_sucos_export()" in flow
    assert "self.next(self.scatter_make_component_reductions)" in flow
    assert "self.pipeline.merge_component_reductions()" in flow
    assert "self.next(self.scatter_make_communities)" in flow
    assert "self.next(self.scatter_make_directed_set_covers)" in flow
    assert "self.pipeline.make_directed_set_covers(self.input)" in flow
    assert "self.next(self.summarize_clusters)" in flow
    assert "self.pipeline.summarize_clusters()" in flow
    assert "self.next(self.finalize_index)" in flow
    assert "self.pipeline.finalize_index()" in flow

    tree = ast.parse(flow)
    flow_class = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "PlinderDataIngestFlow"
    )
    steps = {
        node.name: node
        for node in flow_class.body
        if isinstance(node, ast.FunctionDef)
        and any(
            isinstance(decorator, ast.Name) and decorator.id == "step"
            for decorator in node.decorator_list
        )
    }
    edges = {name: set() for name in steps}
    for name, node in steps.items():
        for call in (child for child in ast.walk(node) if isinstance(child, ast.Call)):
            if not (
                isinstance(call.func, ast.Attribute)
                and call.func.attr == "next"
                and isinstance(call.func.value, ast.Name)
                and call.func.value.id == "self"
            ):
                continue
            for argument in call.args:
                if (
                    isinstance(argument, ast.Attribute)
                    and isinstance(argument.value, ast.Name)
                    and argument.value.id == "self"
                ):
                    edges[name].add(argument.attr)

    assert edges["start"] == {"scatter_make_entries"}
    assert edges["join_collate_entries"] == {"make_dbs"}
    assert edges["make_dbs"] == {"scatter_make_canonical_ligand_archives"}

    reachable = {"start"}
    pending = ["start"]
    while pending:
        for target in edges[pending.pop()]:
            assert target in steps
            if target not in reachable:
                reachable.add(target)
                pending.append(target)
    assert reachable == set(steps)


def test_v3_collation_slurm_uses_local_scratch_and_long_qos_for_global_steps():
    repository = Path(__file__).resolve().parents[3]
    script_path = repository / "scripts" / "slurm" / "collate_v3_shards.sbatch"
    readme_path = repository / "scripts" / "slurm" / "README.md"
    if not script_path.is_file() or not readme_path.is_file():
        pytest.skip("Slurm assets are not installed in wheel-only test layouts")

    script = script_path.read_text()
    documentation = readme_path.read_text()

    assert "${SLURM_TMPDIR:-/scratch/${USER}/plinder-collate-" in script
    assert (
        "sbatch \\\n  --qos=6hours \\\n  "
        "--cpus-per-task=8 --mem=16G \\\n  "
        '--output="${OUTPUT_ROOT}/logs/collate-plan-%j.out"'
    ) in documentation
    assert (
        "sbatch \\\n  --qos=6hours \\\n  " "--cpus-per-task=4 --mem=48G"
    ) in documentation


def test_v3_ingest_configs_use_current_schema_and_stages():
    from plinder.data.pipeline.config import get_config

    repository = Path(__file__).resolve().parents[3]
    config_dir = repository / "flows" / "configs" / "v3"
    if not config_dir.is_dir():
        pytest.skip("ingest configs are not installed in wheel-only test layouts")
    ingest_configs = list(config_dir.glob("*.yaml"))
    assert ingest_configs

    for path in ingest_configs:
        cfg = get_config(config_file=path.as_posix(), cached=False)
        assert set(cfg.flow.run_specific_stages) <= set(tasks.STAGES)
        assert cfg.data.plinder_iteration == "v3"
        if path.name == "make_protein_scores.yaml":
            assert "collate_alignments" in cfg.flow.run_specific_stages
            assert "finalize_alignments" in cfg.flow.run_specific_stages
            assert "export_sucos_shape_pocket_qcov" in cfg.flow.run_specific_stages
            assert "finalize_sucos_export" in cfg.flow.run_specific_stages
            assert "map_batch_alignments" in cfg.flow.run_specific_stages
            assert cfg.flow.map_batch_alignments_batch_size == 1
            assert cfg.foldseek.min_seq_id == 0.0
            assert cfg.mmseqs.min_seq_id == 0.0
            assert tasks.STAGES.index("map_batch_alignments") < tasks.STAGES.index(
                "make_batch_scores"
            )
        if path.name == "make_components.yaml":
            assert "collate_partitions" in cfg.flow.run_specific_stages
            assert "make_component_reductions" in cfg.flow.run_specific_stages
            assert "merge_component_reductions" in cfg.flow.run_specific_stages
            assert "make_communities" in cfg.flow.run_specific_stages
            assert "make_directed_set_covers" in cfg.flow.run_specific_stages
            assert "summarize_clusters" in cfg.flow.run_specific_stages
            assert "finalize_index" in cfg.flow.run_specific_stages
        if path.name == "make_entries_ligands.yaml":
            assert "collate_entries" in cfg.flow.run_specific_stages
            assert "finalize_ligand_archives" in cfg.flow.run_specific_stages
            assert "make_ligands" not in cfg.flow.run_specific_stages


def test_make_canonical_ligand_archives_only_archives_asu_sdfs(tmp_path):
    from plinder.data.pipeline.score import finalize_ligand_archives

    canonical = tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files"
    canonical.mkdir(parents=True)
    (canonical / "A.sdf").write_text("canonical")
    (canonical.parent.parent / "1abc.parquet").touch()
    rotated = tmp_path / "raw_entries" / "ab" / "1abc__1__1.B__1.A" / "ligand_files"
    rotated.mkdir(parents=True)
    (rotated / "1.A.sdf").write_text("rotated")

    tasks.make_canonical_ligand_archives(data_dir=tmp_path, two_char_codes=["ab"])

    archive = tmp_path / "ligand_archives" / "ab.parquet"
    packed = pd.read_parquet(archive)
    assert packed[["pdb_id", "ligand_asym_id"]].to_dict("records") == [
        {"pdb_id": "1abc", "ligand_asym_id": "A"}
    ]
    assert packed["sdf"].tolist() == [b"canonical"]
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame({"entry_pdb_id": ["1abc"], "ligand_asym_id": ["A"]}).to_parquet(
        index / "annotation_table.parquet", index=False
    )
    report = finalize_ligand_archives(tmp_path)
    assert report == {
        "status": "complete",
        "shard_count": 1,
        "ligand_count": 1,
        "compressed_bytes": archive.stat().st_size,
    }


def test_make_sub_dbs_loads_normalized_entry_chain_index(tmp_path, monkeypatch):
    from plinder.core.scores import entries as entry_views

    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (tmp_path / "dbs").mkdir()
    pd.DataFrame({"entry_pdb_id": ["1abc"]}).to_parquet(
        index_dir / "annotation_table.parquet",
        index=False,
    )
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["A"],
        }
    ).to_parquet(index_dir / "entry_chains.parquet", index=False)

    sentinel = {"1abc": object()}

    def fake_entry_views(annotation, *, entry_chains):
        assert annotation["entry_pdb_id"].tolist() == ["1abc"]
        assert entry_chains["chain_asym_id"].tolist() == ["A"]
        return sentinel

    observed = {}
    monkeypatch.setattr(entry_views, "entry_views_from_df", fake_entry_views)
    monkeypatch.setattr(tasks.utils, "get_db_sources", lambda **kwargs: {})
    monkeypatch.setattr(
        tasks.databases,
        "make_sub_dbs",
        lambda db_dir, db_sources, entries, **kwargs: observed.update(
            entries=entries,
            kwargs=kwargs,
        ),
    )

    tasks.make_sub_dbs(data_dir=tmp_path, sub_databases=["holo", "apo", "pred"])

    assert observed["entries"] is sentinel
    assert observed["kwargs"] == {
        "identifiers_by_database": None,
        "tmp_dir": None,
        "threads": 1,
    }


def test_make_holo_sub_dbs_selects_only_protein_receptor_chains(tmp_path, monkeypatch):
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (tmp_path / "dbs").mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "2def"],
            "chain_auth_id": ["X", "Y", "Z"],
            "chain_receptor_type": ["protein", "protein", "dna"],
            "chain_is_holo": [True, False, True],
        }
    ).to_parquet(index_dir / "entry_chains.parquet", index=False)

    observed = {}
    lookup_calls = []
    monkeypatch.setattr(tasks.utils, "get_db_sources", lambda **kwargs: {})
    monkeypatch.setattr(
        tasks.databases,
        "make_sub_dbs",
        lambda db_dir, db_sources, entries, **kwargs: observed.update(
            entries=entries,
            kwargs=kwargs,
        ),
    )
    monkeypatch.setattr(
        tasks,
        "make_alignment_chain_lookup",
        lambda **kwargs: lookup_calls.append(kwargs),
    )

    tasks.make_sub_dbs(data_dir=tmp_path, sub_databases=["holo"])

    assert observed["entries"] is None
    assert observed["kwargs"]["identifiers_by_database"] == {
        "holo_foldseek": {"pdb_00001abc_xyz-enrich_X"},
        "holo_mmseqs": {"1abc_X"},
    }
    assert lookup_calls == [{"data_dir": tmp_path, "scratch_dir": None, "threads": 1}]
