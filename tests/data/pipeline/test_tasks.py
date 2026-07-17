# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import ast
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import pytest
from plinder.data.pipeline import io, tasks


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
        "download_components_cif",
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


def test_final_structure_qc_is_not_a_pipeline_stage():
    assert "structure_qc" not in tasks.STAGES


def test_scoring_finalization_stage_order_and_partitions():
    assert tasks.STAGES.index("make_batch_scores") < tasks.STAGES.index(
        "collate_alignments"
    )
    assert tasks.STAGES.index("collate_alignments") < tasks.STAGES.index(
        "collate_partitions"
    )
    assert tasks.STAGES.index("collate_partitions") < tasks.STAGES.index(
        "make_components_and_communities"
    )
    assert tasks.STAGES.index("make_components_and_communities") < tasks.STAGES.index(
        "finalize_index"
    )
    assert tasks.STAGES.index("finalize_index") < tasks.STAGES.index("make_mmp_index")
    partitions = tasks.scatter_collate_partitions()
    assert len(partitions) == 38
    assert ["0"] in partitions
    assert ["z"] in partitions
    assert ["apo"] in partitions
    assert ["pred"] in partitions


def test_collate_alignments_writes_query_addressable_shards(tmp_path):
    columns = {
        "query_entry": ["1abc", "2abd", "3xyz"],
        "target_entry": ["9zzz", "8yyy", "7xxx"],
        "query_chain_mapped": ["A", "B", "C"],
        "target_chain_mapped": ["D", "E", "F"],
        "source": ["foldseek"] * 3,
        "similarity": [0.8, 0.7, 0.6],
    }
    mapped_dir = tmp_path / "dbs/subdbs/holo_foldseek/mapped_aln"
    mapped_dir.mkdir(parents=True)
    for index, pdb_id in enumerate(columns["query_entry"]):
        pd.DataFrame(
            {column: [values[index]] for column, values in columns.items()}
        ).to_parquet(mapped_dir / f"{pdb_id}.parquet", index=False)

    assert tasks.scatter_collate_alignments(data_dir=tmp_path) == [["ab"], ["xy"]]
    tasks.collate_alignments(data_dir=tmp_path, partition=["ab"])

    shard = pd.read_parquet(
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet"
    )
    assert shard["query_entry"].tolist() == ["1abc", "2abd"]
    assert "3xyz" not in set(shard["query_entry"])


def test_empty_alignment_scatter_has_noop_branch(tmp_path):
    assert tasks.scatter_collate_alignments(data_dir=tmp_path) == [[]]
    tasks.collate_alignments(data_dir=tmp_path, partition=[])


def test_ligand_cluster_scatter_requires_both_node_levels(tmp_path):
    metric = "sucos_shape_pocket_qcov"
    expected = [[(metric, 50)]]
    assert (
        tasks.scatter_make_components_and_communities(
            data_dir=tmp_path,
            metrics=[metric],
            thresholds=[50],
            stop_on_cluster=0,
            skip_existing_clusters=True,
        )
        == expected
    )

    for root in ["clusters", "ligand_clusters"]:
        for cluster, directed in [
            ("components", True),
            ("components", False),
            ("communities", False),
        ]:
            path = (
                tmp_path
                / root
                / f"cluster={cluster}"
                / f"directed={directed}"
                / f"metric={metric}"
                / "threshold=50.parquet"
            )
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        if root == "clusters":
            assert (
                tasks.scatter_make_components_and_communities(
                    data_dir=tmp_path,
                    metrics=[metric],
                    thresholds=[50],
                    stop_on_cluster=0,
                    skip_existing_clusters=True,
                )
                == expected
            )

    assert tasks.scatter_make_components_and_communities(
        data_dir=tmp_path,
        metrics=[metric],
        thresholds=[50],
        stop_on_cluster=0,
        skip_existing_clusters=True,
    ) == [[]]


def test_empty_cached_cluster_chunk_is_a_noop(tmp_path, monkeypatch):
    def fail_if_called(**kwargs):
        pytest.fail(f"cluster computation should not run: {kwargs}")

    monkeypatch.setattr(
        tasks.clusters,
        "make_components_and_communities",
        fail_if_called,
    )

    tasks.make_components_and_communities(
        data_dir=tmp_path,
        metric_threshold=[],
        skip_existing_clusters=True,
    )


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
    assert "self.next(self.scatter_make_canonical_ligand_archives)" in flow
    assert "def make_canonical_ligand_archives" in flow
    assert "self.next(self.scatter_collate_partitions)" in flow
    assert "self.next(self.scatter_collate_alignments)" in flow
    assert "self.pipeline.collate_alignments(self.input)" in flow
    assert "self.next(self.scatter_make_components_and_communities)" in flow
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

    reachable = {"start"}
    pending = ["start"]
    while pending:
        for target in edges[pending.pop()]:
            assert target in steps
            if target not in reachable:
                reachable.add(target)
                pending.append(target)
    assert reachable == set(steps)


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


def test_make_canonical_ligand_archives_only_archives_asu_sdfs(tmp_path):
    canonical = tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files"
    canonical.mkdir(parents=True)
    (canonical / "A.sdf").write_text("canonical")
    (canonical.parent.parent / "1abc.parquet").touch()
    rotated = tmp_path / "raw_entries" / "ab" / "1abc__1__1.B__1.A" / "ligand_files"
    rotated.mkdir(parents=True)
    (rotated / "1.A.sdf").write_text("rotated")

    tasks.make_canonical_ligand_archives(data_dir=tmp_path, two_char_codes=["ab"])

    archive = tmp_path / "ligand_archives" / "ab.zip"
    with ZipFile(archive) as handle:
        assert handle.namelist() == ["1abc/ligand_files/A.sdf"]


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
        lambda db_dir, db_sources, entries: observed.update(entries=entries),
    )

    tasks.make_sub_dbs(data_dir=tmp_path, sub_databases=["holo", "apo", "pred"])

    assert observed["entries"] is sentinel
