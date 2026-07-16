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


def test_final_structure_qc_is_not_a_pipeline_stage():
    assert "structure_qc" not in tasks.STAGES


def test_metaflow_graph_uses_canonical_ligand_archive_stage():
    repository = Path(__file__).resolve().parents[3]
    flow = (repository / "flows" / "data_ingest.py").read_text()

    assert "def scatter_structure_qc" not in flow
    assert "self.pipeline.structure_qc" not in flow
    assert "self.next(self.scatter_make_entries)" in flow
    assert "def join_make_entries" in flow
    assert "self.next(self.scatter_make_canonical_ligand_archives)" in flow
    assert "def make_canonical_ligand_archives" in flow

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
    ingest_configs = list(config_dir.glob("*.yaml"))
    assert ingest_configs

    for path in ingest_configs:
        cfg = get_config(config_file=path.as_posix(), cached=False)
        assert set(cfg.flow.run_specific_stages) <= set(tasks.STAGES)
        assert cfg.data.plinder_iteration == "v3"


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
