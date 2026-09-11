# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path

import pandas as pd
import pytest
from plinder.data.pipeline import config, pipeline


def test_pipeline_noop(tmp_path):
    conf = tmp_path / "test.yaml"
    conf.write_text(
        """\
flow:
    run_specific_stages: download_rcsb_files
    skip_specific_stages: download_rcsb_files
"""
    )
    pipe = pipeline.IngestPipeline(
        config_file=conf.as_posix(), config_args=[], cached=False
    )
    pipe.run()


def test_pipeline_clusters_whole_interfaces_only():
    cfg = config.get_config(cached=False)
    pipe = pipeline.IngestPipeline(conf=cfg)

    assert ("interface", ["interface_qcov"]) in pipe._cluster_entities()


@pytest.mark.parametrize("mode", ["enabled", "metric_removed", "stage_skipped"])
def test_segmented_ligand_ingest_runs_mhfp6_unless_disabled(tmp_path, mode):
    from plinder.data.annotations import get_similarity_scores as scoring

    config_path = (
        Path(__file__).resolve().parents[3]
        / "flows/configs/ingest/make_entries_ligands.yaml"
    )
    if not config_path.is_file():
        pytest.skip("ingest configs are not installed in wheel-only test layouts")
    pipe = pipeline.IngestPipeline(config_file=str(config_path), cached=False)
    pipe.plinder_dir = tmp_path
    if mode == "metric_removed":
        pipe.cfg.flow.cluster_metrics = [
            metric
            for metric in pipe.cfg.flow.cluster_metrics
            if metric != scoring.MHFP6_METRIC
        ]
    elif mode == "stage_skipped":
        pipe.cfg.flow.skip_specific_stages = ["make_mhfp6_scores"]
    else:
        fingerprints = tmp_path / "fingerprints"
        fingerprints.mkdir()
        scoring.write_mhfp6_fingerprints(
            pd.DataFrame(
                {
                    "ligand_smiles_id": [0, 1],
                    "ligand_rdkit_canonical_smiles": ["CCO", "c1ccccc1"],
                }
            ),
            fingerprints,
        )

    chunks = pipe.scatter_make_mhfp6_scores()
    for chunk in chunks:
        pipe.make_mhfp6_scores(chunk)
    outputs = list((tmp_path / scoring.MHFP6_SCORES_DIR).glob("*.parquet"))
    if mode != "enabled":
        assert chunks == [[]]
        assert outputs == []
        return
    assert chunks == [[0, 1]]
    assert len(outputs) == 1
    scores = pd.read_parquet(outputs[0])
    self_scores = scores[scores["query_ligand_id"] == scores["target_ligand_id"]]
    assert set(self_scores["query_ligand_id"]) == {0, 1}
    assert self_scores[scoring.MHFP6_METRIC].tolist() == [100.0, 100.0]
