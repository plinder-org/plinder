# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from plinder.data.pipeline import pipeline


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
