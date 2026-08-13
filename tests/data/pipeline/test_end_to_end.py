# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from plinder.data.pipeline import config, io, pipeline, tasks


def test_end_to_end(mock_alternative_datasets, monkeypatch):
    mock_alternative_datasets("4jvm")  # different from 19hc

    # refresh_bundled_ccd re-syncs biotite's bundled CCD from wwPDB via
    # setup_ccd (a heavy online download); the shipped bt_info snapshot is
    # already usable, so stub it out for the offline end-to-end run.
    monkeypatch.setattr(io, "refresh_bundled_ccd", lambda **kwargs: None)

    stages = ",".join(tasks.STAGES[: tasks.STAGES.index("collate_partitions")])
    import sys

    print(stages, file=sys.stderr, flush=True)
    conf = {
        "context": {
            "two_char_codes": "9h",
        },
        "flow": {
            "run_specific_stages": stages,
        },
        "scorer": {
            "sub_databases": "holo",
        },
    }
    cfg = config.get_config(config=conf)
    pipe = pipeline.IngestPipeline(conf=cfg)
    pipe.run()
