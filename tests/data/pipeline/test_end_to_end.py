# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from plinder.data.pipeline import config, pipeline, tasks


def test_end_to_end(mock_alternative_datasets):
    mock_alternative_datasets("4jvm")  # different from 19hc

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
