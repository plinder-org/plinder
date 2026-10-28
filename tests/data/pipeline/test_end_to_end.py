# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from shutil import copy2

import pyarrow.parquet as pq

from plinder.data.pipeline import config, io, pipeline, tasks


def test_end_to_end(mock_alternative_datasets, monkeypatch, test_dir):
    mock_alternative_datasets("4jvm")  # different from 19hc

    def local_archive(*, kind, data_dir, two_char_code):
        assert two_char_code == "9h"
        destination = data_dir / two_char_code
        destination.mkdir(parents=True, exist_ok=True)
        if kind == "cif":
            source = test_dir / "xx/pdb_000019hc/pdb_000019hc_xyz-enrich.cif.gz"
            destination = destination / source.parent.name
            destination.mkdir(exist_ok=True)
            copy2(source, destination / source.name)
        else:
            assert kind == "val"

    monkeypatch.setattr(io, "rsync_rcsb", local_archive)

    # refresh_bundled_ccd re-syncs biotite's bundled CCD from wwPDB via
    # setup_ccd (a heavy online download); the shipped bt_info snapshot is
    # already usable, so stub it out for the offline end-to-end run.
    monkeypatch.setattr(io, "refresh_bundled_ccd", lambda **kwargs: None)

    # make_ccd_ligand_dbs builds the ~50k-component CCD universe once per CCD
    # release, not per ingest; tests/data/test_ccd_ligand_dbs.py covers it.
    stages = ",".join(
        stage
        for stage in tasks.STAGES[: tasks.STAGES.index("collate_partitions")]
        if stage != "make_ccd_ligand_dbs"
    )
    import sys

    print(stages, file=sys.stderr, flush=True)
    conf = {
        "context": {
            "two_char_codes": "9h",
        },
        "flow": {
            "run_specific_stages": stages,
            "make_entries_cpu": 1,
            "protein_clustering_cpu": 1,
            "make_dbs_cpu": 1,
        },
        "source": {
            "pdb_nextgen_root": "ingest",
            "validation_root": "reports",
        },
        "scorer": {
            "sub_databases": "holo",
        },
    }
    cfg = config.get_config(config=conf)
    pipe = pipeline.IngestPipeline(conf=cfg)
    pipe.run()
    annotations = pq.read_table(pipe.plinder_dir / "index/annotation_table.parquet")
    assert annotations.num_rows > 0
    assert set(annotations["entry_pdb_id"].to_pylist()) == {"19hc"}
