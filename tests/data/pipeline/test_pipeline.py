
from plinder.data.pipeline import config, pipeline

def test_pipeline_noop(tmp_path):

    conf = (tmp_path / "test.yaml")
    conf.write_text("""\
ingest:
    run_specific_stages: download_rcsb_files
    skip_specific_stages: download_rcsb_files
""")
    print("in test", conf.as_posix())
    pipe = pipeline.IngestPipeline(config_file=conf.as_posix(), config_args=[], cached=False)
    pipe.run()
