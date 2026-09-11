import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest
from plinder.data import protein_clusters as clusters

SEQUENCE = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"


def _write_chains(root):
    path = root / "index/entry_chains.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def", "3ghi", "4jkl", "5mno"],
            "chain_asym_id": ["1.A", "B_2", "C", "D", "E"],
            "chain_receptor_type": ["protein", "protein", "protein", "protein", "dna"],
            "chain_sequence": [SEQUENCE, SEQUENCE, SEQUENCE[:25], None, "ATGC"],
            "chain_is_holo": [True, False, False, False, False],
        }
    ).to_parquet(path, index=False)
    return path


def test_sequence_clusters_native_and_reusable(tmp_path, monkeypatch):
    source = _write_chains(tmp_path)
    source_bytes = source.read_bytes()
    kwargs = dict(data_dir=tmp_path, scratch_dir=tmp_path / "scratch", threads=1)
    output = clusters.make_protein_sequence_clusters(**kwargs)
    table = pq.read_table(output)
    assert table.schema.equals(clusters.SEQUENCE_CLUSTER_SCHEMA)
    frame = table.to_pandas().set_index("entry_pdb_id")
    assert set(frame.index) == {"1abc", "2def", "3ghi", "4jkl"}
    assert frame.loc["1abc", "chain_asym_id"] == "1.A"
    assert frame.loc["2def", "chain_asym_id"] == "B_2"
    assert (
        frame.loc["1abc", "representative_entry_pdb_id"]
        == frame.loc["2def", "representative_entry_pdb_id"]
    )
    assert frame.loc["3ghi", "representative_entry_pdb_id"] == "3ghi"
    assert frame.loc["4jkl", "status"] == "missing_sequence"
    assert pd.isna(frame.loc["4jkl", "is_representative"])
    assert pd.isna(frame.loc["4jkl", "representative_chain_asym_id"])
    assert frame["is_representative"].sum() == 2
    parameters = json.loads(table.schema.metadata[clusters.CLUSTER_METADATA_KEY])
    assert parameters["identity"] == 0.4
    assert parameters["coverage"] == 0.8
    assert parameters["coverage_mode"] == 0
    assert parameters["single_step_clustering"] is True
    assert source.read_bytes() == source_bytes
    assert not list((tmp_path / "scratch").iterdir())
    monkeypatch.setattr(
        clusters, "run", lambda _: pytest.fail("recomputed unchanged clusters")
    )
    assert clusters.make_protein_sequence_clusters(**kwargs) == output


@pytest.mark.parametrize("change", ["identity", "coverage", "sequence", "force"])
def test_sequence_cluster_restart_uses_parameters_and_sequences(
    tmp_path, monkeypatch, change
):
    source = _write_chains(tmp_path)
    commands = []

    def fake_run(command):
        commands.append(command)
        Path(command[3] + "_cluster.tsv").write_text(
            "chain_0\tchain_0\nchain_0\tchain_1\nchain_2\tchain_2\n"
        )

    monkeypatch.setattr(clusters, "run", fake_run)
    kwargs = dict(data_dir=tmp_path, scratch_dir=tmp_path / "scratch", threads=2)
    output = clusters.make_protein_sequence_clusters(**kwargs)
    previous = output.read_bytes()
    if change in {"identity", "coverage"}:
        kwargs[change] = 0.9
    elif change == "force":
        kwargs["force_update"] = True
    else:
        frame = pd.read_parquet(source)
        frame.loc[0, "chain_sequence"] = SEQUENCE[::-1]
        frame.to_parquet(source, index=False)

    def fail_run(command):
        raise RuntimeError("external clustering failed")

    monkeypatch.setattr(clusters, "run", fail_run)
    with pytest.raises(RuntimeError, match="external clustering failed"):
        clusters.make_protein_sequence_clusters(**kwargs)
    assert output.read_bytes() == previous
    assert not list((tmp_path / "scratch").iterdir())
    command = commands[0]
    assert command[:2] == ["mmseqs", "easy-cluster"]
    for option, value in {
        "--dbtype": "1",
        "--min-seq-id": "0.4",
        "-c": "0.8",
        "--cov-mode": "0",
        "--cluster-mode": "0",
        "--single-step-clustering": "1",
        "--threads": "2",
    }.items():
        assert command[command.index(option) + 1] == value


@pytest.mark.parametrize(
    "rows",
    [
        "a\ta\n",  # missing member
        "a\ta\na\tb\na\tb\n",  # duplicate member
        "a\ta\na\tb\na\tc\n",  # extra member
        "c\ta\nc\tb\n",  # unknown representative
        "b\ta\na\tb\n",  # representatives do not belong to themselves
    ],
)
def test_cluster_assignments_reject_incomplete_or_inconsistent_output(tmp_path, rows):
    output = tmp_path / "clusters.tsv"
    output.write_text(rows)
    with pytest.raises(ValueError, match="protein cluster"):
        clusters._read_assignments(output, {"a", "b"})


@pytest.mark.parametrize("protein", [True, False])
def test_empty_or_missing_sequences_keep_schema(tmp_path, monkeypatch, protein):
    source = _write_chains(tmp_path)
    frame = pd.read_parquet(source)
    frame = frame.iloc[[3 if protein else 4]]
    frame.to_parquet(source, index=False)
    monkeypatch.setattr(
        clusters, "run", lambda _: pytest.fail("no sequences to cluster")
    )
    output = clusters.make_protein_sequence_clusters(
        data_dir=tmp_path, scratch_dir=tmp_path / "scratch"
    )
    table = pq.read_table(output)
    assert table.schema.equals(clusters.SEQUENCE_CLUSTER_SCHEMA)
    assert table.num_rows == int(protein)


def test_segmented_protein_config_runs_sequence_clustering(tmp_path, monkeypatch):
    from plinder.data.pipeline import pipeline, tasks

    config = (
        Path(__file__).resolve().parents[2]
        / "flows/configs/ingest/make_protein_scores.yaml"
    )
    if not config.is_file():
        pytest.skip("ingest configs are not installed in wheel-only test layouts")
    pipe = pipeline.IngestPipeline(
        config_file=str(config), config_args=[], cached=False
    )
    pipe.plinder_dir = tmp_path
    calls = []
    monkeypatch.setattr(
        tasks, "make_protein_sequence_clusters", lambda **kw: calls.append(kw)
    )
    pipe.make_protein_sequence_clusters()
    assert len(calls) == 1
    assert calls[0]["data_dir"] == tmp_path
    assert calls[0]["identity"] == 0.4
    assert calls[0]["coverage"] == 0.8
    assert tasks.STAGES.index("collate_entries") < tasks.STAGES.index(
        "make_protein_sequence_clusters"
    )
    pipe.cfg.flow.skip_specific_stages = ["make_protein_sequence_clusters"]
    pipe.make_protein_sequence_clusters()
    assert len(calls) == 1
