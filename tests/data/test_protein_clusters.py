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


def _write_structure_sources(root):
    import gzip

    import biotite.structure as struc
    import numpy as np
    from biotite.structure.io import pdbx
    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )
    from plinder.data.pipeline.ingest import resolve_entry_paths

    fixture = (
        Path(__file__).resolve().parents[1]
        / "test_data/system_instance_dataframe/plinder_final_dir_structure/apo/7OS1__1__1.A.cif"
    )
    atoms = get_structure_with_altloc(read_mmcif_file(fixture))
    atoms = atoms[atoms.chain_id == atoms.chain_id[0]]
    starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    atoms = atoms[: starts[100]].copy()
    _write_chains(root)
    source_root = root / "cifs"
    paths = []
    for pdb_id, asym_id in [
        ("1abc", "1.A"),
        ("2def", "B_2"),
        ("3ghi", "C"),
        ("4jkl", "D"),
    ]:
        selected = atoms.copy()
        if pdb_id == "3ghi":
            selected = selected[: starts[25]]
        elif pdb_id == "4jkl":
            selected = selected[selected.atom_name == "N"]
        selected.chain_id[:] = asym_id
        if pdb_id == "1abc":
            extra = atoms.copy()
            extra.chain_id[:] = "Z"
            selected = selected + extra
            second_model = selected.copy()
            second_model.coord[:, 0] *= 5
            selected = struc.stack([selected, second_model])
        cif = pdbx.CIFFile()
        pdbx.set_structure(cif, selected)
        # Exercise label-asym selection independently of author chain names.
        cif.block["atom_site"]["auth_asym_id"] = pdbx.CIFColumn(
            np.full(cif.block["atom_site"].row_count, "author_chain")
        )
        path, _ = resolve_entry_paths(
            pdb_id, cif_root=source_root, validation_root=source_root
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as handle:
            cif.write(handle)
        paths.append(path)
    return source_root, paths, atoms


def test_structure_input_keeps_full_first_model_and_selected_asym_only(tmp_path):
    import numpy as np
    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )

    _, paths, original = _write_structure_sources(tmp_path)
    input_dir = tmp_path / "prepared"
    input_dir.mkdir()
    chains = pd.DataFrame({"chain_asym_id": ["1.A"], "member": ["chain_0_A"]})
    assert clusters._prepare_structure_entry(paths[0], chains, input_dir) == {
        "chain_0_A"
    }
    assert [path.name for path in input_dir.iterdir()] == ["chain_0.cif"]
    prepared = get_structure_with_altloc(read_mmcif_file(input_dir / "chain_0.cif"))
    assert set(prepared.chain_id) == {"A"}
    assert len(prepared) == len(original)
    np.testing.assert_array_equal(prepared.atom_name, original.atom_name)
    np.testing.assert_allclose(prepared.coord, original.coord, atol=1e-4)


def test_structure_clusters_native_and_reusable(tmp_path, monkeypatch):
    source_root, _, _ = _write_structure_sources(tmp_path)
    kwargs = dict(
        data_dir=tmp_path,
        cif_root=source_root,
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )
    output = clusters.make_protein_structure_clusters(**kwargs)
    assert output == tmp_path / "protein_clusters/structure.parquet"
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
    assert frame.loc["4jkl", "status"] == "insufficient_coordinates"
    assert pd.isna(frame.loc["4jkl", "representative_chain_asym_id"])
    assert pd.isna(frame.loc["4jkl", "is_representative"])
    assert frame["is_representative"].sum() == 2
    parameters = json.loads(table.schema.metadata[clusters.CLUSTER_METADATA_KEY])
    assert parameters["lddt"] == 0.7
    assert parameters["coverage"] == 0.8
    assert parameters["coverage_mode"] == 0
    assert parameters["single_step_clustering"] is True
    assert not list((tmp_path / "scratch").iterdir())
    monkeypatch.setattr(
        clusters,
        "_prepare_structure_entry",
        lambda *_: pytest.fail("reread unchanged CIFs"),
    )
    monkeypatch.setattr(
        clusters, "run", lambda _: pytest.fail("recomputed unchanged clusters")
    )
    assert clusters.make_protein_structure_clusters(**kwargs) == output


@pytest.mark.parametrize("change", ["lddt", "coverage", "source", "force"])
def test_structure_cluster_restart_checks_parameters_and_sources(
    tmp_path, monkeypatch, change
):
    import os

    source_root, paths, _ = _write_structure_sources(tmp_path)
    calls = []

    def prepare(path, frame, input_dir):
        return set(frame.loc[frame["entry_pdb_id"] != "4jkl", "member"])

    def fake_run(command):
        calls.append(command)
        Path(command[3] + "_cluster.tsv").write_text(
            "chain_0_A\tchain_0_A\nchain_0_A\tchain_1_A\nchain_2_A\tchain_2_A\n"
        )

    monkeypatch.setattr(clusters, "_prepare_structure_entry", prepare)
    monkeypatch.setattr(clusters, "run", fake_run)
    kwargs = dict(
        data_dir=tmp_path,
        cif_root=source_root,
        scratch_dir=tmp_path / "scratch",
        threads=2,
    )
    output = clusters.make_protein_structure_clusters(**kwargs)
    before = output.read_bytes()
    if change in {"lddt", "coverage"}:
        kwargs[change] = 0.9
    elif change == "force":
        kwargs["force_update"] = True
    else:
        stat = paths[0].stat()
        os.utime(paths[0], ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))

    def fail_run(command):
        raise RuntimeError("external clustering failed")

    monkeypatch.setattr(clusters, "run", fail_run)
    with pytest.raises(RuntimeError, match="external clustering failed"):
        clusters.make_protein_structure_clusters(**kwargs)
    assert output.read_bytes() == before
    assert not list((tmp_path / "scratch").iterdir())
    assert calls[0][:2] == ["foldseek", "easy-cluster"]
    for option, value in {
        "--lddt-threshold": "0.7",
        "-c": "0.8",
        "--cov-mode": "0",
        "--cluster-mode": "0",
        "--single-step-clustering": "1",
        "--alignment-type": "2",
        "--threads": "2",
    }.items():
        assert calls[0][calls[0].index(option) + 1] == value


@pytest.mark.parametrize("protein", [True, False])
def test_structure_empty_or_unresolved_input(tmp_path, monkeypatch, protein):
    source_root, _, _ = _write_structure_sources(tmp_path)
    source = tmp_path / "index/entry_chains.parquet"
    frame = pd.read_parquet(source).iloc[[3 if protein else 4]]
    frame.to_parquet(source, index=False)
    monkeypatch.setattr(
        clusters, "run", lambda _: pytest.fail("no coordinates to cluster")
    )
    output = clusters.make_protein_structure_clusters(
        data_dir=tmp_path, cif_root=source_root, scratch_dir=tmp_path / "scratch"
    )
    assert pq.read_table(output).num_rows == int(protein)


def test_structure_missing_source_stops_stage(tmp_path):
    source_root, paths, _ = _write_structure_sources(tmp_path)
    paths[0].unlink()
    with pytest.raises(FileNotFoundError):
        clusters.make_protein_structure_clusters(
            data_dir=tmp_path, cif_root=source_root, scratch_dir=tmp_path / "scratch"
        )
    assert not (tmp_path / "protein_clusters/structure.parquet").exists()


def test_segmented_protein_config_runs_structure_clustering(tmp_path, monkeypatch):
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
    pipe.cfg.source.pdb_nextgen_root = str(tmp_path / "source")
    calls = []
    monkeypatch.setattr(
        tasks, "make_protein_structure_clusters", lambda **kw: calls.append(kw)
    )
    pipe.make_protein_structure_clusters()
    assert len(calls) == 1
    assert calls[0]["cif_root"] == tmp_path / "source"
    assert calls[0]["lddt"] == 0.7
    assert calls[0]["coverage"] == 0.8
    pipe.cfg.flow.skip_specific_stages = ["make_protein_structure_clusters"]
    pipe.make_protein_structure_clusters()
    assert len(calls) == 1
