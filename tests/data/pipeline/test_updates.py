# Copyright (c) 2026, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

import gzip
import json
import sys

import pandas as pd
import pytest
from plinder.data.pipeline import updates
from plinder.data.pipeline.ingest import load_manifest


def _json_gz(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        json.dump(value, handle)


def _cif(root, pdb_id, revision):
    path = (
        root
        / "data/entries/divided"
        / pdb_id[1:3]
        / f"pdb_0000{pdb_id}"
        / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(
            f"data_{pdb_id}\nloop_\n"
            "_pdbx_audit_revision_history.ordinal\n"
            "_pdbx_audit_revision_history.data_content_type\n"
            "_pdbx_audit_revision_history.major_revision\n"
            "_pdbx_audit_revision_history.minor_revision\n"
            f"1 'Structure model' {revision[0]} {revision[1]}\n"
            "2 'Structure factors' 99 99\n#\n"
        )
    return path


def _sources(rows):
    return pd.DataFrame(rows, columns=["entry_pdb_id", *updates.REVISION_COLUMNS])


@pytest.fixture
def inputs(tmp_path):
    release = tmp_path / "release"
    (release / "index").mkdir(parents=True)
    _sources(
        [
            ("1abc", 1, 0),
            ("2abc", 1, 0),
            ("3abc", 1, 0),
            ("4abc", 1, 0),
        ]
    ).to_parquet(release / "index/entry_sources.parquet", index=False)
    root = tmp_path / "nextgen"
    revisions = {"1abc": (1, 0), "2abc": (1, 1), "3abc": (2, 0), "5abc": (1, 0)}
    _json_gz(
        root / "holdings/current_file_holdings.json.gz",
        {
            f"pdb_0000{key}": {"mmcif": [f"/pdb_nextgen/{key}.cif.gz"]}
            for key in revisions
        },
    )
    _json_gz(
        root / "holdings/released_structures_last_modified_dates.json.gz",
        {f"pdb_0000{key}": "2026-06-01T00:00:00+00:00" for key in revisions},
    )
    for key, revision in revisions.items():
        _cif(root, key, revision)
    obsolete = tmp_path / "obsolete.dat"
    obsolete.write_text(
        " LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS\n"
        "OBSLTE    01-JUN-26 4ABC     5ABC\n"
    )
    return {"data_dir": release, "nextgen_root": root, "obsolete_path": obsolete}


def _bytes(root):
    return {
        str(path.relative_to(root)): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file()
    }


def test_plan_classifies_entries_and_both_search_directions(inputs, tmp_path):
    before = _bytes(inputs["data_dir"])
    output = tmp_path / "plan"
    plan = updates.plan_update(**inputs, output_dir=output, threads=2)
    assert plan["status"] == "ready"
    assert plan["counts"] == {"added": 1, "revised": 2, "obsolete": 1, "unchanged": 1}
    assert load_manifest(output / "ingest.txt") == ["2abc", "3abc", "5abc"]
    assert load_manifest(output / "remove.txt") == ["4abc"]
    assert load_manifest(output / "invalidate.txt") == ["2abc", "3abc", "4abc"]
    assert load_manifest(output / "active.txt") == ["1abc", "2abc", "3abc", "5abc"]
    assert load_manifest(output / "unchanged.txt") == ["1abc"]
    pairs = set()
    for work in plan["work"]["searches"]:
        pairs.update(
            (q, t)
            for q in load_manifest(output / work["queries"])
            for t in load_manifest(output / work["targets"])
        )
    active = set(load_manifest(output / "active.txt"))
    changed = set(load_manifest(output / "ingest.txt"))
    assert pairs == {
        (q, t) for q in active for t in active if q in changed or t in changed
    }
    report = pd.read_parquet(output / "entries.parquet").set_index("pdb_id")
    assert list(report.loc["4abc", "replacement_pdb_ids"]) == ["5abc"]
    assert pd.read_csv(output / "entries.tsv", sep="\t").shape[0] == 5
    assert json.loads((output / "plan.json").read_text()) == plan
    assert before == _bytes(inputs["data_dir"])


def test_obsolete_preserves_multiple_successors_and_withdrawals(tmp_path):
    path = tmp_path / "obsolete.dat"
    path.write_text(
        " LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS\n"
        "OBSLTE    01-JUN-26 1ABC     2ABC 3ABC\n"
        "OBSLTE    01-JUN-26 1ABC     4ABC\n"
        "OBSLTE    01-JUN-26 5ABC\n"
    )
    assert updates.read_obsolete(path) == {
        "1abc": ["2abc", "3abc", "4abc"],
        "5abc": [],
    }


@pytest.mark.parametrize("line", ["OBSLTE", "OBSLTE 01-JUN-26", "OBSLTE 01-JUN-26 BAD"])
def test_malformed_obsolete_record_fails(tmp_path, line):
    path = tmp_path / "obsolete.dat"
    path.write_text(line)
    with pytest.raises(ValueError):
        updates.read_obsolete(path)


@pytest.mark.parametrize(
    "problem", ["absent", "conflict", "older", "missing_cif", "missing_revision"]
)
def test_source_problems_block_plan_without_removal_lists(inputs, tmp_path, problem):
    if problem == "absent":
        inputs["obsolete_path"].write_text("")
    elif problem == "conflict":
        with inputs["obsolete_path"].open("a") as handle:
            handle.write("OBSLTE    01-JUN-26 1ABC\n")
    elif problem == "older":
        _cif(inputs["nextgen_root"], "1abc", (0, 9))
    elif problem == "missing_cif":
        _cif(inputs["nextgen_root"], "1abc", (1, 0)).unlink()
    else:
        path = _cif(inputs["nextgen_root"], "1abc", (1, 0))
        with gzip.open(path, "wt") as handle:
            handle.write("data_1abc\n_entry.id 1abc\n")
    output = tmp_path / "blocked"
    plan = updates.plan_update(**inputs, output_dir=output)
    assert plan["status"] == "blocked"
    assert plan["counts"]["blocked"] == 1
    assert not list(output.glob("*.txt"))
    assert (
        pd.read_parquet(output / "entries.parquet")
        .query("action == 'blocked'")
        .reason.str.len()
        .gt(0)
        .all()
    )


def test_null_release_revision_blocks(inputs, tmp_path):
    path = inputs["data_dir"] / "index/entry_sources.parquet"
    sources = pd.read_parquet(path)
    sources.loc[0, updates.REVISION_COLUMNS[0]] = None
    sources.to_parquet(path, index=False)
    plan = updates.plan_update(**inputs, output_dir=tmp_path / "plan")
    assert plan["status"] == "blocked"


def test_snapshot_cache_reuses_only_unchanged_successful_reads(
    inputs, tmp_path, monkeypatch
):
    root = inputs["nextgen_root"]
    snapshot = updates.read_snapshot(root)
    cache = tmp_path / "snapshot.parquet"
    snapshot.loc[snapshot.entry_pdb_id == "2abc", "error"] = "previous read failed"
    snapshot.loc[snapshot.entry_pdb_id == "3abc", "last_modified"] = "older timestamp"
    snapshot.to_parquet(cache, index=False)
    original = updates._read_revision
    calls = []

    def read(path):
        calls.append(path.parent.name)
        return original(path)

    monkeypatch.setattr(updates, "_read_revision", read)
    new = updates.read_snapshot(root, previous_snapshot=cache, threads=1)
    assert calls == ["pdb_00002abc", "pdb_00003abc"]
    assert new.error.isna().all()


def test_inventory_cache_does_not_mark_unapplied_changes_complete(inputs, tmp_path):
    first = tmp_path / "first"
    plan = updates.plan_update(**inputs, output_dir=first)
    second = tmp_path / "second"
    repeated = updates.plan_update(
        **inputs, output_dir=second, previous_snapshot=first / "snapshot.parquet"
    )
    assert repeated == plan
    pd.testing.assert_frame_equal(
        pd.read_parquet(first / "entries.parquet"),
        pd.read_parquet(second / "entries.parquet"),
    )


def test_no_changes_after_release_matches_snapshot(inputs, tmp_path):
    snapshot = updates.read_snapshot(inputs["nextgen_root"])
    snapshot[["entry_pdb_id", *updates.REVISION_COLUMNS]].to_parquet(
        inputs["data_dir"] / "index/entry_sources.parquet",
        index=False,
    )
    plan = updates.plan_update(**inputs, output_dir=tmp_path / "plan")
    assert plan["status"] == "no_changes"
    assert plan["counts"] == {"unchanged": 4}
    assert plan["work"]["refresh"] == []
    assert load_manifest(tmp_path / "plan/ingest.txt") == []


def test_removal_only_plan(inputs, tmp_path):
    plan = updates.plan_update(**inputs, output_dir=tmp_path / "plan", pdb_ids=["4ABC"])
    assert plan["status"] == "ready"
    assert plan["scope"] == ["4abc"]
    assert plan["counts"] == {"obsolete": 1}
    assert load_manifest(tmp_path / "plan/ingest.txt") == []
    assert load_manifest(tmp_path / "plan/remove.txt") == ["4abc"]


def test_duplicate_source_ids_are_rejected(inputs, tmp_path):
    path = inputs["data_dir"] / "index/entry_sources.parquet"
    sources = pd.read_parquet(path)
    pd.concat([sources, sources.iloc[:1]]).to_parquet(path, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        updates.plan_update(**inputs, output_dir=tmp_path / "plan")
    assert not (tmp_path / "plan").exists()


def test_plan_pins_resolved_snapshot_root(inputs, tmp_path):
    alias = tmp_path / "latest"
    alias.symlink_to(inputs["nextgen_root"], target_is_directory=True)
    inputs["nextgen_root"] = alias
    plan = updates.plan_update(**inputs, output_dir=tmp_path / "plan")
    assert plan["nextgen_root"] == str(alias.resolve())
    assert str(alias) not in plan["inputs"]["holdings"]["path"]


def test_do_not_overwrite_release_or_existing_plan(inputs, tmp_path):
    with pytest.raises(ValueError, match="outside"):
        updates.plan_update(**inputs, output_dir=inputs["data_dir"] / "plan")
    output = tmp_path / "plan"
    updates.plan_update(**inputs, output_dir=output)
    before = _bytes(output)
    with pytest.raises(FileExistsError):
        updates.plan_update(**inputs, output_dir=output)
    assert before == _bytes(output)


def test_cli_blocked_exit_and_report(inputs, tmp_path, monkeypatch, capsys):
    inputs["obsolete_path"].write_text("")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "updates",
            str(inputs["data_dir"]),
            str(tmp_path / "plan"),
            "--nextgen-root",
            str(inputs["nextgen_root"]),
            "--obsolete",
            str(inputs["obsolete_path"]),
            "--threads",
            "1",
        ],
    )
    with pytest.raises(SystemExit) as exc:
        updates.main()
    assert exc.value.code == 1
    assert json.loads(capsys.readouterr().out)["status"] == "blocked"


def test_changed_catalogue_during_read_does_not_write_plan(
    inputs, tmp_path, monkeypatch
):
    original = updates.read_snapshot

    def read(*args, **kwargs):
        frame = original(*args, **kwargs)
        inputs["obsolete_path"].write_text("")
        return frame

    monkeypatch.setattr(updates, "read_snapshot", read)
    with pytest.raises(RuntimeError, match="catalogues changed"):
        updates.plan_update(**inputs, output_dir=tmp_path / "plan")
    assert not (tmp_path / "plan").exists()


def test_missing_timestamp_always_reads_revision(inputs, tmp_path, monkeypatch):
    root = inputs["nextgen_root"]
    cached = updates.read_snapshot(root)
    path = tmp_path / "cache.parquet"
    cached.to_parquet(path, index=False)
    _json_gz(root / "holdings/released_structures_last_modified_dates.json.gz", {})
    original = updates._read_revision
    calls = []

    def read(path):
        calls.append(path)
        return original(path)

    monkeypatch.setattr(updates, "_read_revision", read)
    snapshot = updates.read_snapshot(root, previous_snapshot=path)
    assert len(calls) == 4
    assert snapshot.last_modified.isna().all()


def test_empty_scope_is_not_interpreted_as_all_entries(inputs, tmp_path):
    plan = updates.plan_update(**inputs, output_dir=tmp_path / "empty", pdb_ids=[])
    assert plan["status"] == "no_changes"
    assert plan["counts"] == {}
    assert pd.read_parquet(tmp_path / "empty/entries.parquet").empty
