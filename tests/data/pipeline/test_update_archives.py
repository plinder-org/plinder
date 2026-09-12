import json
from pathlib import Path
from shutil import rmtree

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from plinder.core.utils.files import file_sha256, write_json_atomic
from plinder.data.pipeline import collate, tasks, update_archives
from plinder.data.pipeline.score import finalize_ligand_archives
from plinder.data.pipeline.update_entries import _file_stats


def _write_poses(root, poses):
    for pdb_id, asym, sdf in poses:
        entry = root / "raw_entries" / pdb_id[1:3] / pdb_id
        (entry / "ligand_files").mkdir(parents=True, exist_ok=True)
        (entry.parent / f"{pdb_id}.parquet").touch()
        (entry / "ligand_files" / f"{asym}.sdf").write_bytes(sdf)


def _write_tables(root, poses):
    for name, path in collate.entry_table_paths(root).items():
        path.parent.mkdir(parents=True, exist_ok=True)
        table = pa.table({"entry_pdb_id": [pdb_id for pdb_id, _, _ in poses]})
        if name == "annotation":
            table = table.append_column(
                "ligand_asym_id",
                pa.array([asym for _, asym, _ in poses], type=pa.string()),
            )
        pq.write_table(table, path)
    write_json_atomic(root / "index/collation.json", {"status": "complete"})


@pytest.fixture(params=[True, False], ids=["with_poses", "ligand_free"])
def archive_update_case(tmp_path, request):
    base, workspace = tmp_path / "base", tmp_path / "updated"
    original = [
        (pdb_id, "L", f"old {pdb_id}".encode())
        for pdb_id in ("1abc", "2abd", "3def", "4ghi")
    ]
    _write_tables(base, original)
    _write_poses(base, original)
    tasks.make_canonical_ligand_archives(
        data_dir=base, two_char_codes=["ab", "de", "gh"]
    )
    finalize_ligand_archives(base)
    added = (
        [("1abc", "M", b"revised pose"), ("5jkl", "N", b"added pose")]
        if request.param
        else []
    )
    expected = [row for row in original if row[0] in {"2abd", "4ghi"}] + added
    _write_tables(workspace, expected)
    _write_poses(workspace / ".incoming", added)
    prepared = {
        **collate.entry_table_paths(workspace),
        "collation": workspace / "index/collation.json",
    }
    write_json_atomic(
        workspace / "entry_update.json",
        {
            "status": collate.REPAIR_REQUIRED_STATUS,
            "base_release": str(base),
            "inputs": {"base_tables": _file_stats(collate.entry_table_paths(base))},
            "outputs": _file_stats(prepared),
            "ingested_pdb_ids": ["1abc", "5jkl"],
            "obsolete_pdb_ids": ["3def"],
        },
    )
    # Incremental archive updates must work without the old per-entry files.
    rmtree(base / "raw_entries")
    return base, workspace, expected


def _archive_bytes(root):
    return {
        path.name: path.read_bytes() for path in (root / "ligand_archives").iterdir()
    }


def test_archive_update_matches_full_repack_and_reuses_unchanged_shards(
    archive_update_case, tmp_path
):
    base, workspace, expected = archive_update_case
    before = _archive_bytes(base)
    report = update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    assert report["status"] == "complete"
    assert report["updated_pdb_ids"] == ["1abc", "3def", "5jkl"]
    assert report["ligand_count"] == len(expected)
    assert _archive_bytes(base) == before
    assert not (workspace / "ligand_archives/de.parquet").exists()
    assert (workspace / "ligand_archives/gh.parquet").stat().st_ino == (
        base / "ligand_archives/gh.parquet"
    ).stat().st_ino
    assert (workspace / "ligand_archives/ab.parquet").stat().st_ino != (
        base / "ligand_archives/ab.parquet"
    ).stat().st_ino
    full = tmp_path / "full"
    _write_poses(full, expected)
    tasks.make_canonical_ligand_archives(
        data_dir=full, two_char_codes=sorted({row[0][1:3] for row in expected})
    )
    for path in (full / "ligand_archives").glob("*.parquet"):
        actual = pq.read_table(workspace / "ligand_archives" / path.name)
        assert actual.to_pylist() == pq.read_table(path).to_pylist()
    assert (
        update_archives.update_ligand_archives(workspace, memory_limit="1GB") == report
    )
    assert _archive_bytes(base) == before
    assert (
        json.loads((workspace / "entry_update.json").read_text())["status"]
        == collate.REPAIR_REQUIRED_STATUS
    )


def test_failed_archive_install_preserves_previous_output_and_can_retry(
    archive_update_case, monkeypatch
):
    base, workspace, _ = archive_update_case
    update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    before, base_before = _archive_bytes(workspace), _archive_bytes(base)
    rename = Path.rename

    def fail_install(path, target):
        if path.name == ".ligand_archives.pending":
            raise OSError("install failed")
        return rename(path, target)

    with monkeypatch.context() as patch:
        patch.setattr(Path, "rename", fail_install)
        with pytest.raises(OSError, match="install failed"):
            update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    assert _archive_bytes(workspace) == before
    assert _archive_bytes(base) == base_before
    update_archives.update_ligand_archives(workspace, memory_limit="1GB")


def test_archive_update_recovers_interrupted_directory_swap(archive_update_case):
    _, workspace, _ = archive_update_case
    first = update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    (workspace / "ligand_archives").rename(workspace / ".ligand_archives.previous")
    assert (
        update_archives.update_ligand_archives(workspace, memory_limit="1GB") == first
    )
    assert not (workspace / ".ligand_archives.previous").exists()


def test_missing_new_pose_leaves_previous_archives_untouched(archive_update_case):
    base, workspace, expected = archive_update_case
    if not any(row[0] == "1abc" for row in expected):
        return
    update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    before, base_before = _archive_bytes(workspace), _archive_bytes(base)
    (workspace / ".incoming/raw_entries/ab/1abc/ligand_files/M.sdf").unlink()
    with pytest.raises(ValueError, match="missing_index_ligands=1"):
        update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    assert _archive_bytes(workspace) == before
    assert _archive_bytes(base) == base_before


def test_archive_update_rejects_changed_entry_tables(archive_update_case):
    _, workspace, _ = archive_update_case
    path = workspace / "index/entry_metadata.parquet"
    path.write_bytes(b"changed")
    with pytest.raises(ValueError, match="entry tables changed"):
        update_archives.update_ligand_archives(workspace)


def test_archive_update_rejects_incomplete_base_archives(archive_update_case):
    base, workspace, _ = archive_update_case
    (base / "ligand_archives/ab.parquet").unlink()
    with pytest.raises(ValueError, match="completion manifest"):
        update_archives.update_ligand_archives(workspace)


def test_archive_update_detects_input_changes_during_build(
    archive_update_case, monkeypatch
):
    _, workspace, _ = archive_update_case
    update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    before = _archive_bytes(workspace)
    validate = update_archives.finalize_ligand_archives

    def change_inputs(*args, **kwargs):
        result = validate(*args, **kwargs)
        write_json_atomic(workspace / "index/collation.json", {"status": "changed"})
        return result

    monkeypatch.setattr(update_archives, "finalize_ligand_archives", change_inputs)
    with pytest.raises(RuntimeError, match="inputs changed"):
        update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    assert _archive_bytes(workspace) == before


def test_archive_update_can_remove_every_ligand(archive_update_case):
    base, workspace, _ = archive_update_case
    rmtree(workspace / ".incoming", ignore_errors=True)
    annotation = workspace / "index/annotation_table.parquet"
    pq.write_table(
        pa.table(
            {
                "entry_pdb_id": pa.array([], type=pa.string()),
                "ligand_asym_id": pa.array([], type=pa.string()),
            }
        ),
        annotation,
    )
    marker = workspace / "entry_update.json"
    report = json.loads(marker.read_text())
    report["ingested_pdb_ids"] = []
    report["obsolete_pdb_ids"] = ["1abc", "2abd", "3def", "4ghi"]
    report["outputs"] = _file_stats(
        {
            **collate.entry_table_paths(workspace),
            "collation": workspace / "index/collation.json",
        }
    )
    write_json_atomic(marker, report)
    before = file_sha256(base / "ligand_archives/manifest.json")
    archives = update_archives.update_ligand_archives(workspace, memory_limit="1GB")
    assert archives["shard_count"] == archives["ligand_count"] == 0
    assert list((workspace / "ligand_archives").iterdir()) == [
        workspace / "ligand_archives/manifest.json"
    ]
    assert file_sha256(base / "ligand_archives/manifest.json") == before
