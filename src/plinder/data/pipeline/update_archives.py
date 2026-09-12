"""Carry canonical ligand archives into a prepared entry-update workspace."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from shutil import rmtree
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.core.utils.files import file_sha256, link_or_copy_file, write_json_atomic
from plinder.data.pipeline import collate, tasks
from plinder.data.pipeline.score import finalize_ligand_archives
from plinder.data.pipeline.update_entries import _file_stats


def update_ligand_archives(
    workspace: Path, *, threads: int = 4, memory_limit: str = "8GB"
) -> dict[str, Any]:
    """Replace changed entries' canonical SDFs while keeping the base read-only.

    Unchanged shards are hard-linked, or copied across filesystems. Treat release
    files as immutable: subsequent updates replace whole files, never edit them
    in place. Reruns rebuild affected shards and validate the complete archive set.
    Scores, clusters and other downstream work remain pending.
    """
    workspace = workspace.resolve(strict=True)
    entry_marker = workspace / "entry_update.json"
    entries = json.loads(entry_marker.read_text())
    if entries["status"] != collate.REPAIR_REQUIRED_STATUS:
        raise ValueError(
            "finish preparing entry tables before updating ligand archives"
        )
    base = Path(entries["base_release"]).resolve(strict=True)
    if workspace == base or base in workspace.parents:
        raise ValueError("update ligand archives outside the existing release")
    prepared = {
        **collate.entry_table_paths(workspace),
        "collation": workspace / "index" / collate.FINAL_MARKER_NAME,
    }
    base_tables = collate.entry_table_paths(base)
    if (
        _file_stats(prepared) != entries["outputs"]
        or _file_stats(base_tables) != entries["inputs"]["base_tables"]
    ):
        raise ValueError("entry tables changed; prepare a new update workspace")

    source_manifest = base / "ligand_archives/manifest.json"
    source_report = json.loads(source_manifest.read_text())
    old_archives = {
        path.stem: path for path in sorted((base / "ligand_archives").glob("*.parquet"))
    }
    if (
        source_report["status"] != "complete"
        or source_report["shard_count"] != len(old_archives)
        or source_report["ligand_count"]
        != sum(pq.read_metadata(path).num_rows for path in old_archives.values())
        or source_report["compressed_bytes"]
        != sum(path.stat().st_size for path in old_archives.values())
    ):
        raise ValueError("base ligand archives do not match their completion manifest")
    changed = sorted(
        set(entries["ingested_pdb_ids"]) | set(entries["obsolete_pdb_ids"])
    )
    changed_codes = {pdb_id[1:3] for pdb_id in changed}
    incoming = workspace / ".incoming"
    sdf_sources = {
        str(path): path
        for pdb_id in entries["ingested_pdb_ids"]
        for path in (
            incoming / "raw_entries" / pdb_id[1:3] / pdb_id / "ligand_files"
        ).glob("*.sdf")
    }
    sources = {
        str(path): path
        for path in [
            entry_marker,
            source_manifest,
            *old_archives.values(),
            *prepared.values(),
            *base_tables.values(),
        ]
    }
    sources.update(sdf_sources)
    before = _file_stats(sources)
    target = workspace / "ligand_archives"
    stage = workspace / ".ligand_archives.pending"
    backup = workspace / ".ligand_archives.previous"
    if backup.exists():
        if not target.exists():
            backup.rename(target)
        else:
            rmtree(backup)
    if stage.exists():
        rmtree(stage)
    stage.mkdir()
    try:
        tasks.make_canonical_ligand_archives(
            data_dir=incoming, two_char_codes=sorted(changed_codes)
        )
        expected_shards = []
        with duckdb.connect() as connection:
            collate._configure_duckdb(
                connection, threads=threads, memory_limit=memory_limit, scratch_dir=None
            )
            connection.register(
                "changed_entries",
                pa.table({"pdb_id": pa.array(changed, type=pa.string())}),
            )
            for code in sorted(set(old_archives) | changed_codes):
                output = stage / f"{code}.parquet"
                if code not in changed_codes:
                    link_or_copy_file(old_archives[code].resolve(strict=True), output)
                    expected_shards.append(code)
                    continue
                queries = []
                if code in old_archives:
                    connection.read_parquet(str(old_archives[code])).create_view(
                        "old_ligands", replace=True
                    )
                    queries.append(
                        "SELECT old_ligands.* FROM old_ligands ANTI JOIN changed_entries USING (pdb_id)"
                    )
                new_archive = incoming / "ligand_archives" / f"{code}.parquet"
                if new_archive.is_file():
                    connection.read_parquet(str(new_archive)).create_view(
                        "new_ligands", replace=True
                    )
                    queries.append("SELECT * FROM new_ligands")
                if not queries:
                    continue
                query = " UNION ALL BY NAME ".join(queries)
                if not collate._fetch_scalar(
                    connection, f"SELECT count(*) FROM ({query})"
                ):
                    continue
                collate._copy_query(
                    connection,
                    f"SELECT * FROM ({query}) ORDER BY pdb_id, ligand_asym_id",
                    output,
                    row_group_size=2_048,
                )
                expected_shards.append(code)
        report = finalize_ligand_archives(
            workspace,
            archive_dir=stage,
            expected_shards=expected_shards,
            threads=threads,
            memory_limit=memory_limit,
        )
        if before != _file_stats(sources):
            raise RuntimeError("archive inputs changed while preparing the update")
        report.update(
            updated_pdb_ids=changed,
            entry_update_sha256=file_sha256(entry_marker),
            base_manifest_sha256=file_sha256(source_manifest),
        )
        write_json_atomic(stage / "manifest.json", report)
        if target.exists():
            target.rename(backup)
        stage.rename(target)
    except BaseException:
        if stage.exists():
            rmtree(stage)
        if backup.exists() and not target.exists():
            backup.rename(target)
        raise
    if backup.exists():
        rmtree(backup)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-limit", default="8GB")
    args = parser.parse_args()
    report = update_ligand_archives(
        args.workspace, threads=args.threads, memory_limit=args.memory_limit
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
