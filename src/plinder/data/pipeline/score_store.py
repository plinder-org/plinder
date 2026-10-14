# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Repackage existing scores for selective downloads, without rescoring.

Run after score finalization and before preparing a public release. Output is
atomic per artifact, immutable once installed, and always separate from inputs.
"""

from __future__ import annotations

import argparse
import base64
import json
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.core.release import PlinderRelease
from plinder.core.scores.store import SCORE_KEYS


def _literal(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _bucket(column: str, *, numeric: bool) -> str:
    quoted = _identifier(column)
    if numeric:
        return f"'n' || CAST(CAST(floor({quoted} / 1024.0) AS BIGINT) AS VARCHAR)"
    return (
        f"CASE WHEN regexp_matches({quoted}, '^[a-zA-Z0-9]{{4}}(__|$)') "
        f"THEN 'p' || lower(substr({quoted}, 2, 2)) "
        f"ELSE 'h' || substr(sha256(split_part({quoted}, '__', 1)), 1, 2) END"
    )


def _partition_files(root: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for path in sorted(root.rglob("*.parquet")):
        bucket = path.parent.name.removeprefix("_bucket=")
        result.setdefault(bucket, []).append(path.relative_to(root.parent).as_posix())
    return result


@contextmanager
def _connection(
    threads: int, memory_limit: str, scratch: Path
) -> Iterator[duckdb.DuckDBPyConnection]:
    """Give each bounded unit its own engine and release its state on exit."""
    with duckdb.connect() as con:
        con.execute(f"SET threads = {threads}")
        con.execute(f"SET memory_limit = {_literal(memory_limit)}")
        con.execute(f"SET temp_directory = {_literal(scratch)}")
        con.execute("SET preserve_insertion_order = false")
        con.execute("SET partitioned_write_max_open_files = 8")
        yield con


def _sort_bucket(
    bucket_dir: Path,
    sorted_file: Path,
    columns: str,
    order: str,
    *,
    threads: int,
    memory_limit: str,
    scratch: Path,
) -> tuple[int, int, int]:
    """Sort/check one bucket without retaining previous buckets' sort state."""
    raw = f"read_parquet({_literal(bucket_dir / '*.parquet')}, hive_partitioning=false)"
    sorted_file.parent.mkdir(parents=True, exist_ok=True)
    with _connection(threads, memory_limit, scratch) as con:
        con.execute(
            f"COPY (SELECT {columns} FROM {raw} ORDER BY {order}) "
            f"TO {_literal(sorted_file)} (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 32768)"
        )
        fingerprint = (
            f"count(*), sum(CAST(hash({columns}) AS HUGEINT)), bit_xor(hash({columns}))"
        )
        before = con.execute(f"SELECT {fingerprint} FROM {raw}").fetchone()
        after = con.execute(
            f"SELECT {fingerprint} FROM read_parquet({_literal(sorted_file)})"
        ).fetchone()
        if before != after:
            raise ValueError(f"score-store fingerprint mismatch: {bucket_dir.name}")
        return int(after[0]), int(after[1]), int(after[2])


def _target_index(
    staging: Path,
    data: dict[str, list[str]],
    column: str,
    *,
    numeric: bool,
    threads: int,
    memory_limit: str,
    scratch: Path,
) -> None:
    """Deduplicate one query bucket at a time, then sort eight coarse groups.

    No release-wide DISTINCT, sort, or thousands-of-writers operation is held
    in memory. The temporary coarse partition has just eight writers; final
    partitions in separate coarse groups never overlap.
    """
    seeds = staging / f"index-parts-{column}"
    seeds.mkdir(parents=True, exist_ok=True)
    target = staging / f"target-{column}"
    target.mkdir(parents=True, exist_ok=True)
    quoted = _identifier(column)
    for number, (shard, relative_paths) in enumerate(sorted(data.items()), start=1):
        inputs = "[" + ", ".join(_literal(staging / p) for p in relative_paths) + "]"
        with _connection(threads, memory_limit, scratch) as con:
            con.execute(
                f"COPY (SELECT *, hash(target_bucket) % 8 AS _coarse FROM ("
                f"SELECT DISTINCT {_literal(column)} AS column, "
                f"CAST({quoted} AS VARCHAR) AS value, {_literal(shard)} AS shard, "
                f"{_bucket(column, numeric=numeric)} AS target_bucket "
                f"FROM read_parquet({inputs}, hive_partitioning=false))) "
                f"TO {_literal(seeds / shard)} (FORMAT parquet, PARTITION_BY (_coarse), "
                "COMPRESSION zstd, ROW_GROUP_SIZE 32768)"
            )
        if number % 100 == 0 or number == len(data):
            print(f"{column}: indexed {number}/{len(data)} query buckets", flush=True)
    for coarse in range(8):
        seed_files = sorted(seeds.glob(f"*/_coarse={coarse}/*.parquet"))
        if not seed_files:
            continue
        inputs = "[" + ", ".join(_literal(p) for p in seed_files) + "]"
        with _connection(threads, memory_limit, scratch) as con:
            con.execute(
                'COPY (SELECT "column", value, shard, target_bucket FROM '
                f"read_parquet({inputs}, hive_partitioning=false) "
                "ORDER BY target_bucket, value, shard) "
                f"TO {_literal(target / f'group-{coarse}')} "
                "(FORMAT parquet, PARTITION_BY (target_bucket), COMPRESSION zstd, ROW_GROUP_SIZE 32768)"
            )
        for path in seed_files:
            path.unlink()
        print(f"{column}: finalized coarse group {coarse + 1}/8", flush=True)
    if seeds.exists():
        shutil.rmtree(seeds)


def export_score_store(
    source: Path,
    destination: Path,
    artifact: str,
    *,
    threads: int = 2,
    memory_limit: str = "2GB",
    scratch_dir: Path | None = None,
) -> dict[str, Any]:
    """Create sorted query shards and exact target-to-shard routing indices.

    Original row multiplicity, columns, score direction, and Arrow types are
    retained. Installation refuses overwrite, including a partially staged
    destination, so an existing release cannot accidentally be modified.
    """
    query_keys, target_keys = SCORE_KEYS[artifact]
    release = PlinderRelease(source)
    input_path = release.fetch(artifact).resolve()
    output = destination.resolve() / "score_store" / artifact
    if output.exists():
        raise FileExistsError(output)
    if (
        input_path == output
        or input_path in output.parents
        or output in input_path.parents
    ):
        raise ValueError("score-store input and output must be separate")
    if threads < 1:
        raise ValueError("threads must be positive")
    dataset = pq.ParquetDataset(input_path)
    schema = dataset.schema
    # Arrow schemas record dictionary types, not their values. Keep the Hive
    # categories even when a routed query selects no rows (or no shards).
    partition_dictionaries = {}
    if dataset.partitioning is not None:
        for field, dictionary in zip(
            dataset.partitioning.schema, dataset.partitioning.dictionaries
        ):
            if dictionary is None or not pa.types.is_dictionary(field.type):
                continue
            table = pa.table({field.name: dictionary})
            sink = pa.BufferOutputStream()
            with pa.ipc.new_stream(sink, table.schema) as writer:
                writer.write_table(table)
            partition_dictionaries[field.name] = base64.b64encode(
                sink.getvalue().to_pybytes()
            ).decode()
    columns = ", ".join(_identifier(name) for name in schema.names)
    numeric = artifact == "ligand_scores"
    inputs = (
        str(input_path / "**/*.parquet") if input_path.is_dir() else str(input_path)
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{artifact}-", dir=output.parent))
    scratch = scratch_dir or staging / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    try:
        with _connection(threads, memory_limit, scratch) as con:
            con.execute(
                f"CREATE VIEW scores AS SELECT {columns} FROM "
                f"read_parquet({_literal(inputs)}, hive_partitioning=true)"
            )
            count = int(con.execute("SELECT count(*) FROM scores").fetchone()[0])
            print(f"{artifact}: exporting {count:,} rows", flush=True)
            order = ", ".join(_identifier(c) for c in (*query_keys, *target_keys))
            nulls = " OR ".join(
                f"{_identifier(c)} IS NULL" for c in (*query_keys, *target_keys)
            )
            if con.execute(
                f"SELECT EXISTS(SELECT 1 FROM scores WHERE {nulls})"
            ).fetchone()[0]:
                raise ValueError("score routing keys cannot be null")
            for alias in query_keys[1:]:
                mismatch = f"{_bucket(alias, numeric=numeric)} != {_bucket(query_keys[0], numeric=numeric)}"
                if con.execute(
                    f"SELECT EXISTS(SELECT 1 FROM scores WHERE {mismatch})"
                ).fetchone()[0]:
                    raise ValueError(
                        f"inconsistent query entry in routing alias: {alias}"
                    )
            if count:
                # Partition first: sorting a whole multi-billion-row export
                # would spill vastly more than sorting one bounded bucket.
                con.execute(
                    f"COPY (SELECT {columns}, {_bucket(query_keys[0], numeric=numeric)} "
                    f"AS _bucket FROM scores) "
                    f"TO {_literal(staging / 'raw')} (FORMAT parquet, PARTITION_BY (_bucket), "
                    "COMPRESSION zstd, ROW_GROUP_SIZE 32768)"
                )
        fingerprint_sum = 0
        fingerprint_xor = 0
        actual = 0
        bucket_dirs = sorted((staging / "raw").glob("_bucket=*"))
        for bucket_number, bucket_dir in enumerate(bucket_dirs, start=1):
            sorted_file = staging / "data" / bucket_dir.name / "data_0.parquet"
            rows, hash_sum, hash_xor = _sort_bucket(
                bucket_dir,
                sorted_file,
                columns,
                order,
                threads=threads,
                memory_limit=memory_limit,
                scratch=scratch,
            )
            actual += rows
            fingerprint_sum += hash_sum
            fingerprint_xor ^= hash_xor
            shutil.rmtree(bucket_dir)
            if bucket_number % 100 == 0 or bucket_number == len(bucket_dirs):
                print(
                    f"{artifact}: sorted {bucket_number}/{len(bucket_dirs)} buckets",
                    flush=True,
                )
        if (staging / "raw").exists():
            (staging / "raw").rmdir()
        if actual != count:
            raise ValueError(f"score-store row-count mismatch: {actual} != {count}")
        data = _partition_files(staging / "data")
        for column in target_keys if count else ():
            _target_index(
                staging,
                data,
                column,
                numeric=numeric,
                threads=threads,
                memory_limit=memory_limit,
                scratch=scratch,
            )
        targets: dict[str, list[str]] = {}
        for column in target_keys:
            for path in sorted((staging / f"target-{column}").rglob("*.parquet")):
                bucket = path.parent.name.removeprefix("target_bucket=")
                targets.setdefault(bucket, []).append(
                    path.relative_to(staging).as_posix()
                )
        manifest: dict[str, Any] = {
            "version": 1,
            "artifact": artifact,
            "schema": base64.b64encode(schema.serialize().to_pybytes()).decode(),
            "partition_dictionaries": partition_dictionaries,
            "rows": count,
            "fingerprint": {
                "hash_sum": str(fingerprint_sum),
                "hash_xor": str(fingerprint_xor),
            },
            "data": data,
            "targets": targets,
            "bytes": sum(p.stat().st_size for p in staging.rglob("*.parquet")),
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, sort_keys=True, indent=2) + "\n"
        )
        if scratch == staging / "scratch":
            shutil.rmtree(scratch)
        staging.rename(output)
        print(f"{artifact}: installed {len(data)} buckets at {output}", flush=True)
        return manifest
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--artifact", action="append", choices=tuple(SCORE_KEYS))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-limit", default="2GB")
    parser.add_argument("--scratch-dir", type=Path)
    args = parser.parse_args()
    for artifact in args.artifact or SCORE_KEYS:
        export_score_store(
            args.source,
            args.destination,
            artifact,
            threads=args.threads,
            memory_limit=args.memory_limit,
            scratch_dir=args.scratch_dir,
        )


if __name__ == "__main__":
    main()
