"""Precompute chain scores for proteins outside ligand pockets and interfaces."""

from __future__ import annotations

import argparse
import shutil
from collections.abc import Iterator
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from plinder.core.scores.custom import resolve_search_database
from plinder.core.utils.schemas import PROTEIN_SIMILARITY_EXPORT_SCHEMA
from plinder.data import databases
from plinder.data.annotations.get_similarity_scores import (
    get_sequence_similarity_helper,
    run_alignment,
)
from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig


def _chain_ids(data_dir: Path, backend: str) -> dict[str, tuple[str, str]]:
    chains = pd.read_parquet(
        data_dir / "index/alignment_chain_lookup.parquet",
        columns=["entry_pdb_id", "chain_auth_id", "chain_asym_id"],
    )
    prefix = "pdb_0000" if backend == "foldseek" else ""
    infix = "_xyz-enrich_" if backend == "foldseek" else "_"
    identifiers = (
        prefix
        + chains["entry_pdb_id"].astype(str)
        + infix
        + chains["chain_auth_id"].astype(str)
    )
    if identifiers.duplicated().any():
        raise ValueError("release chain lookup has ambiguous author-chain IDs")
    return dict(
        zip(
            identifiers,
            zip(chains["entry_pdb_id"], chains["chain_asym_id"], strict=True),
            strict=True,
        )
    )


def _percent(values: pd.Series, *, allow_missing: bool = False) -> pa.Array:
    numeric = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    missing = np.isnan(numeric)
    if (not allow_missing and missing.any()) or np.isinf(numeric).any():
        raise ValueError("protein alignment contains a non-finite score")
    scaled = np.floor(np.nan_to_num(numeric, nan=0) * 100 + 0.5)
    return pa.array(scaled.clip(0, 100).astype("uint8"), mask=missing)


def _alignment_frames(scanner: ds.Scanner) -> Iterator[pd.DataFrame]:
    """Combine small search row groups before converting and writing scores."""
    pending: list[pa.RecordBatch] = []
    count = 0
    for batch in scanner.to_batches():
        if not batch.num_rows:
            continue
        pending.append(batch)
        count += batch.num_rows
        if count >= 100_000:
            yield pa.Table.from_batches(pending).to_pandas()
            pending.clear()
            count = 0
    if pending:
        yield pa.Table.from_batches(pending).to_pandas()


def _score_alignment_batches(
    raw_path: Path,
    output_path: Path,
    *,
    backend: str,
    chain_ids: dict[str, tuple[str, str]],
) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(".tmp.parquet")
    count = 0

    @lru_cache(maxsize=100_000)
    def sequence_similarity(query: str, target: str) -> float:
        return get_sequence_similarity_helper(query, target)

    with pq.ParquetWriter(
        temporary, PROTEIN_SIMILARITY_EXPORT_SCHEMA, compression="zstd"
    ) as writer:
        scanner = ds.dataset(raw_path, format="parquet").scanner(
            columns=[
                "query",
                "target",
                "qcov",
                "tcov",
                "fident",
                "qaln",
                "taln",
                *(["lddt"] if backend == "foldseek" else []),
            ],
            batch_size=100_000,
        )
        for rows in _alignment_frames(scanner):
            query_ids = rows["query"]
            target_ids = rows["target"]
            if backend == "foldseek":
                model_prefix = r"_xyz-enrich_MODEL_\d+_"
                query_ids = query_ids.str.replace(
                    model_prefix, "_xyz-enrich_", regex=True
                )
                target_ids = target_ids.str.replace(
                    model_prefix, "_xyz-enrich_", regex=True
                )
            queries = query_ids.map(chain_ids)
            targets = target_ids.map(chain_ids)
            if queries.isna().any() or targets.isna().any():
                unknown = {
                    "query": rows.loc[queries.isna(), "query"].head(5).tolist(),
                    "target": rows.loc[targets.isna(), "target"].head(5).tolist(),
                }
                raise ValueError(f"monomer alignment has unknown chains: {unknown}")
            query_entry, query_chain = zip(*queries, strict=True)
            target_entry, target_chain = zip(*targets, strict=True)
            similarity = pd.Series(
                [
                    sequence_similarity(str(q).upper(), str(t).upper())
                    for q, t in zip(rows["qaln"], rows["taln"], strict=True)
                ]
            )
            size = len(rows)
            table = pa.Table.from_arrays(
                [
                    pa.array(query_entry, type=pa.string()),
                    pa.array(target_entry, type=pa.string()),
                    pa.array(query_chain, type=pa.string()),
                    pa.array(target_chain, type=pa.string()),
                    pa.array([backend] * size, type=pa.string()),
                    _percent(rows["qcov"]),
                    _percent(rows["tcov"]),
                    _percent(rows["fident"]),
                    _percent(similarity),
                    (
                        _percent(rows["lddt"], allow_missing=True)
                        if backend == "foldseek"
                        else pa.nulls(size, type=pa.uint8())
                    ),
                ],
                schema=PROTEIN_SIMILARITY_EXPORT_SCHEMA,
            )
            writer.write_table(table)
            count += size
    temporary.replace(output_path)
    return count


def score_monomer_pairs(
    data_dir: Path,
    *,
    backend: str,
    target_kind: str,
    scratch_dir: Path,
    output_path: Path,
    threads: int,
    pilot_queries: int | None = None,
    query_batch_index: int | None = None,
    query_batch_count: int | None = None,
) -> int:
    """Score unrepresented protein query chains against one target universe."""
    if backend not in {"foldseek", "mmseqs"}:
        raise ValueError(f"unsupported backend: {backend}")
    if target_kind not in {"holo", "monomer"}:
        raise ValueError(f"unsupported target universe: {target_kind}")
    if threads < 1:
        raise ValueError("threads must be positive")
    data_dir = Path(data_dir)
    scratch_dir = Path(scratch_dir)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    query = data_dir / "dbs/subdbs" / f"monomer_{backend}" / f"monomer_{backend}"
    if pilot_queries is not None or query_batch_index is not None:
        identifiers = sorted(databases.database_identifiers(query))
        if pilot_queries is not None:
            if pilot_queries < 1 or query_batch_index is not None:
                raise ValueError("pilot_queries must be positive and unbatched")
            identifiers = identifiers[:pilot_queries]
        else:
            if (
                query_batch_count is None
                or not 0 <= query_batch_index < query_batch_count
            ):
                raise ValueError("invalid monomer query batch")
            start = len(identifiers) * query_batch_index // query_batch_count
            stop = len(identifiers) * (query_batch_index + 1) // query_batch_count
            identifiers = identifiers[start:stop]
        print(f"Searching {len(identifiers):,} monomer query chains", flush=True)
        subset_root = scratch_dir / "query_subset"
        subset_root.mkdir(exist_ok=True)
        missing = databases.make_sub_db(set(identifiers), query, subset_root, backend)
        if missing:
            raise ValueError(f"missing pilot query chains: {sorted(missing)[:10]}")
        query = subset_root / subset_root.name

    target = resolve_search_database(
        backend, data_dir=data_dir, monomer=target_kind == "monomer"
    )
    config = FoldseekConfig() if backend == "foldseek" else MMSeqsConfig()
    raw = scratch_dir / "raw_alignment"
    run_alignment(
        aln_type=backend,
        query_db=query,
        target_db=target.conversion_target,
        search_target_db=target.search_target,
        cluster_alignment_db=target.cluster_alignments,
        search_db=scratch_dir / "results",
        aln_file=raw,
        alignment_config=config,
        tmp_dir=scratch_dir / "search_tmp",
        threads=threads,
    )
    try:
        return _score_alignment_batches(
            raw.with_suffix(".parquet"),
            Path(output_path),
            backend=backend,
            chain_ids=_chain_ids(data_dir, backend),
        )
    finally:
        shutil.rmtree(raw.with_suffix(".parquet"), ignore_errors=True)


def finalize_monomer_similarity_scores(
    data_dir: Path, *, query_batch_count: int = 8
) -> int:
    """Install the complete monomer-pair table after every batch succeeds."""
    exports = Path(data_dir) / "exports"
    staging = exports / ".monomer_similarity_scores_staging"
    installing = exports / ".monomer_similarity_scores_installing"
    installed = exports / "monomer_similarity_scores"
    if installed.exists():
        raise FileExistsError(installed)
    total_rows = 0
    for backend in ("mmseqs", "foldseek"):
        for target_kind in ("holo", "monomer"):
            for batch_index in range(query_batch_count):
                path = (
                    staging
                    / f"alignment_type={backend}"
                    / f"target_kind={target_kind}"
                    / f"part-{batch_index}.parquet"
                )
                if not path.is_file():
                    raise FileNotFoundError(path)
                file = pq.ParquetFile(path)
                if not file.schema_arrow.equals(PROTEIN_SIMILARITY_EXPORT_SCHEMA):
                    raise ValueError(f"invalid monomer score schema: {path}")
                total_rows += file.metadata.num_rows
    if installing.exists():
        shutil.rmtree(installing)
    for backend in ("mmseqs", "foldseek"):
        backend_root = installing / f"alignment_type={backend}"
        backend_root.mkdir(parents=True)
        for target_kind in ("holo", "monomer"):
            target_root = (
                staging / f"alignment_type={backend}" / f"target_kind={target_kind}"
            )
            for batch_index in range(query_batch_count):
                (
                    backend_root / f"monomer_{target_kind}_{batch_index}.parquet"
                ).hardlink_to(target_root / f"part-{batch_index}.parquet")
    installing.rename(installed)
    shutil.rmtree(staging, ignore_errors=True)
    return total_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("--backend", choices=("foldseek", "mmseqs"), required=True)
    parser.add_argument("--target-kind", choices=("holo", "monomer"), required=True)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--pilot-queries", type=int)
    parser.add_argument("--query-batch-index", type=int)
    parser.add_argument("--query-batch-count", type=int)
    args = parser.parse_args()
    count = score_monomer_pairs(
        args.data_dir,
        backend=args.backend,
        target_kind=args.target_kind,
        scratch_dir=args.scratch_dir,
        output_path=args.output,
        threads=args.threads,
        pilot_queries=args.pilot_queries,
        query_batch_index=args.query_batch_index,
        query_batch_count=args.query_batch_count,
    )
    print(f"Wrote {count:,} {args.backend}/{args.target_kind} pair rows", flush=True)


if __name__ == "__main__":
    main()
