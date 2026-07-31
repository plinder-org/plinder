# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Shared CLI for sharded V3 protein search and derived scoring."""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
from collections.abc import Iterable
from pathlib import Path
from shutil import copyfile, rmtree
from textwrap import dedent
from time import perf_counter
from typing import Any, cast

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from plinder.core.scores.metrics import DEFAULT_CLUSTER_METRICS
from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger
from plinder.data import clusters, databases
from plinder.data.pipeline import config, tasks

LOG = setup_logger(__name__)

MANIFEST_RELATIVE = Path("manifests/protein_scoring_queries.parquet")
PLAN_RELATIVE = Path("manifests/protein_scoring_plan.json")
FOLDSEEK_INPUT_RELATIVE = Path("manifests/foldseek_createdb_inputs.tsv")
SCORE_WORK_RELATIVE = Path("manifests/protein_scoring_work.parquet")
LIGAND_3D_CANDIDATE_MANIFEST_RELATIVE = Path("manifests/ligand_3d_candidates.parquet")
LIGAND_3D_WORK_RELATIVE = Path("manifests/ligand_3d_work.parquet")
LIGAND_3D_RETRY_WORK_RELATIVE = Path("manifests/ligand_3d_retry_work.parquet")
LIGAND_3D_PAIR_VALIDATION_RELATIVE = Path("scores/ligand_3d_pair_validation.json")
LIGAND_ARCHIVE_MANIFEST_RELATIVE = Path("ligand_archives/manifest.json")
DROPPED_QUERY_RELATIVE = Path("manifests/dropped_queries.parquet")
SCORE_REPAIR_RELATIVE = Path("manifests/score_repair.parquet")
SCORE_REPAIR_LIGAND_3D_RELATIVE = Path("manifests/score_repair_ligand_3d.parquet")
SCORE_REPAIR_DROP_DIR_RELATIVE = Path("manifests/score_repair_query_drops")
SCORE_REPAIR_DROP_RELATIVE = Path("manifests/score_repair_incomplete_queries.parquet")
INTERFACE_SCORE_WORK_RELATIVE = Path("manifests/interface_scoring_work.parquet")
INTERFACE_SCORE_PLAN_RELATIVE = Path("manifests/interface_scoring_plan.json")
INTERFACE_SCORE_ROOT_RELATIVE = Path("interface_scores")
INTERFACE_QCOV_EXPORT_RELATIVE = Path("exports/all_interface_qcov.parquet")
DEFAULT_CLUSTER_THRESHOLDS = (30, 50, 70, 90, 100)


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(exist_ok=True, parents=True)
    temporary = path.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(exist_ok=True, parents=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def _source_signature(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def dropped_query_ids(data_dir: Path) -> set[str]:
    """Return PDB IDs intentionally retained only as scoring targets."""
    path = data_dir / DROPPED_QUERY_RELATIVE
    if not path.is_file():
        return set()
    frame = pd.read_parquet(path, columns=["pdb_id"])
    return set(frame["pdb_id"].dropna().astype(str))


def active_scoring_query_ids(data_dir: Path) -> set[str]:
    """Return planned derived-score queries after intentional query drops."""
    work = pd.read_parquet(data_dir / SCORE_WORK_RELATIVE, columns=["pdb_id"])
    return set(work["pdb_id"].dropna().astype(str)).difference(
        dropped_query_ids(data_dir)
    )


def _load_pdb_id_manifest(path: Path) -> list[str]:
    """Load a PDB-ID column from Parquet or a newline-delimited ingest manifest."""
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix == ".parquet":
        frame = pd.read_parquet(path, columns=["pdb_id"])
        values = frame["pdb_id"].dropna().astype(str).str.lower().tolist()
    else:
        from plinder.data.pipeline.ingest import load_manifest

        values = load_manifest(path)
    if len(values) != len(set(values)):
        raise ValueError(f"PDB manifest contains duplicate IDs: {path}")
    return cast(list[str], values)


def _query_eligible_entry_ids(
    annotation: pd.DataFrame,
    entry_chains: pd.DataFrame,
    *,
    max_query_protein_chains: int,
    max_query_proper_ligand_chains: int,
) -> set[str]:
    """Return entries with at least one query-eligible holo system.

    The ligand limit applies only to proper ligand instances.  Keeping this
    calculation at system granularity avoids excluding a useful query because
    its system also contains ions, waters, artifacts, or other non-proper
    ligands.
    """
    required = {
        "entry_pdb_id",
        "system_id",
        "system_type",
        "system_protein_chains_asym_id",
        "ligand_id",
        "ligand_is_proper",
    }
    missing = required.difference(annotation.columns)
    if missing:
        raise ValueError(f"annotation is missing query eligibility columns: {missing}")
    holo = annotation.loc[
        annotation["system_type"].astype(str).eq("holo"), list(required)
    ].copy()
    if holo.empty:
        return set()
    keys = ["entry_pdb_id", "system_id"]
    chain_required = {"entry_pdb_id", "chain_asym_id", "chain_receptor_type"}
    missing_chains = chain_required.difference(entry_chains.columns)
    if missing_chains:
        raise ValueError(
            f"entry chains are missing eligibility columns: {missing_chains}"
        )
    protein_chains = {
        (str(row.entry_pdb_id), str(row.chain_asym_id))
        for row in entry_chains.loc[
            entry_chains["chain_receptor_type"].fillna("").astype(str).eq("protein")
        ].itertuples(index=False)
    }
    system_chains = holo.drop_duplicates(keys).set_index(keys)[
        "system_protein_chains_asym_id"
    ]
    system_protein_counts = pd.Series(
        {
            (entry_id, system_id): len(
                {
                    str(instance_chain)
                    for instance_chain in (
                        chains
                        if isinstance(chains, Iterable)
                        and not isinstance(chains, (str, bytes))
                        else []
                    )
                    if (str(entry_id), str(instance_chain).split(".", 1)[-1])
                    in protein_chains
                }
            )
            for (entry_id, system_id), chains in system_chains.items()
        },
        name="protein_chains",
    )
    system_protein_counts.index = system_protein_counts.index.set_names(keys)
    systems = system_protein_counts.to_frame()
    proper = holo["ligand_is_proper"].fillna(False).astype(bool)
    proper_counts = (
        holo.loc[proper]
        .groupby(keys, observed=True, sort=False)["ligand_id"]
        .nunique()
        .rename("proper_ligand_chains")
    )
    systems = systems.join(proper_counts, how="left").fillna(
        {"proper_ligand_chains": 0}
    )
    eligible = systems[
        systems["protein_chains"].gt(0)
        & systems["protein_chains"].le(max_query_protein_chains)
        & systems["proper_ligand_chains"].gt(0)
        & systems["proper_ligand_chains"].le(max_query_proper_ligand_chains)
    ]
    return set(eligible.index.get_level_values("entry_pdb_id").astype(str))


def _assign_repair_batches(
    records: list[dict[str, Any]],
    *,
    batch_size: int,
    index_offset: int,
) -> int:
    """Balance repair work while keeping every batch below its row limit."""
    batch_count = math.ceil(len(records) / batch_size)
    heap = [(0, 0, index_offset + index) for index in range(batch_count)]
    heapq.heapify(heap)
    for record in sorted(
        records,
        key=lambda item: (-int(item["estimated_work"]), str(item["pdb_id"])),
    ):
        load, count, index = heapq.heappop(heap)
        record["repair_batch_index"] = index
        count += 1
        if count < batch_size:
            heapq.heappush(heap, (load + int(record["estimated_work"]), count, index))
    return batch_count


def plan_score_repair(
    data_dir: Path,
    *,
    affected_manifest: Path,
    additional_full_query_manifest: Path | None = None,
    output_path: Path | None = None,
    batch_size: int = 10,
    target_batch_size: int | None = None,
    threads: int = 4,
    scratch_dir: Path | None = None,
    memory_limit: str = "32GB",
) -> dict[str, Any]:
    """Plan full-query and target-only rescoring for affected PDB entries.

    ``additional_full_query_manifest`` can promote queries whose scoreable
    system set changed independently of the affected target entries. These
    entries are rescored as queries but are not treated as changed targets.
    """
    if target_batch_size is None:
        target_batch_size = batch_size
    if min(batch_size, target_batch_size, threads) < 1:
        raise ValueError("repair batch size and threads must be positive")
    affected = set(_load_pdb_id_manifest(affected_manifest))
    if not affected:
        raise ValueError("affected-entry repair manifest is empty")
    protein_queries = set(
        pd.read_parquet(data_dir / MANIFEST_RELATIVE, columns=["pdb_id"])[
            "pdb_id"
        ].astype(str)
    ).difference(dropped_query_ids(data_dir))
    prior_plan = _load_plan(data_dir)
    max_query_protein_chains = int(prior_plan["score_max_query_protein_chains"])
    max_query_proper_ligand_chains = int(
        prior_plan["score_max_query_proper_ligand_chains"]
    )
    annotation = pd.read_parquet(
        data_dir / "index/annotation_table.parquet",
        columns=[
            "entry_pdb_id",
            "system_id",
            "system_type",
            "system_protein_chains_asym_id",
            "ligand_id",
            "ligand_is_proper",
        ],
    )
    entry_chains = pd.read_parquet(
        data_dir / "index/entry_chains.parquet",
        columns=["entry_pdb_id", "chain_asym_id", "chain_receptor_type"],
    )
    eligible_queries = _query_eligible_entry_ids(
        annotation,
        entry_chains,
        max_query_protein_chains=max_query_protein_chains,
        max_query_proper_ligand_chains=max_query_proper_ligand_chains,
    )
    active = protein_queries.intersection(eligible_queries)
    planned_active = active_scoring_query_ids(data_dir)
    if active != planned_active:
        raise RuntimeError(
            "the protein-scoring work manifest does not match the repaired index; "
            "rerun plan-score-batches before planning score repair"
        )
    requested_additional_full = (
        set(_load_pdb_id_manifest(additional_full_query_manifest))
        if additional_full_query_manifest is not None
        else set()
    )
    additional_full = requested_additional_full.intersection(active)
    full_queries = affected.intersection(active).union(additional_full)
    existing_affected_queries = {
        pdb_id
        for pdb_id in affected
        if (
            data_dir / "dbs" / "subdbs" / "search_db=holo" / f"{pdb_id}.parquet"
        ).is_file()
    }
    inactive_queries = existing_affected_queries.difference(active)
    alignment_paths = sorted(
        (data_dir / "alignments" / "search_db=holo").glob(
            "alignment_type=*/shard=*.parquet"
        )
    )
    if not alignment_paths:
        raise FileNotFoundError("no mapped holo alignment shards found")

    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    connection.sql(f"SET memory_limit='{memory_limit}'")
    if scratch_dir is not None:
        scratch_dir.mkdir(exist_ok=True, parents=True)
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
    connection.register(
        "affected_targets", pd.DataFrame({"target_entry": sorted(affected)})
    )
    connection.register("active_queries", pd.DataFrame({"query_entry": sorted(active)}))
    connection.register(
        "affected_queries", pd.DataFrame({"affected_query": sorted(full_queries)})
    )
    paths_sql = ", ".join(f"'{path.as_posix()}'" for path in alignment_paths)
    hits = connection.sql(
        f"""
        SELECT
            alignments.query_entry,
            list_sort(
                list(DISTINCT alignments.target_entry)
                FILTER (WHERE affected_targets.target_entry IS NOT NULL)
            ) AS target_entries,
            count(*) FILTER (
                WHERE affected_targets.target_entry IS NOT NULL
            )::BIGINT AS target_alignment_rows,
            count(*) FILTER (
                WHERE affected_queries.affected_query IS NOT NULL
            )::BIGINT AS full_alignment_rows
        FROM read_parquet([{paths_sql}], union_by_name=true) AS alignments
        INNER JOIN active_queries USING (query_entry)
        LEFT JOIN affected_targets USING (target_entry)
        LEFT JOIN affected_queries
          ON alignments.query_entry = affected_queries.affected_query
        WHERE affected_targets.target_entry IS NOT NULL
           OR affected_queries.affected_query IS NOT NULL
        GROUP BY alignments.query_entry
        """
    ).df()
    connection.close()

    targets_by_query: dict[str, list[str]] = {}
    work_by_query: dict[str, int] = {}
    if not hits.empty:
        for row in hits.itertuples(index=False):
            query_id = str(row.query_entry)
            target_entries = row.target_entries
            targets_by_query[query_id] = (
                list(map(str, target_entries))
                if isinstance(target_entries, Iterable)
                and not isinstance(target_entries, (str, bytes))
                else []
            )
            work_by_query[query_id] = (
                int(row.full_alignment_rows)
                if query_id in full_queries
                else int(row.target_alignment_rows)
            )

    all_queries = sorted(
        set(targets_by_query).union(full_queries).union(inactive_queries)
    )
    if not all_queries:
        raise ValueError("no active score queries are affected by the repair")
    score_work_path = data_dir / SCORE_WORK_RELATIVE
    full_work = {pdb_id: work_by_query.get(pdb_id, 1) for pdb_id in full_queries}
    if score_work_path.is_file():
        score_work = pd.read_parquet(
            score_work_path, columns=["pdb_id", "estimated_work"]
        )
        for row in score_work.itertuples(index=False):
            pdb_id = str(row.pdb_id)
            if pdb_id in full_queries:
                full_work[pdb_id] = max(
                    full_work.get(pdb_id, 1), max(1, int(row.estimated_work))
                )
    records: list[dict[str, Any]] = [
        {
            "pdb_id": pdb_id,
            "repair_mode": (
                "full"
                if pdb_id in full_queries
                else "drop"
                if pdb_id in inactive_queries
                else "targets"
            ),
            "target_pdb_ids": (
                []
                if pdb_id in full_queries or pdb_id in inactive_queries
                else targets_by_query[pdb_id]
            ),
            "affected_target_count": (
                0
                if pdb_id in full_queries or pdb_id in inactive_queries
                else len(targets_by_query[pdb_id])
            ),
            "estimated_work": (
                full_work.get(pdb_id, 1)
                if pdb_id in full_queries
                else 1
                if pdb_id in inactive_queries
                else max(1, work_by_query[pdb_id])
            ),
            "shard": pdb_id[1:3],
        }
        for pdb_id in all_queries
    ]

    full_records = [record for record in records if record["repair_mode"] == "full"]
    target_records = [record for record in records if record["repair_mode"] != "full"]
    full_batch_count = _assign_repair_batches(
        full_records,
        batch_size=batch_size,
        index_offset=0,
    )
    target_batch_count = _assign_repair_batches(
        target_records,
        batch_size=target_batch_size,
        index_offset=full_batch_count,
    )
    batch_count = full_batch_count + target_batch_count

    frame = pd.DataFrame(records).sort_values(
        ["repair_batch_index", "estimated_work", "pdb_id"],
        ascending=[True, False, True],
        ignore_index=True,
    )
    output = output_path or data_dir / SCORE_REPAIR_RELATIVE
    _atomic_parquet(frame, output)
    summary = {
        "status": "complete",
        "affected_manifest": _source_signature(affected_manifest),
        "additional_full_query_manifest": (
            _source_signature(additional_full_query_manifest)
            if additional_full_query_manifest is not None
            else None
        ),
        "output": _source_signature(output),
        "affected_entry_count": len(affected),
        "additional_full_query_count": len(additional_full),
        "ignored_additional_full_query_count": len(
            requested_additional_full.difference(active)
        ),
        "full_query_count": len(full_queries),
        "dropped_query_count": len(inactive_queries),
        "target_only_query_count": (
            len(all_queries) - len(full_queries) - len(inactive_queries)
        ),
        "query_count": len(all_queries),
        "query_shard_count": int(frame["shard"].nunique()) if not frame.empty else 0,
        "batch_size": max(batch_size, target_batch_size),
        "full_query_batch_size": batch_size,
        "full_query_batch_count": full_batch_count,
        "target_query_batch_size": target_batch_size,
        "target_query_batch_count": target_batch_count,
        "batch_count": batch_count,
    }
    _atomic_json(summary, output.with_suffix(".json"))
    return summary


def _score_repair_batch(
    manifest_path: Path, *, batch_index: int, batch_size: int
) -> list[dict[str, Any]]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("repair batch index and size must be positive")
    frame = pd.read_parquet(manifest_path)
    selected = frame[frame["repair_batch_index"].eq(batch_index)]
    if len(selected) > batch_size:
        raise ValueError(
            f"repair batch {batch_index} contains {len(selected)} rows; "
            f"batch size is {batch_size}"
        )
    return cast(list[dict[str, Any]], selected.to_dict("records"))


def _score_repair_query_paths(data_dir: Path, pdb_id: str) -> tuple[Path, Path]:
    return (
        data_dir / "dbs/subdbs/search_db=holo" / f"{pdb_id}.parquet",
        data_dir
        / "scores/ligand_3d_candidates/search_db=holo"
        / f"shard={pdb_id[1:3]}"
        / f"{pdb_id}.parquet",
    )


def _score_repair_query_is_current(
    data_dir: Path, pdb_id: str, *, repair_started_ns: int
) -> bool:
    score_path, candidate_path = _score_repair_query_paths(data_dir, pdb_id)
    try:
        return (
            score_path.stat().st_mtime_ns > repair_started_ns
            and candidate_path.stat().st_mtime_ns > repair_started_ns
        )
    except OSError:
        return False


def _score_repair_run_id(repair_manifest: Path) -> str:
    """Return a stable identifier for one immutable repair manifest."""
    signature = json.dumps(
        _source_signature(repair_manifest), sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(signature.encode()).hexdigest()[:16]


def _score_repair_marker_dir(data_dir: Path, repair_manifest: Path) -> Path:
    """Keep resource-limit markers isolated by repair run."""
    return (
        data_dir
        / SCORE_REPAIR_DROP_DIR_RELATIVE
        / _score_repair_run_id(repair_manifest)
    )


def _score_repair_marked_targets(
    data_dir: Path,
    repair_manifest: Path,
) -> set[str]:
    """Load and validate target-query drop markers for one repair run."""
    marker_dir = _score_repair_marker_dir(data_dir, repair_manifest)
    repair_run_id = _score_repair_run_id(repair_manifest)
    marked_targets: set[str] = set()
    if not marker_dir.is_dir():
        return marked_targets
    for path in sorted(marker_dir.glob("*.json")):
        if path.name.endswith(".tmp.json"):
            continue
        payload = json.loads(path.read_text())
        pdb_id = str(payload.get("pdb_id", ""))
        if (
            payload.get("repair_run_id") != repair_run_id
            or payload.get("repair_mode") != "targets"
            or not pdb_id
            or path.stem != pdb_id
        ):
            raise ValueError(f"invalid target repair drop marker: {path}")
        marked_targets.add(pdb_id)
    return marked_targets


def repair_target_score_remainders(
    data_dir: Path,
    *,
    repair_manifest: Path,
    batch_index: int,
    batch_size: int,
    scorer_cfg: Any,
    scratch_dir: Path,
    threads: int,
) -> dict[str, Any]:
    """Drop a timed-out target query and resume the rest of its batch."""
    repairs = _score_repair_batch(
        repair_manifest, batch_index=batch_index, batch_size=batch_size
    )
    repair_started_ns = repair_manifest.stat().st_mtime_ns
    marked_targets = _score_repair_marked_targets(data_dir, repair_manifest)
    for repair in repairs:
        if str(repair["repair_mode"]) == "drop":
            for path in _score_repair_query_paths(data_dir, str(repair["pdb_id"])):
                path.unlink(missing_ok=True)
    incomplete = [
        repair
        for repair in repairs
        if str(repair["repair_mode"]) != "drop"
        and str(repair["pdb_id"]) not in marked_targets
        and not _score_repair_query_is_current(
            data_dir,
            str(repair["pdb_id"]),
            repair_started_ns=repair_started_ns,
        )
    ]
    if not incomplete:
        return {
            "status": "complete",
            "batch_index": batch_index,
            "dropped_query": None,
            "resumed_query_count": 0,
        }
    invalid_modes = sorted(
        {
            str(repair["repair_mode"])
            for repair in incomplete
            if str(repair["repair_mode"]) != "targets"
        }
    )
    if invalid_modes:
        raise ValueError(
            "target remainder repair received non-target modes: " f"{invalid_modes}"
        )

    stalled = incomplete[0]
    stalled_pdb_id = str(stalled["pdb_id"])
    for path in _score_repair_query_paths(data_dir, stalled_pdb_id):
        path.unlink(missing_ok=True)
    marker = (
        _score_repair_marker_dir(data_dir, repair_manifest) / f"{stalled_pdb_id}.json"
    )
    _atomic_json(
        {
            "pdb_id": stalled_pdb_id,
            "repair_run_id": _score_repair_run_id(repair_manifest),
            "repair_batch_index": batch_index,
            "repair_mode": "targets",
            "reason": "resource_limit",
        },
        marker,
    )
    remainder = incomplete[1:]
    if remainder:
        tasks.repair_batch_scores(
            data_dir=data_dir,
            repairs=remainder,
            scorer_cfg=scorer_cfg,
            scratch_dir=scratch_dir,
            threads=threads,
        )
    return {
        "status": "complete",
        "batch_index": batch_index,
        "dropped_query": stalled_pdb_id,
        "resumed_query_count": len(remainder),
    }


def finalize_score_repair_queries(
    data_dir: Path,
    *,
    repair_manifest: Path,
    max_new_drops: int,
) -> dict[str, Any]:
    """Validate repaired queries and record only resource-limited tail queries."""
    if max_new_drops < 0:
        raise ValueError("max_new_drops must be non-negative")
    repairs = pd.read_parquet(repair_manifest)
    required = {"pdb_id", "repair_mode"}
    missing = sorted(required.difference(repairs.columns))
    if missing:
        raise ValueError(f"repair manifest is missing columns {missing}")
    repair_started_ns = repair_manifest.stat().st_mtime_ns
    incomplete_by_mode: dict[str, set[str]] = {}
    for row in repairs.itertuples(index=False):
        mode = str(row.repair_mode)
        pdb_id = str(row.pdb_id)
        if mode == "drop" or _score_repair_query_is_current(
            data_dir, pdb_id, repair_started_ns=repair_started_ns
        ):
            continue
        incomplete_by_mode.setdefault(mode, set()).add(pdb_id)

    marked_targets = _score_repair_marked_targets(data_dir, repair_manifest)
    incomplete_targets = incomplete_by_mode.get("targets", set())
    unexpected_targets = sorted(incomplete_targets.difference(marked_targets))
    if unexpected_targets:
        raise ValueError(
            "target score repairs remain incomplete without resource-limit markers: "
            f"{unexpected_targets[:20]}"
        )
    stale_markers = sorted(marked_targets.difference(incomplete_targets))
    if stale_markers:
        raise ValueError(
            f"target repair drop markers unexpectedly have current outputs: "
            f"{stale_markers[:20]}"
        )
    invalid_modes = sorted(set(incomplete_by_mode).difference({"full", "targets"}))
    if invalid_modes:
        raise ValueError(f"unsupported incomplete repair modes: {invalid_modes}")

    new_drops = marked_targets | incomplete_by_mode.get("full", set())
    if len(new_drops) > max_new_drops:
        raise ValueError(
            f"score repair would drop {len(new_drops)} queries, exceeding "
            f"the configured limit {max_new_drops}"
        )
    for pdb_id in new_drops:
        for path in _score_repair_query_paths(data_dir, pdb_id):
            path.unlink(missing_ok=True)

    existing_columns = ["pdb_id", "stage", "reason", "details"]
    existing = pd.DataFrame(columns=existing_columns)
    preserved_drops: set[str] = set()
    existing_drop_path = data_dir / DROPPED_QUERY_RELATIVE
    if existing_drop_path.is_file():
        existing = pd.read_parquet(existing_drop_path)
        missing = sorted(set(existing_columns).difference(existing.columns))
        if missing:
            raise ValueError(f"existing dropped-query manifest is missing {missing}")
        preserved_drops = set(
            existing.loc[existing["stage"].eq("derived_scoring"), "pdb_id"].astype(str)
        )
    requested_drops = sorted(preserved_drops | new_drops)
    source = data_dir / SCORE_REPAIR_DROP_RELATIVE
    _atomic_parquet(pd.DataFrame({"pdb_id": requested_drops}), source)
    existing_ids = set(existing["pdb_id"].astype(str))
    newly_recorded = sorted(new_drops.difference(existing_ids))
    details = json.dumps({"source_manifest": source.name}, sort_keys=True)
    additions = pd.DataFrame(
        [
            {
                "pdb_id": pdb_id,
                "stage": "derived_scoring",
                "reason": "resource_limit",
                "details": details,
            }
            for pdb_id in newly_recorded
        ],
        columns=existing_columns,
    )
    dropped = pd.concat([existing[existing_columns], additions], ignore_index=True)
    if dropped["pdb_id"].astype(str).duplicated().any():
        raise ValueError("a query cannot be dropped at multiple stages")
    dropped = dropped.sort_values("pdb_id", ignore_index=True)
    plan_path = data_dir / PLAN_RELATIVE
    if not plan_path.is_file():
        raise FileNotFoundError(plan_path)
    plan = json.loads(plan_path.read_text())
    _atomic_parquet(dropped, existing_drop_path)
    stage_counts = {
        str(stage): int(count)
        for stage, count in dropped.groupby("stage", dropna=False).size().items()
    }
    plan["dropped_queries"] = _source_signature(existing_drop_path)
    plan["dropped_query_counts"] = {
        **stage_counts,
        "total": len(dropped),
    }
    _atomic_json(plan, plan_path)
    dropped_report = {
        "new": len(newly_recorded),
        "derived_scoring": stage_counts.get("derived_scoring", 0),
        "total": len(dropped),
    }
    report = {
        "status": "complete",
        "new_target_query_drops": len(marked_targets),
        "new_full_query_drops": len(incomplete_by_mode.get("full", set())),
        "new_query_drops": len(new_drops),
        "preserved_query_drops": len(preserved_drops),
        "dropped_queries": dropped_report,
    }
    _atomic_json(report, source.with_suffix(".json"))
    return report


def _score_repair_shard_batch(
    manifest_path: Path, *, batch_index: int, batch_size: int
) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("repair shard batch index and size must be positive")
    shards = sorted(
        set(pd.read_parquet(manifest_path, columns=["shard"])["shard"].astype(str))
    )
    start = batch_index * batch_size
    return shards[start : start + batch_size]


def plan_score_repair_ligand_3d(
    data_dir: Path,
    *,
    repair_manifest: Path,
    batch_size: int,
    threads: int,
    scratch_dir: Path,
    memory_limit: str = "8GB",
) -> dict[str, Any]:
    """Plan only canonical ligand pairs absent from the existing score cache."""
    if min(batch_size, threads) < 1:
        raise ValueError("repair ligand batch size and threads must be positive")
    shards = sorted(
        set(pd.read_parquet(repair_manifest, columns=["shard"])["shard"].astype(str))
    )
    candidate_paths = [
        data_dir / "scores/ligand_3d_pair_candidate_shards" / f"shard={shard}.parquet"
        for shard in shards
    ]
    cached_paths = [
        data_dir / "scores/ligand_3d_by_query" / f"{shard}.parquet" for shard in shards
    ]
    missing_inputs = [
        path for path in [*candidate_paths, *cached_paths] if not path.is_file()
    ]
    if missing_inputs:
        raise FileNotFoundError(
            f"missing ligand 3D repair inputs: {missing_inputs[:10]}"
        )
    scratch_dir.mkdir(exist_ok=True, parents=True)
    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
    connection.sql(f"SET memory_limit='{memory_limit}'")
    candidate_sql = ", ".join(f"'{path.as_posix()}'" for path in candidate_paths)
    cached_sql = ", ".join(f"'{path.as_posix()}'" for path in cached_paths)
    annotation = data_dir / "index/annotation_table.parquet"
    temporary = scratch_dir / "score-repair-ligand-3d.parquet"
    temporary.unlink(missing_ok=True)
    connection.sql(
        f"""
        COPY (
            WITH candidates AS (
                SELECT DISTINCT
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id
                FROM read_parquet([{candidate_sql}])
                WHERE pocket_qcov > 0
            ), cached AS (
                SELECT DISTINCT
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id
                FROM read_parquet([{cached_sql}])
            ), missing AS (
                SELECT candidates.*
                FROM candidates ANTI JOIN cached USING (
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id
                )
            ), ligand_sizes AS (
                SELECT
                    entry_pdb_id AS entry,
                    ligand_asym_id AS asym_id,
                    max(greatest(coalesce(ligand_num_heavy_atoms, 1), 1))::BIGINT
                        AS heavy_atoms
                FROM read_parquet('{annotation.as_posix()}')
                WHERE system_type = 'holo'
                  AND coalesce(ligand_is_proper, false)
                  AND coalesce(ligand_is_3d_score_able, false)
                GROUP BY entry_pdb_id, ligand_asym_id
            ), sized AS (
                SELECT
                    missing.*,
                    (query.heavy_atoms * target.heavy_atoms)::BIGINT
                        AS estimated_work
                FROM missing
                INNER JOIN ligand_sizes AS query
                  ON missing.query_entry = query.entry
                 AND missing.query_ligand_asym_id = query.asym_id
                INNER JOIN ligand_sizes AS target
                  ON missing.target_entry = target.entry
                 AND missing.target_ligand_asym_id = target.asym_id
            )
            SELECT
                *,
                substr(query_entry, 2, 2) AS shard,
                floor((row_number() OVER (
                    ORDER BY
                        query_entry,
                        query_ligand_asym_id,
                        estimated_work DESC,
                        target_entry,
                        target_ligand_asym_id
                ) - 1) / {batch_size})::INTEGER AS repair_pair_batch_index
            FROM sized
        ) TO '{temporary.as_posix()}' (
            FORMAT PARQUET,
            COMPRESSION ZSTD,
            ROW_GROUP_SIZE {max(2_048, batch_size)}
        )
        """
    )
    counts = connection.sql(
        f"""
        SELECT
            count(*)::BIGINT,
            coalesce(max(repair_pair_batch_index) + 1, 0)::BIGINT,
            count(DISTINCT shard)::BIGINT
        FROM read_parquet('{temporary.as_posix()}')
        """
    ).fetchone()
    connection.close()
    if counts is None:
        raise RuntimeError("failed to count missing repair ligand pairs")
    output = data_dir / SCORE_REPAIR_LIGAND_3D_RELATIVE
    output.parent.mkdir(exist_ok=True, parents=True)
    install = output.with_suffix(output.suffix + ".tmp")
    copyfile(temporary, install)
    install.replace(output)
    temporary.unlink(missing_ok=True)
    pair_count, batch_count, shard_count = map(int, counts)
    report = {
        "status": "complete",
        "repair_manifest": _source_signature(repair_manifest),
        "output": _source_signature(output),
        "pair_count": pair_count,
        "batch_size": batch_size,
        "batch_count": batch_count,
        "shard_count": shard_count,
    }
    _atomic_json(report, output.with_suffix(".json"))
    return report


def _score_repair_ligand_3d_batch(
    data_dir: Path, *, batch_index: int, batch_size: int
) -> pd.DataFrame:
    work = pd.read_parquet(
        data_dir / SCORE_REPAIR_LIGAND_3D_RELATIVE,
        filters=[("repair_pair_batch_index", "==", batch_index)],
    )
    if len(work) > batch_size:
        raise ValueError(
            f"repair ligand batch {batch_index} contains {len(work)} pairs; "
            f"batch size is {batch_size}"
        )
    return work


def merge_score_repair_ligand_3d(
    data_dir: Path,
    *,
    shards: list[str],
    scratch_dir: Path,
    threads: int,
) -> dict[str, Any]:
    """Atomically append newly scored repair pairs to canonical query caches."""
    if threads < 1:
        raise ValueError("repair merge threads must be positive")
    import duckdb

    work_path = data_dir / SCORE_REPAIR_LIGAND_3D_RELATIVE
    repair_dir = data_dir / "scores/ligand_3d_pair_repairs"
    scratch_dir.mkdir(exist_ok=True, parents=True)
    updated: list[str] = []
    for shard in shards:
        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
        batch_rows = connection.sql(
            f"""
            SELECT DISTINCT repair_pair_batch_index
            FROM read_parquet('{work_path.as_posix()}')
            WHERE shard = '{shard}'
            ORDER BY repair_pair_batch_index
            """
        ).fetchall()
        if not batch_rows:
            connection.close()
            continue
        repair_paths = [repair_dir / f"{int(row[0])}.parquet" for row in batch_rows]
        absent = [path for path in repair_paths if not path.is_file()]
        if absent:
            connection.close()
            raise FileNotFoundError(
                f"missing ligand 3D repair batches for shard {shard}: {absent[:10]}"
            )
        cached = data_dir / "scores/ligand_3d_by_query" / f"{shard}.parquet"
        paths_sql = ", ".join(f"'{path.as_posix()}'" for path in repair_paths)
        keys = ", ".join(
            [
                "query_entry",
                "query_ligand_asym_id",
                "target_entry",
                "target_ligand_asym_id",
            ]
        )
        validation = connection.sql(
            f"""
            WITH expected AS (
                SELECT {keys}
                FROM read_parquet('{work_path.as_posix()}')
                WHERE shard = '{shard}'
            ), observed AS (
                SELECT scored.*
                FROM read_parquet([{paths_sql}]) AS scored
                INNER JOIN expected USING ({keys})
            )
            SELECT
                (SELECT count(*) FROM expected),
                (SELECT count(*) FROM observed),
                (SELECT count(*) - count(DISTINCT ({keys})) FROM observed),
                (SELECT count(*) FROM observed
                 INNER JOIN read_parquet('{cached.as_posix()}') USING ({keys}))
            """
        ).fetchone()
        if (
            validation is None
            or validation[0] != validation[1]
            or validation[2]
            or validation[3]
        ):
            connection.close()
            raise ValueError(
                f"invalid ligand 3D repair coverage for shard {shard}: {validation}"
            )
        temporary = scratch_dir / f"{shard}.parquet"
        temporary.unlink(missing_ok=True)
        connection.sql(
            f"""
            COPY (
                WITH expected AS (
                    SELECT {keys}
                    FROM read_parquet('{work_path.as_posix()}')
                    WHERE shard = '{shard}'
                ), repaired AS (
                    SELECT scored.*
                    FROM read_parquet([{paths_sql}]) AS scored
                    INNER JOIN expected USING ({keys})
                )
                SELECT * FROM read_parquet('{cached.as_posix()}')
                UNION ALL
                SELECT * FROM repaired
                ORDER BY {keys}
            ) TO '{temporary.as_posix()}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
        connection.close()
        install = cached.with_suffix(cached.suffix + ".tmp")
        copyfile(temporary, install)
        install.replace(cached)
        temporary.unlink(missing_ok=True)
        updated.append(shard)
    return {"status": "complete", "updated_shards": updated}


def record_dropped_queries(
    data_dir: Path,
    *,
    score_query_manifest: Path,
) -> dict[str, Any]:
    """Record mapping and resource-limited score queries in one manifest."""
    if not score_query_manifest.is_file():
        raise FileNotFoundError(score_query_manifest)
    score_drops = pd.read_parquet(score_query_manifest)
    if "pdb_id" not in score_drops.columns:
        raise ValueError(
            f"dropped-query source is missing pdb_id: {score_query_manifest}"
        )
    requested_score_ids = set(score_drops["pdb_id"].dropna().astype(str))
    if score_drops["pdb_id"].dropna().astype(str).duplicated().any():
        raise ValueError(
            f"dropped-query source contains duplicate PDB IDs: {score_query_manifest}"
        )

    alignment_manifest = data_dir / "alignments" / "manifest.json"
    if not alignment_manifest.is_file():
        raise FileNotFoundError(alignment_manifest)
    alignment_payload = json.loads(alignment_manifest.read_text())
    mapping_details = alignment_payload.get("skipped_queries")
    if not isinstance(mapping_details, dict):
        raise ValueError("alignment manifest has no skipped-query mapping")
    mapping_ids = set(map(str, mapping_details))

    output = data_dir / DROPPED_QUERY_RELATIVE
    if output.is_file():
        existing = pd.read_parquet(output)
        required = {"pdb_id", "stage", "reason", "details"}
        missing = sorted(required.difference(existing.columns))
        if missing:
            raise ValueError(f"existing dropped-query manifest is missing {missing}")
    score_ids = set()
    for pdb_id in requested_score_ids.difference(mapping_ids):
        score_path = (
            data_dir / "dbs" / "subdbs" / "search_db=holo" / f"{pdb_id}.parquet"
        )
        candidate_path = (
            data_dir
            / "scores"
            / "ligand_3d_candidates"
            / "search_db=holo"
            / f"shard={pdb_id[1:3]}"
            / f"{pdb_id}.parquet"
        )
        try:
            metadata = pq.read_schema(score_path).metadata or {}
        except (OSError, ValueError):
            metadata = {}
        complete = (
            metadata.get(b"plinder.ligand_3d") in {b"deferred", b"complete"}
            and candidate_path.is_file()
        )
        if not complete:
            score_ids.add(pdb_id)

    planned_ids = set(
        pd.read_parquet(data_dir / SCORE_WORK_RELATIVE, columns=["pdb_id"])[
            "pdb_id"
        ].astype(str)
    )
    unexpected = sorted((mapping_ids | score_ids).difference(planned_ids))
    if unexpected:
        raise ValueError(
            f"dropped queries are absent from the scoring plan: {unexpected[:10]}"
        )

    rows = [
        {
            "pdb_id": pdb_id,
            "stage": "alignment_mapping",
            "reason": str(details.get("reason", "alignment_mapping_skipped")),
            "details": json.dumps(details, sort_keys=True),
        }
        for pdb_id, details in sorted(mapping_details.items())
    ]
    score_detail = json.dumps(
        {"source_manifest": score_query_manifest.name}, sort_keys=True
    )
    rows.extend(
        {
            "pdb_id": pdb_id,
            "stage": "derived_scoring",
            "reason": "resource_limit",
            "details": score_detail,
        }
        for pdb_id in sorted(score_ids)
    )
    frame = pd.DataFrame(rows, columns=["pdb_id", "stage", "reason", "details"])
    if frame["pdb_id"].duplicated().any():
        raise ValueError("a query cannot be dropped at multiple stages")
    plan = _load_plan(data_dir)
    _atomic_parquet(frame, output)

    plan["dropped_queries"] = _source_signature(output)
    plan["dropped_query_counts"] = {
        "alignment_mapping": len(mapping_ids),
        "derived_scoring": len(score_ids),
        "total": len(frame),
    }
    _atomic_json(plan, data_dir / PLAN_RELATIVE)
    return {
        "status": "complete",
        "alignment_mapping": len(mapping_ids),
        "requested_derived_scoring": len(requested_score_ids),
        "completed_since_retry_manifest": len(requested_score_ids)
        - len(mapping_ids.intersection(requested_score_ids))
        - len(score_ids),
        "derived_scoring": len(score_ids),
        "total": len(frame),
        "path": str(output),
    }


def plan_protein_scoring(
    data_dir: Path,
    *,
    pdb_ids: list[str] | None = None,
    two_char_codes: list[str] | None = None,
    max_seqs: int = 10_000,
) -> dict[str, Any]:
    """Freeze the protein-query universe used by every array stage."""
    if max_seqs < 1:
        raise ValueError("max_seqs must be positive")
    chunks = tasks.scatter_protein_scoring(
        data_dir=data_dir,
        batch_size=1_000_000,
        two_char_codes=two_char_codes or [],
        pdb_ids=pdb_ids or [],
    )
    query_ids = [pdb_id for chunk in chunks for pdb_id in chunk]
    if not query_ids:
        raise ValueError("no protein-containing V3 entries were selected")
    chain_path = data_dir / "index" / "entry_chains.parquet"
    chains = tasks._protein_scoring_chains(data_dir)
    chains = chains[chains["entry_pdb_id"].isin(query_ids)]
    counts = chains.groupby("entry_pdb_id")["chain_asym_id"].nunique().to_dict()
    manifest = pd.DataFrame(
        {
            "pdb_id": query_ids,
            "shard": [pdb_id[1:3] for pdb_id in query_ids],
            "protein_chain_count": [counts[pdb_id] for pdb_id in query_ids],
        }
    ).sort_values(["shard", "pdb_id"], ignore_index=True)
    manifest_path = data_dir / MANIFEST_RELATIVE
    _atomic_parquet(manifest, manifest_path)
    payload: dict[str, Any] = {
        "query_count": len(manifest),
        "protein_chain_count": int(manifest["protein_chain_count"].sum()),
        "shard_count": int(manifest["shard"].nunique()),
        "max_seqs": max_seqs,
        "search_databases": ["holo"],
        "alignment_types": ["foldseek", "mmseqs"],
        "target_clustering": {
            "identity": 1.0,
            "coverage": 1.0,
            "coverage_mode": 0,
            "expand_to_chain_level": True,
        },
        "entry_chains": _source_signature(chain_path),
        "manifest": _source_signature(manifest_path),
    }
    interface_path = data_dir / "index" / "interface_annotation_table.parquet"
    if interface_path.is_file():
        payload["interface_annotation"] = _source_signature(interface_path)
    _atomic_json(payload, data_dir / PLAN_RELATIVE)
    return payload


def _load_plan(
    data_dir: Path,
    *,
    recheck_source: bool = False,
    recheck_score_inputs: bool = True,
) -> dict[str, Any]:
    plan_path = data_dir / PLAN_RELATIVE
    if not plan_path.is_file():
        raise FileNotFoundError(f"missing protein scoring plan: {plan_path}")
    plan: dict[str, Any] = json.loads(plan_path.read_text())
    manifest_path = data_dir / MANIFEST_RELATIVE
    if _source_signature(manifest_path) != plan["manifest"]:
        raise ValueError("protein scoring query manifest changed after planning")
    dropped_signature = plan.get("dropped_queries")
    if dropped_signature is not None:
        dropped_path = data_dir / DROPPED_QUERY_RELATIVE
        if (
            not dropped_path.is_file()
            or _source_signature(dropped_path) != dropped_signature
        ):
            raise ValueError("dropped-query manifest changed after planning")
    if recheck_source:
        for key, description in [
            ("entry_chains", "entry chain index"),
            ("interface_annotation", "interface annotation index"),
        ]:
            signature = plan.get(key)
            if not isinstance(signature, dict):
                continue
            source = Path(str(signature["path"]))
            if _source_signature(source) != signature:
                raise ValueError(
                    f"{description} changed after protein scoring was planned"
                )
    score_annotation = plan.get("score_annotation")
    if recheck_score_inputs and isinstance(score_annotation, dict):
        source = Path(str(score_annotation["path"]))
        if _source_signature(source) != score_annotation:
            raise ValueError(
                "annotation index changed after score batches were planned"
            )
    return plan


def _query_batch(data_dir: Path, batch_index: int, batch_size: int) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    _load_plan(data_dir)
    queries = pd.read_parquet(data_dir / MANIFEST_RELATIVE, columns=["pdb_id"])
    start = batch_index * batch_size
    values = queries["pdb_id"].iloc[start : start + batch_size].tolist()
    return [str(value) for value in values]


def _score_batch(data_dir: Path, batch_index: int, batch_size: int) -> list[str]:
    """Load one cost-balanced score batch."""
    work_path = data_dir / SCORE_WORK_RELATIVE
    if not work_path.is_file():
        raise FileNotFoundError(
            f"missing score work plan; run plan-score-batches first: {work_path}"
        )
    plan = _load_plan(data_dir)
    if batch_size != int(plan.get("score_batch_size", -1)):
        raise ValueError(
            f"score batch size is {plan.get('score_batch_size')}, got {batch_size}"
        )
    if _source_signature(work_path) != plan.get("score_work"):
        raise ValueError("protein scoring work manifest changed after planning")
    work = pd.read_parquet(
        work_path,
        columns=["pdb_id"],
        filters=[("score_batch_index", "==", batch_index)],
    )
    dropped = dropped_query_ids(data_dir)
    if dropped:
        work = work[~work["pdb_id"].astype(str).isin(dropped)]
    return [str(value) for value in work["pdb_id"]]


def _score_manifest_batch(
    data_dir: Path,
    manifest_path: Path,
    batch_index: int,
    batch_size: int,
) -> list[str]:
    """Load a fixed batch from an explicit PDB retry manifest."""
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = pd.read_parquet(manifest_path)
    if "pdb_id" not in manifest.columns:
        raise ValueError(f"score retry manifest is missing pdb_id: {manifest_path}")
    pdb_ids = manifest["pdb_id"].dropna().astype(str)
    if pdb_ids.duplicated().any():
        raise ValueError(
            f"score retry manifest contains duplicate PDB IDs: {manifest_path}"
        )
    dropped = dropped_query_ids(data_dir)
    if dropped:
        pdb_ids = pdb_ids[~pdb_ids.isin(dropped)]
    if "retry_batch_index" in manifest.columns:
        selected = manifest.loc[
            manifest["retry_batch_index"].eq(batch_index), "pdb_id"
        ].dropna()
        if len(selected) > batch_size:
            raise ValueError(
                f"retry batch {batch_index} contains {len(selected)} PDB IDs, "
                f"exceeding batch size {batch_size}"
            )
        return [str(value) for value in selected]
    start = batch_index * batch_size
    return [str(value) for value in pdb_ids.iloc[start : start + batch_size]]


def _score_shard_batch(data_dir: Path, batch_index: int, batch_size: int) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    _load_plan(data_dir)
    shards = sorted({pdb_id[1:3] for pdb_id in active_scoring_query_ids(data_dir)})
    start = batch_index * batch_size
    return shards[start : start + batch_size]


def plan_score_batches(
    data_dir: Path,
    *,
    batch_size: int,
    threads: int,
    scratch_dir: Path,
    max_query_protein_chains: int = 30,
    max_query_proper_ligand_chains: int = 30,
) -> dict[str, Any]:
    """Estimate mapped-hit work and greedily balance fixed-size score batches."""
    if (
        min(
            batch_size,
            threads,
            max_query_protein_chains,
            max_query_proper_ligand_chains,
        )
        < 1
    ):
        raise ValueError("batch size, threads, and chain limits must be positive")
    plan = _load_plan(
        data_dir,
        recheck_source=True,
        recheck_score_inputs=False,
    )
    annotation = data_dir / "index" / "annotation_table.parquet"
    entry_chains = data_dir / "index" / "entry_chains.parquet"
    alignment_paths = sorted(
        (data_dir / "alignments" / "search_db=holo").glob(
            "alignment_type=*/shard=*.parquet"
        )
    )
    if not alignment_paths:
        raise FileNotFoundError("no mapped alignment release shards found")

    alignment_paths_sql = ", ".join(f"'{path.as_posix()}'" for path in alignment_paths)
    hit_work_sql = dedent(
        f"""
        SELECT
            hit.query_entry,
            count(*)::BIGINT AS alignment_rows,
            sum(coalesce(target.proper_ligand_rows, 0))::DOUBLE
                AS target_ligand_hit_weight,
            sum(coalesce(target.scoreable_canonical_ligands, 0))::DOUBLE
                AS target_canonical_hit_weight
        FROM read_parquet([{alignment_paths_sql}], union_by_name = true) AS hit
        LEFT JOIN target_ligand_counts AS target
          ON hit.target_entry = target.entry_pdb_id
        GROUP BY hit.query_entry
        """
    )

    import duckdb

    scratch_dir.mkdir(exist_ok=True, parents=True)
    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
    work = connection.sql(
        dedent(
            f"""
            WITH holo_annotation AS (
                SELECT *
                FROM read_parquet('{annotation.as_posix()}')
                WHERE system_type = 'holo'
            ), system_ligand_counts AS (
                SELECT
                    entry_pdb_id,
                    system_id,
                    count(DISTINCT ligand_id) FILTER (
                        WHERE coalesce(ligand_is_proper, false)
                    )::BIGINT AS proper_ligand_chains
                FROM holo_annotation
                GROUP BY entry_pdb_id, system_id
            ), system_receptor_chains AS (
                SELECT DISTINCT
                    entry_pdb_id,
                    system_id,
                    unnest(system_protein_chains_asym_id)::VARCHAR
                        AS instance_chain
                FROM holo_annotation
            ), system_protein_counts AS (
                SELECT
                    receptor.entry_pdb_id,
                    receptor.system_id,
                    count(*) FILTER (
                        WHERE chains.chain_receptor_type = 'protein'
                    )::BIGINT AS protein_chains
                FROM system_receptor_chains AS receptor
                LEFT JOIN read_parquet('{entry_chains.as_posix()}') AS chains
                  ON receptor.entry_pdb_id = chains.entry_pdb_id
                 AND split_part(receptor.instance_chain, '.', 2)
                     = chains.chain_asym_id
                GROUP BY receptor.entry_pdb_id, receptor.system_id
            ), system_counts AS (
                SELECT
                    ligands.entry_pdb_id,
                    ligands.system_id,
                    coalesce(proteins.protein_chains, 0)::BIGINT
                        AS protein_chains,
                    ligands.proper_ligand_chains
                FROM system_ligand_counts AS ligands
                LEFT JOIN system_protein_counts AS proteins
                  USING (entry_pdb_id, system_id)
            ), eligible_systems AS (
                SELECT entry_pdb_id, system_id
                FROM system_counts
                WHERE protein_chains BETWEEN 1 AND {max_query_protein_chains}
                  AND proper_ligand_chains > 0
                  AND proper_ligand_chains <= {max_query_proper_ligand_chains}
            ), target_systems AS (
                SELECT entry_pdb_id, system_id
                FROM system_counts
                WHERE protein_chains > 0
                  AND proper_ligand_chains > 0
            ), target_ligand_counts AS (
                SELECT
                    annotation.entry_pdb_id,
                    count(*) FILTER (
                        WHERE coalesce(annotation.ligand_is_proper, false)
                    )::DOUBLE AS proper_ligand_rows,
                    count(DISTINCT annotation.ligand_asym_id) FILTER (
                        WHERE coalesce(annotation.ligand_is_proper, false)
                          AND coalesce(annotation.ligand_is_3d_score_able, false)
                    )::DOUBLE AS scoreable_canonical_ligands
                FROM holo_annotation AS annotation
                INNER JOIN target_systems USING (entry_pdb_id, system_id)
                GROUP BY annotation.entry_pdb_id
            ), query_ligand_counts AS (
                SELECT
                    annotation.entry_pdb_id,
                    count(*) FILTER (
                        WHERE coalesce(annotation.ligand_is_proper, false)
                    )::DOUBLE AS proper_ligand_rows,
                    count(DISTINCT annotation.ligand_asym_id) FILTER (
                        WHERE coalesce(annotation.ligand_is_proper, false)
                          AND coalesce(annotation.ligand_is_3d_score_able, false)
                    )::DOUBLE AS scoreable_canonical_ligands
                FROM holo_annotation AS annotation
                INNER JOIN eligible_systems USING (entry_pdb_id, system_id)
                GROUP BY annotation.entry_pdb_id
            ), eligible_queries AS (
                SELECT DISTINCT entry_pdb_id
                FROM eligible_systems
            ), hit_work_parts AS (
                {hit_work_sql}
            ), hit_work AS (
                SELECT
                    query_entry,
                    sum(alignment_rows)::BIGINT AS alignment_rows,
                    sum(target_ligand_hit_weight)::DOUBLE
                        AS target_ligand_hit_weight,
                    sum(target_canonical_hit_weight)::DOUBLE
                        AS target_canonical_hit_weight
                FROM hit_work_parts
                GROUP BY query_entry
            )
            SELECT
                query.pdb_id,
                query.shard,
                coalesce(hit.alignment_rows, 0)::BIGINT AS alignment_rows,
                coalesce(query_ligands.proper_ligand_rows, 0)
                    * coalesce(hit.target_ligand_hit_weight, 0)
                    AS ligand_pair_work,
                coalesce(query_ligands.scoreable_canonical_ligands, 0)
                    * coalesce(hit.target_canonical_hit_weight, 0)
                    AS canonical_3d_work
            FROM read_parquet(
                '{(data_dir / MANIFEST_RELATIVE).as_posix()}'
            ) AS query
            INNER JOIN eligible_queries AS eligible
              ON query.pdb_id = eligible.entry_pdb_id
            LEFT JOIN query_ligand_counts AS query_ligands
              ON query.pdb_id = query_ligands.entry_pdb_id
            LEFT JOIN hit_work AS hit
              ON query.pdb_id = hit.query_entry
            """
        )
    ).df()
    connection.close()
    dropped = dropped_query_ids(data_dir)
    if dropped:
        work = work[~work["pdb_id"].astype(str).isin(dropped)].copy()
    if work.empty:
        raise ValueError("no scoreable systems satisfy the configured chain limits")
    work["estimated_work"] = work["ligand_pair_work"] + 2.0 * work["canonical_3d_work"]

    assignment: dict[str, int] = {}
    loads: list[float] = []
    counts: list[int] = []
    for _, shard_work in work.groupby("shard", sort=True):
        shard_batch_count = math.ceil(len(shard_work) / batch_size)
        offset = len(loads)
        loads.extend([0.0] * shard_batch_count)
        counts.extend([0] * shard_batch_count)
        for row in shard_work.sort_values(
            ["estimated_work", "pdb_id"], ascending=[False, True]
        ).itertuples(index=False):
            eligible = [
                offset + local_index
                for local_index in range(shard_batch_count)
                if counts[offset + local_index] < batch_size
            ]
            selected = min(
                eligible, key=lambda index: (loads[index], counts[index], index)
            )
            assignment[str(row.pdb_id)] = selected
            loads[selected] += float(row.estimated_work)
            counts[selected] += 1
    batch_count = len(loads)
    work["score_batch_index"] = work["pdb_id"].map(assignment).astype("int32")
    work = work.sort_values(
        ["score_batch_index", "estimated_work", "pdb_id"],
        ascending=[True, False, True],
        ignore_index=True,
    )
    output = data_dir / SCORE_WORK_RELATIVE
    _atomic_parquet(work, output)
    plan.pop("score_max_protein_chains", None)
    plan.pop("score_max_ligand_chains", None)
    plan.update(
        {
            "score_batch_size": batch_size,
            "score_batch_count": batch_count,
            "score_max_query_protein_chains": max_query_protein_chains,
            "score_max_query_proper_ligand_chains": (max_query_proper_ligand_chains),
            "score_annotation": _source_signature(annotation),
            "score_batch_estimated_work_min": min(loads),
            "score_batch_estimated_work_max": max(loads),
            "score_work": _source_signature(output),
        }
    )
    _atomic_json(plan, data_dir / PLAN_RELATIVE)
    return {
        "query_count": len(work),
        "batch_size": batch_size,
        "batch_count": batch_count,
        "estimated_work_min": min(loads),
        "estimated_work_max": max(loads),
    }


def _ligand_3d_candidate_files(data_dir: Path) -> list[Path]:
    return sorted(
        (data_dir / "scores" / "ligand_3d_candidate_shards").glob("shard=*.parquet")
    )


def plan_ligand_3d_batches(
    data_dir: Path,
    *,
    batch_size: int,
    threads: int,
    scratch_dir: Path,
    memory_limit: str = "8GB",
) -> dict[str, Any]:
    """Freeze every scoreable canonical pair with positive pocket coverage."""
    started = perf_counter()
    if batch_size < 1 or threads < 1:
        raise ValueError("batch size and threads must be positive")
    plan = _load_plan(data_dir)
    protein_work = pd.read_parquet(data_dir / SCORE_WORK_RELATIVE, columns=["pdb_id"])
    expected_pdb_ids = set(protein_work["pdb_id"].astype(str)).difference(
        dropped_query_ids(data_dir)
    )
    expected_shards = {pdb_id[1:3] for pdb_id in expected_pdb_ids}
    candidate_files = _ligand_3d_candidate_files(data_dir)
    observed_shards = {path.stem.removeprefix("shard=") for path in candidate_files}
    missing_shards = expected_shards - observed_shards
    extra_shards = observed_shards - expected_shards
    if missing_shards or extra_shards:
        raise ValueError(
            "ligand 3D candidate shards do not cover planned score queries: "
            f"missing={sorted(missing_shards)[:10]}, "
            f"extra={sorted(extra_shards)[:10]}"
        )
    candidate_signatures = []
    pair_candidate_paths: list[Path] = []
    manifested_pdb_ids: set[str] = set()
    for path in candidate_files:
        shard = path.stem.removeprefix("shard=")
        manifest_path = path.with_suffix(".json")
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid ligand 3D candidate shard manifest: {manifest_path}"
            ) from exc
        schema = pq.read_schema(path)
        missing = sorted(
            set(schemas.LIGAND_3D_CANDIDATE_SCHEMA.names).difference(schema.names)
        )
        if missing:
            raise ValueError(f"candidate file {path} is missing columns {missing}")
        stat = path.stat()
        output_signature = {
            "path": str(path.resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "rows": pq.ParquetFile(path).metadata.num_rows,
        }
        inputs = manifest.get("inputs")
        if (
            manifest.get("shard") != shard
            or manifest.get("output") != output_signature
            or not isinstance(inputs, list)
        ):
            raise ValueError(f"stale ligand 3D candidate shard: {path}")
        pair_output = manifest.get("pair_output")
        pair_path = (
            data_dir
            / "scores"
            / "ligand_3d_pair_candidate_shards"
            / f"shard={shard}.parquet"
        )
        if not isinstance(pair_output, dict) or not pair_path.is_file():
            raise ValueError(f"missing compact ligand 3D candidate shard: {pair_path}")
        pair_schema = pq.read_schema(pair_path)
        pair_missing = sorted(
            set(schemas.LIGAND_3D_PAIR_CANDIDATE_SCHEMA.names).difference(
                pair_schema.names
            )
        )
        if pair_missing:
            raise ValueError(
                f"compact candidate file {pair_path} is missing {pair_missing}"
            )
        pair_stat = pair_path.stat()
        pair_output_signature = {
            "path": str(pair_path.resolve()),
            "size": pair_stat.st_size,
            "mtime_ns": pair_stat.st_mtime_ns,
            "rows": pq.ParquetFile(pair_path).metadata.num_rows,
        }
        if pair_output != pair_output_signature:
            raise ValueError(f"stale compact ligand 3D candidate shard: {pair_path}")
        pair_candidate_paths.append(pair_path)
        for source in inputs:
            if not isinstance(source, dict):
                raise ValueError(f"invalid candidate input in {manifest_path}")
            source_pdb_id = str(source.get("pdb_id", ""))
            if not source_pdb_id:
                raise ValueError(f"invalid candidate input in {manifest_path}")
            manifested_pdb_ids.add(source_pdb_id)
        candidate_signatures.append(
            {
                "shard": shard,
                **output_signature,
                "pair_path": pair_output_signature["path"],
                "pair_size": pair_output_signature["size"],
                "pair_mtime_ns": pair_output_signature["mtime_ns"],
                "pair_rows": pair_output_signature["rows"],
            }
        )
    LOG.info(
        "ligand 3D planning: validated %d candidate shards in %.1fs",
        len(candidate_files),
        perf_counter() - started,
    )
    if manifested_pdb_ids != expected_pdb_ids:
        raise ValueError(
            "collated ligand 3D candidates do not cover their query shards: "
            f"missing={sorted(expected_pdb_ids - manifested_pdb_ids)[:10]}, "
            f"extra={sorted(manifested_pdb_ids - expected_pdb_ids)[:10]}"
        )

    candidate_manifest = data_dir / LIGAND_3D_CANDIDATE_MANIFEST_RELATIVE
    work_path = data_dir / LIGAND_3D_WORK_RELATIVE
    candidate_frame = pd.DataFrame(candidate_signatures).sort_values(
        "shard", ignore_index=True
    )
    plan_is_current = False
    if candidate_manifest.is_file() and work_path.is_file():
        try:
            plan_is_current = (
                int(plan.get("ligand_3d_batch_size", -1)) == batch_size
                and _source_signature(candidate_manifest)
                == plan.get("ligand_3d_candidate_manifest")
                and _source_signature(work_path) == plan.get("ligand_3d_work")
                and pd.read_parquet(candidate_manifest).equals(candidate_frame)
            )
        except (OSError, TypeError, ValueError):
            plan_is_current = False

    if plan_is_current:
        return {
            "pair_count": int(plan["ligand_3d_pair_count"]),
            "batch_size": batch_size,
            "batch_count": int(plan["ligand_3d_batch_count"]),
            "estimated_work_min": float(plan["ligand_3d_estimated_work_min"]),
            "estimated_work_max": float(plan["ligand_3d_estimated_work_max"]),
            "cached": True,
        }

    scratch_dir.mkdir(exist_ok=True, parents=True)
    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
    connection.sql(f"SET memory_limit='{memory_limit}'")
    connection.sql("SET preserve_insertion_order=false")
    annotation = data_dir / "index" / "annotation_table.parquet"
    pair_paths_sql = ", ".join(f"'{path.as_posix()}'" for path in pair_candidate_paths)
    planned_path = scratch_dir / "ligand-3d-work.parquet"
    planned_path.unlink(missing_ok=True)
    connection.sql(
        dedent(
            f"""
            COPY (
                WITH ligand_sizes AS (
                    SELECT
                        entry_pdb_id AS entry,
                        ligand_asym_id AS asym_id,
                        max(greatest(coalesce(ligand_num_heavy_atoms, 1), 1))::BIGINT
                            AS heavy_atoms
                    FROM read_parquet('{annotation.as_posix()}')
                    WHERE system_type = 'holo'
                      AND coalesce(ligand_is_proper, false)
                      AND coalesce(ligand_is_3d_score_able, false)
                    GROUP BY entry_pdb_id, ligand_asym_id
                ), scoreable_pairs AS (
                    SELECT candidate_pairs.*,
                           (query.heavy_atoms * target.heavy_atoms)::BIGINT
                               AS estimated_work
                    FROM read_parquet([{pair_paths_sql}]) AS candidate_pairs
                    INNER JOIN ligand_sizes AS query
                      ON candidate_pairs.query_entry = query.entry
                     AND candidate_pairs.query_ligand_asym_id = query.asym_id
                    INNER JOIN ligand_sizes AS target
                      ON candidate_pairs.target_entry = target.entry
                     AND candidate_pairs.target_ligand_asym_id = target.asym_id
                    WHERE candidate_pairs.pocket_qcov > 0
                )
                SELECT
                    query_entry,
                    query_ligand_asym_id,
                    target_entry,
                    target_ligand_asym_id,
                    estimated_work,
                    floor((row_number() OVER (
                        ORDER BY
                            query_entry,
                            query_ligand_asym_id,
                            estimated_work DESC,
                            target_entry,
                            target_ligand_asym_id
                    ) - 1) / {batch_size})::INTEGER AS ligand_3d_batch_index
                FROM scoreable_pairs
            ) TO '{planned_path.as_posix()}' (
                FORMAT PARQUET,
                COMPRESSION ZSTD,
                ROW_GROUP_SIZE {max(2_048, batch_size)}
            )
            """
        )
    )

    counts = connection.sql(
        f"""
        SELECT
            count(*)::BIGINT AS pair_count,
            count(DISTINCT substr(query_entry, 2, 2))::BIGINT AS query_shards,
            coalesce(max(ligand_3d_batch_index) + 1, 0)::BIGINT AS batch_count
        FROM read_parquet('{planned_path.as_posix()}')
        """
    ).fetchone()
    load_bounds = connection.sql(
        f"""
        SELECT
            coalesce(min(batch_work), 0)::DOUBLE,
            coalesce(max(batch_work), 0)::DOUBLE,
            coalesce(max(batch_rows), 0)::BIGINT
        FROM (
            SELECT
                ligand_3d_batch_index,
                sum(estimated_work)::DOUBLE AS batch_work,
                count(*)::BIGINT AS batch_rows
            FROM read_parquet('{planned_path.as_posix()}')
            GROUP BY ligand_3d_batch_index
        )
        """
    ).fetchone()
    connection.close()
    if counts is None or load_bounds is None:
        raise RuntimeError("failed to validate ligand 3D work plan")
    pair_count, query_shard_count, batch_count = (int(value) for value in counts)
    estimated_work_min, estimated_work_max, maximum_batch_rows = load_bounds
    if maximum_batch_rows > batch_size:
        raise ValueError(
            f"ligand 3D batch exceeds size {batch_size}: {maximum_batch_rows}"
        )
    work_path.parent.mkdir(exist_ok=True, parents=True)
    install = work_path.with_suffix(work_path.suffix + ".tmp")
    copyfile(planned_path, install)
    install.replace(work_path)
    planned_path.unlink(missing_ok=True)
    _atomic_parquet(candidate_frame, candidate_manifest)
    plan.update(
        {
            "ligand_3d_batch_size": batch_size,
            "ligand_3d_batch_count": batch_count,
            "ligand_3d_pair_count": pair_count,
            "ligand_3d_plan_complete": True,
            "ligand_3d_query_shard_count": query_shard_count,
            "ligand_3d_estimated_work_min": estimated_work_min,
            "ligand_3d_estimated_work_max": estimated_work_max,
            "ligand_3d_candidate_manifest": _source_signature(candidate_manifest),
            "ligand_3d_work": _source_signature(work_path),
        }
    )
    for key in list(plan):
        if key.startswith("ligand_3d_retry_"):
            plan.pop(key)
    _atomic_json(plan, data_dir / PLAN_RELATIVE)
    LOG.info(
        "ligand 3D planning: wrote %d pairs in %d batches (total %.1fs)",
        pair_count,
        batch_count,
        perf_counter() - started,
    )
    return {
        "pair_count": pair_count,
        "batch_size": batch_size,
        "batch_count": batch_count,
        "estimated_work_min": estimated_work_min,
        "estimated_work_max": estimated_work_max,
        "cached": False,
    }


def _ligand_3d_batch(data_dir: Path, batch_index: int, batch_size: int) -> pd.DataFrame:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    plan = _load_plan(data_dir)
    if batch_size != int(plan.get("ligand_3d_batch_size", -1)):
        raise ValueError(
            f"ligand 3D batch size is {plan.get('ligand_3d_batch_size')}, "
            f"got {batch_size}"
        )
    work_path = data_dir / LIGAND_3D_WORK_RELATIVE
    if _source_signature(work_path) != plan.get("ligand_3d_work"):
        raise ValueError("ligand 3D work manifest changed after planning")
    return pd.read_parquet(
        work_path,
        filters=[("ligand_3d_batch_index", "==", batch_index)],
    )


def plan_ligand_3d_retries(
    data_dir: Path,
    *,
    batch_size: int,
) -> dict[str, Any]:
    """Split only missing large 3D batches into small resumable retry units."""
    if batch_size < 1:
        raise ValueError("ligand 3D retry batch size must be positive")
    plan = _load_plan(data_dir)
    if not plan.get("ligand_3d_plan_complete", False):
        raise ValueError("complete ligand 3D candidate planning before retry planning")
    work_path = data_dir / LIGAND_3D_WORK_RELATIVE
    if _source_signature(work_path) != plan.get("ligand_3d_work"):
        raise ValueError("ligand 3D work manifest changed after planning")
    batch_count = int(plan["ligand_3d_batch_count"])
    output_dir = data_dir / "scores" / "ligand_3d_pairs"
    missing_original = [
        index
        for index in range(batch_count)
        if not (output_dir / f"{index}.parquet").is_file()
    ]
    if missing_original:
        work = pd.read_parquet(
            work_path,
            filters=[("ligand_3d_batch_index", "in", missing_original)],
        )
    else:
        work = pd.DataFrame()
    parts: list[pd.DataFrame] = []
    retry_index = 0
    if not work.empty:
        order = [
            "query_entry",
            "query_ligand_asym_id",
            "estimated_work",
            "target_entry",
            "target_ligand_asym_id",
        ]
        ascending = [True, True, False, True, True]
        for original_index, group in work.groupby("ligand_3d_batch_index", sort=True):
            group = group.sort_values(order, ascending=ascending, ignore_index=True)
            for start in range(0, len(group), batch_size):
                part = group.iloc[start : start + batch_size].copy()
                part["original_ligand_3d_batch_index"] = int(original_index)
                part["retry_batch_index"] = retry_index
                parts.append(part)
                retry_index += 1
    retry_work = (
        pd.concat(parts, ignore_index=True)
        if parts
        else pd.DataFrame(
            columns=[
                "query_entry",
                "query_ligand_asym_id",
                "target_entry",
                "target_ligand_asym_id",
                "estimated_work",
                "ligand_3d_batch_index",
                "original_ligand_3d_batch_index",
                "retry_batch_index",
            ]
        )
    )
    retry_path = data_dir / LIGAND_3D_RETRY_WORK_RELATIVE
    _atomic_parquet(retry_work, retry_path)
    plan.update(
        {
            "ligand_3d_retry_batch_size": batch_size,
            "ligand_3d_retry_batch_count": retry_index,
            "ligand_3d_retry_original_batches": missing_original,
            "ligand_3d_retry_work": _source_signature(retry_path),
        }
    )
    _atomic_json(plan, data_dir / PLAN_RELATIVE)
    return {
        "status": "planned",
        "missing_original_batch_count": len(missing_original),
        "retry_pair_count": len(retry_work),
        "retry_batch_size": batch_size,
        "retry_batch_count": retry_index,
    }


def _ligand_3d_retry_batch(
    data_dir: Path, batch_index: int, batch_size: int
) -> pd.DataFrame:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch index must be non-negative and batch size positive")
    plan = _load_plan(data_dir)
    if batch_size != int(plan.get("ligand_3d_retry_batch_size", -1)):
        raise ValueError(
            f"ligand 3D retry batch size is "
            f"{plan.get('ligand_3d_retry_batch_size')}, got {batch_size}"
        )
    retry_path = data_dir / LIGAND_3D_RETRY_WORK_RELATIVE
    if _source_signature(retry_path) != plan.get("ligand_3d_retry_work"):
        raise ValueError("ligand 3D retry work changed after planning")
    return pd.read_parquet(
        retry_path,
        filters=[("retry_batch_index", "==", batch_index)],
    )


def finalize_ligand_3d_retries(
    data_dir: Path,
) -> dict[str, Any]:
    """Reassemble complete fine-grained retries into their original batches."""
    plan = _load_plan(data_dir)
    retry_path = data_dir / LIGAND_3D_RETRY_WORK_RELATIVE
    if _source_signature(retry_path) != plan.get("ligand_3d_retry_work"):
        raise ValueError("ligand 3D retry work changed after planning")
    retry_work = pd.read_parquet(retry_path)
    pair_columns = [
        "query_entry",
        "query_ligand_asym_id",
        "target_entry",
        "target_ligand_asym_id",
    ]
    retry_output_dir = data_dir / "scores" / "ligand_3d_pair_retries"
    valid_outputs: dict[int, pd.DataFrame] = {}
    incomplete: list[int] = []
    for retry_index, expected in retry_work.groupby("retry_batch_index", sort=True):
        retry_index = int(retry_index)
        path = retry_output_dir / f"{retry_index}.parquet"
        try:
            schema = pq.read_schema(path)
            observed = pd.read_parquet(path)
            expected_pairs = pd.MultiIndex.from_frame(expected[pair_columns])
            observed_pairs = pd.MultiIndex.from_frame(observed[pair_columns])
            valid = set(schemas.LIGAND_3D_SCORE_SCHEMA.names).issubset(
                schema.names
            ) and observed_pairs.equals(expected_pairs)
        except (OSError, ValueError):
            valid = False
            observed = pd.DataFrame()
        if valid:
            valid_outputs[retry_index] = observed
        else:
            incomplete.append(retry_index)
    if incomplete:
        raise ValueError(f"ligand 3D retries remain incomplete: {incomplete[:20]}")

    canonical_dir = data_dir / "scores" / "ligand_3d_pairs"
    canonical_dir.mkdir(exist_ok=True, parents=True)
    for original_index, expected_original in retry_work.groupby(
        "original_ligand_3d_batch_index", sort=True
    ):
        scored_parts: list[pd.DataFrame] = []
        for retry_index, _expected in expected_original.groupby(
            "retry_batch_index", sort=True
        ):
            retry_index = int(retry_index)
            scored_parts.append(valid_outputs[retry_index])
        combined = pd.concat(scored_parts, ignore_index=True)
        expected_pairs = pd.MultiIndex.from_frame(expected_original[pair_columns])
        combined_pairs = pd.MultiIndex.from_frame(combined[pair_columns])
        if len(combined) != len(expected_original) or set(combined_pairs) != set(
            expected_pairs
        ):
            raise ValueError(
                f"retry collation changed canonical pairs for batch {original_index}"
            )
        combined = expected_original[pair_columns].merge(
            combined,
            on=pair_columns,
            how="left",
            validate="one_to_one",
        )
        output = canonical_dir / f"{int(original_index)}.parquet"
        temporary = output.with_suffix(".tmp.parquet")
        combined.to_parquet(
            temporary,
            schema=schemas.LIGAND_3D_SCORE_SCHEMA,
            index=False,
        )
        temporary.replace(output)

    plan["ligand_3d_retry_complete"] = True
    _atomic_json(plan, data_dir / PLAN_RELATIVE)
    return {
        "status": "complete",
        "retry_batch_count": int(plan["ligand_3d_retry_batch_count"]),
        "completed_retry_batch_count": len(valid_outputs),
    }


def _ligand_3d_shard_batch(
    data_dir: Path, batch_index: int, batch_size: int
) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    _load_plan(data_dir)
    shards = sorted({pdb_id[1:3] for pdb_id in active_scoring_query_ids(data_dir)})
    start = batch_index * batch_size
    return shards[start : start + batch_size]


def finalize_ligand_archives(data_dir: Path) -> dict[str, Any]:
    """Validate packed canonical SDF coverage against the collated index."""
    archives = sorted((data_dir / "ligand_archives").glob("*.parquet"))
    expected_shards = sorted(
        path.name for path in (data_dir / "raw_entries").iterdir() if path.is_dir()
    )
    observed_shards = [path.stem for path in archives]
    if observed_shards != expected_shards:
        raise ValueError(
            "canonical ligand archive shards do not match raw entries: "
            f"missing={sorted(set(expected_shards) - set(observed_shards))[:10]}, "
            f"extra={sorted(set(observed_shards) - set(expected_shards))[:10]}"
        )

    import duckdb

    paths_sql = ", ".join(f"'{path.as_posix()}'" for path in archives)
    annotation = data_dir / "index" / "annotation_table.parquet"
    connection = duckdb.connect()
    archive_counts = connection.sql(
        dedent(
            f"""
            SELECT
                count(*)::BIGINT,
                count(DISTINCT (pdb_id, ligand_asym_id))::BIGINT,
                count(*) FILTER (WHERE octet_length(sdf) = 0)::BIGINT
            FROM read_parquet([{paths_sql}], union_by_name=true)
            """
        )
    ).fetchone()
    if archive_counts is None:
        connection.close()
        raise RuntimeError("failed to validate canonical ligand archives")
    row_count, unique_count, empty_count = archive_counts
    missing_result = connection.sql(
        dedent(
            f"""
            WITH expected AS (
                SELECT DISTINCT
                    entry_pdb_id AS pdb_id,
                    ligand_asym_id
                FROM read_parquet('{annotation.as_posix()}')
                WHERE ligand_asym_id IS NOT NULL
            ), packed AS (
                SELECT pdb_id, ligand_asym_id
                FROM read_parquet([{paths_sql}], union_by_name=true)
            )
            SELECT count(*)::BIGINT
            FROM expected
            ANTI JOIN packed USING (pdb_id, ligand_asym_id)
            """
        )
    ).fetchone()
    if missing_result is None:
        connection.close()
        raise RuntimeError("failed to validate canonical ligand archive coverage")
    missing_count = missing_result[0]
    connection.close()
    if row_count != unique_count:
        raise ValueError(
            f"packed canonical ligand keys are not unique: {row_count} rows, "
            f"{unique_count} unique keys"
        )
    if empty_count or missing_count:
        raise ValueError(
            "invalid packed canonical ligands: "
            f"empty_sdfs={empty_count}, missing_index_ligands={missing_count}"
        )
    report = {
        "status": "complete",
        "shard_count": len(archives),
        "ligand_count": row_count,
        "compressed_bytes": sum(path.stat().st_size for path in archives),
    }
    _atomic_json(report, data_dir / LIGAND_ARCHIVE_MANIFEST_RELATIVE)
    return report


def make_foldseek_input_manifest(data_dir: Path, cif_root: Path) -> Path:
    """List only planned protein-entry CIFs without recursively walking PDB."""
    _load_plan(data_dir, recheck_source=True)
    pdb_ids = pd.read_parquet(data_dir / MANIFEST_RELATIVE, columns=["pdb_id"])[
        "pdb_id"
    ].astype(str)
    output = data_dir / FOLDSEEK_INPUT_RELATIVE
    output.parent.mkdir(exist_ok=True, parents=True)
    temporary = output.with_suffix(".tmp.tsv")
    with temporary.open("w") as handle:
        for pdb_id in pdb_ids:
            cif_path = (
                cif_root
                / pdb_id[1:3]
                / f"pdb_0000{pdb_id}"
                / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
            )
            handle.write(f"{cif_path}\n")
    temporary.replace(output)
    return output


def _shard_batch(data_dir: Path, batch_index: int, batch_size: int) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    _load_plan(data_dir)
    shards = sorted(
        pd.read_parquet(data_dir / MANIFEST_RELATIVE, columns=["shard"])[
            "shard"
        ].unique()
    )
    start = batch_index * batch_size
    return [str(value) for value in shards[start : start + batch_size]]


def _entry_shard_batch(data_dir: Path, batch_index: int, batch_size: int) -> list[str]:
    """Select raw-entry shards, including entries without protein receptors."""
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    shards = sorted(
        path.name for path in (data_dir / "raw_entries").iterdir() if path.is_dir()
    )
    start = batch_index * batch_size
    return shards[start : start + batch_size]


def _score_partition_batch(batch_index: int, batch_size: int) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    partitions = [chunk[0] for chunk in tasks.scatter_collate_partitions()]
    start = batch_index * batch_size
    return partitions[start : start + batch_size]


def _cluster_parameters(
    *, metrics: list[str] | None, thresholds: list[int] | None
) -> tuple[list[str], list[int]]:
    selected_metrics = list(dict.fromkeys(metrics or DEFAULT_CLUSTER_METRICS))
    selected_thresholds = sorted(
        set(thresholds or DEFAULT_CLUSTER_THRESHOLDS), reverse=True
    )
    if not selected_metrics:
        raise ValueError("at least one clustering metric is required")
    unsupported = sorted(set(selected_metrics).difference(DEFAULT_CLUSTER_METRICS))
    if unsupported:
        if set(unsupported).intersection({"shape", "color", "sucos_shape"}):
            raise ValueError(
                "raw 3D ligand diagnostics cannot be clustered; use "
                "sucos_shape_pocket_qcov"
            )
        raise ValueError(f"unsupported release clustering metrics: {unsupported}")
    if not selected_thresholds or not all(
        0 <= threshold <= 100 for threshold in selected_thresholds
    ):
        raise ValueError("clustering thresholds must be in [0, 100]")
    return selected_metrics, selected_thresholds


def plan_clustering(
    data_dir: Path,
    *,
    metrics: list[str] | None = None,
    thresholds: list[int] | None = None,
    source_batch_size: int = 20,
    community_batch_size: int = 1,
    symmetric_bucket_count: int = clusters.SYMMETRIC_EDGE_BUCKET_COUNT,
) -> dict[str, Any]:
    """Plan reciprocal-minimum components and centroid clusters."""
    if source_batch_size < 1 or community_batch_size < 1:
        raise ValueError("clustering batch sizes must be positive")
    selected_metrics, selected_thresholds = _cluster_parameters(
        metrics=metrics,
        thresholds=thresholds,
    )
    if any(metric != "tanimoto_similarity_ecfp4_1024" for metric in selected_metrics):
        clusters.prepare_component_node_universe(data_dir)
    symmetric_plan = clusters.prepare_symmetric_edge_plan(
        data_dir=data_dir,
        metrics=selected_metrics,
        source_batch_size=source_batch_size,
        bucket_count=symmetric_bucket_count,
    )
    cluster_plan_path = data_dir / "ligand_clusters" / "clustering_plan.json"
    cluster_selection = {
        "version": 1,
        "metrics": selected_metrics,
        "thresholds": selected_thresholds,
        "symmetric_plan_hash": symmetric_plan["plan_hash"],
    }
    prior_selection: dict[str, Any] | None = None
    if cluster_plan_path.is_file():
        try:
            prior_selection = json.loads(cluster_plan_path.read_text())
        except (OSError, ValueError, json.JSONDecodeError):
            pass
    if prior_selection != cluster_selection:
        for obsolete in [
            data_dir / "ligand_clusters/cluster=components",
            data_dir / "ligand_clusters/cluster=communities",
        ]:
            if obsolete.exists():
                rmtree(obsolete)
        sampling_root = data_dir / "ligand_sampling" / "directed_set_cover"
        for output_dir in sampling_root.glob("metric=*"):
            rmtree(output_dir)
        for diagnostics in [
            data_dir / "ligand_clusters/stats.parquet",
            data_dir / "ligand_clusters/stats.json",
        ]:
            diagnostics.unlink(missing_ok=True)
        _atomic_json(cluster_selection, cluster_plan_path)
    else:
        # Old directed artifacts are never valid under the current selection.
        for obsolete in [
            data_dir / "ligand_clusters/cluster=components/directed=True",
            data_dir / "ligand_clusters/cluster=communities/directed=True",
        ]:
            if obsolete.exists():
                rmtree(obsolete)
    symmetric_shard_count = len(selected_metrics) * symmetric_bucket_count
    community_task_count = len(selected_metrics) * len(selected_thresholds)
    return {
        "status": "planned",
        "metrics": selected_metrics,
        "thresholds": selected_thresholds,
        "source_batch_size": source_batch_size,
        "symmetric_bucket_count": symmetric_bucket_count,
        "symmetric_fragment_batch_count": len(symmetric_plan["batches"]),
        "symmetric_edge_shard_count": symmetric_shard_count,
        "component_reduction_batch_count": symmetric_shard_count,
        "community_task_count": community_task_count,
        "community_batch_size": community_batch_size,
        "community_batch_count": math.ceil(community_task_count / community_batch_size),
        "directed_cover_task_count": community_task_count,
        "directed_cover_batch_count": math.ceil(
            community_task_count / community_batch_size
        ),
    }


def summarize_clustering_artifacts(
    data_dir: Path,
    *,
    metrics: list[str] | None = None,
    thresholds: list[int] | None = None,
) -> dict[str, Any]:
    """Validate published ligand clusters and write compact diagnostic statistics."""
    selected_metrics, selected_thresholds = _cluster_parameters(
        metrics=metrics,
        thresholds=thresholds,
    )
    output_dir = data_dir / "ligand_clusters"
    sampling_dir = data_dir / "ligand_sampling" / "directed_set_cover"
    artifacts: list[tuple[str, int, str, bool, Path]] = []
    for metric in selected_metrics:
        for threshold in selected_thresholds:
            for cluster in ["components", "communities"]:
                artifacts.append(
                    (
                        metric,
                        threshold,
                        cluster,
                        False,
                        output_dir
                        / f"cluster={cluster}"
                        / "directed=False"
                        / f"metric={metric}"
                        / f"threshold={threshold}.parquet",
                    )
                )
            artifacts.append(
                (
                    metric,
                    threshold,
                    "directed_set_cover",
                    True,
                    sampling_dir
                    / f"metric={metric}"
                    / f"threshold={threshold}.parquet",
                )
            )
    expected_paths = {path for *_, path in artifacts}
    rows: list[dict[str, Any]] = []
    issues: list[str] = []
    started = perf_counter()
    for index, (metric, threshold, cluster, directed, path) in enumerate(
        artifacts, start=1
    ):
        if not path.is_file():
            issues.append(f"missing cluster artifact: {path}")
            continue
        frame = pd.read_parquet(path, columns=["ligand_id", "label"])
        duplicate_nodes = int(frame["ligand_id"].duplicated().sum())
        null_labels = int(frame["label"].isna().sum())
        sizes = frame["label"].dropna().astype(str).value_counts()
        rows.append(
            {
                "metric": metric,
                "threshold": threshold,
                "cluster": cluster,
                "directed": directed,
                "node_count": len(frame),
                "cluster_count": len(sizes),
                "singleton_cluster_count": int((sizes == 1).sum()),
                "largest_cluster_size": int(sizes.max()) if len(sizes) else 0,
                "median_cluster_size": float(sizes.median()) if len(sizes) else 0.0,
                "p95_cluster_size": (
                    float(sizes.quantile(0.95)) if len(sizes) else 0.0
                ),
                "duplicate_node_count": duplicate_nodes,
                "null_label_count": null_labels,
                "file_size_bytes": path.stat().st_size,
            }
        )
        if duplicate_nodes:
            issues.append(f"{path} contains {duplicate_nodes} duplicate ligand IDs")
        if null_labels:
            issues.append(f"{path} contains {null_labels} null labels")
        if index % 10 == 0 or index == len(artifacts):
            elapsed = perf_counter() - started
            rate = index / elapsed
            LOG.info(
                "cluster statistics progress: "
                f"validated={index}/{len(artifacts)} rate={rate:.2f}/s "
                f"eta_seconds={(len(artifacts) - index) / rate:.1f}"
            )

    stats = pd.DataFrame(rows)
    if len(stats) == len(artifacts):
        for metric in selected_metrics:
            metric_stats = stats[stats["metric"].eq(metric)]
            for cluster, directed in [
                ("components", False),
                ("communities", False),
            ]:
                subset = metric_stats[
                    metric_stats["cluster"].eq(cluster)
                    & metric_stats["directed"].eq(directed)
                ].set_index("threshold")
                node_counts = set(subset["node_count"].astype(int))
                if len(node_counts) != 1:
                    issues.append(
                        f"{metric} {cluster} directed={directed} has inconsistent "
                        f"node counts: {sorted(node_counts)}"
                    )
                if cluster != "components":
                    continue
                previous_count: int | None = None
                for threshold in selected_thresholds:
                    cluster_count = int(subset.loc[threshold, "cluster_count"])
                    if previous_count is not None and cluster_count > previous_count:
                        issues.append(
                            f"{metric} {cluster} directed={directed} gains clusters "
                            f"when lowering threshold to {threshold}: "
                            f"{previous_count} -> {cluster_count}"
                        )
                    previous_count = cluster_count
            for threshold in selected_thresholds:
                at_threshold = metric_stats[metric_stats["threshold"].eq(threshold)]
                threshold_node_counts = set(at_threshold["node_count"].astype(int))
                if len(threshold_node_counts) != 1:
                    issues.append(
                        f"{metric} threshold={threshold} has inconsistent node "
                        f"counts across cluster types: {sorted(threshold_node_counts)}"
                    )
                counts = at_threshold.set_index("cluster")["cluster_count"]
                if int(counts["communities"]) < int(counts["components"]):
                    issues.append(
                        f"{metric} threshold={threshold} has fewer centroid "
                        "clusters than connected components"
                    )
    published_paths = set(
        output_dir.glob("cluster=*/directed=*/metric=*/threshold=*.parquet")
    ) | set(sampling_dir.glob("metric=*/threshold=*.parquet"))
    unexpected_artifacts = sorted(
        path.relative_to(data_dir).as_posix()
        for path in published_paths.difference(expected_paths)
    )
    if unexpected_artifacts:
        issues.append(
            f"unexpected published cluster artifacts: {unexpected_artifacts[:20]}"
        )

    selected_metric_set = set(selected_metrics)
    published_metrics = {
        path.parent.name.removeprefix("metric=") for path in published_paths
    }
    unexpected_metrics = sorted(published_metrics.difference(selected_metric_set))
    if unexpected_metrics:
        issues.append(f"unexpected published cluster metrics: {unexpected_metrics}")

    stats_path = output_dir / "stats.parquet"
    _atomic_parquet(stats, stats_path)
    report = {
        "status": "complete" if not issues else "invalid",
        "metric_count": len(selected_metrics),
        "thresholds": selected_thresholds,
        "artifact_count": len(stats),
        "expected_artifact_count": len(artifacts),
        "unexpected_artifacts": unexpected_artifacts,
        "unexpected_metrics": unexpected_metrics,
        "issues": issues,
        "stats_path": stats_path.relative_to(data_dir).as_posix(),
        "elapsed_seconds": perf_counter() - started,
    }
    _atomic_json(report, output_dir / "stats.json")
    if issues:
        raise ValueError(f"invalid clustering artifacts: {issues[:10]}")
    return report


def _component_source_batch(
    data_dir: Path,
    *,
    metrics: list[str],
    batch_index: int,
    batch_size: int,
) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    batches = tasks.scatter_component_reduction_sources(
        data_dir=data_dir,
        metrics=metrics,
        batch_size=batch_size,
    )
    return batches[batch_index] if batch_index < len(batches) else []


def _symmetric_fragment_batch(
    data_dir: Path,
    *,
    batch_index: int,
    batch_size: int,
) -> list[dict[str, Any]]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    plan = clusters.load_symmetric_edge_plan(data_dir)
    start = batch_index * batch_size
    return cast(list[dict[str, Any]], plan["batches"][start : start + batch_size])


def _symmetric_edge_shard_batch(
    data_dir: Path,
    *,
    batch_index: int,
    batch_size: int,
) -> list[tuple[str, int]]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    plan = clusters.load_symmetric_edge_plan(data_dir)
    work = [
        (str(metric), bucket)
        for metric in plan["metrics"]
        for bucket in range(int(plan["bucket_count"]))
    ]
    start = batch_index * batch_size
    return work[start : start + batch_size]


def _community_batch(
    *,
    metrics: list[str],
    thresholds: list[int],
    batch_index: int,
    batch_size: int,
) -> list[tuple[str, int]]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    work = [(metric, threshold) for metric in metrics for threshold in thresholds]
    start = batch_index * batch_size
    return work[start : start + batch_size]


def _scoring_config(
    data_dir: Path,
    max_seqs: int,
    *,
    max_query_protein_chains: int = 30,
    max_query_proper_ligand_chains: int = 30,
) -> Any:
    return config.get_config(
        config={
            "scorer": {
                "sub_databases": "holo",
                "max_query_protein_chains": max_query_protein_chains,
                "max_query_proper_ligand_chains": (max_query_proper_ligand_chains),
            },
            "foldseek": {"max_seqs": max_seqs, "min_seq_id": 0.0},
            "mmseqs": {"max_seqs": max_seqs, "min_seq_id": 0.0},
        },
        cached=False,
    )


def _database_query_ids(database: Path, alignment_type: str) -> set[str]:
    identifiers = databases.database_identifiers(database)
    if alignment_type == "foldseek":
        return {identifier.replace("pdb_0000", "", 1)[:4] for identifier in identifiers}
    return {identifier.split("_", maxsplit=1)[0] for identifier in identifiers}


def finalize_alignment_artifacts(data_dir: Path) -> dict[str, Any]:
    """Validate compact release shards from their atomic mapping manifests."""
    plan = _load_plan(data_dir, recheck_source=True)
    expected_queries = set(
        pd.read_parquet(data_dir / MANIFEST_RELATIVE, columns=["pdb_id"])[
            "pdb_id"
        ].astype(str)
    )
    counts: dict[str, int] = {}
    missing: dict[str, list[str]] = {}
    cluster_manifests: dict[str, dict[str, Any]] = {}
    skipped_query_details: dict[str, dict[str, Any]] = {}
    backend_queries: dict[str, set[str]] = {}
    for alignment_type in plan["alignment_types"]:
        root = data_dir / "dbs" / "subdbs" / f"holo_{alignment_type}"
        cluster_manifest_path = root / "exact_cluster.json"
        if not cluster_manifest_path.is_file():
            raise FileNotFoundError(
                f"missing exact-cluster manifest: {cluster_manifest_path}"
            )
        cluster_manifest = json.loads(cluster_manifest_path.read_text())
        if cluster_manifest.get("alignment_type") != alignment_type:
            raise ValueError(f"invalid exact-cluster manifest: {cluster_manifest_path}")
        clustering = plan["target_clustering"]
        for key in ["identity", "coverage", "coverage_mode"]:
            if cluster_manifest.get(key) != clustering[key]:
                raise ValueError(
                    f"{alignment_type} target clustering does not match plan: "
                    f"{key}={cluster_manifest.get(key)!r}"
                )
        source_index = cluster_manifest.get("source_index")
        if not isinstance(source_index, dict):
            raise ValueError(
                f"invalid exact-cluster source signature: {cluster_manifest_path}"
            )
        source_path = root / str(source_index["name"])
        source_stat = source_path.stat()
        observed_source = {
            "name": source_path.name,
            "size": source_stat.st_size,
            "mtime_ns": source_stat.st_mtime_ns,
        }
        if observed_source != source_index:
            raise ValueError(
                f"{alignment_type} target DB changed after exact clustering"
            )
        for artifact_key in ["search_target", "conversion_target"]:
            artifact = root / str(cluster_manifest[artifact_key])
            if not artifact.with_suffix(".dbtype").is_file():
                raise FileNotFoundError(
                    f"missing {alignment_type} {artifact_key}: {artifact}"
                )
        if databases._has_external_database_links(root):
            raise ValueError(
                f"{alignment_type} target database contains external or broken links"
            )
        cluster_manifests[alignment_type] = cluster_manifest
        backend_queries[alignment_type] = _database_query_ids(
            root / root.name, alignment_type
        ).intersection(expected_queries)

    manifest_root = data_dir / "alignments" / "manifests"
    manifest_paths = sorted(manifest_root.glob("shard=*.json"))
    expected_shards = {pdb_id[1:3] for pdb_id in expected_queries}
    observed_manifest_shards = {
        path.stem.removeprefix("shard=") for path in manifest_paths
    }
    absent_manifests = sorted(expected_shards.difference(observed_manifest_shards))
    if absent_manifests:
        missing["mapping_manifests"] = absent_manifests

    lookup_signature = tasks._completed_alignment_chain_lookup(data_dir)
    raw_queries_by_backend: dict[str, set[str]] = {
        alignment_type: set() for alignment_type in plan["alignment_types"]
    }
    mapped_queries_by_backend: dict[str, set[str]] = {
        alignment_type: set() for alignment_type in plan["alignment_types"]
    }
    release_shards_by_backend: dict[str, set[str]] = {
        alignment_type: set() for alignment_type in plan["alignment_types"]
    }
    invalid_manifests: list[str] = []
    for manifest_path in manifest_paths:
        shard = manifest_path.stem.removeprefix("shard=")
        try:
            shard_manifest = json.loads(manifest_path.read_text())
        except (OSError, TypeError, ValueError):
            invalid_manifests.append(shard)
            continue
        inputs = shard_manifest.get("inputs")
        outputs = shard_manifest.get("outputs")
        shard_skipped = shard_manifest.get("skipped_queries")
        if (
            shard_manifest.get("shard") != shard
            or shard_manifest.get("alignment_chain_lookup") != lookup_signature
            or not isinstance(inputs, dict)
            or not isinstance(outputs, dict)
            or not isinstance(shard_skipped, dict)
        ):
            invalid_manifests.append(shard)
            continue
        skipped_query_details.update(
            {str(key): value for key, value in shard_skipped.items()}
        )
        for alignment_type in plan["alignment_types"]:
            signatures = inputs.get(alignment_type)
            if not isinstance(signatures, list):
                invalid_manifests.append(shard)
                continue
            raw_queries = {
                Path(str(signature["name"])).stem for signature in signatures
            }
            raw_queries_by_backend[alignment_type].update(raw_queries)
            mapped_queries_by_backend[alignment_type].update(
                raw_queries.difference(shard_skipped)
            )
            expected_output = outputs.get(alignment_type)
            if not signatures:
                if expected_output is not None:
                    invalid_manifests.append(shard)
                continue
            release = tasks._alignment_release_path(
                data_dir=data_dir,
                search_db="holo",
                alignment_type=alignment_type,
                shard=shard,
            )
            if not isinstance(expected_output, dict) or not release.is_file():
                invalid_manifests.append(shard)
                continue
            stat = release.stat()
            if expected_output != {
                "name": release.name,
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }:
                invalid_manifests.append(shard)
                continue
            columns = set(pq.read_schema(release).names)
            if not schemas.mapped_alignment_schema_is_current(
                columns, alignment_type=alignment_type
            ):
                invalid_manifests.append(shard)
                continue
            release_shards_by_backend[alignment_type].add(shard)
    if invalid_manifests:
        missing["invalid_mapping_manifests"] = sorted(set(invalid_manifests))[:100]

    observed_raw_all: set[str] = set()
    observed_mapped_all: set[str] = set()
    for alignment_type in plan["alignment_types"]:
        raw_queries = raw_queries_by_backend[alignment_type]
        mapped_queries = mapped_queries_by_backend[alignment_type]
        skipped_queries = raw_queries.intersection(skipped_query_details)
        observed_raw_all.update(raw_queries)
        observed_mapped_all.update(mapped_queries)
        counts[f"{alignment_type}_aln"] = len(raw_queries)
        counts[f"{alignment_type}_mapped_aln"] = len(mapped_queries)
        counts[f"{alignment_type}_skipped_queries"] = len(skipped_queries)
        counts[f"{alignment_type}_release_shards"] = len(
            release_shards_by_backend[alignment_type]
        )
        absent_raw = sorted(backend_queries[alignment_type].difference(raw_queries))
        if absent_raw:
            missing[f"{alignment_type}_aln"] = absent_raw[:100]
        unexpected_raw = sorted(raw_queries.difference(expected_queries))
        if unexpected_raw:
            missing[f"{alignment_type}_aln_unexpected"] = unexpected_raw[:100]
        raw_without_mapping = sorted(
            raw_queries.difference(mapped_queries).difference(skipped_queries)
        )
        if raw_without_mapping:
            missing[f"{alignment_type}_mapping"] = raw_without_mapping[:100]
    absent_raw_all = sorted(expected_queries.difference(observed_raw_all))
    if absent_raw_all:
        missing["all_backends_aln"] = absent_raw_all[:100]
    absent_mapped_all = sorted(
        expected_queries.difference(skipped_query_details).difference(
            observed_mapped_all
        )
    )
    if absent_mapped_all:
        missing["all_backends_mapped_aln"] = absent_mapped_all[:100]
    if missing:
        raise ValueError(f"protein alignment artifacts are incomplete: {missing}")
    report = {
        **{
            key: plan[key]
            for key in [
                "query_count",
                "protein_chain_count",
                "shard_count",
                "max_seqs",
                "search_databases",
                "alignment_types",
                "target_clustering",
            ]
        },
        "status": "complete",
        "artifact_counts": counts,
        "exact_cluster_manifests": cluster_manifests,
        "skipped_queries": skipped_query_details,
    }
    _atomic_json(report, data_dir / "alignments" / "manifest.json")
    return report


def _ligand_3d_pair_validation_inputs(
    *, data_dir: Path, work_path: Path, pair_paths: list[Path]
) -> dict[str, Any]:
    """Fingerprint inputs covered by the expensive pair-score audit."""
    digest = hashlib.sha256()
    pair_root = data_dir / "scores" / "ligand_3d_pairs"
    for path in pair_paths:
        stat = path.stat()
        relative = path.relative_to(pair_root).as_posix()
        digest.update(f"{relative}\0{stat.st_size}\0{stat.st_mtime_ns}\n".encode())
    return {
        "work": _source_signature(work_path),
        "pair_file_count": len(pair_paths),
        "pair_files_sha256": digest.hexdigest(),
        "required_schema": schemas.LIGAND_3D_SCORE_SCHEMA.names,
    }


def finalize_ligand_3d_artifacts(data_dir: Path) -> dict[str, Any]:
    """Validate canonical pair coverage and merged score completion."""
    started = perf_counter()
    LOG.info("score finalization: validating plan and candidate manifest")
    plan = _load_plan(data_dir)
    if not plan.get("ligand_3d_plan_complete", False):
        raise ValueError("ligand 3D planning did not complete")
    work_path = data_dir / LIGAND_3D_WORK_RELATIVE
    if _source_signature(work_path) != plan.get("ligand_3d_work"):
        raise ValueError("ligand 3D work manifest changed after planning")
    candidate_manifest_path = data_dir / LIGAND_3D_CANDIDATE_MANIFEST_RELATIVE
    if _source_signature(candidate_manifest_path) != plan.get(
        "ligand_3d_candidate_manifest"
    ):
        raise ValueError("ligand 3D candidate manifest changed after planning")
    candidate_manifest = pd.read_parquet(candidate_manifest_path)
    candidate_count = len(candidate_manifest)
    phase_started = perf_counter()
    for index, row in enumerate(candidate_manifest.itertuples(index=False), start=1):
        path = Path(str(row.path))
        if not path.is_file():
            raise FileNotFoundError(path)
        stat = path.stat()
        if stat.st_size != int(row.size) or stat.st_mtime_ns != int(row.mtime_ns):
            raise ValueError(f"ligand 3D candidate changed after planning: {path}")
        if index % 250 == 0 or index == candidate_count:
            LOG.info(
                "score finalization candidate manifests: "
                f"validated={index}/{candidate_count} "
                f"elapsed_seconds={perf_counter() - phase_started:.1f}"
            )

    pair_dir = data_dir / "scores" / "ligand_3d_pairs"
    pair_paths = sorted(pair_dir.glob("*.parquet"))
    expected_batches = set(range(int(plan["ligand_3d_batch_count"])))
    observed_batches = {int(path.stem) for path in pair_paths if path.stem.isdigit()}
    if observed_batches != expected_batches:
        raise ValueError(
            "ligand 3D pair score batches are incomplete: "
            f"missing={sorted(expected_batches - observed_batches)[:10]}, "
            f"extra={sorted(observed_batches - expected_batches)[:10]}"
        )
    expected_pdb_ids = active_scoring_query_ids(data_dir)
    expected_query_shards = {pdb_id[1:3] for pdb_id in expected_pdb_ids}
    final_score_dir = data_dir / "scores" / "search_db=holo"
    final_score_paths = sorted(final_score_dir.glob("*.parquet"))
    observed_final_shards = {path.stem for path in final_score_paths}
    if observed_final_shards != expected_query_shards:
        raise ValueError(
            "packed ligand 3D score shards are incomplete: "
            f"missing={sorted(expected_query_shards - observed_final_shards)[:10]}, "
            f"extra={sorted(observed_final_shards - expected_query_shards)[:10]}"
        )
    LOG.info(
        f"score finalization: validating {len(final_score_paths)} packed score schemas"
    )
    phase_started = perf_counter()
    for index, path in enumerate(final_score_paths, start=1):
        schema = pq.read_schema(path)
        missing = sorted(
            set(schemas.PROTEIN_SIMILARITY_SCHEMA.names).difference(schema.names)
        )
        if missing:
            raise ValueError(f"packed score shard {path} is missing {missing}")
        if index % 100 == 0 or index == len(final_score_paths):
            LOG.info(
                "score finalization packed schemas: "
                f"validated={index}/{len(final_score_paths)} "
                f"elapsed_seconds={perf_counter() - phase_started:.1f}"
            )

    pair_validation_inputs = _ligand_3d_pair_validation_inputs(
        data_dir=data_dir,
        work_path=work_path,
        pair_paths=pair_paths,
    )
    validation_path = data_dir / LIGAND_3D_PAIR_VALIDATION_RELATIVE
    cached_validation: dict[str, Any] = {}
    if validation_path.is_file():
        try:
            cached_validation = json.loads(validation_path.read_text())
        except (OSError, TypeError, ValueError):
            cached_validation = {}
    cached_coverage = cached_validation.get("coverage")
    coverage_columns = [
        "expected_rows",
        "observed_rows",
        "missing_rows",
        "extra_rows",
        "duplicate_rows",
    ]
    if (
        cached_validation.get("status") == "complete"
        and cached_validation.get("inputs") == pair_validation_inputs
        and isinstance(cached_coverage, dict)
        and set(coverage_columns).issubset(cached_coverage)
    ):
        coverage = pd.Series(
            {column: int(cached_coverage[column]) for column in coverage_columns}
        )
        LOG.info(
            "score finalization: reusing signature-matched canonical-pair "
            f"validation counts={coverage.to_dict()}"
        )
    else:
        LOG.info(f"score finalization: validating {len(pair_paths)} pair-score schemas")
        phase_started = perf_counter()
        for index, path in enumerate(pair_paths, start=1):
            schema = pq.read_schema(path)
            missing = sorted(
                set(schemas.LIGAND_3D_SCORE_SCHEMA.names).difference(schema.names)
            )
            if missing:
                raise ValueError(f"ligand 3D score batch {path} is missing {missing}")
            if index % 1_000 == 0 or index == len(pair_paths):
                elapsed = perf_counter() - phase_started
                rate = index / elapsed
                LOG.info(
                    "score finalization pair schemas: "
                    f"validated={index}/{len(pair_paths)} rate={rate:.1f}/s "
                    f"eta_seconds={(len(pair_paths) - index) / rate:.1f}"
                )

        import duckdb

        connection = duckdb.connect()
        pair_keys = ", ".join(
            [
                "query_entry",
                "query_ligand_asym_id",
                "target_entry",
                "target_ligand_asym_id",
            ]
        )
        pair_paths_sql = ", ".join(f"'{path.as_posix()}'" for path in pair_paths)
        LOG.info("score finalization: starting canonical-pair coverage query")
        phase_started = perf_counter()
        coverage = (
            connection.sql(
                dedent(
                    f"""
                WITH expected AS (
                    SELECT {pair_keys}
                    FROM read_parquet('{work_path.as_posix()}')
                ), observed AS (
                    SELECT {pair_keys}
                    FROM read_parquet([{pair_paths_sql}])
                )
                SELECT
                    (SELECT count(*) FROM expected) AS expected_rows,
                    (SELECT count(*) FROM observed) AS observed_rows,
                    (SELECT count(*) FROM (
                        SELECT * FROM expected EXCEPT SELECT * FROM observed
                    )) AS missing_rows,
                    (SELECT count(*) FROM (
                        SELECT * FROM observed EXCEPT SELECT * FROM expected
                    )) AS extra_rows,
                    (SELECT count(*) - count(DISTINCT ({pair_keys})) FROM observed)
                        AS duplicate_rows
                """
                )
            )
            .df()
            .iloc[0]
        )
        connection.close()
        LOG.info(
            "score finalization: canonical-pair coverage query complete "
            f"elapsed_seconds={perf_counter() - phase_started:.1f} "
            f"counts={coverage.to_dict()}"
        )
        if any(
            int(coverage[column])
            for column in ["missing_rows", "extra_rows", "duplicate_rows"]
        ):
            raise ValueError(f"ligand 3D pair coverage mismatch: {coverage.to_dict()}")
        _atomic_json(
            {
                "status": "complete",
                "inputs": pair_validation_inputs,
                "coverage": {
                    column: int(coverage[column]) for column in coverage.index
                },
            },
            validation_path,
        )
    report = {
        "status": "complete",
        "canonical_pair_count": int(coverage["expected_rows"]),
        "pair_batch_count": len(pair_paths),
        "query_shard_count": len(final_score_paths),
        "query_count": len(expected_pdb_ids),
    }
    _atomic_json(report, data_dir / "scores" / "ligand_3d_manifest.json")
    LOG.info(
        "score finalization complete: "
        f"elapsed_seconds={perf_counter() - started:.1f} report={report}"
    )
    return report


def finalize_score_repair_artifacts(
    data_dir: Path,
    *,
    repair_manifest: Path,
) -> dict[str, Any]:
    """Validate ligand-pair caches and score shards rebuilt by a repair run.

    A repair may intentionally refresh candidate shards after the original
    all-vs-all ligand plan was finalized.  The original finalizer must keep
    rejecting that mutation; this repair-specific finalizer instead verifies
    that every repaired shard was rebuilt from the current candidate and pair
    inputs.
    """
    repairs = pd.read_parquet(repair_manifest)
    if "shard" not in repairs.columns:
        raise ValueError(f"repair manifest is missing shard: {repair_manifest}")
    shards = sorted(set(repairs["shard"].dropna().astype(str)))
    if not shards:
        raise ValueError(f"repair manifest contains no shards: {repair_manifest}")

    repair_started_ns = repair_manifest.stat().st_mtime_ns
    candidate_dir = data_dir / "scores" / "ligand_3d_candidate_shards"
    pair_candidate_dir = data_dir / "scores" / "ligand_3d_pair_candidate_shards"
    pair_dir = data_dir / "scores" / "ligand_3d_by_query"
    score_dir = data_dir / "scores" / "search_db=holo"
    candidate_rows = 0
    pair_rows = 0
    score_rows = 0
    started = perf_counter()
    for index, shard in enumerate(shards, start=1):
        candidate = candidate_dir / f"shard={shard}.parquet"
        pair_candidate = pair_candidate_dir / f"shard={shard}.parquet"
        pair = pair_dir / f"{shard}.parquet"
        score = score_dir / f"{shard}.parquet"
        manifest_path = candidate.with_suffix(".json")
        required_paths = [candidate, pair_candidate, pair, score, manifest_path]
        missing_paths = [path for path in required_paths if not path.is_file()]
        if missing_paths:
            raise FileNotFoundError(
                f"missing repaired score artifacts for shard {shard}: {missing_paths}"
            )

        try:
            candidate_manifest = json.loads(manifest_path.read_text())
        except (OSError, TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid ligand 3D candidate manifest: {manifest_path}"
            ) from exc
        candidate_stat = candidate.stat()
        pair_candidate_stat = pair_candidate.stat()
        candidate_row_count = pq.ParquetFile(candidate).metadata.num_rows
        pair_candidate_row_count = pq.ParquetFile(pair_candidate).metadata.num_rows
        expected_candidate = {
            "path": str(candidate.resolve()),
            "size": candidate_stat.st_size,
            "mtime_ns": candidate_stat.st_mtime_ns,
            "rows": candidate_row_count,
        }
        expected_pair_candidate = {
            "path": str(pair_candidate.resolve()),
            "size": pair_candidate_stat.st_size,
            "mtime_ns": pair_candidate_stat.st_mtime_ns,
            "rows": pair_candidate_row_count,
        }
        if (
            candidate_manifest.get("shard") != shard
            or candidate_manifest.get("output") != expected_candidate
            or candidate_manifest.get("pair_output") != expected_pair_candidate
        ):
            raise ValueError(
                f"stale ligand 3D candidate manifest for repaired shard {shard}"
            )

        candidate_missing = sorted(
            set(schemas.LIGAND_3D_CANDIDATE_SCHEMA.names).difference(
                pq.read_schema(candidate).names
            )
        )
        pair_candidate_missing = sorted(
            set(schemas.LIGAND_3D_PAIR_CANDIDATE_SCHEMA.names).difference(
                pq.read_schema(pair_candidate).names
            )
        )
        pair_missing = sorted(
            set(schemas.LIGAND_3D_SCORE_SCHEMA.names).difference(
                pq.read_schema(pair).names
            )
        )
        score_missing = sorted(
            set(schemas.PROTEIN_SIMILARITY_SCHEMA.names).difference(
                pq.read_schema(score).names
            )
        )
        missing_columns = {
            name: columns
            for name, columns in [
                ("candidate", candidate_missing),
                ("pair_candidate", pair_candidate_missing),
                ("pair", pair_missing),
                ("score", score_missing),
            ]
            if columns
        }
        if missing_columns:
            raise ValueError(
                f"repaired shard {shard} has incomplete schemas: {missing_columns}"
            )

        score_stat = score.stat()
        newest_input_ns = max(
            repair_started_ns,
            candidate_stat.st_mtime_ns,
            pair_candidate_stat.st_mtime_ns,
            pair.stat().st_mtime_ns,
        )
        if score_stat.st_mtime_ns < newest_input_ns:
            raise ValueError(
                f"repaired score shard {score} predates its current inputs"
            )
        candidate_rows += candidate_row_count
        pair_rows += pq.ParquetFile(pair).metadata.num_rows
        score_rows += pq.ParquetFile(score).metadata.num_rows
        if index % 100 == 0 or index == len(shards):
            elapsed = perf_counter() - started
            rate = index / elapsed
            LOG.info(
                "score repair finalization: "
                f"validated={index}/{len(shards)} rate={rate:.1f}/s "
                f"eta_seconds={(len(shards) - index) / rate:.1f}"
            )

    report = {
        "status": "complete",
        "repair_manifest": _source_signature(repair_manifest),
        "shard_count": len(shards),
        "candidate_rows": candidate_rows,
        "pair_rows": pair_rows,
        "score_rows": score_rows,
    }
    _atomic_json(report, data_dir / "scores" / "score_repair_manifest.json")
    return report


def _interface_score_alignment_manifest(data_dir: Path) -> tuple[Path, dict[str, Any]]:
    path = data_dir / "alignments" / "manifest.json"
    try:
        payload = cast(dict[str, Any], json.loads(path.read_text()))
    except (OSError, TypeError, ValueError) as exc:
        raise FileNotFoundError(
            f"missing or invalid finalized alignment manifest: {path}"
        ) from exc
    if payload.get("status") != "complete":
        raise ValueError(f"protein alignments are not finalized: {path}")
    if not isinstance(payload.get("skipped_queries"), dict):
        raise ValueError(f"alignment manifest has no skipped-query map: {path}")
    return path, payload


def plan_interface_scoring(
    data_dir: Path,
    *,
    batch_size: int = 1,
) -> dict[str, Any]:
    """Freeze query-sharded all-vs-all interface scoring work.

    Every interface whose entry has a completed mapped protein query is kept
    as a query. Mapping-budget drops remain valid targets but are omitted as
    queries, matching the mapped-alignment release contract.
    """
    if batch_size < 1:
        raise ValueError("interface score batch size must be positive")
    protein_plan = _load_plan(
        data_dir,
        recheck_source=True,
        recheck_score_inputs=False,
    )
    alignment_manifest_path, alignment_manifest = _interface_score_alignment_manifest(
        data_dir
    )
    interface_path = data_dir / "index" / "interface_annotation_table.parquet"
    if not interface_path.is_file():
        raise FileNotFoundError(interface_path)
    interface_schema = pq.read_schema(interface_path)
    required = {"entry_pdb_id", "system_id"}
    missing = required.difference(interface_schema.names)
    if missing:
        raise ValueError(f"interface annotation is missing columns: {sorted(missing)}")

    planned_queries = set(
        pd.read_parquet(data_dir / MANIFEST_RELATIVE, columns=["pdb_id"])[
            "pdb_id"
        ].astype(str)
    )
    skipped_queries = {str(value) for value in alignment_manifest["skipped_queries"]}
    query_entries = planned_queries.difference(skipped_queries)
    interfaces = pd.read_parquet(
        interface_path,
        columns=["entry_pdb_id", "system_id"],
    )
    interfaces["entry_pdb_id"] = interfaces["entry_pdb_id"].astype(str)
    interfaces = interfaces[interfaces["entry_pdb_id"].isin(planned_queries)].copy()
    if interfaces["system_id"].duplicated().any():
        duplicate = str(
            interfaces.loc[interfaces["system_id"].duplicated(), "system_id"].iloc[0]
        )
        raise ValueError(f"interface system IDs are not unique: {duplicate}")
    query_interfaces = interfaces[interfaces["entry_pdb_id"].isin(query_entries)].copy()
    query_interfaces["shard"] = query_interfaces["entry_pdb_id"].str[1:3]

    records: list[dict[str, Any]] = []
    for shard, shard_interfaces in query_interfaces.groupby("shard", sort=True):
        alignment_rows = 0
        backend_count = 0
        for alignment_type in ["foldseek", "mmseqs"]:
            path = tasks._alignment_release_path(
                data_dir=data_dir,
                search_db="holo",
                alignment_type=alignment_type,
                shard=str(shard),
            )
            if path.is_file():
                backend_count += 1
                alignment_rows += int(pq.ParquetFile(path).metadata.num_rows)
        if backend_count == 0:
            raise FileNotFoundError(
                f"no mapped alignment backend is available for interface shard {shard}"
            )
        records.append(
            {
                "shard": str(shard),
                "query_entry_count": int(shard_interfaces["entry_pdb_id"].nunique()),
                "query_interface_count": len(shard_interfaces),
                "alignment_rows": alignment_rows,
            }
        )
    work = pd.DataFrame.from_records(
        records,
        columns=[
            "shard",
            "query_entry_count",
            "query_interface_count",
            "alignment_rows",
        ],
    )
    work_path = data_dir / INTERFACE_SCORE_WORK_RELATIVE
    _atomic_parquet(work, work_path)
    payload = {
        "status": "planned",
        "batch_size": batch_size,
        "batch_count": math.ceil(len(work) / batch_size),
        "shard_count": len(work),
        "query_entry_count": int(query_interfaces["entry_pdb_id"].nunique()),
        "query_interface_count": len(query_interfaces),
        "target_interface_count": len(interfaces),
        "skipped_query_count": len(skipped_queries.intersection(planned_queries)),
        "protein_manifest": protein_plan["manifest"],
        "alignment_manifest": _source_signature(alignment_manifest_path),
        "interface_annotation": _source_signature(interface_path),
        "work": _source_signature(work_path),
    }
    _atomic_json(payload, data_dir / INTERFACE_SCORE_PLAN_RELATIVE)

    # A new frozen plan replaces the previous query universe. Remove only
    # query-sharded interface outputs that the new plan no longer names so a
    # resumed finalizer cannot mix generations or fail on obsolete shards.
    expected_shards = set(work["shard"].astype(str))
    score_root = data_dir / INTERFACE_SCORE_ROOT_RELATIVE
    for stale in score_root.glob("shard=*.parquet"):
        if stale.stem.removeprefix("shard=") not in expected_shards:
            stale.unlink()
            stale.with_suffix(".json").unlink(missing_ok=True)
    return payload


def _load_interface_score_plan(data_dir: Path) -> dict[str, Any]:
    path = data_dir / INTERFACE_SCORE_PLAN_RELATIVE
    try:
        payload = cast(dict[str, Any], json.loads(path.read_text()))
    except (OSError, TypeError, ValueError) as exc:
        raise FileNotFoundError(
            f"missing or invalid interface score plan: {path}"
        ) from exc
    checks = {
        "protein_manifest": data_dir / MANIFEST_RELATIVE,
        "alignment_manifest": data_dir / "alignments" / "manifest.json",
        "interface_annotation": (
            data_dir / "index" / "interface_annotation_table.parquet"
        ),
        "work": data_dir / INTERFACE_SCORE_WORK_RELATIVE,
    }
    for key, source in checks.items():
        if not source.is_file() or payload.get(key) != _source_signature(source):
            raise ValueError(f"interface scoring input changed after planning: {key}")
    if payload.get("status") != "planned":
        raise ValueError(
            f"invalid interface score plan status: {payload.get('status')}"
        )
    return payload


def _interface_score_shard_batch(
    data_dir: Path,
    *,
    batch_index: int,
    batch_size: int,
) -> list[str]:
    if batch_index < 0 or batch_size < 1:
        raise ValueError("batch_index must be non-negative and batch_size positive")
    plan = _load_interface_score_plan(data_dir)
    if batch_size != int(plan["batch_size"]):
        raise ValueError(
            f"interface score batch size is {plan['batch_size']}, got {batch_size}"
        )
    work = pd.read_parquet(data_dir / INTERFACE_SCORE_WORK_RELATIVE, columns=["shard"])
    start = batch_index * batch_size
    return cast(
        list[str],
        work.iloc[start : start + batch_size]["shard"].astype(str).tolist(),
    )


def _interface_score_inputs(data_dir: Path, shard: str) -> dict[str, Any]:
    interface_path = data_dir / "index" / "interface_annotation_table.parquet"
    mapping_manifest = tasks._alignment_mapping_manifest_path(
        data_dir=data_dir,
        shard=shard,
    )
    if not mapping_manifest.is_file():
        raise FileNotFoundError(mapping_manifest)
    alignments: dict[str, dict[str, int | str]] = {}
    for alignment_type in ["foldseek", "mmseqs"]:
        path = tasks._alignment_release_path(
            data_dir=data_dir,
            search_db="holo",
            alignment_type=alignment_type,
            shard=shard,
        )
        if path.is_file():
            alignments[alignment_type] = _source_signature(path)
    if not alignments:
        raise FileNotFoundError(
            f"no mapped alignment backend is available for interface shard {shard}"
        )
    return {
        "interface_annotation": _source_signature(interface_path),
        "alignment_manifest": _source_signature(
            data_dir / "alignments" / "manifest.json"
        ),
        "mapping_manifest": _source_signature(mapping_manifest),
        "alignments": alignments,
    }


def _interface_score_output_is_current(
    *,
    output: Path,
    manifest: Path,
    inputs: dict[str, Any],
) -> dict[str, Any] | None:
    if not output.is_file() or not manifest.is_file():
        return None
    try:
        payload = cast(dict[str, Any], json.loads(manifest.read_text()))
        if (
            payload.get("inputs") == inputs
            and payload.get("output") == _source_signature(output)
            and pq.read_schema(output).equals(schemas.INTERFACE_SCORE_SHARD_SCHEMA)
            and int(payload["rows"]) == pq.ParquetFile(output).metadata.num_rows
        ):
            return payload
    except (OSError, TypeError, ValueError):
        pass
    return None


def score_interface_qcov_shards(
    data_dir: Path,
    *,
    shards: list[str],
    scratch_dir: Path,
    threads: int = 4,
    memory_limit: str = "32GB",
    force_update: bool = False,
) -> dict[str, Any]:
    """Score all directed interface pairs for complete query shards."""
    if threads < 1:
        raise ValueError("interface scoring threads must be positive")
    _load_interface_score_plan(data_dir)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    output_root = data_dir / INTERFACE_SCORE_ROOT_RELATIVE
    output_root.mkdir(exist_ok=True, parents=True)
    interface_path = data_dir / "index" / "interface_annotation_table.parquet"
    reports: list[dict[str, Any]] = []
    started = perf_counter()
    for index, shard in enumerate(shards, start=1):
        inputs = _interface_score_inputs(data_dir, shard)
        output = output_root / f"shard={shard}.parquet"
        manifest = output.with_suffix(".json")
        current = _interface_score_output_is_current(
            output=output,
            manifest=manifest,
            inputs=inputs,
        )
        if not force_update and current is not None:
            reports.append(current)
            LOG.info(
                "interface score progress: shards=%d/%d shard=%s rows=%s "
                "elapsed_seconds=%.1f cached=true",
                index,
                len(shards),
                shard,
                current["rows"],
                perf_counter() - started,
            )
            continue

        backend_selects: list[str] = []
        alignments = cast(dict[str, dict[str, Any]], inputs["alignments"])
        for alignment_type, signature in alignments.items():
            path = Path(str(signature["path"]))
            backend_selects.append(
                dedent(
                    f"""
                    SELECT
                        '{alignment_type}' AS backend,
                        row_number() OVER ()::BIGINT AS alignment_id,
                        query_entry,
                        target_entry,
                        query_chain_mapped,
                        target_chain_mapped,
                        query_selected_residue_numbers,
                        target_selected_residue_numbers
                    FROM read_parquet('{path.as_posix()}')
                    """
                ).strip()
            )
        mapped_sql = "\nUNION ALL\n".join(backend_selects)
        local_root = scratch_dir / shard
        local_root.mkdir(exist_ok=True, parents=True)
        local_output = local_root / output.name
        local_output.unlink(missing_ok=True)
        import duckdb

        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{local_root.as_posix()}'")
        connection.sql(f"SET memory_limit='{memory_limit}'")
        connection.sql("SET preserve_insertion_order=false")
        connection.sql(
            dedent(
                f"""
                COPY (
                    WITH interface_sides AS (
                        SELECT
                            entry_pdb_id,
                            system_id,
                            1::TINYINT AS side,
                            interface_chain_1 AS instance_chain,
                            split_part(interface_chain_1, '.', 2) AS asym_id,
                            interface_chain_1_residue_numbers AS residue_numbers,
                            len(interface_chain_1_residue_numbers)::DOUBLE
                                AS residue_count
                        FROM read_parquet('{interface_path.as_posix()}')
                        UNION ALL
                        SELECT
                            entry_pdb_id,
                            system_id,
                            2::TINYINT AS side,
                            interface_chain_2 AS instance_chain,
                            split_part(interface_chain_2, '.', 2) AS asym_id,
                            interface_chain_2_residue_numbers AS residue_numbers,
                            len(interface_chain_2_residue_numbers)::DOUBLE
                                AS residue_count
                        FROM read_parquet('{interface_path.as_posix()}')
                    ), mapped AS (
                        {mapped_sql}
                    ), mapped_interface_residues AS (
                        SELECT
                            mapped.backend,
                            mapped.alignment_id,
                            query_side.system_id AS iface1,
                            target_side.system_id AS iface2,
                            query_side.side AS query_side,
                            target_side.side AS target_side,
                            query_side.instance_chain AS query_instance_chain,
                            target_side.instance_chain AS target_instance_chain,
                            query_side.residue_numbers AS query_interface_residues,
                            target_side.residue_numbers AS target_interface_residues,
                            query_side.residue_count AS query_residue_count,
                            unnest(mapped.query_selected_residue_numbers)
                                AS query_residue,
                            unnest(mapped.target_selected_residue_numbers)
                                AS target_residue
                        FROM mapped
                        INNER JOIN interface_sides AS query_side
                          ON mapped.query_entry = query_side.entry_pdb_id
                         AND mapped.query_chain_mapped = query_side.asym_id
                        INNER JOIN interface_sides AS target_side
                          ON mapped.target_entry = target_side.entry_pdb_id
                         AND mapped.target_chain_mapped = target_side.asym_id
                        WHERE query_side.system_id != target_side.system_id
                    ), coverage_per_alignment AS (
                        SELECT
                            backend,
                            alignment_id,
                            iface1,
                            iface2,
                            query_side,
                            target_side,
                            query_instance_chain,
                            target_instance_chain,
                            count(DISTINCT query_residue)::DOUBLE
                                / query_residue_count AS qcov
                        FROM mapped_interface_residues
                        WHERE list_contains(
                                  query_interface_residues, query_residue
                              )
                          AND list_contains(
                                  target_interface_residues, target_residue
                              )
                        GROUP BY
                            backend,
                            alignment_id,
                            iface1,
                            iface2,
                            query_side,
                            target_side,
                            query_instance_chain,
                            target_instance_chain,
                            query_residue_count
                    ), side_coverage AS (
                        SELECT
                            backend,
                            iface1,
                            iface2,
                            query_side,
                            target_side,
                            query_instance_chain,
                            target_instance_chain,
                            max(qcov) AS qcov
                        FROM coverage_per_alignment
                        GROUP BY
                            backend,
                            iface1,
                            iface2,
                            query_side,
                            target_side,
                            query_instance_chain,
                            target_instance_chain
                    ), assignments AS (
                        SELECT
                            first.backend,
                            first.iface1,
                            first.iface2,
                            first.qcov AS iface1_qcov,
                            second.qcov AS iface2_qcov,
                            first.qcov * second.qcov AS final_score,
                            first.query_instance_chain || ':'
                                || first.target_instance_chain || ';'
                                || second.query_instance_chain || ':'
                                || second.target_instance_chain AS mapping,
                            0::TINYINT AS assignment_order
                        FROM side_coverage AS first
                        INNER JOIN side_coverage AS second
                          ON first.backend = second.backend
                         AND first.iface1 = second.iface1
                         AND first.iface2 = second.iface2
                        WHERE first.query_side = 1
                          AND first.target_side = 1
                          AND second.query_side = 2
                          AND second.target_side = 2
                        UNION ALL
                        SELECT
                            first.backend,
                            first.iface1,
                            first.iface2,
                            first.qcov AS iface1_qcov,
                            second.qcov AS iface2_qcov,
                            first.qcov * second.qcov AS final_score,
                            first.query_instance_chain || ':'
                                || first.target_instance_chain || ';'
                                || second.query_instance_chain || ':'
                                || second.target_instance_chain AS mapping,
                            1::TINYINT AS assignment_order
                        FROM side_coverage AS first
                        INNER JOIN side_coverage AS second
                          ON first.backend = second.backend
                         AND first.iface1 = second.iface1
                         AND first.iface2 = second.iface2
                        WHERE first.query_side = 1
                          AND first.target_side = 2
                          AND second.query_side = 2
                          AND second.target_side = 1
                    ), backend_best AS (
                        SELECT * EXCLUDE (assignment_rank, assignment_order)
                        FROM (
                            SELECT
                                *,
                                row_number() OVER (
                                    PARTITION BY backend, iface1, iface2
                                    ORDER BY final_score DESC,
                                             assignment_order,
                                             mapping
                                ) AS assignment_rank
                            FROM assignments
                        )
                        WHERE assignment_rank = 1
                    ), pair_best AS (
                        SELECT iface1, iface2, max(final_score) AS final_score
                        FROM backend_best
                        GROUP BY iface1, iface2
                    ), winning_backends AS (
                        SELECT backend_best.*
                        FROM backend_best
                        INNER JOIN pair_best USING (iface1, iface2)
                        WHERE abs(
                            backend_best.final_score - pair_best.final_score
                        ) < 1e-12
                    )
                    SELECT
                        iface1::VARCHAR AS query_system,
                        iface2::VARCHAR AS target_system,
                        first(
                            mapping ORDER BY
                            CASE WHEN backend = 'foldseek' THEN 0 ELSE 1 END,
                            mapping
                        )::VARCHAR AS mapping,
                        CASE
                            WHEN count(DISTINCT backend) = 2 THEN 'both'
                            ELSE min(backend)
                        END::VARCHAR AS source,
                        first(
                            iface1_qcov ORDER BY
                            CASE WHEN backend = 'foldseek' THEN 0 ELSE 1 END,
                            mapping
                        )::FLOAT AS iface1_qcov,
                        first(
                            iface2_qcov ORDER BY
                            CASE WHEN backend = 'foldseek' THEN 0 ELSE 1 END,
                            mapping
                        )::FLOAT AS iface2_qcov,
                        CAST(
                            least(100, greatest(0, round(max(final_score) * 100)))
                            AS TINYINT
                        ) AS similarity
                    FROM winning_backends
                    GROUP BY iface1, iface2
                    HAVING max(final_score) > 0
                    ORDER BY query_system, target_system
                ) TO '{local_output.as_posix()}' (
                    FORMAT PARQUET,
                    COMPRESSION ZSTD,
                    ROW_GROUP_SIZE 500000
                )
                """
            )
        )
        connection.close()
        observed_schema = pq.read_schema(local_output)
        if not observed_schema.equals(schemas.INTERFACE_SCORE_SHARD_SCHEMA):
            local_output.unlink(missing_ok=True)
            raise ValueError(
                f"interface score shard has unexpected schema: {observed_schema}"
            )
        rows = int(pq.ParquetFile(local_output).metadata.num_rows)
        install = output.with_suffix(output.suffix + ".tmp")
        copyfile(local_output, install)
        install.replace(output)
        local_output.unlink(missing_ok=True)
        payload = {
            "status": "complete",
            "shard": shard,
            "inputs": inputs,
            "output": _source_signature(output),
            "rows": rows,
        }
        _atomic_json(payload, manifest)
        reports.append(payload)
        LOG.info(
            "interface score progress: shards=%d/%d shard=%s rows=%d "
            "elapsed_seconds=%.1f",
            index,
            len(shards),
            shard,
            rows,
            perf_counter() - started,
        )
    return {
        "status": "complete",
        "shard_count": len(reports),
        "row_count": sum(int(report["rows"]) for report in reports),
        "shards": [str(report["shard"]) for report in reports],
    }


def finalize_interface_qcov_scores(
    data_dir: Path,
    *,
    output: Path | None = None,
    scratch_dir: Path,
    threads: int = 8,
    memory_limit: str = "32GB",
) -> dict[str, Any]:
    """Validate interface score shards and publish the compact release table."""
    if threads < 1:
        raise ValueError("interface score finalization threads must be positive")
    plan = _load_interface_score_plan(data_dir)
    output = (output or data_dir / INTERFACE_QCOV_EXPORT_RELATIVE).resolve()
    score_root = data_dir / INTERFACE_SCORE_ROOT_RELATIVE
    expected_shards = (
        pd.read_parquet(
            data_dir / INTERFACE_SCORE_WORK_RELATIVE,
            columns=["shard"],
        )["shard"]
        .astype(str)
        .tolist()
    )
    observed_paths = sorted(score_root.glob("shard=*.parquet"))
    observed_shards = [path.stem.removeprefix("shard=") for path in observed_paths]
    if observed_shards != expected_shards:
        raise ValueError(
            "interface score shards are incomplete: "
            f"missing={sorted(set(expected_shards) - set(observed_shards))[:10]}, "
            f"extra={sorted(set(observed_shards) - set(expected_shards))[:10]}"
        )

    sources: list[dict[str, Any]] = []
    for index, (shard, path) in enumerate(
        zip(expected_shards, observed_paths), start=1
    ):
        manifest = path.with_suffix(".json")
        inputs = _interface_score_inputs(data_dir, shard)
        payload = _interface_score_output_is_current(
            output=path,
            manifest=manifest,
            inputs=inputs,
        )
        if payload is None or payload.get("shard") != shard:
            raise ValueError(f"stale interface score shard: {path}")
        sources.append(
            {
                "shard": shard,
                **_source_signature(path),
                "rows": int(payload["rows"]),
            }
        )
        if index % 100 == 0 or index == len(observed_paths):
            LOG.info(
                "interface score finalization validation: shards=%d/%d",
                index,
                len(observed_paths),
            )

    release_manifest = output.with_suffix(output.suffix + ".json")
    if output.is_file() and release_manifest.is_file():
        try:
            current = cast(dict[str, Any], json.loads(release_manifest.read_text()))
            if (
                current.get("interface_plan")
                == _source_signature(data_dir / INTERFACE_SCORE_PLAN_RELATIVE)
                and current.get("sources") == sources
                and current.get("output") == _source_signature(output)
                and pq.read_schema(output).equals(schemas.INTERFACE_QCOV_EXPORT_SCHEMA)
            ):
                return current
        except (OSError, TypeError, ValueError):
            pass

    scratch_dir.mkdir(exist_ok=True, parents=True)
    local_output = scratch_dir / output.name
    local_output.unlink(missing_ok=True)
    if observed_paths:
        import duckdb

        paths_sql = ", ".join(f"'{path.as_posix()}'" for path in observed_paths)
        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
        connection.sql(f"SET memory_limit='{memory_limit}'")
        connection.sql("SET preserve_insertion_order=false")
        result = connection.sql(
            dedent(
                f"""
                SELECT
                    count(*)::BIGINT AS rows,
                    count(*) FILTER (
                        WHERE query_system = target_system
                    )::BIGINT AS self_rows
                FROM read_parquet([{paths_sql}])
                """
            )
        ).fetchone()
        if result is None:
            connection.close()
            raise RuntimeError("failed to validate interface score shards")
        rows, self_rows = (int(value) for value in result)
        expected_rows = sum(int(source["rows"]) for source in sources)
        if rows != expected_rows or self_rows:
            connection.close()
            raise ValueError(
                "invalid interface score rows: "
                f"rows={rows}/{expected_rows}, self_rows={self_rows}"
            )
        connection.sql(
            dedent(
                f"""
                COPY (
                    SELECT
                        query_system,
                        target_system,
                        iface1_qcov,
                        iface2_qcov,
                        similarity
                    FROM read_parquet([{paths_sql}])
                ) TO '{local_output.as_posix()}' (
                    FORMAT PARQUET,
                    COMPRESSION ZSTD,
                    ROW_GROUP_SIZE 500000
                )
                """
            )
        )
        connection.close()
    else:
        rows = 0
        pq.write_table(
            pa.Table.from_pylist([], schema=schemas.INTERFACE_QCOV_EXPORT_SCHEMA),
            local_output,
            compression="zstd",
        )
    if not pq.read_schema(local_output).equals(schemas.INTERFACE_QCOV_EXPORT_SCHEMA):
        local_output.unlink(missing_ok=True)
        raise ValueError("compact interface score export has an unexpected schema")
    if pq.ParquetFile(local_output).metadata.num_rows != rows:
        local_output.unlink(missing_ok=True)
        raise ValueError("compact interface score export row count changed")
    output.parent.mkdir(exist_ok=True, parents=True)
    install = output.with_suffix(output.suffix + ".tmp")
    copyfile(local_output, install)
    install.replace(output)
    local_output.unlink(missing_ok=True)
    report = {
        "status": "complete",
        "metric": "interface_qcov",
        "direction": "query_system_to_target_system",
        "interface_plan": _source_signature(data_dir / INTERFACE_SCORE_PLAN_RELATIVE),
        "sources": sources,
        "shard_count": len(sources),
        "row_count": rows,
        "query_interface_count": int(plan["query_interface_count"]),
        "target_interface_count": int(plan["target_interface_count"]),
        "output": _source_signature(output),
    }
    _atomic_json(report, release_manifest)
    return report


def export_sucos_shape_pocket_qcov_batch(
    data_dir: Path,
    *,
    output_dir: Path,
    batch_index: int | None = None,
    batch_size: int | None = None,
    shards: list[str] | None = None,
    scratch_dir: Path,
    threads: int = 4,
    memory_limit: str = "16GB",
) -> dict[str, Any]:
    """Export complete SuCOS-pocket scores as resumable query shards."""
    if threads < 1:
        raise ValueError("export threads must be positive")
    if shards is None:
        if batch_index is None or batch_size is None:
            raise ValueError(
                "batch index and size are required when shards are omitted"
            )
        shards = _ligand_3d_shard_batch(data_dir, batch_index, batch_size)
    elif batch_index is not None or batch_size is not None:
        raise ValueError("pass explicit shards or a batch index and size, not both")
    output_dir = output_dir.resolve()
    output_dir.mkdir(exist_ok=True, parents=True)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    import duckdb

    pair_keys = ", ".join(
        [
            "query_entry",
            "query_ligand_asym_id",
            "target_entry",
            "target_ligand_asym_id",
        ]
    )
    expected_columns = [
        "query_system",
        "query_ligand_id",
        "target_system",
        "target_ligand_id",
        "similarity",
    ]
    reports: list[dict[str, Any]] = []
    exported_rows = 0
    started = perf_counter()
    for index, shard in enumerate(shards, start=1):
        candidate_path = (
            data_dir
            / "scores"
            / "ligand_3d_candidate_shards"
            / f"shard={shard}.parquet"
        )
        pair_path = data_dir / "scores" / "ligand_3d_by_query" / f"{shard}.parquet"
        if not candidate_path.is_file() or not pair_path.is_file():
            raise FileNotFoundError(
                f"missing SuCOS export input for shard {shard}: "
                f"candidate={candidate_path.is_file()} pairs={pair_path.is_file()}"
            )
        inputs = {
            "candidate": _source_signature(candidate_path),
            "pairs": _source_signature(pair_path),
        }
        output = output_dir / f"{shard}.parquet"
        manifest = output.with_suffix(".json")
        if output.is_file() and manifest.is_file():
            try:
                payload = json.loads(manifest.read_text())
                if (
                    payload.get("inputs") == inputs
                    and payload.get("output") == _source_signature(output)
                    and pq.read_schema(output).names == expected_columns
                ):
                    reports.append(payload)
                    exported_rows += int(payload["rows"])
                    LOG.info(
                        "SuCOS export progress: "
                        f"shards={index}/{len(shards)} rows={exported_rows} "
                        f"elapsed_seconds={perf_counter() - started:.1f} cached=true"
                    )
                    continue
            except (OSError, TypeError, ValueError):
                pass

        local_output = scratch_dir / output.name
        local_output.unlink(missing_ok=True)
        connection = duckdb.connect()
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
        connection.sql(f"SET memory_limit='{memory_limit}'")
        connection.sql("SET preserve_insertion_order=false")
        connection.sql(
            f"""
            COPY (
                SELECT
                    candidates.query_system,
                    candidates.query_ligand_id,
                    candidates.target_system,
                    candidates.target_ligand_id,
                    round(pairs.sucos_shape * candidates.pocket_qcov * 100)
                        ::TINYINT AS similarity
                FROM read_parquet('{candidate_path.as_posix()}') AS candidates
                INNER JOIN read_parquet('{pair_path.as_posix()}') AS pairs
                USING ({pair_keys})
                WHERE pairs.sucos_shape IS NOT NULL
            ) TO '{local_output.as_posix()}' (
                FORMAT PARQUET,
                COMPRESSION ZSTD,
                ROW_GROUP_SIZE 500000
            )
            """
        )
        connection.close()
        install = output.with_suffix(output.suffix + ".tmp")
        copyfile(local_output, install)
        install.replace(output)
        local_output.unlink(missing_ok=True)
        metadata = pq.ParquetFile(output).metadata
        payload = {
            "status": "complete",
            "shard": shard,
            "inputs": inputs,
            "output": _source_signature(output),
            "rows": int(metadata.num_rows),
        }
        _atomic_json(payload, manifest)
        reports.append(payload)
        exported_rows += int(payload["rows"])
        elapsed = perf_counter() - started
        LOG.info(
            "SuCOS export progress: "
            f"shards={index}/{len(shards)} rows={exported_rows} "
            f"elapsed_seconds={elapsed:.1f}"
        )
    return {
        "status": "complete",
        "selected_count": len(shards),
        "row_count": exported_rows,
        "outputs": [str(output_dir / f"{shard}.parquet") for shard in shards],
    }


def finalize_sucos_shape_pocket_qcov_export(
    data_dir: Path,
    *,
    source_dir: Path,
    output: Path,
    scratch_dir: Path,
    threads: int = 8,
    memory_limit: str = "32GB",
) -> dict[str, Any]:
    """Validate SuCOS query shards and publish one release Parquet."""
    if threads < 1:
        raise ValueError("export threads must be positive")
    source_dir = source_dir.resolve()
    output = output.resolve()
    if output.parent == source_dir or output.is_relative_to(source_dir):
        raise ValueError("final SuCOS export must be outside its shard directory")
    expected_shards = sorted(
        {pdb_id[1:3] for pdb_id in active_scoring_query_ids(data_dir)}
    )
    paths = sorted(source_dir.glob("*.parquet"))
    observed_shards = [path.stem for path in paths]
    if observed_shards != expected_shards:
        raise ValueError(
            "SuCOS export shards are incomplete: "
            f"missing={sorted(set(expected_shards) - set(observed_shards))[:10]}, "
            f"extra={sorted(set(observed_shards) - set(expected_shards))[:10]}"
        )

    sources: list[dict[str, Any]] = []
    for index, path in enumerate(paths, start=1):
        manifest = path.with_suffix(".json")
        try:
            payload = json.loads(manifest.read_text())
        except (OSError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid SuCOS export manifest: {manifest}") from exc
        shard = path.stem
        candidate_path = (
            data_dir
            / "scores"
            / "ligand_3d_candidate_shards"
            / f"shard={shard}.parquet"
        )
        pair_path = data_dir / "scores" / "ligand_3d_by_query" / f"{shard}.parquet"
        expected_inputs = {
            "candidate": _source_signature(candidate_path),
            "pairs": _source_signature(pair_path),
        }
        if (
            payload.get("shard") != shard
            or payload.get("inputs") != expected_inputs
            or payload.get("output") != _source_signature(path)
        ):
            raise ValueError(f"stale SuCOS export shard: {path}")
        sources.append(
            {
                "shard": shard,
                **_source_signature(path),
                "rows": int(payload["rows"]),
            }
        )
        if index % 100 == 0 or index == len(paths):
            LOG.info(
                "SuCOS finalization validation progress: "
                f"shards={index}/{len(paths)}"
            )

    release_manifest = output.with_suffix(output.suffix + ".json")
    if output.is_file() and release_manifest.is_file():
        try:
            existing = cast(dict[str, Any], json.loads(release_manifest.read_text()))
            if existing.get("sources") == sources and existing.get(
                "output"
            ) == _source_signature(output):
                return existing
        except (OSError, TypeError, ValueError):
            pass

    scratch_dir.mkdir(exist_ok=True, parents=True)
    local_output = scratch_dir / output.name
    local_output.unlink(missing_ok=True)
    paths_sql = ", ".join(f"'{path.as_posix()}'" for path in paths)
    import duckdb

    connection = duckdb.connect()
    connection.sql(f"SET threads={threads}")
    connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
    connection.sql(f"SET memory_limit='{memory_limit}'")
    connection.sql("SET preserve_insertion_order=false")
    LOG.info(
        "SuCOS finalization: concatenating "
        f"shards={len(paths)} rows={sum(int(source['rows']) for source in sources)}"
    )
    connection.sql(
        f"""
        COPY (
            SELECT *
            FROM read_parquet([{paths_sql}])
        ) TO '{local_output.as_posix()}' (
            FORMAT PARQUET,
            COMPRESSION ZSTD,
            ROW_GROUP_SIZE 500000
        )
        """
    )
    connection.close()
    metadata = pq.ParquetFile(local_output).metadata
    expected_rows = sum(int(source["rows"]) for source in sources)
    if metadata.num_rows != expected_rows:
        local_output.unlink(missing_ok=True)
        raise ValueError(
            f"final SuCOS export has {metadata.num_rows} rows; expected {expected_rows}"
        )
    output.parent.mkdir(exist_ok=True, parents=True)
    install = output.with_suffix(output.suffix + ".tmp")
    copyfile(local_output, install)
    install.replace(output)
    local_output.unlink(missing_ok=True)
    report = {
        "status": "complete",
        "metric": "sucos_shape_pocket_qcov",
        "sources": sources,
        "row_count": expected_rows,
        "output": _source_signature(output),
    }
    _atomic_json(report, release_manifest)
    return report


def _add_cluster_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument(
        "--metric",
        action="append",
        dest="metrics",
        help="cluster only this metric; repeat to select multiple metrics",
    )
    command.add_argument(
        "--threshold",
        action="append",
        type=int,
        dest="thresholds",
        help="cluster at this 0-100 threshold; repeat to select multiple values",
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan")
    plan.add_argument("data_dir", type=Path)
    plan.add_argument("--max-seqs", type=int, default=10_000)
    plan.add_argument("--pdb-id", action="append", default=[])
    plan.add_argument("--two-char-code", action="append", default=[])

    databases = subparsers.add_parser("create-dbs")
    databases.add_argument("data_dir", type=Path)
    databases.add_argument("--cif-root", type=Path, required=True)
    databases.add_argument("--seqres-path", type=Path, required=True)
    databases.add_argument("--threads", type=int, default=1)
    databases.add_argument("--scratch-dir", type=Path)
    databases.add_argument("--force", action="store_true")

    sub_dbs = subparsers.add_parser("make-sub-dbs")
    sub_dbs.add_argument("data_dir", type=Path)
    sub_dbs.add_argument("--threads", type=int, default=1)
    sub_dbs.add_argument("--scratch-dir", type=Path, required=True)
    lookup_refresh = subparsers.add_parser("refresh-alignment-lookup")
    lookup_refresh.add_argument("data_dir", type=Path)
    lookup_refresh.add_argument("--threads", type=int, default=1)
    lookup_refresh.add_argument("--scratch-dir", type=Path, required=True)

    score_plan = subparsers.add_parser("plan-score-batches")
    score_plan.add_argument("data_dir", type=Path)
    score_plan.add_argument("--batch-size", type=int, default=50)
    score_plan.add_argument("--threads", type=int, default=4)
    score_plan.add_argument("--scratch-dir", type=Path, required=True)
    score_plan.add_argument("--max-query-protein-chains", type=int, default=30)
    score_plan.add_argument("--max-query-proper-ligand-chains", type=int, default=30)

    ligand_3d_plan = subparsers.add_parser("plan-ligand-3d")
    ligand_3d_plan.add_argument("data_dir", type=Path)
    ligand_3d_plan.add_argument("--batch-size", type=int, default=30_000)
    ligand_3d_plan.add_argument("--threads", type=int, default=4)
    ligand_3d_plan.add_argument("--scratch-dir", type=Path, required=True)
    ligand_3d_plan.add_argument("--memory-limit", default="8GB")

    ligand_3d_retry_plan = subparsers.add_parser("plan-ligand-3d-retries")
    ligand_3d_retry_plan.add_argument("data_dir", type=Path)
    ligand_3d_retry_plan.add_argument("--batch-size", type=int, default=500)
    interface_plan = subparsers.add_parser("plan-interface-scores")
    interface_plan.add_argument("data_dir", type=Path)
    interface_plan.add_argument("--batch-size", type=int, default=1)
    repair_plan = subparsers.add_parser("plan-score-repair")
    repair_plan.add_argument("data_dir", type=Path)
    repair_plan.add_argument("--affected-manifest", type=Path, required=True)
    repair_plan.add_argument("--additional-full-query-manifest", type=Path)
    repair_plan.add_argument("--output", type=Path)
    repair_plan.add_argument("--batch-size", type=int, default=10)
    repair_plan.add_argument("--target-batch-size", type=int, default=100)
    repair_plan.add_argument("--threads", type=int, default=4)
    repair_plan.add_argument("--scratch-dir", type=Path)
    repair_plan.add_argument("--memory-limit", default="32GB")
    repair_ligand_plan = subparsers.add_parser("plan-score-repair-ligand-3d")
    repair_ligand_plan.add_argument("data_dir", type=Path)
    repair_ligand_plan.add_argument("--repair-manifest", type=Path, required=True)
    repair_ligand_plan.add_argument("--batch-size", type=int, default=30_000)
    repair_ligand_plan.add_argument("--threads", type=int, default=4)
    repair_ligand_plan.add_argument("--scratch-dir", type=Path, required=True)
    repair_ligand_plan.add_argument("--memory-limit", default="8GB")
    repair_finalizer = subparsers.add_parser("finalize-score-repair-queries")
    repair_finalizer.add_argument("data_dir", type=Path)
    repair_finalizer.add_argument("--repair-manifest", type=Path, required=True)
    repair_finalizer.add_argument("--max-new-drops", type=int, required=True)
    cluster_plan = subparsers.add_parser("plan-clusters")
    cluster_plan.add_argument("data_dir", type=Path)
    cluster_plan.add_argument("--source-batch-size", type=int, default=20)
    cluster_plan.add_argument("--community-batch-size", type=int, default=1)
    cluster_plan.add_argument(
        "--symmetric-bucket-count",
        type=int,
        default=clusters.SYMMETRIC_EDGE_BUCKET_COUNT,
    )
    _add_cluster_arguments(cluster_plan)

    cluster_stats = subparsers.add_parser("cluster-stats")
    cluster_stats.add_argument("data_dir", type=Path)
    _add_cluster_arguments(cluster_stats)
    index_finalizer = subparsers.add_parser("finalize-index")
    index_finalizer.add_argument("data_dir", type=Path)
    dropped_queries = subparsers.add_parser("drop-score-queries")
    dropped_queries.add_argument("data_dir", type=Path)
    dropped_queries.add_argument("--pdb-manifest", type=Path, required=True)

    for name in [
        "search",
        "map",
        "pack-ligands",
        "score",
        "score-pdbs",
        "repair-score-pdbs",
        "repair-score-target-remainders",
        "repair-pack-ligands",
        "repair-candidates",
        "repair-ligand-3d",
        "repair-merge-ligand-3d",
        "repair-score-shards",
        "repair-sucos-shards",
        "collate-ligand-3d-candidates",
        "score-ligand-3d",
        "score-ligand-3d-retry",
        "collate-ligand-3d",
        "merge-ligand-3d",
        "export-sucos-shards",
        "score-interface-shards",
        "collate-alignments",
        "collate-score-partitions",
        "symmetric-edge-fragments",
        "symmetric-edge-shards",
        "component-reductions",
        "communities",
        "directed-covers",
    ]:
        command = subparsers.add_parser(name)
        command.add_argument("data_dir", type=Path)
        command.add_argument("--batch-index", type=int, required=True)
        command.add_argument("--batch-size", type=int, required=True)
        command.add_argument("--threads", type=int, default=1)
        command.add_argument("--scratch-dir", type=Path, required=True)
        command.add_argument("--force", action="store_true")
        if name == "score-pdbs":
            command.add_argument("--pdb-manifest", type=Path, required=True)
        if name.startswith("repair-"):
            command.add_argument("--repair-manifest", type=Path, required=True)
        if name == "search":
            command.add_argument(
                "--alignment-type",
                action="append",
                choices=["foldseek", "mmseqs"],
                required=True,
            )
        if name in {"component-reductions", "communities", "directed-covers"}:
            _add_cluster_arguments(command)
        if name in {"export-sucos-shards", "repair-sucos-shards"}:
            command.add_argument("--output-dir", type=Path, required=True)
            command.add_argument("--memory-limit", default="16GB")
        if name == "score-interface-shards":
            command.add_argument("--memory-limit", default="32GB")

    finalizer = subparsers.add_parser("finalize-alignments")
    finalizer.add_argument("data_dir", type=Path)
    ligand_finalizer = subparsers.add_parser("finalize-ligands")
    ligand_finalizer.add_argument("data_dir", type=Path)
    score_finalizer = subparsers.add_parser("finalize-scores")
    score_finalizer.add_argument("data_dir", type=Path)
    repair_score_finalizer = subparsers.add_parser("validate-score-repair")
    repair_score_finalizer.add_argument("data_dir", type=Path)
    repair_score_finalizer.add_argument("--repair-manifest", type=Path, required=True)
    retry_finalizer = subparsers.add_parser("finalize-ligand-3d-retries")
    retry_finalizer.add_argument("data_dir", type=Path)
    sucos_finalizer = subparsers.add_parser("finalize-sucos-export")
    sucos_finalizer.add_argument("data_dir", type=Path)
    sucos_finalizer.add_argument("--source-dir", type=Path, required=True)
    sucos_finalizer.add_argument("--output", type=Path, required=True)
    sucos_finalizer.add_argument("--scratch-dir", type=Path, required=True)
    sucos_finalizer.add_argument("--threads", type=int, default=8)
    sucos_finalizer.add_argument("--memory-limit", default="32GB")
    interface_finalizer = subparsers.add_parser("finalize-interface-scores")
    interface_finalizer.add_argument("data_dir", type=Path)
    interface_finalizer.add_argument("--output", type=Path)
    interface_finalizer.add_argument("--scratch-dir", type=Path, required=True)
    interface_finalizer.add_argument("--threads", type=int, default=8)
    interface_finalizer.add_argument("--memory-limit", default="32GB")
    component_merger = subparsers.add_parser("merge-components")
    component_merger.add_argument("data_dir", type=Path)
    _add_cluster_arguments(component_merger)
    return parser


def main() -> None:
    args = _parser().parse_args()
    data_dir = args.data_dir.resolve()
    if args.command == "plan":
        result = plan_protein_scoring(
            data_dir,
            pdb_ids=args.pdb_id,
            two_char_codes=args.two_char_code,
            max_seqs=args.max_seqs,
        )
    elif args.command == "create-dbs":
        foldseek_input = make_foldseek_input_manifest(
            data_dir,
            args.cif_root.resolve(),
        )
        tasks.make_dbs(
            data_dir=data_dir,
            sub_databases=["holo"],
            cpu=args.threads,
            cif_root=foldseek_input,
            seqres_path=args.seqres_path.resolve(),
            create=True,
            index=False,
            force_update=args.force,
            build_dir=(
                args.scratch_dir.resolve() if args.scratch_dir is not None else None
            ),
        )
        result = {"status": "complete"}
    elif args.command == "make-sub-dbs":
        tasks.make_sub_dbs(
            data_dir=data_dir,
            sub_databases=["holo"],
            cpu=args.threads,
            scratch_dir=args.scratch_dir.resolve(),
        )
        result = {"status": "complete"}
    elif args.command == "refresh-alignment-lookup":
        lookup = data_dir / tasks.ALIGNMENT_CHAIN_LOOKUP_RELATIVE
        before = _source_signature(lookup) if lookup.is_file() else None
        tasks.make_alignment_chain_lookup(
            data_dir=data_dir,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            force_update=True,
        )
        after = _source_signature(lookup)
        if before is not None and after != before:
            raise RuntimeError(
                "alignment-chain lookup content changed; mapped alignment shards "
                "must be regenerated before score repair"
            )
        result = {"status": "complete", "lookup_unchanged": before == after}
    elif args.command == "plan-score-batches":
        result = plan_score_batches(
            data_dir,
            batch_size=args.batch_size,
            threads=args.threads,
            scratch_dir=args.scratch_dir.resolve(),
            max_query_protein_chains=args.max_query_protein_chains,
            max_query_proper_ligand_chains=(args.max_query_proper_ligand_chains),
        )
    elif args.command == "plan-ligand-3d":
        result = plan_ligand_3d_batches(
            data_dir,
            batch_size=args.batch_size,
            threads=args.threads,
            scratch_dir=args.scratch_dir.resolve(),
            memory_limit=args.memory_limit,
        )
    elif args.command == "plan-ligand-3d-retries":
        result = plan_ligand_3d_retries(data_dir, batch_size=args.batch_size)
    elif args.command == "plan-interface-scores":
        result = plan_interface_scoring(data_dir, batch_size=args.batch_size)
    elif args.command == "plan-score-repair":
        result = plan_score_repair(
            data_dir,
            affected_manifest=args.affected_manifest.resolve(),
            additional_full_query_manifest=(
                args.additional_full_query_manifest.resolve()
                if args.additional_full_query_manifest is not None
                else None
            ),
            output_path=args.output.resolve() if args.output is not None else None,
            batch_size=args.batch_size,
            target_batch_size=args.target_batch_size,
            threads=args.threads,
            scratch_dir=(
                args.scratch_dir.resolve() if args.scratch_dir is not None else None
            ),
            memory_limit=args.memory_limit,
        )
    elif args.command == "plan-score-repair-ligand-3d":
        result = plan_score_repair_ligand_3d(
            data_dir,
            repair_manifest=args.repair_manifest.resolve(),
            batch_size=args.batch_size,
            threads=args.threads,
            scratch_dir=args.scratch_dir.resolve(),
            memory_limit=args.memory_limit,
        )
    elif args.command == "finalize-score-repair-queries":
        result = finalize_score_repair_queries(
            data_dir,
            repair_manifest=args.repair_manifest.resolve(),
            max_new_drops=args.max_new_drops,
        )
    elif args.command == "plan-clusters":
        result = plan_clustering(
            data_dir,
            metrics=args.metrics,
            thresholds=args.thresholds,
            source_batch_size=args.source_batch_size,
            community_batch_size=args.community_batch_size,
            symmetric_bucket_count=args.symmetric_bucket_count,
        )
    elif args.command == "cluster-stats":
        result = summarize_clustering_artifacts(
            data_dir,
            metrics=args.metrics,
            thresholds=args.thresholds,
        )
    elif args.command == "finalize-index":
        tasks.finalize_index(data_dir=data_dir)
        result = {"status": "complete"}
    elif args.command == "drop-score-queries":
        result = record_dropped_queries(
            data_dir,
            score_query_manifest=args.pdb_manifest.resolve(),
        )
    elif args.command == "finalize-alignments":
        result = finalize_alignment_artifacts(data_dir)
    elif args.command == "finalize-ligands":
        result = finalize_ligand_archives(data_dir)
    elif args.command == "finalize-scores":
        result = finalize_ligand_3d_artifacts(data_dir)
    elif args.command == "validate-score-repair":
        result = finalize_score_repair_artifacts(
            data_dir,
            repair_manifest=args.repair_manifest.resolve(),
        )
    elif args.command == "finalize-ligand-3d-retries":
        result = finalize_ligand_3d_retries(data_dir)
    elif args.command == "finalize-sucos-export":
        result = finalize_sucos_shape_pocket_qcov_export(
            data_dir,
            source_dir=args.source_dir,
            output=args.output,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            memory_limit=args.memory_limit,
        )
    elif args.command == "finalize-interface-scores":
        result = finalize_interface_qcov_scores(
            data_dir,
            output=args.output.resolve() if args.output is not None else None,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            memory_limit=args.memory_limit,
        )
    elif args.command == "score-interface-shards":
        shards = _interface_score_shard_batch(
            data_dir,
            batch_index=args.batch_index,
            batch_size=args.batch_size,
        )
        result = score_interface_qcov_shards(
            data_dir,
            shards=shards,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            memory_limit=args.memory_limit,
            force_update=args.force,
        )
    elif args.command in {"export-sucos-shards", "repair-sucos-shards"}:
        explicit_shards = None
        batch_index = args.batch_index
        batch_size = args.batch_size
        if args.command == "repair-sucos-shards":
            explicit_shards = _score_repair_shard_batch(
                args.repair_manifest.resolve(),
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
            batch_index = None
            batch_size = None
        result = export_sucos_shape_pocket_qcov_batch(
            data_dir,
            output_dir=args.output_dir,
            batch_index=batch_index,
            batch_size=batch_size,
            shards=explicit_shards,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            memory_limit=args.memory_limit,
        )
    elif args.command == "symmetric-edge-fragments":
        batches = _symmetric_fragment_batch(
            data_dir,
            batch_index=args.batch_index,
            batch_size=args.batch_size,
        )
        tasks.make_symmetric_edge_fragments(
            data_dir=data_dir,
            batches=batches,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            force_update=args.force,
        )
        result = {
            "status": "complete",
            "batch_keys": [batch["key"] for batch in batches],
        }
    elif args.command == "symmetric-edge-shards":
        metric_buckets = _symmetric_edge_shard_batch(
            data_dir,
            batch_index=args.batch_index,
            batch_size=args.batch_size,
        )
        tasks.make_symmetric_edge_shards(
            data_dir=data_dir,
            metric_buckets=metric_buckets,
            scratch_dir=args.scratch_dir.resolve(),
            threads=args.threads,
            force_update=args.force,
        )
        result = {
            "status": "complete",
            "metric_buckets": metric_buckets,
        }
    elif args.command == "component-reductions":
        metrics, thresholds = _cluster_parameters(
            metrics=args.metrics,
            thresholds=args.thresholds,
        )
        sources = _component_source_batch(
            data_dir,
            metrics=metrics,
            batch_index=args.batch_index,
            batch_size=args.batch_size,
        )
        tasks.make_component_reductions(
            data_dir=data_dir,
            source_paths=sources,
            metrics=metrics,
            thresholds=thresholds,
            scratch_dir=args.scratch_dir.resolve(),
            force_update=args.force,
            metric_workers=args.threads,
        )
        result = {"status": "complete", "sources": sources}
    elif args.command == "merge-components":
        metrics, thresholds = _cluster_parameters(
            metrics=args.metrics,
            thresholds=args.thresholds,
        )
        tasks.merge_component_reductions(
            data_dir=data_dir,
            metrics=metrics,
            thresholds=thresholds,
        )
        result = {"status": "complete", "metrics": metrics}
    elif args.command in {"communities", "directed-covers"}:
        metrics, thresholds = _cluster_parameters(
            metrics=args.metrics,
            thresholds=args.thresholds,
        )
        work = _community_batch(
            metrics=metrics,
            thresholds=thresholds,
            batch_index=args.batch_index,
            batch_size=args.batch_size,
        )
        for metric_threshold in work:
            if args.command == "communities":
                tasks.make_communities(
                    data_dir=data_dir,
                    metric_threshold=[metric_threshold],
                    skip_existing_clusters=not args.force,
                    scratch_dir=args.scratch_dir.resolve(),
                    threads=args.threads,
                )
            else:
                tasks.make_directed_set_covers(
                    data_dir=data_dir,
                    metric_threshold=[metric_threshold],
                    skip_existing=not args.force,
                    scratch_dir=args.scratch_dir.resolve(),
                    threads=args.threads,
                )
        result = {"status": "complete", "metric_thresholds": work}
    elif args.command == "pack-ligands":
        shards = _entry_shard_batch(data_dir, args.batch_index, args.batch_size)
        tasks.make_canonical_ligand_archives(
            data_dir=data_dir,
            two_char_codes=shards,
            scratch_dir=args.scratch_dir.resolve(),
        )
        result = {"status": "complete", "shards": shards}
    else:
        plan = _load_plan(data_dir)
        cfg = _scoring_config(
            data_dir,
            int(plan["max_seqs"]),
            max_query_protein_chains=int(plan["score_max_query_protein_chains"]),
            max_query_proper_ligand_chains=int(
                plan["score_max_query_proper_ligand_chains"]
            ),
        )
        scratch_dir = args.scratch_dir.resolve()
        scratch_dir.mkdir(exist_ok=True, parents=True)
        if args.command == "collate-alignments":
            shards = _shard_batch(data_dir, args.batch_index, args.batch_size)
            for shard in shards:
                tasks.collate_alignments(
                    data_dir=data_dir,
                    partition=[shard],
                    scratch_dir=scratch_dir,
                    threads=args.threads,
                )
            result = {"status": "complete", "shards": shards}
        elif args.command == "repair-candidates":
            shards = _score_repair_shard_batch(
                args.repair_manifest.resolve(),
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
            outputs = tasks.collate_ligand_3d_candidates(
                data_dir=data_dir,
                shards=shards,
                scratch_dir=scratch_dir,
                threads=args.threads,
            )
            result = {
                "status": "complete",
                "shards": shards,
                "output_count": len(outputs),
            }
        elif args.command == "repair-pack-ligands":
            shards = _score_repair_shard_batch(
                args.repair_manifest.resolve(),
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
            tasks.make_canonical_ligand_archives(
                data_dir=data_dir,
                two_char_codes=shards,
                scratch_dir=scratch_dir,
            )
            result = {"status": "complete", "shards": shards}
        elif args.command == "repair-score-shards":
            shards = _score_repair_shard_batch(
                args.repair_manifest.resolve(),
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
            outputs = tasks.merge_ligand_3d_scores(
                data_dir=data_dir,
                shards=shards,
                scorer_cfg=cfg.scorer,
                force_update=True,
                scratch_dir=scratch_dir,
                threads=args.threads,
                reuse_cached_pairs=True,
            )
            result = {
                "status": "complete",
                "shards": shards,
                "output_count": len(outputs),
            }
        elif args.command == "repair-merge-ligand-3d":
            shards = _score_repair_shard_batch(
                args.repair_manifest.resolve(),
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
            result = merge_score_repair_ligand_3d(
                data_dir,
                shards=shards,
                scratch_dir=scratch_dir,
                threads=args.threads,
            )
        elif args.command == "repair-ligand-3d":
            pairs = _score_repair_ligand_3d_batch(
                data_dir,
                batch_index=args.batch_index,
                batch_size=args.batch_size,
            )
            output = None
            if not pairs.empty:
                output = tasks.make_ligand_3d_scores(
                    data_dir=data_dir,
                    pairs=pairs,
                    batch_index=args.batch_index,
                    scorer_cfg=cfg.scorer,
                    force_update=args.force,
                    scratch_dir=scratch_dir,
                    threads=args.threads,
                    output_path=(
                        data_dir
                        / "scores/ligand_3d_pair_repairs"
                        / f"{args.batch_index}.parquet"
                    ),
                )
            result = {
                "status": "complete",
                "pair_count": len(pairs),
                "output": str(output) if output is not None else None,
            }
        elif args.command == "collate-score-partitions":
            partitions = _score_partition_batch(args.batch_index, args.batch_size)
            for partition in partitions:
                tasks.collate_partitions(
                    data_dir=data_dir,
                    partition=[partition],
                    scratch_dir=scratch_dir,
                    threads=args.threads,
                )
            result = {"status": "complete", "partitions": partitions}
        elif args.command == "collate-ligand-3d-candidates":
            shards = _score_shard_batch(data_dir, args.batch_index, args.batch_size)
            outputs = tasks.collate_ligand_3d_candidates(
                data_dir=data_dir,
                shards=shards,
                scratch_dir=scratch_dir,
                threads=args.threads,
            )
            result = {
                "status": "complete",
                "shards": shards,
                "output_count": len(outputs),
            }
        elif args.command == "map":
            shards = _shard_batch(data_dir, args.batch_index, args.batch_size)
            tasks.map_batch_alignments(
                data_dir=data_dir,
                shards=shards,
                scorer_cfg=cfg.scorer,
                force_update=args.force,
                scratch_dir=scratch_dir,
            )
            result = {"status": "complete", "shards": shards}
        elif args.command == "score-ligand-3d":
            pairs = _ligand_3d_batch(data_dir, args.batch_index, args.batch_size)
            output = tasks.make_ligand_3d_scores(
                data_dir=data_dir,
                pairs=pairs,
                batch_index=args.batch_index,
                scorer_cfg=cfg.scorer,
                force_update=args.force,
                scratch_dir=scratch_dir,
                threads=args.threads,
            )
            result = {
                "status": "complete",
                "pair_count": len(pairs),
                "output": str(output),
            }
        elif args.command == "score-ligand-3d-retry":
            pairs = _ligand_3d_retry_batch(data_dir, args.batch_index, args.batch_size)
            output = tasks.make_ligand_3d_scores(
                data_dir=data_dir,
                pairs=pairs,
                batch_index=args.batch_index,
                scorer_cfg=cfg.scorer,
                force_update=args.force,
                scratch_dir=scratch_dir,
                threads=args.threads,
                output_path=(
                    data_dir
                    / "scores"
                    / "ligand_3d_pair_retries"
                    / f"{args.batch_index}.parquet"
                ),
            )
            result = {
                "status": "complete",
                "pair_count": len(pairs),
                "output": str(output),
            }
        elif args.command == "collate-ligand-3d":
            shards = _ligand_3d_shard_batch(data_dir, args.batch_index, args.batch_size)
            outputs = tasks.collate_ligand_3d_scores(
                data_dir=data_dir,
                shards=shards,
                scratch_dir=scratch_dir,
                threads=args.threads,
            )
            result = {
                "status": "complete",
                "shards": shards,
                "output_count": len(outputs),
            }
        elif args.command == "merge-ligand-3d":
            shards = _ligand_3d_shard_batch(data_dir, args.batch_index, args.batch_size)
            outputs = tasks.merge_ligand_3d_scores(
                data_dir=data_dir,
                shards=shards,
                scorer_cfg=cfg.scorer,
                force_update=args.force,
                scratch_dir=scratch_dir,
                threads=args.threads,
            )
            result = {
                "status": "complete",
                "shards": shards,
                "output_count": len(outputs),
            }
        else:
            if args.command == "repair-score-target-remainders":
                result = repair_target_score_remainders(
                    data_dir,
                    repair_manifest=args.repair_manifest.resolve(),
                    batch_index=args.batch_index,
                    batch_size=args.batch_size,
                    scorer_cfg=cfg.scorer,
                    scratch_dir=scratch_dir,
                    threads=args.threads,
                )
                print(json.dumps(result, indent=2, sort_keys=True))
                return
            if args.command == "repair-score-pdbs":
                repairs = _score_repair_batch(
                    args.repair_manifest.resolve(),
                    batch_index=args.batch_index,
                    batch_size=args.batch_size,
                )
                tasks.repair_batch_scores(
                    data_dir=data_dir,
                    repairs=repairs,
                    scorer_cfg=cfg.scorer,
                    scratch_dir=scratch_dir,
                    threads=args.threads,
                )
                result = {
                    "status": "complete",
                    "pdb_ids": [str(repair["pdb_id"]) for repair in repairs],
                }
                print(json.dumps(result, indent=2, sort_keys=True))
                return
            if args.command == "score":
                pdb_ids = _score_batch(data_dir, args.batch_index, args.batch_size)
            elif args.command == "score-pdbs":
                pdb_ids = _score_manifest_batch(
                    data_dir,
                    args.pdb_manifest.resolve(),
                    args.batch_index,
                    args.batch_size,
                )
            else:
                pdb_ids = _query_batch(data_dir, args.batch_index, args.batch_size)
            if args.command == "search":
                tasks.run_batch_searches(
                    data_dir=data_dir,
                    pdb_ids=pdb_ids,
                    scorer_cfg=cfg.scorer,
                    foldseek_cfg=cfg.foldseek,
                    mmseqs_cfg=cfg.mmseqs,
                    cpu=args.threads,
                    scratch_dir=scratch_dir,
                    alignment_types=args.alignment_type,
                    force_update=args.force,
                )
            elif args.command in {"score", "score-pdbs"}:
                tasks.make_batch_scores(
                    data_dir=data_dir,
                    pdb_ids=pdb_ids,
                    scorer_cfg=cfg.scorer,
                    force_update=args.force,
                    scratch_dir=scratch_dir,
                    threads=args.threads,
                    defer_ligand_3d=True,
                )
            result = {"status": "complete", "pdb_ids": pdb_ids}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
