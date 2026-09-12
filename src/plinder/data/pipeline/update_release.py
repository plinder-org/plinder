# Copyright (c) 2026, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Finish a planned PDB update in a new release directory."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from shutil import copyfile, rmtree
from typing import Any, Iterable, cast

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from omegaconf import DictConfig, OmegaConf

from plinder.core.release import RELEASE_PATHS
from plinder.core.utils import schemas
from plinder.core.utils.files import (
    file_sha256,
    link_or_copy_file,
    read_json_cache,
    write_json_atomic,
)
from plinder.data import clusters, databases, protein_clusters
from plinder.data.annotations import get_similarity_scores
from plinder.data.pipeline import config, score, tasks, utils
from plinder.data.pipeline.update_archives import update_ligand_archives
from plinder.data.pipeline.update_entries import apply_entry_update
from plinder.data.pipeline.updates import load_update_plan

_BASE_ARTIFACT_ROOTS = (
    "alignments",
    "ccd_dbs",
    "dbs",
    "exports",
    "fingerprints",
    "index",
    "interface_clusters",
    "interface_sampling",
    "interface_scores",
    "ligand_clusters",
    "ligand_sampling",
    "ligand_scores",
    "manifests",
    "mhfp6_scores",
    "protein_clusters",
    "scores",
    "search_databases",
)
_TRANSIENT_BASE_DIRECTORIES = {".staging", "ligand_3d_pair_repairs"}
_TRANSIENT_BASE_SUFFIXES = (".installing", ".previous", ".tmp", ".tmp.parquet")

_STAGES = (
    "entries_and_archives",
    "search_inputs",
    "alignments",
    "scores",
    "ligand_chemistry",
    "similarity_covers",
    "final_tables",
)


def _reuse_base_artifacts(base: Path, workspace: Path) -> None:
    """Hard-link unchanged downstream files without replacing prepared files."""
    for root_name in _BASE_ARTIFACT_ROOTS:
        source_root = base / root_name
        if not source_root.exists():
            continue
        for current, directories, filenames in os.walk(source_root, topdown=True):
            source_dir = Path(current)
            relative = source_dir.relative_to(base)
            target_dir = workspace / relative
            target_dir.mkdir(parents=True, exist_ok=True)
            for name in list(directories):
                if name in _TRANSIENT_BASE_DIRECTORIES or name.endswith(
                    _TRANSIENT_BASE_SUFFIXES
                ):
                    directories.remove(name)
                    continue
                source = source_dir / name
                target = target_dir / name
                if source.is_symlink():
                    if not target.exists() and not target.is_symlink():
                        link_or_copy_file(source, target)
                    directories.remove(name)
                else:
                    target.mkdir(exist_ok=True)
            for name in filenames:
                if name.endswith(_TRANSIENT_BASE_SUFFIXES):
                    continue
                source = source_dir / name
                target = target_dir / name
                if not target.exists() and not target.is_symlink():
                    link_or_copy_file(source, target)


def _chains_by_entry(chains: pd.DataFrame) -> dict[str, set[str]]:
    return {
        str(entry_id): set(rows["chain_auth_id"].dropna().astype(str))
        for entry_id, rows in chains.groupby("entry_pdb_id", sort=False)
    }


def _scoring_chains(data_dir: Path, search_db: str) -> pd.DataFrame:
    if search_db == "apo":
        return tasks._apo_scoring_chains(data_dir)
    chains = tasks._protein_scoring_chains(data_dir)
    auth_ids = pd.read_parquet(
        data_dir / "index" / "entry_chains.parquet",
        columns=["entry_pdb_id", "chain_asym_id", "chain_auth_id"],
    )
    return chains.merge(
        auth_ids,
        on=["entry_pdb_id", "chain_asym_id"],
        how="left",
        validate="one_to_one",
    )


def _query_chains(data_dir: Path, search_db: str) -> pd.DataFrame:
    return (
        tasks._linked_apo_query_chains(data_dir)
        if search_db == "apo"
        else _scoring_chains(data_dir, search_db)
    )


def _database_identifiers(chains: pd.DataFrame, alignment_type: str) -> set[str]:
    chains = chains.loc[chains["chain_auth_id"].notna()]
    if alignment_type == "foldseek":
        return {
            f"pdb_0000{row.entry_pdb_id}_xyz-enrich_{row.chain_auth_id}"
            for row in chains.itertuples(index=False)
        }
    return {
        f"{row.entry_pdb_id}_{row.chain_auth_id}"
        for row in chains.itertuples(index=False)
    }


def _queries_with_old_target_hits(
    data_dir: Path, *, search_db: str, affected: set[str]
) -> set[str]:
    if not affected:
        return set()
    paths = [
        path
        for alignment_type in ("foldseek", "mmseqs")
        for path in sorted(
            (
                data_dir / "dbs" / "subdbs" / f"{search_db}_{alignment_type}" / "aln"
            ).glob("*.parquet")
        )
    ]
    if not paths:
        return set()
    connection = duckdb.connect()
    try:
        connection.register(
            "affected_targets",
            pd.DataFrame({"target_pdb_id": sorted(affected)}),
        )
        sources = ", ".join(f"'{path.as_posix()}'" for path in paths)
        hits = connection.sql(
            f"""
            SELECT DISTINCT regexp_extract(filename, '/([^/]+)\\.parquet$', 1)
                AS query_entry
            FROM read_parquet(
                [{sources}], union_by_name=true, filename=true
            ) AS alignments
            INNER JOIN affected_targets USING (target_pdb_id)
            """
        ).df()
    finally:
        connection.close()
    return set(hits["query_entry"].dropna().astype(str))


def _changed_target_database(
    data_dir: Path,
    *,
    search_db: str,
    affected: set[str],
    scratch_dir: Path,
    threads: int,
) -> Path | None:
    chains = _scoring_chains(data_dir, search_db)
    chains = chains.loc[chains["entry_pdb_id"].astype(str).isin(affected)]
    identifiers = {
        f"{search_db}_{alignment_type}": _database_identifiers(chains, alignment_type)
        for alignment_type in ("foldseek", "mmseqs")
    }
    if not any(identifiers.values()):
        return None
    target_dir = scratch_dir / f"changed-{search_db}-targets"
    if target_dir.exists():
        rmtree(target_dir)
    databases.make_sub_dbs(
        target_dir,
        utils.get_db_sources(data_dir=data_dir, sub_databases=[search_db]),
        identifiers_by_database=identifiers,
        tmp_dir=scratch_dir / f"changed-{search_db}-target-work",
        threads=threads,
    )
    return target_dir


def _search_changed_targets(
    data_dir: Path,
    *,
    search_db: str,
    query_ids: list[str],
    target_database_dir: Path | None,
    scorer_cfg: DictConfig,
    foldseek_cfg: DictConfig,
    mmseqs_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    batch_size: int,
) -> set[str]:
    if target_database_dir is None or not query_ids:
        return set()
    query_auth_ids = _chains_by_entry(_query_chains(data_dir, search_db))
    hits: set[str] = set()
    for start in range(0, len(query_ids), batch_size):
        batch = query_ids[start : start + batch_size]
        batch_root = scratch_dir / f"reverse-{search_db}-{start // batch_size:05d}"
        scorer, _, work = utils.get_scorer(
            data_dir=data_dir,
            pdb_ids=batch,
            scorer_cfg=scorer_cfg,
            load_entries=False,
            foldseek_cfg=foldseek_cfg,
            mmseqs_cfg=mmseqs_cfg,
            scratch_dir=batch_root / "queries",
        )
        try:
            hits.update(
                scorer.run_alignments(
                    entry_ids=batch,
                    search_db=search_db,
                    output_folder=work,
                    threads=threads,
                    query_chain_auth_ids=query_auth_ids,
                    target_database_dir=target_database_dir,
                    result_database_dir=batch_root / "results",
                    write_empty_results=False,
                )
            )
        finally:
            rmtree(batch_root, ignore_errors=True)
    return hits


def _remove_inactive_queries(
    data_dir: Path, *, search_db: str, pdb_ids: Iterable[str]
) -> None:
    for pdb_id in pdb_ids:
        for alignment_type in ("foldseek", "mmseqs"):
            (
                data_dir
                / "dbs"
                / "subdbs"
                / f"{search_db}_{alignment_type}"
                / "aln"
                / f"{pdb_id}.parquet"
            ).unlink(missing_ok=True)
        (
            data_dir / "dbs" / "subdbs" / f"search_db={search_db}" / f"{pdb_id}.parquet"
        ).unlink(missing_ok=True)


def _refresh_alignment_shards(
    data_dir: Path,
    *,
    search_db: str,
    shards: set[str],
    scorer_cfg: DictConfig,
    scratch_dir: Path,
) -> None:
    for shard in sorted(shards):
        raw_exists = any(
            any(
                path.stem[1:3] == shard
                for path in (
                    data_dir
                    / "dbs"
                    / "subdbs"
                    / f"{search_db}_{alignment_type}"
                    / "aln"
                ).glob("*.parquet")
            )
            for alignment_type in ("foldseek", "mmseqs")
        )
        if raw_exists:
            tasks.map_batch_alignments(
                data_dir=data_dir,
                shards=[shard],
                scorer_cfg=scorer_cfg,
                force_update=True,
                scratch_dir=scratch_dir,
                search_db=search_db,
            )
            continue
        for alignment_type in ("foldseek", "mmseqs"):
            (
                data_dir
                / "alignments"
                / f"search_db={search_db}"
                / f"alignment_type={alignment_type}"
                / f"shard={shard}.parquet"
            ).unlink(missing_ok=True)
        manifest_root = data_dir / "alignments" / "manifests"
        if search_db != "holo":
            manifest_root = manifest_root / f"search_db={search_db}"
        (manifest_root / f"shard={shard}.json").unlink(missing_ok=True)


def _rebase_unchanged_alignment_manifests(
    data_dir: Path, *, search_db: str, repaired_shards: set[str]
) -> None:
    """Point untouched shard reports at the updated chain lookup."""
    lookup = tasks._completed_alignment_chain_lookup(data_dir)
    if lookup is None:
        raise ValueError("alignment chain lookup is incomplete")
    manifest_root = data_dir / "alignments" / "manifests"
    if search_db != "holo":
        manifest_root = manifest_root / f"search_db={search_db}"
    for path in sorted(manifest_root.glob("shard=*.json")):
        shard = path.stem.removeprefix("shard=")
        if shard in repaired_shards:
            continue
        payload = json.loads(path.read_text())
        inputs = tasks._alignment_input_signatures(
            data_dir=data_dir, search_db=search_db, shard=shard
        )
        if payload.get("inputs") != inputs:
            raise ValueError(f"untouched raw alignment inputs changed: {path}")
        outputs = payload.get("outputs")
        if not isinstance(outputs, dict):
            raise ValueError(f"invalid alignment report: {path}")
        for alignment_type, signatures in inputs.items():
            output = tasks._alignment_release_path(
                data_dir=data_dir,
                search_db=search_db,
                alignment_type=alignment_type,
                shard=shard,
            )
            if not signatures:
                if outputs.get(alignment_type) is not None:
                    raise ValueError(f"unexpected mapped alignment output: {output}")
                continue
            if not output.is_file() or not utils._mapped_alignment_file_is_current(
                output, alignment_type=alignment_type
            ):
                raise ValueError(f"invalid mapped alignment output: {output}")
            stat = output.stat()
            if outputs.get(alignment_type) != {
                "name": output.name,
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }:
                raise ValueError(f"mapped alignment report changed: {output}")
        payload["alignment_chain_lookup"] = lookup
        write_json_atomic(path, payload)


def _plan_alignment_repairs(
    data_dir: Path,
    *,
    search_db: str,
    affected: set[str],
    scorer_cfg: DictConfig,
    foldseek_cfg: DictConfig,
    mmseqs_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    batch_size: int,
) -> set[str]:
    """Find queries whose capped search can change in either direction."""
    manifest_path = (
        data_dir / score.LINKED_APO_QUERY_MANIFEST_RELATIVE
        if search_db == "apo"
        else data_dir / score.MANIFEST_RELATIVE
    )
    active = set(
        pd.read_parquet(manifest_path, columns=["pdb_id"])["pdb_id"].astype(str)
    )
    if search_db == "pred":
        return active.intersection(affected)
    old_hits = _queries_with_old_target_hits(
        data_dir, search_db=search_db, affected=affected
    )
    changed_targets = _changed_target_database(
        data_dir,
        search_db=search_db,
        affected=affected,
        scratch_dir=scratch_dir,
        threads=threads,
    )
    reverse_hits = _search_changed_targets(
        data_dir,
        search_db=search_db,
        query_ids=sorted(active.difference(affected)),
        target_database_dir=changed_targets,
        scorer_cfg=scorer_cfg,
        foldseek_cfg=foldseek_cfg,
        mmseqs_cfg=mmseqs_cfg,
        scratch_dir=scratch_dir,
        threads=threads,
        batch_size=batch_size,
    )
    return active.intersection(affected | old_hits | reverse_hits)


def _load_or_plan_alignment_repairs(
    data_dir: Path,
    *,
    search_databases: list[str],
    affected: set[str],
    scorer_cfg: DictConfig,
    foldseek_cfg: DictConfig,
    mmseqs_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    batch_size: int,
    state_path: Path,
    report: dict[str, Any],
) -> dict[str, set[str]]:
    """Persist the full repair set before any raw alignment is replaced."""
    stored = report.get("alignment_repair_queries")
    if stored is None:
        planned = {
            search_db: sorted(
                _plan_alignment_repairs(
                    data_dir,
                    search_db=search_db,
                    affected=affected,
                    scorer_cfg=scorer_cfg,
                    foldseek_cfg=foldseek_cfg,
                    mmseqs_cfg=mmseqs_cfg,
                    scratch_dir=scratch_dir / search_db,
                    threads=threads,
                    batch_size=batch_size,
                )
            )
            for search_db in search_databases
        }
        report["alignment_repair_queries"] = planned
        write_json_atomic(state_path, report)
        stored = planned
    if not isinstance(stored, dict) or set(stored) != set(search_databases):
        raise ValueError("weekly alignment repair plan does not match configuration")
    return {
        search_db: set(map(str, stored[search_db])) for search_db in search_databases
    }


def repair_alignments(
    data_dir: Path,
    *,
    search_db: str,
    affected: set[str],
    full_queries: set[str],
    scorer_cfg: DictConfig,
    foldseek_cfg: DictConfig,
    mmseqs_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    batch_size: int,
) -> set[str]:
    """Replace the previously planned alignment queries and their mapped shards."""
    manifest_path = (
        data_dir / score.LINKED_APO_QUERY_MANIFEST_RELATIVE
        if search_db == "apo"
        else data_dir / score.MANIFEST_RELATIVE
    )
    active = set(
        pd.read_parquet(manifest_path, columns=["pdb_id"])["pdb_id"].astype(str)
    )
    full_queries = active.intersection(full_queries)
    inactive = affected.difference(active)
    _remove_inactive_queries(data_dir, search_db=search_db, pdb_ids=inactive)
    selected_cfg = cast(
        DictConfig, OmegaConf.merge(scorer_cfg, {"sub_databases": [search_db]})
    )
    for start in range(0, len(full_queries), batch_size):
        tasks.run_batch_searches(
            data_dir=data_dir,
            pdb_ids=sorted(full_queries)[start : start + batch_size],
            scorer_cfg=selected_cfg,
            foldseek_cfg=foldseek_cfg,
            mmseqs_cfg=mmseqs_cfg,
            cpu=threads,
            scratch_dir=scratch_dir / f"full-{search_db}-{start // batch_size:05d}",
            force_update=True,
        )
    touched_shards = {pdb_id[1:3] for pdb_id in full_queries | inactive}
    _refresh_alignment_shards(
        data_dir,
        search_db=search_db,
        shards=touched_shards,
        scorer_cfg=selected_cfg,
        scratch_dir=scratch_dir / f"map-{search_db}",
    )
    _rebase_unchanged_alignment_manifests(
        data_dir, search_db=search_db, repaired_shards=touched_shards
    )
    return full_queries


def _write_pdb_manifest(path: Path, pdb_ids: Iterable[str]) -> Path:
    path.parent.mkdir(exist_ok=True, parents=True)
    frame = pd.DataFrame({"pdb_id": sorted(set(map(str, pdb_ids)))})
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)
    return path


def _remove_affected_shape_scores(
    data_dir: Path,
    *,
    shards: Iterable[str],
    affected: set[str],
    scratch_dir: Path,
    threads: int,
) -> int:
    """Remove cached ligand-pair scores whose coordinates may have changed."""
    if not affected:
        return 0
    scratch_dir.mkdir(parents=True, exist_ok=True)
    affected_entries = pd.DataFrame({"pdb_id": sorted(affected)})
    removed = 0
    with duckdb.connect() as connection:
        connection.sql(f"SET threads={threads}")
        connection.sql(f"SET temp_directory='{scratch_dir.as_posix()}'")
        connection.register("affected_entries", affected_entries)
        for shard in sorted(set(map(str, shards))):
            cached = data_dir / "scores/ligand_3d_by_query" / f"{shard}.parquet"
            if not cached.is_file():
                continue
            temporary = scratch_dir / f"{shard}.parquet"
            temporary.unlink(missing_ok=True)
            previous_rows = pq.ParquetFile(cached).metadata.num_rows
            connection.execute(
                f"""
                COPY (
                    SELECT scores.*
                    FROM read_parquet('{cached.as_posix()}') AS scores
                    ANTI JOIN affected_entries AS query_entries
                      ON scores.query_entry = query_entries.pdb_id
                    ANTI JOIN affected_entries AS target_entries
                      ON scores.target_entry = target_entries.pdb_id
                    ORDER BY
                        query_entry,
                        query_ligand_asym_id,
                        target_entry,
                        target_ligand_asym_id
                ) TO '{temporary.as_posix()}' (
                    FORMAT PARQUET, COMPRESSION ZSTD
                )
                """
            )
            if not pq.read_schema(temporary).equals(schemas.LIGAND_3D_SCORE_SCHEMA):
                temporary.unlink(missing_ok=True)
                raise ValueError(
                    f"shape-score cache has an unexpected schema: {cached}"
                )
            current_rows = pq.ParquetFile(temporary).metadata.num_rows
            install = cached.with_suffix(cached.suffix + ".tmp")
            copyfile(temporary, install)
            install.replace(cached)
            temporary.unlink(missing_ok=True)
            removed += previous_rows - current_rows
    return removed


def repair_holo_scores(
    data_dir: Path,
    *,
    base_data_dir: Path,
    affected: set[str],
    full_alignment_queries: set[str],
    scorer_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    memory_limit: str,
    score_batch_size: int,
    ligand_batch_size: int,
) -> dict[str, Any]:
    """Repair ligand-pocket scores and cached ligand-pair shape scores."""
    rmtree(data_dir / "scores/ligand_3d_pair_repairs", ignore_errors=True)
    holo_cfg = cast(
        DictConfig, OmegaConf.merge(scorer_cfg, {"sub_databases": ["holo"]})
    )
    score.plan_score_batches(
        data_dir,
        batch_size=score_batch_size,
        threads=threads,
        scratch_dir=scratch_dir / "plan-scores",
        max_query_protein_chains=int(scorer_cfg.max_query_protein_chains),
        max_query_proper_ligand_chains=int(scorer_cfg.max_query_proper_ligand_chains),
        reuse_mapped_alignments=True,
    )
    affected_manifest = _write_pdb_manifest(
        data_dir / "manifests" / "weekly_affected_entries.parquet", affected
    )
    full_manifest = _write_pdb_manifest(
        data_dir / "manifests" / "weekly_full_score_queries.parquet",
        full_alignment_queries,
    )
    try:
        plan = score.plan_score_repair(
            data_dir,
            affected_manifest=affected_manifest,
            additional_full_query_manifest=full_manifest,
            batch_size=score_batch_size,
            target_batch_size=score_batch_size,
            threads=threads,
            scratch_dir=scratch_dir / "plan-repair",
            memory_limit=memory_limit,
        )
    except ValueError as exc:
        if str(exc) != "no active score queries are affected by the repair":
            raise
        return {"status": "complete", "query_count": 0, "pair_count": 0}

    repair_manifest = data_dir / score.SCORE_REPAIR_RELATIVE
    repairs = pd.read_parquet(repair_manifest)
    for batch_index, batch in repairs.groupby("repair_batch_index", sort=True):
        tasks.repair_batch_scores(
            data_dir=data_dir,
            repairs=batch.to_dict("records"),
            scorer_cfg=holo_cfg,
            scratch_dir=scratch_dir / "queries" / str(batch_index),
            threads=threads,
        )
    score.finalize_score_repair_queries(
        data_dir,
        repair_manifest=repair_manifest,
        max_new_drops=0,
    )

    replacement_queries = set(repairs["pdb_id"].astype(str))
    active_replacements = replacement_queries.intersection(
        score.published_scoring_query_ids(data_dir)
    )
    shards = sorted(set(repairs["shard"].astype(str)))
    existing_shards: list[str] = []
    new_shards: list[str] = []
    for shard in shards:
        base_files = [
            base_data_dir
            / "scores"
            / "ligand_3d_candidate_shards"
            / f"shard={shard}.parquet",
            base_data_dir
            / "scores"
            / "ligand_pair_score_shards"
            / f"shard={shard}.parquet",
            base_data_dir / "scores" / "search_db=holo" / f"{shard}.parquet",
        ]
        if all(path.is_file() for path in base_files):
            existing_shards.append(shard)
        elif any(path.is_file() for path in base_files):
            raise FileNotFoundError(
                f"base release has an incomplete packed score shard {shard}: "
                f"{[str(path) for path in base_files if not path.is_file()]}"
            )
        else:
            new_shards.append(shard)
    if existing_shards:
        tasks.collate_ligand_3d_candidates(
            data_dir=data_dir,
            shards=existing_shards,
            scratch_dir=scratch_dir / "candidate-shards",
            threads=threads,
            replacement_query_ids=replacement_queries,
            source_query_ids=active_replacements,
        )
    if new_shards:
        tasks.collate_ligand_3d_candidates(
            data_dir=data_dir,
            shards=new_shards,
            scratch_dir=scratch_dir / "candidate-shards",
            threads=threads,
        )
    pair_cache_root = data_dir / "scores" / "ligand_3d_by_query"
    pair_cache_root.mkdir(exist_ok=True, parents=True)
    for shard in new_shards:
        pair_cache = pair_cache_root / f"{shard}.parquet"
        if not pair_cache.is_file():
            pq.write_table(
                pa.Table.from_pylist([], schema=schemas.LIGAND_3D_SCORE_SCHEMA),
                pair_cache,
            )
    removed_shape_scores = _remove_affected_shape_scores(
        data_dir,
        shards=shards,
        affected=affected,
        scratch_dir=scratch_dir / "invalidate-ligand-3d",
        threads=threads,
    )
    ligand_plan = score.plan_score_repair_ligand_3d(
        data_dir,
        repair_manifest=repair_manifest,
        batch_size=ligand_batch_size,
        threads=threads,
        scratch_dir=scratch_dir / "plan-ligand-3d",
        memory_limit=memory_limit,
    )
    repair_pair_dir = data_dir / "scores" / "ligand_3d_pair_repairs"
    for batch_index in range(int(ligand_plan["batch_count"])):
        pairs = score._score_repair_ligand_3d_batch(
            data_dir,
            batch_index=batch_index,
            batch_size=ligand_batch_size,
        )
        tasks.make_ligand_3d_scores(
            data_dir=data_dir,
            pairs=pairs,
            batch_index=batch_index,
            scorer_cfg=holo_cfg,
            force_update=False,
            scratch_dir=scratch_dir / "ligand-3d" / str(batch_index),
            threads=threads,
            output_path=repair_pair_dir / f"{batch_index}.parquet",
        )
    score.merge_score_repair_ligand_3d(
        data_dir,
        shards=shards,
        scratch_dir=scratch_dir / "merge-ligand-3d",
        threads=threads,
    )
    if existing_shards:
        tasks.merge_ligand_3d_scores(
            data_dir=data_dir,
            shards=existing_shards,
            scorer_cfg=holo_cfg,
            force_update=True,
            scratch_dir=scratch_dir / "merge-score-shards",
            threads=threads,
            reuse_cached_pairs=True,
            replacement_query_ids=replacement_queries,
        )
    if new_shards:
        tasks.merge_ligand_3d_scores(
            data_dir=data_dir,
            shards=new_shards,
            scorer_cfg=holo_cfg,
            force_update=True,
            scratch_dir=scratch_dir / "merge-score-shards",
            threads=threads,
            reuse_cached_pairs=True,
        )
    repair_report = score.finalize_score_repair_artifacts(
        data_dir,
        repair_manifest=repair_manifest,
    )
    active_shards = {
        pdb_id[1:3] for pdb_id in score.published_scoring_query_ids(data_dir)
    }
    for shard in set(shards).difference(active_shards):
        for path in [
            data_dir / "scores/ligand_3d_candidate_shards" / f"shard={shard}.parquet",
            data_dir / "scores/ligand_3d_candidate_shards" / f"shard={shard}.json",
            data_dir / "scores/ligand_pair_score_shards" / f"shard={shard}.parquet",
            data_dir
            / "scores/ligand_3d_pair_candidate_shards"
            / f"shard={shard}.parquet",
            data_dir / "scores/ligand_3d_by_query" / f"{shard}.parquet",
            data_dir / "scores/search_db=holo" / f"{shard}.parquet",
        ]:
            path.unlink(missing_ok=True)
    return {
        **plan,
        "validation": repair_report,
        "pair_count": int(ligand_plan["pair_count"]),
        "removed_shape_scores": removed_shape_scores,
        "repaired_shards": sorted(active_shards.intersection(shards)),
    }


def repair_non_holo_scores(
    data_dir: Path,
    *,
    search_db: str,
    affected: set[str],
    full_alignment_queries: set[str],
    scorer_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    memory_limit: str,
) -> dict[str, Any]:
    """Replace protein scores for an enabled apo or predicted target database."""
    if search_db not in {"apo", "pred"}:
        raise ValueError(f"expected apo or pred score database, got {search_db}")
    query_manifest = (
        data_dir / score.LINKED_APO_QUERY_MANIFEST_RELATIVE
        if search_db == "apo"
        else data_dir / score.MANIFEST_RELATIVE
    )
    active = set(
        pd.read_parquet(query_manifest, columns=["pdb_id"])["pdb_id"].astype(str)
    )
    queries = sorted(full_alignment_queries.intersection(active))
    selected_cfg = cast(
        DictConfig, OmegaConf.merge(scorer_cfg, {"sub_databases": [search_db]})
    )
    if queries:
        tasks.make_batch_scores(
            data_dir=data_dir,
            pdb_ids=queries,
            scorer_cfg=selected_cfg,
            force_update=True,
            scratch_dir=scratch_dir / f"{search_db}-scores",
            threads=threads,
            defer_ligand_3d=False,
        )
    _remove_inactive_queries(
        data_dir,
        search_db=search_db,
        pdb_ids=affected.difference(active),
    )
    source_dir = data_dir / "dbs" / "subdbs" / f"search_db={search_db}"
    output = data_dir / "scores" / f"search_db={search_db}" / f"{search_db}.parquet"
    if any(source_dir.glob("*.parquet")):
        tasks.collate_partitions(
            data_dir=data_dir,
            partition=[search_db],
            scratch_dir=scratch_dir / f"{search_db}-partition",
            threads=threads,
            memory_limit=memory_limit,
        )
    else:
        output.unlink(missing_ok=True)
    return {
        "status": "complete",
        "query_count": len(queries),
        "score_path": str(output),
    }


def repair_apo_scores(
    data_dir: Path,
    *,
    affected: set[str],
    full_alignment_queries: set[str],
    scorer_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    memory_limit: str,
) -> dict[str, Any]:
    """Recompute the lightweight protein scores used for linked apo chains."""
    report = repair_non_holo_scores(
        data_dir,
        search_db="apo",
        affected=affected,
        full_alignment_queries=full_alignment_queries,
        scorer_cfg=scorer_cfg,
        scratch_dir=scratch_dir,
        threads=threads,
        memory_limit=memory_limit,
    )
    output = tasks.make_linked_apo_structures(
        data_dir=data_dir,
        scratch_dir=scratch_dir / "linked-apo",
        threads=threads,
        memory_limit=memory_limit,
    )
    return {
        **report,
        "linked_apo_structures": str(output),
    }


def repair_interface_scores(
    data_dir: Path,
    *,
    affected: set[str],
    full_alignment_queries: set[str],
    scratch_dir: Path,
    threads: int,
    memory_limit: str,
) -> dict[str, Any]:
    """Replace interface rows for changed searches and retain other query rows."""
    plan = score.plan_interface_scoring(data_dir, batch_size=1)
    work = pd.read_parquet(data_dir / score.INTERFACE_SCORE_WORK_RELATIVE)
    representatives = pd.read_parquet(
        data_dir / tasks.INTERFACE_REPRESENTATIVES_RELATIVE,
        columns=["entry_pdb_id"],
    )
    planned_queries = set(
        pd.read_parquet(data_dir / score.MANIFEST_RELATIVE, columns=["pdb_id"])[
            "pdb_id"
        ].astype(str)
    )
    alignment_report = json.loads(
        (data_dir / "alignments" / "manifest.json").read_text()
    )
    skipped_queries = set(map(str, alignment_report.get("skipped_queries", {})))
    current_queries = (
        set(representatives["entry_pdb_id"].astype(str))
        & planned_queries - skipped_queries
    )
    replacements = affected | full_alignment_queries
    active_replacements = replacements.intersection(current_queries)
    repair_root = data_dir / "interface_score_repairs" / "weekly"
    repair_paths: dict[str, Path] = {}
    replacements_by_shard: dict[str, set[str]] = {}
    for pdb_id in active_replacements:
        replacements_by_shard.setdefault(pdb_id[1:3], set()).add(pdb_id)
    for shard, entries in sorted(replacements_by_shard.items()):
        output = repair_root / f"shard={shard}.parquet"
        score.score_interface_qcov_shards(
            data_dir,
            shards=[str(shard)],
            scratch_dir=scratch_dir / "interface-parts" / str(shard),
            threads=threads,
            memory_limit=memory_limit,
            force_update=True,
            query_entries_by_shard={str(shard): entries},
            output_paths_by_shard={str(shard): output},
        )
        repair_paths[str(shard)] = output

    score_root = data_dir / score.INTERFACE_SCORE_ROOT_RELATIVE
    active_shards = set(work["shard"].astype(str))
    for path in score_root.glob("shard=*.parquet"):
        if path.stem.removeprefix("shard=") not in active_shards:
            path.unlink()
            path.with_suffix(".json").unlink(missing_ok=True)
    reports: list[dict[str, Any]] = []
    for shard in work["shard"].astype(str):
        output = score_root / f"shard={shard}.parquet"
        repair_path = repair_paths.get(shard)
        removed = sorted(pdb_id for pdb_id in replacements if pdb_id[1:3] == shard)
        local_root = scratch_dir / "interface-merge" / shard
        local_root.mkdir(parents=True, exist_ok=True)
        local_output = local_root / output.name
        local_output.unlink(missing_ok=True)
        connection = duckdb.connect()
        try:
            connection.execute(f"SET threads={threads}")
            connection.execute(f"SET temp_directory='{local_root.as_posix()}'")
            connection.execute(f"SET memory_limit='{memory_limit}'")
            selects = []
            if output.is_file():
                if removed:
                    connection.register(
                        "replacement_queries",
                        pd.DataFrame({"entry_pdb_id": removed}),
                    )
                    selects.append(
                        f"""
                        SELECT existing.*
                        FROM read_parquet('{output.as_posix()}') AS existing
                        ANTI JOIN replacement_queries
                          ON split_part(existing.query_system, '__', 1)
                           = replacement_queries.entry_pdb_id
                        """
                    )
                else:
                    selects.append(f"SELECT * FROM read_parquet('{output.as_posix()}')")
            if repair_path is not None:
                selects.append(
                    f"SELECT * FROM read_parquet('{repair_path.as_posix()}')"
                )
            if not selects:
                raise FileNotFoundError(f"interface shard has no source rows: {shard}")
            combined = " UNION ALL BY NAME ".join(selects)
            connection.execute(
                f"""
                COPY (
                    SELECT
                        query_system::VARCHAR AS query_system,
                        target_system::VARCHAR AS target_system,
                        mapping::VARCHAR AS mapping,
                        source::VARCHAR AS source,
                        metric::VARCHAR AS metric,
                        iface1_qcov::FLOAT AS iface1_qcov,
                        iface2_qcov::FLOAT AS iface2_qcov,
                        similarity::TINYINT AS similarity
                    FROM ({combined})
                    ORDER BY query_system, target_system, metric
                ) TO '{local_output.as_posix()}' (
                    FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500000
                )
                """
            )
            duplicate = connection.execute(
                f"""
                SELECT count(*) FROM (
                    SELECT query_system, target_system, metric
                    FROM read_parquet('{local_output.as_posix()}')
                    GROUP BY ALL HAVING count(*) > 1
                )
                """
            ).fetchone()
        finally:
            connection.close()
        if duplicate is None or int(duplicate[0]):
            raise ValueError(f"interface repair produced duplicate rows: {shard}")
        if not pq.read_schema(local_output).equals(
            schemas.INTERFACE_SCORE_SHARD_SCHEMA
        ):
            raise ValueError(f"interface repair produced the wrong schema: {shard}")
        output.parent.mkdir(parents=True, exist_ok=True)
        install = output.with_suffix(output.suffix + ".tmp")
        copyfile(local_output, install)
        install.replace(output)
        local_output.unlink(missing_ok=True)
        payload = {
            "status": "complete",
            "shard": shard,
            "inputs": score._interface_score_inputs(data_dir, shard),
            "output": score._source_signature(output),
            "rows": int(pq.ParquetFile(output).metadata.num_rows),
        }
        write_json_atomic(output.with_suffix(".json"), payload)
        reports.append(payload)
    export = score.finalize_interface_similarity_scores(
        data_dir,
        scratch_dir=scratch_dir / "interface-export",
        threads=threads,
        memory_limit=memory_limit,
    )
    return {
        **plan,
        "repaired_query_count": len(active_replacements),
        "score_rows": sum(int(report["rows"]) for report in reports),
        "export_rows": int(export["row_count"]),
    }


def _chemical_score_query_ids(paths: list[Path], metric: str) -> set[int]:
    if not paths:
        return set()
    sources = ", ".join(f"'{path.as_posix()}'" for path in paths)
    with duckdb.connect() as connection:
        rows = connection.sql(
            f"""
            SELECT DISTINCT query_ligand_id
            FROM read_parquet([{sources}], union_by_name=true)
            WHERE query_ligand_id = target_ligand_id AND {metric} = 100
            """
        ).fetchall()
    return {int(row[0]) for row in rows}


def _refresh_chemical_score_shards(
    *,
    data_dir: Path,
    directory: str,
    metric: str,
    batch_size: int,
    number_id_col: str,
    minimum_similarity: float,
    scratch_dir: Path,
    prior_ligands: pd.DataFrame,
    prior_manifest: dict[str, str | float],
    current_manifest: dict[str, str | float],
) -> dict[str, int | str]:
    """Keep old chemical edges and calculate both directions for new nodes."""
    if batch_size < 1:
        raise ValueError("chemical score batch size must be positive")
    root = data_dir / directory
    stage = data_dir / f".{directory}.pending"
    backup = data_dir / f".{directory}.previous"
    if backup.exists():
        if not root.exists():
            backup.rename(root)
        else:
            rmtree(backup)
    rmtree(stage, ignore_errors=True)
    stage.mkdir(parents=True)
    smiles_column = "ligand_rdkit_canonical_smiles"
    fingerprint_column = "fingerprint" if directory == "ligand_scores" else "mhfp6"
    fingerprint_file = (
        "ligands_per_smiles.parquet"
        if directory == "ligand_scores"
        else get_similarity_scores.MHFP6_FINGERPRINT_FILE
    )
    current_ligands = pd.read_parquet(
        data_dir / "fingerprints" / fingerprint_file,
        columns=[number_id_col, smiles_column, fingerprint_column],
    )
    for frame, label in ((prior_ligands, "prior"), (current_ligands, "current")):
        missing = sorted(
            {number_id_col, smiles_column, fingerprint_column}.difference(frame.columns)
        )
        if missing:
            raise ValueError(f"{label} ligand fingerprints are missing {missing}")
        if (
            frame[number_id_col].duplicated().any()
            or frame[smiles_column].duplicated().any()
        ):
            raise ValueError(f"{label} ligand fingerprint IDs are not unique")
    active_ids = set(map(int, current_ligands[number_id_col]))
    id_mapping = prior_ligands.merge(
        current_ligands,
        on=[smiles_column, fingerprint_column],
        how="inner",
        suffixes=("_old", "_new"),
        validate="one_to_one",
    ).rename(
        columns={
            f"{number_id_col}_old": "old_id",
            f"{number_id_col}_new": "new_id",
        }
    )[["old_id", "new_id"]]
    id_mapping = id_mapping.astype({"old_id": "int32", "new_id": "int32"})
    prior_paths = sorted(root.glob("*.parquet"))
    if (
        read_json_cache(root / get_similarity_scores.LIGAND_SCORE_MANIFEST)
        != prior_manifest
    ):
        prior_paths = []
    prior_queries = _chemical_score_query_ids(prior_paths, metric)
    complete_mapping = id_mapping.loc[id_mapping["old_id"].isin(prior_queries)]
    mapped_queries = set(map(int, complete_mapping["new_id"]))
    missing_queries = active_ids.difference(mapped_queries)
    removed_ids = set(map(int, prior_ligands[number_id_col])).difference(
        set(map(int, id_mapping["old_id"]))
    )
    identity_mapping = (
        not removed_ids
        and len(id_mapping) == len(prior_ligands)
        and prior_queries == set(map(int, prior_ligands[number_id_col]))
        and id_mapping["old_id"].equals(id_mapping["new_id"])
    )

    if prior_paths and identity_mapping:
        for index, source in enumerate(prior_paths):
            link_or_copy_file(source, stage / f"retained-{index:05d}.parquet")
    elif prior_paths:
        scratch_dir.mkdir(parents=True, exist_ok=True)
        retained = scratch_dir / "retained.parquet"
        retained.unlink(missing_ok=True)
        sources = ", ".join(f"'{path.as_posix()}'" for path in prior_paths)
        retained_targets = id_mapping.loc[~id_mapping["new_id"].isin(missing_queries)]
        with duckdb.connect() as connection:
            connection.register("query_mapping", complete_mapping)
            connection.register("target_mapping", retained_targets)
            connection.execute(
                f"""
                COPY (
                    SELECT
                        query_mapping.new_id::INTEGER AS query_ligand_id,
                        target_mapping.new_id::INTEGER AS target_ligand_id,
                        scores.{metric}::FLOAT AS {metric}
                    FROM read_parquet([{sources}], union_by_name=true) AS scores
                    INNER JOIN query_mapping
                      ON scores.query_ligand_id = query_mapping.old_id
                    INNER JOIN target_mapping
                      ON scores.target_ligand_id = target_mapping.old_id
                ) TO '{retained.as_posix()}' (
                    FORMAT PARQUET, COMPRESSION ZSTD
                )
                """
            )
        link_or_copy_file(retained, stage / "retained.parquet")

    added_paths: list[Path] = []
    for index, start in enumerate(range(0, len(missing_queries), batch_size)):
        ligand_ids = sorted(missing_queries)[start : start + batch_size]
        output = stage / f"added-{index:05d}.parquet"
        if directory == "ligand_scores":
            get_similarity_scores.ligand_scores(
                ligand_ids=ligand_ids,
                data_dir=data_dir,
                output_path=output,
                number_id_col=number_id_col,
                minimum_similarity=minimum_similarity,
            )
        else:
            get_similarity_scores.mhfp6_ligand_scores(
                ligand_ids=ligand_ids,
                data_dir=data_dir,
                output_path=output,
                number_id_col=number_id_col,
                minimum_similarity=minimum_similarity,
            )
        added_paths.append(output)

    if added_paths and active_ids.difference(missing_queries):
        scratch_dir.mkdir(parents=True, exist_ok=True)
        reverse = scratch_dir / "reverse.parquet"
        reverse.unlink(missing_ok=True)
        sources = ", ".join(f"'{path.as_posix()}'" for path in added_paths)
        new = pd.DataFrame({"ligand_id": sorted(missing_queries)})
        with duckdb.connect() as connection:
            connection.register("new_ligands", new)
            connection.execute(
                f"""
                COPY (
                    SELECT
                        scores.target_ligand_id::INTEGER AS query_ligand_id,
                        scores.query_ligand_id::INTEGER AS target_ligand_id,
                        scores.{metric}::FLOAT AS {metric}
                    FROM read_parquet([{sources}]) AS scores
                    ANTI JOIN new_ligands
                      ON scores.target_ligand_id = new_ligands.ligand_id
                ) TO '{reverse.as_posix()}' (
                    FORMAT PARQUET, COMPRESSION ZSTD
                )
                """
            )
        link_or_copy_file(reverse, stage / "reverse.parquet")

    installed_paths = sorted(stage.glob("*.parquet"))
    get_similarity_scores._check_self_score_coverage(
        installed_paths,
        metric=metric,
        expected_query_ids=active_ids,
        label=metric,
    )
    write_json_atomic(
        stage / get_similarity_scores.LIGAND_SCORE_MANIFEST,
        current_manifest,
    )
    if root.exists():
        root.rename(backup)
    try:
        stage.rename(root)
    except BaseException:
        if backup.exists() and not root.exists():
            backup.rename(root)
        raise
    rmtree(backup, ignore_errors=True)
    return {
        "status": "complete",
        "active_ligands": len(active_ids),
        "new_query_scores": len(missing_queries),
        "removed_ligands": len(removed_ids),
        "shards": len(installed_paths),
    }


def refresh_ligand_chemistry(
    data_dir: Path,
    *,
    cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
) -> dict[str, Any]:
    """Update unique-ligand chemistry, pair similarities, and MMP tables."""
    fingerprint_path = data_dir / "fingerprints/ligands_per_smiles.parquet"
    prior_ligands = pd.read_parquet(
        fingerprint_path,
        columns=[
            "ligand_smiles_id",
            "ligand_rdkit_canonical_smiles",
            "fingerprint",
        ],
    )
    number_id_col = str(cfg.ligand.number_id_col)
    minimum_similarity = float(cfg.ligand.minimum_similarity)
    prior_tanimoto_manifest = get_similarity_scores.ligand_score_manifest_payload(
        fingerprint_path=fingerprint_path,
        fingerprint_col="fingerprint",
        metric="tanimoto_similarity_ecfp4_1024",
        minimum_similarity=minimum_similarity,
        number_id_col=number_id_col,
    )
    use_mhfp6 = get_similarity_scores.MHFP6_METRIC in cfg.flow.cluster_metrics
    if use_mhfp6:
        prior_mhfp6_ligands = pd.read_parquet(
            data_dir / "fingerprints" / get_similarity_scores.MHFP6_FINGERPRINT_FILE,
            columns=[
                "ligand_smiles_id",
                "ligand_rdkit_canonical_smiles",
                "mhfp6",
            ],
        )
        prior_mhfp6_manifest = get_similarity_scores.ligand_score_manifest_payload(
            fingerprint_path=(
                data_dir / "fingerprints" / get_similarity_scores.MHFP6_FINGERPRINT_FILE
            ),
            fingerprint_col="mhfp6",
            metric=get_similarity_scores.MHFP6_METRIC,
            minimum_similarity=minimum_similarity,
            number_id_col=number_id_col,
        )
    tasks.compute_ligand_fingerprints(
        data_dir=data_dir,
        cofactor_similarity_threshold=float(cfg.ligand.cofactor_similarity_threshold),
        retain_score_shards=True,
    )
    ligand_scores = _refresh_chemical_score_shards(
        data_dir=data_dir,
        directory="ligand_scores",
        metric="tanimoto_similarity_ecfp4_1024",
        batch_size=int(cfg.flow.make_ligands_batch_size),
        number_id_col=number_id_col,
        minimum_similarity=minimum_similarity,
        scratch_dir=scratch_dir / "tanimoto",
        prior_ligands=prior_ligands,
        prior_manifest=prior_tanimoto_manifest,
        current_manifest=get_similarity_scores.ligand_score_manifest_payload(
            fingerprint_path=fingerprint_path,
            fingerprint_col="fingerprint",
            metric="tanimoto_similarity_ecfp4_1024",
            minimum_similarity=minimum_similarity,
            number_id_col=number_id_col,
        ),
    )
    mhfp6_scores = None
    if use_mhfp6:
        mhfp6_scores = _refresh_chemical_score_shards(
            data_dir=data_dir,
            directory=get_similarity_scores.MHFP6_SCORES_DIR,
            metric=get_similarity_scores.MHFP6_METRIC,
            batch_size=int(cfg.flow.make_ligands_batch_size),
            number_id_col=number_id_col,
            minimum_similarity=minimum_similarity,
            scratch_dir=scratch_dir / "mhfp6",
            prior_ligands=prior_mhfp6_ligands,
            prior_manifest=prior_mhfp6_manifest,
            current_manifest=get_similarity_scores.ligand_score_manifest_payload(
                fingerprint_path=(
                    data_dir
                    / "fingerprints"
                    / get_similarity_scores.MHFP6_FINGERPRINT_FILE
                ),
                fingerprint_col="mhfp6",
                metric=get_similarity_scores.MHFP6_METRIC,
                minimum_similarity=minimum_similarity,
                number_id_col=number_id_col,
            ),
        )
    tasks.annotate_ligand_similarity(
        data_dir=data_dir,
        minimum_similarity=minimum_similarity,
        number_id_col=number_id_col,
    )
    mmp = tasks.make_ligand_mmp_pairs(
        data_dir=data_dir,
        scratch_dir=scratch_dir / "mmp",
        threads=threads,
    )
    ccd = None
    if "make_ccd_ligand_dbs" not in cfg.flow.skip_specific_stages:
        ccd = tasks.make_ccd_ligand_dbs(
            data_dir=data_dir,
            scratch_dir=scratch_dir / "ccd",
            threads=threads,
            minimum_similarity=float(cfg.ligand.minimum_similarity),
        )
    return {
        "status": "complete",
        "tanimoto": ligand_scores,
        "mhfp6": mhfp6_scores,
        "mmp": str(mmp),
        "ccd_matches": str(ccd) if ccd is not None else None,
    }


def refresh_ligand_score_export(
    data_dir: Path,
    *,
    scratch_dir: Path,
    threads: int,
    memory_limit: str,
) -> dict[str, Any]:
    """Update the compact public ligand-pair score table."""
    shards = sorted(
        {pdb_id[1:3] for pdb_id in score.published_scoring_query_ids(data_dir)}
    )
    shard_dir = data_dir / "exports" / "ligand_similarity_score_shards"
    active_shards = set(shards)
    for path in shard_dir.glob("*.parquet"):
        if path.stem not in active_shards:
            path.unlink()
            path.with_suffix(".json").unlink(missing_ok=True)
    score.export_ligand_similarity_scores_batch(
        data_dir,
        output_dir=shard_dir,
        shards=shards,
        scratch_dir=scratch_dir / "ligand-score-shards",
        threads=threads,
        memory_limit=memory_limit,
    )
    return score.finalize_ligand_similarity_scores(
        data_dir,
        source_dir=shard_dir,
        output=data_dir / RELEASE_PATHS["ligand_similarity_scores"],
        scratch_dir=scratch_dir / "ligand-score-export",
        threads=threads,
        memory_limit=memory_limit,
    )


def rebuild_similarity_covers(
    data_dir: Path,
    *,
    cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
) -> dict[str, Any]:
    """Recompute release-wide ligand and interface representative covers."""
    reports: dict[str, Any] = {}
    entities: list[tuple[clusters.ClusterEntity, list[str]]] = [
        ("ligand", list(cfg.flow.cluster_metrics)),
        ("interface", list(score.INTERFACE_CLUSTER_METRICS)),
    ]
    thresholds = list(map(int, cfg.flow.cluster_thresholds))
    for entity_type, metrics in entities:
        plan = score.plan_clustering(
            data_dir,
            metrics=metrics,
            thresholds=thresholds,
            source_batch_size=int(cfg.flow.symmetric_edge_source_batch_size),
            cover_batch_size=1,
            symmetric_bucket_count=int(cfg.flow.symmetric_edge_bucket_count),
            entity_type=entity_type,
        )
        symmetric = clusters.load_symmetric_edge_plan(data_dir, entity_type=entity_type)
        for batch in symmetric["batches"]:
            tasks.make_symmetric_edge_fragments(
                data_dir=data_dir,
                batches=[batch],
                scratch_dir=scratch_dir / entity_type / "edge-fragments",
                threads=threads,
                force_update=True,
                entity_type=entity_type,
            )
        for metric in metrics:
            for bucket in range(int(symmetric["bucket_count"])):
                tasks.make_symmetric_edge_shards(
                    data_dir=data_dir,
                    metric_buckets=[(metric, bucket)],
                    scratch_dir=scratch_dir / entity_type / "edge-shards",
                    threads=threads,
                    force_update=True,
                    entity_type=entity_type,
                )
        for sources in tasks.scatter_component_reduction_sources(
            data_dir=data_dir,
            metrics=metrics,
            batch_size=int(cfg.flow.component_reduction_source_batch_size),
            entity_type=entity_type,
        ):
            tasks.make_component_reductions(
                data_dir=data_dir,
                source_paths=sources,
                metrics=metrics,
                thresholds=thresholds,
                scratch_dir=scratch_dir / entity_type / "components",
                force_update=True,
                metric_workers=min(
                    threads, int(cfg.flow.component_reduction_metric_workers)
                ),
                entity_type=entity_type,
            )
        tasks.merge_component_reductions(
            data_dir=data_dir,
            metrics=metrics,
            thresholds=thresholds,
            entity_type=entity_type,
        )
        for work in tasks.scatter_make_set_covers(
            data_dir=data_dir,
            metrics=metrics,
            thresholds=thresholds,
            stop_on_cluster=0,
            skip_existing_clusters=False,
            entity_type=entity_type,
        ):
            tasks.make_set_covers(
                data_dir=data_dir,
                metric_threshold=work,
                skip_existing_clusters=False,
                scratch_dir=scratch_dir / entity_type / "set-covers",
                threads=threads,
                entity_type=entity_type,
            )
        for work in tasks.scatter_make_directed_set_covers(
            data_dir=data_dir,
            metrics=metrics,
            thresholds=thresholds,
            stop_on_cluster=0,
            skip_existing=False,
            entity_type=entity_type,
        ):
            tasks.make_directed_set_covers(
                data_dir=data_dir,
                metric_threshold=work,
                skip_existing=False,
                scratch_dir=scratch_dir / entity_type / "directed-covers",
                threads=threads,
                entity_type=entity_type,
            )
        reports[entity_type] = {
            **plan,
            "summary": score.summarize_clustering_artifacts(
                data_dir,
                metrics=metrics,
                thresholds=thresholds,
                entity_type=entity_type,
            ),
        }
    return {"status": "complete", **reports}


def _load_configuration(path: Path) -> DictConfig:
    return config.get_config(config=OmegaConf.load(path), config_args=[], cached=False)


def _weekly_inputs(
    plan_dir: Path, config_path: Path, validation_root: Path
) -> dict[str, Any]:
    return {
        "plan": file_sha256(plan_dir / "plan.json"),
        "entries": file_sha256(plan_dir / "entries.parquet"),
        "config": file_sha256(config_path),
        "validation_root": str(validation_root.resolve(strict=True)),
    }


def _read_weekly_report(path: Path, *, inputs: dict[str, Any]) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    report: dict[str, Any] = json.loads(path.read_text())
    if report.get("inputs") != inputs:
        raise ValueError("weekly update inputs changed; use a new workspace")
    stages = report.get("completed_stages")
    if not isinstance(stages, list) or stages != list(_STAGES[: len(stages)]):
        raise ValueError(f"invalid weekly update stage report: {path}")
    return report


def _finish_stage(
    path: Path,
    report: dict[str, Any],
    stage: str,
    details: dict[str, Any],
) -> None:
    stages = cast(list[str], report["completed_stages"])
    expected = _STAGES[len(stages)]
    if stage != expected:
        raise RuntimeError(f"expected weekly stage {expected}, got {stage}")
    report[stage] = details
    stages.append(stage)
    report["status"] = "complete" if stage == _STAGES[-1] else "running"
    write_json_atomic(path, report)


def apply_release_update(
    plan_dir: Path,
    workspace: Path,
    *,
    validation_root: Path,
    config_path: Path,
    scratch_dir: Path,
    threads: int = 8,
    memory_limit: str = "32GB",
    search_batch_size: int = 5_000,
) -> dict[str, Any]:
    """Apply one planned PDB update and finish every derived release table."""
    if min(threads, search_batch_size) < 1:
        raise ValueError("threads and search batch size must be positive")
    plan_dir = plan_dir.resolve(strict=True)
    config_path = config_path.resolve(strict=True)
    validation_root = validation_root.resolve(strict=True)
    scratch_dir = scratch_dir.resolve()
    scratch_dir.mkdir(parents=True, exist_ok=True)
    inputs = _weekly_inputs(plan_dir, config_path, validation_root)
    state_path = workspace.resolve() / "weekly_update.json"
    report = _read_weekly_report(state_path, inputs=inputs)
    if report is not None and report["status"] in {"complete", "finalizing"}:
        marker_path = workspace / "index" / "collation.json"
        marker = json.loads(marker_path.read_text())
        if marker.get("status") != "complete":
            if report["status"] == "complete":
                raise ValueError("completed weekly update has an incomplete index")
        else:
            if report["status"] == "finalizing":
                _finish_stage(
                    state_path,
                    report,
                    "final_tables",
                    {"collation": str(marker_path)},
                )
            return report

    cfg = _load_configuration(config_path)
    search_databases = list(dict.fromkeys(map(str, cfg.scorer.sub_databases)))
    pdb_search_databases = [
        search_db for search_db in search_databases if search_db in {"holo", "apo"}
    ]
    plan, entries = load_update_plan(plan_dir)
    base = Path(plan["data_dir"]).resolve(strict=True)
    affected = set(
        entries.loc[
            entries.action.isin(["added", "revised", "obsolete"]), "pdb_id"
        ].astype(str)
    )
    if report is None:
        entry_report = apply_entry_update(
            plan_dir,
            workspace,
            validation_root=validation_root,
            annotation_cfg=dict(cfg.annotation),
            entry_cfg=dict(cfg.entry),
            interface_cfg=dict(cfg.interface),
            threads=threads,
            memory_limit=memory_limit,
            scratch_dir=scratch_dir / "entries",
        )
        archive_report = update_ligand_archives(
            workspace, threads=threads, memory_limit=memory_limit
        )
        workspace = workspace.resolve(strict=True)
        _reuse_base_artifacts(base, workspace)
        report = {
            "status": "running",
            "inputs": inputs,
            "base_release": str(base),
            "affected_pdb_ids": sorted(affected),
            "completed_stages": [],
        }
        _finish_stage(
            state_path,
            report,
            "entries_and_archives",
            {"entries": entry_report, "ligand_archives": archive_report},
        )
    else:
        workspace = workspace.resolve(strict=True)

    completed = set(report["completed_stages"])
    if "search_inputs" not in completed:
        lookup = tasks.make_alignment_chain_lookup(
            data_dir=workspace,
            scratch_dir=scratch_dir / "scoring-inputs",
            threads=threads,
            memory_limit=memory_limit,
        )
        protein_plan = score.plan_protein_scoring(
            workspace, max_seqs=int(cfg.foldseek.max_seqs)
        )
        sequence_clusters = protein_clusters.make_protein_sequence_clusters(
            data_dir=workspace,
            scratch_dir=scratch_dir / "protein-sequence-clusters",
            threads=threads,
            identity=float(cfg.flow.protein_sequence_cluster_identity),
            coverage=float(cfg.flow.protein_cluster_coverage),
        )
        structure_clusters = protein_clusters.make_protein_structure_clusters(
            data_dir=workspace,
            cif_root=Path(plan["nextgen_root"]) / "data/entries/divided",
            scratch_dir=scratch_dir / "protein-structure-clusters",
            threads=threads,
            lddt=float(cfg.flow.protein_structure_cluster_lddt),
            coverage=float(cfg.flow.protein_cluster_coverage),
        )
        if pdb_search_databases:
            foldseek_input = score.make_foldseek_input_manifest(
                workspace, Path(plan["nextgen_root"]) / "data/entries/divided"
            )
            mmseqs_input = score.make_mmseqs_input_fasta(workspace)
            tasks.make_dbs(
                data_dir=workspace,
                sub_databases=pdb_search_databases,
                cpu=threads,
                cif_root=foldseek_input,
                seqres_path=mmseqs_input,
                scratch_dir=scratch_dir / "databases",
                build_dir=scratch_dir / "database-build",
                index=False,
            )
            tasks.make_sub_dbs(
                data_dir=workspace,
                sub_databases=pdb_search_databases,
                cpu=threads,
                scratch_dir=scratch_dir / "exact-search-dbs",
            )
        search_details: dict[str, Any] = {
            "alignment_chain_lookup": str(lookup),
            "protein_plan": protein_plan,
            "protein_sequence_clusters": str(sequence_clusters),
            "protein_structure_clusters": str(structure_clusters),
            "search_databases": search_databases,
        }
        if "apo" in search_databases:
            search_details["apo_plan"] = score.plan_linked_apo_scoring(
                workspace, max_seqs=int(cfg.foldseek.max_seqs)
            )
        _finish_stage(
            state_path,
            report,
            "search_inputs",
            search_details,
        )

    completed = set(report["completed_stages"])
    if "alignments" not in completed:
        planned_queries = _load_or_plan_alignment_repairs(
            workspace,
            search_databases=search_databases,
            affected=affected,
            scorer_cfg=cfg.scorer,
            foldseek_cfg=cfg.foldseek,
            mmseqs_cfg=cfg.mmseqs,
            scratch_dir=scratch_dir / "alignment-planning",
            threads=threads,
            batch_size=search_batch_size,
            state_path=state_path,
            report=report,
        )
        full_queries = {
            search_db: sorted(
                repair_alignments(
                    workspace,
                    search_db=search_db,
                    affected=affected,
                    full_queries=planned_queries[search_db],
                    scorer_cfg=cfg.scorer,
                    foldseek_cfg=cfg.foldseek,
                    mmseqs_cfg=cfg.mmseqs,
                    scratch_dir=scratch_dir / "alignment-repair" / search_db,
                    threads=threads,
                    batch_size=search_batch_size,
                )
            )
            for search_db in search_databases
        }
        alignment_details: dict[str, Any] = {"full_queries": full_queries}
        if "holo" in search_databases:
            alignment_details["holo"] = score.finalize_alignment_artifacts(workspace)
        _finish_stage(
            state_path,
            report,
            "alignments",
            alignment_details,
        )
    full_query_sets = {
        name: set(map(str, values))
        for name, values in report["alignments"]["full_queries"].items()
    }

    completed = set(report["completed_stages"])
    if "scores" not in completed:
        score_details: dict[str, Any] = {}
        if "holo" in search_databases:
            score_details["ligand"] = repair_holo_scores(
                workspace,
                base_data_dir=base,
                affected=affected,
                full_alignment_queries=full_query_sets["holo"],
                scorer_cfg=cfg.scorer,
                scratch_dir=scratch_dir / "holo-score-repair",
                threads=threads,
                memory_limit=memory_limit,
                score_batch_size=int(cfg.flow.make_batch_scores_batch_size),
                ligand_batch_size=int(cfg.flow.make_ligand_3d_scores_batch_size),
            )
            score_details["interface"] = repair_interface_scores(
                workspace,
                affected=affected,
                full_alignment_queries=full_query_sets["holo"],
                scratch_dir=scratch_dir / "interface-score-repair",
                threads=threads,
                memory_limit=memory_limit,
            )
            score_details["ligand_export"] = refresh_ligand_score_export(
                workspace,
                scratch_dir=scratch_dir / "ligand-score-export",
                threads=threads,
                memory_limit=memory_limit,
            )
        if "apo" in search_databases:
            score_details["apo"] = repair_apo_scores(
                workspace,
                affected=affected,
                full_alignment_queries=full_query_sets["apo"],
                scorer_cfg=cfg.scorer,
                scratch_dir=scratch_dir / "apo-score-repair",
                threads=threads,
                memory_limit=memory_limit,
            )
        if "pred" in search_databases:
            score_details["pred"] = repair_non_holo_scores(
                workspace,
                search_db="pred",
                affected=affected,
                full_alignment_queries=full_query_sets["pred"],
                scorer_cfg=cfg.scorer,
                scratch_dir=scratch_dir / "pred-score-repair",
                threads=threads,
                memory_limit=memory_limit,
            )
        _finish_stage(
            state_path,
            report,
            "scores",
            score_details,
        )

    completed = set(report["completed_stages"])
    if "ligand_chemistry" not in completed:
        chemistry = refresh_ligand_chemistry(
            workspace,
            cfg=cfg,
            scratch_dir=scratch_dir / "ligand-chemistry",
            threads=threads,
        )
        _finish_stage(state_path, report, "ligand_chemistry", chemistry)

    completed = set(report["completed_stages"])
    if "similarity_covers" not in completed:
        covers = rebuild_similarity_covers(
            workspace,
            cfg=cfg,
            scratch_dir=scratch_dir / "similarity-covers",
            threads=threads,
        )
        _finish_stage(state_path, report, "similarity_covers", covers)

    if "final_tables" not in set(report["completed_stages"]):
        report["status"] = "finalizing"
        write_json_atomic(state_path, report)
        tasks.finalize_index(data_dir=workspace)
        marker_path = workspace / "index" / "collation.json"
        marker = json.loads(marker_path.read_text())
        if marker.get("status") != "complete":
            raise ValueError("final table assembly did not complete the update marker")
        _finish_stage(
            state_path,
            report,
            "final_tables",
            {"collation": str(marker_path)},
        )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan_dir", type=Path)
    parser.add_argument("workspace", type=Path)
    parser.add_argument("--validation-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-limit", default="32GB")
    parser.add_argument("--search-batch-size", type=int, default=5_000)
    args = parser.parse_args()
    report = apply_release_update(
        args.plan_dir,
        args.workspace,
        validation_root=args.validation_root,
        config_path=args.config,
        scratch_dir=args.scratch_dir,
        threads=args.threads,
        memory_limit=args.memory_limit,
        search_batch_size=args.search_batch_size,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
