# Copyright (c) 2026, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Finish a planned PDB update in a new release directory."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import replace
from multiprocessing import get_context
from pathlib import Path
from shutil import copyfile, rmtree
from typing import Any, Iterable, cast

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
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
from plinder.data import clusters, databases
from plinder.data.annotations import get_similarity_scores
from plinder.data.pipeline import config, score, tasks, utils
from plinder.data.pipeline.update_archives import update_ligand_archives
from plinder.data.pipeline.update_entries import apply_entry_update
from plinder.data.pipeline.updates import load_update_plan
from plinder.data.pipeline.weekly_clusters import (
    extend_interface_clusters,
    extend_ligand_clusters,
    extend_protein_clusters,
)

LOG = logging.getLogger(__name__)

_BASE_ARTIFACT_ROOTS = (
    "alignments",
    "alignment_cigars",
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
_TRANSIENT_BASE_DIRECTORIES = {".staging", "ligand_3d_pair_repairs", "ligand_3d_pairs"}
_TRANSIENT_BASE_SUFFIXES = (".installing", ".previous", ".tmp", ".tmp.parquet")

_STAGES = (
    "entries_and_archives",
    "search_inputs",
    "alignments",
    "scores",
    "ligand_chemistry",
    "cluster_assignments",
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
                if (
                    name in _TRANSIENT_BASE_DIRECTORIES
                    or (root_name == "dbs" and name == "aln")
                    or (
                        relative.as_posix() == "dbs/subdbs"
                        and name in {"search_db=holo", "search_db=apo"}
                    )
                    or name.endswith(_TRANSIENT_BASE_SUFFIXES)
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
    if search_db == "interface_apo":
        return tasks._interface_apo_scoring_chains(data_dir)
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


def _query_manifest(data_dir: Path, search_db: str) -> Path:
    if search_db == "apo":
        return data_dir / score.LINKED_APO_QUERY_MANIFEST_RELATIVE
    if search_db == "interface_apo":
        return data_dir / score.INTERFACE_APO_QUERY_MANIFEST_RELATIVE
    return data_dir / score.MANIFEST_RELATIVE


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


def _prepare_weekly_search_overlay(
    *,
    base: Path,
    data_dir: Path,
    affected: set[str],
    nextgen_root: Path,
    search_databases: list[str],
    scratch_dir: Path,
    threads: int,
) -> set[str]:
    """Keep the full search DB fixed and index cumulative changed chains."""
    shadowed_path = Path("manifests/weekly_shadowed_entries.parquet")
    previous = base / shadowed_path
    shadowed = set(affected)
    if previous.is_file():
        shadowed.update(
            map(str, pq.read_table(previous, columns=["pdb_id"])["pdb_id"].to_pylist())
        )
    chains = pd.read_parquet(
        data_dir / "index/entry_chains.parquet",
        columns=[
            "entry_pdb_id",
            "chain_auth_id",
            "chain_receptor_type",
            "chain_sequence",
        ],
        filters=[("entry_pdb_id", "in", sorted(shadowed))],
    )
    proteins = chains.loc[
        chains.chain_receptor_type.eq("protein") & chains.chain_auth_id.notna()
    ].copy()
    active = set(proteins.entry_pdb_id.astype(str))
    overlay = data_dir / "dbs/weekly_delta"
    if active:
        staging = scratch_dir / "weekly_delta_build"
        rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True)
        foldseek_input = staging / "foldseek_inputs.tsv"
        with foldseek_input.open("w") as handle:
            for pdb_id in sorted(active):
                source = (
                    nextgen_root
                    / "data/entries/divided"
                    / pdb_id[1:3]
                    / f"pdb_0000{pdb_id}"
                    / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
                )
                if not source.is_file():
                    raise FileNotFoundError(source)
                handle.write(f"{source}\n")
        proteins["identifier"] = (
            proteins.entry_pdb_id.astype(str) + "_" + proteins.chain_auth_id.astype(str)
        )
        conflicts = proteins.groupby("identifier").chain_sequence.nunique()
        if conflicts.gt(1).any():
            raise ValueError("changed chains have conflicting auth-chain sequences")
        sequences = proteins.dropna(subset=["chain_sequence"]).drop_duplicates(
            "identifier"
        )
        sequences = sequences.loc[sequences.chain_sequence.ne("")]
        if not sequences.chain_sequence.astype(str).str.fullmatch("[A-Za-z]+").all():
            raise ValueError("changed protein sequences contain invalid residues")
        mmseqs_input = staging / "mmseqs_inputs.fasta"
        with mmseqs_input.open("w") as handle:
            for row in sequences.itertuples(index=False):
                handle.write(f">{row.identifier}\n{row.chain_sequence}\n")
        sources_by_backend = {"foldseek": foldseek_input}
        if not sequences.empty:
            sources_by_backend["mmseqs"] = mmseqs_input
        for backend, source in sources_by_backend.items():
            databases.create_db(source, staging / backend, backend, threads=threads)

        selected: dict[str, set[str]] = {}
        sources: dict[str, Path] = {}
        for search_db in search_databases:
            if search_db == "pred":
                continue
            target_chains = _scoring_chains(data_dir, search_db)
            target_chains = target_chains.loc[
                target_chains.entry_pdb_id.astype(str).isin(active)
            ]
            for backend in ("foldseek", "mmseqs"):
                if backend not in sources_by_backend:
                    continue
                identifiers = _database_identifiers(target_chains, backend)
                if identifiers:
                    key = f"{search_db}_{backend}"
                    selected[key] = identifiers
                    sources[key] = staging / backend / backend
        if sources:
            databases.make_sub_dbs(
                staging / "subdbs",
                sources,
                identifiers_by_database=selected,
                tmp_dir=scratch_dir / "exact-search-dbs",
                threads=threads,
            )
        databases.install_database_directory(staging, overlay)
    elif overlay.exists():
        rmtree(overlay)

    target = data_dir / shadowed_path
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".parquet.tmp")
    pq.write_table(
        pa.table({"pdb_id": pa.array(sorted(shadowed), type=pa.string())}),
        temporary,
    )
    temporary.replace(target)
    return shadowed


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
                data_dir
                / "alignment_cigars"
                / f"search_db={search_db}"
                / f"alignment_type={alignment_type}"
            ).glob("shard=*.parquet")
        )
    ]
    if not paths:
        return set()
    connection = duckdb.connect()
    try:
        connection.register(
            "affected_targets",
            pd.DataFrame({"target_entry": sorted(affected)}),
        )
        sources = ", ".join(
            f"'{path.as_posix().replace(chr(39), chr(39) * 2)}'" for path in paths
        )
        hits = connection.sql(
            f"""
            SELECT DISTINCT query_entry
            FROM read_parquet([{sources}], union_by_name=true) AS alignments
            INNER JOIN affected_targets USING (target_entry)
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
    selected = {
        key: ids
        for key, ids in identifiers.items()
        if ids
        and (
            data_dir
            / "dbs/weekly_delta"
            / key.rsplit("_", 1)[1]
            / f"{key.rsplit('_', 1)[1]}.dbtype"
        ).is_file()
    }
    if not selected:
        return None
    databases.make_sub_dbs(
        target_dir,
        {
            key: data_dir
            / "dbs/weekly_delta"
            / key.rsplit("_", 1)[1]
            / key.rsplit("_", 1)[1]
            for key in selected
        },
        identifiers_by_database=selected,
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
) -> set[str]:
    if target_database_dir is None or not query_ids:
        return set()
    scorer, _, work = utils.get_scorer(
        data_dir=data_dir,
        pdb_ids=[],
        scorer_cfg=scorer_cfg,
        load_entries=False,
        foldseek_cfg=foldseek_cfg,
        mmseqs_cfg=mmseqs_cfg,
        scratch_dir=scratch_dir / "queries",
    )
    hits: set[str] = set()
    active = set(query_ids)
    try:
        for aln_type in ("mmseqs", "foldseek"):
            if not (
                target_database_dir / f"{search_db}_{aln_type}/exact_cluster.json"
            ).is_file():
                continue
            changed_db, _, _ = databases.exact_search_database_paths(
                target_database_dir, search_db, aln_type
            )
            probe_dir = work / search_db / aln_type
            probe_dir.mkdir(parents=True, exist_ok=True)
            # Candidate discovery is permissive; the usual forward search
            # applies the requested cutoff and produces the release scores.
            probe_config = replace(
                scorer.get_config(search_db, aln_type),
                evalue=2.0,
                max_seqs=2_000_000,
                coverage=0.0,
                min_seq_id=0.0,
            )
            full_databases = [scorer.source_to_full_db_file[f"holo_{aln_type}"]]
            overlay = data_dir / "dbs/weekly_delta" / aln_type / aln_type
            if overlay.with_suffix(".dbtype").is_file():
                full_databases.append(overlay)
            for index, full_db in enumerate(full_databases):
                target_identifiers = get_similarity_scores.run_alignment(
                    aln_type=aln_type,
                    query_db=changed_db,
                    target_db=full_db,
                    search_target_db=full_db,
                    search_db=probe_dir / f"search-{index}",
                    aln_file=probe_dir / f"hits-{index}.tsv",
                    alignment_config=probe_config,
                    tmp_dir=probe_dir / f"tmp-{index}",
                    threads=threads,
                    query_ids_only=True,
                    id_column="target",
                )
                hits.update(
                    identifier.removeprefix("pdb_0000")[:4]
                    for identifier in target_identifiers or ()
                    if identifier.removeprefix("pdb_0000")[:4] in active
                )
    finally:
        rmtree(work, ignore_errors=True)
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
    replacement_query_ids: set[str],
    replacement_target_ids: set[str] | None = None,
    scorer_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
) -> None:
    def refresh(shard: str) -> None:
        tasks.map_batch_alignments(
            data_dir=data_dir,
            shards=[shard],
            scorer_cfg=scorer_cfg,
            force_update=True,
            scratch_dir=scratch_dir,
            search_db=search_db,
            replacement_query_ids={
                pdb_id for pdb_id in replacement_query_ids if pdb_id[1:3] == shard
            },
            replacement_target_ids=replacement_target_ids,
        )
        if any(
            tasks._alignment_release_path(
                data_dir=data_dir,
                search_db=search_db,
                alignment_type=alignment_type,
                shard=shard,
            ).is_file()
            for alignment_type in ("foldseek", "mmseqs")
        ):
            return
        manifest_root = data_dir / "alignments" / "manifests"
        if search_db != "holo":
            manifest_root = manifest_root / f"search_db={search_db}"
        (manifest_root / f"shard={shard}.json").unlink(missing_ok=True)

    ordered = sorted(shards)
    workers = min(4, threads, len(ordered))
    if workers == 1:
        for shard in ordered:
            refresh(shard)
    elif workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(refresh, ordered))


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
        inputs = payload.get("inputs")
        if not isinstance(inputs, dict):
            raise ValueError(f"invalid alignment report inputs: {path}")
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
            if not output.is_file() or not (
                schemas.release_alignment_mapping_schema_is_current(
                    set(pq.read_schema(output).names),
                    alignment_type=alignment_type,
                    include_scores=search_db != "holo",
                )
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
) -> set[str]:
    """Find queries whose capped search can change in either direction."""
    manifest_path = _query_manifest(data_dir, search_db)
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
    targeted_queries: set[str] | None = None,
    scorer_cfg: DictConfig,
    foldseek_cfg: DictConfig,
    mmseqs_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
    batch_size: int,
) -> set[str]:
    """Fully search changed queries and add changed-target hits to old queries."""
    manifest_path = _query_manifest(data_dir, search_db)
    active = set(
        pd.read_parquet(manifest_path, columns=["pdb_id"])["pdb_id"].astype(str)
    )
    full_queries = active.intersection(full_queries)
    targeted_queries = active.intersection(targeted_queries or set()) - full_queries
    inactive = affected.difference(active)
    _remove_inactive_queries(data_dir, search_db=search_db, pdb_ids=inactive)
    configured_search_db = "apo" if search_db == "interface_apo" else search_db
    selected_cfg = cast(
        DictConfig,
        OmegaConf.merge(scorer_cfg, {"sub_databases": [configured_search_db]}),
    )
    for alignment_type in ("foldseek", "mmseqs"):
        raw_dir = data_dir / "dbs/subdbs" / f"{search_db}_{alignment_type}" / "aln"
        for pdb_id in full_queries | targeted_queries:
            (raw_dir / f"{pdb_id}.parquet").unlink(missing_ok=True)
    if targeted_queries:
        changed_targets = _changed_target_database(
            data_dir,
            search_db=search_db,
            affected=affected,
            scratch_dir=scratch_dir / "changed-targets",
            threads=threads,
        )
        if changed_targets is None:
            # Obsolete targets have no database to search. Re-search their old
            # query partners to remove stale hits.
            full_queries.update(targeted_queries)
            targeted_queries.clear()
        else:
            # Arrow limits the partitioned raw output to 1,024 query PDBs
            # per write batch.
            targeted_batch_size = min(batch_size, 500)
            available_backends = [
                backend
                for backend in ("foldseek", "mmseqs")
                if (
                    changed_targets / f"{search_db}_{backend}/exact_cluster.json"
                ).is_file()
            ]
            shadowed_path = data_dir / "manifests/weekly_shadowed_entries.parquet"
            shadowed = set(
                map(
                    str,
                    pq.read_table(shadowed_path, columns=["pdb_id"])[
                        "pdb_id"
                    ].to_pylist(),
                )
            )
            for source_label, query_ids in (
                ("base", targeted_queries - shadowed),
                ("overlay", targeted_queries & shadowed),
            ):
                ordered_targets = sorted(query_ids)
                for start in range(0, len(ordered_targets), targeted_batch_size):
                    tasks.run_batch_searches(
                        data_dir=data_dir,
                        pdb_ids=ordered_targets[start : start + targeted_batch_size],
                        scorer_cfg=selected_cfg,
                        foldseek_cfg=foldseek_cfg,
                        mmseqs_cfg=mmseqs_cfg,
                        cpu=threads,
                        scratch_dir=(
                            scratch_dir / f"targeted-{search_db}-{source_label}-"
                            f"{start // targeted_batch_size:05d}"
                        ),
                        search_databases=[search_db],
                        alignment_types=available_backends,
                        force_update=True,
                        target_database_dir=changed_targets,
                        query_database_dir=(
                            data_dir / "dbs/weekly_delta"
                            if source_label == "overlay"
                            else None
                        ),
                    )
            targeted_shards = {pdb_id[1:3] for pdb_id in targeted_queries}
            _refresh_alignment_shards(
                data_dir,
                search_db=search_db,
                shards=targeted_shards,
                replacement_query_ids=targeted_queries,
                replacement_target_ids=affected,
                scorer_cfg=selected_cfg,
                scratch_dir=scratch_dir / f"map-targeted-{search_db}",
                threads=threads,
            )
            for alignment_type in ("foldseek", "mmseqs"):
                raw_dir = (
                    data_dir / "dbs/subdbs" / f"{search_db}_{alignment_type}" / "aln"
                )
                for pdb_id in targeted_queries:
                    (raw_dir / f"{pdb_id}.parquet").unlink(missing_ok=True)
    shadowed_path = data_dir / "manifests/weekly_shadowed_entries.parquet"
    shadowed = set(
        map(str, pq.read_table(shadowed_path, columns=["pdb_id"])["pdb_id"].to_pylist())
    )
    for source_label, query_ids in (
        ("base", full_queries - shadowed),
        ("overlay", full_queries & shadowed),
    ):
        ordered = sorted(query_ids)
        for start in range(0, len(ordered), batch_size):
            tasks.run_batch_searches(
                data_dir=data_dir,
                pdb_ids=ordered[start : start + batch_size],
                scorer_cfg=selected_cfg,
                foldseek_cfg=foldseek_cfg,
                mmseqs_cfg=mmseqs_cfg,
                cpu=threads,
                scratch_dir=(
                    scratch_dir
                    / f"full-{search_db}-{source_label}-{start // batch_size:05d}"
                ),
                search_databases=[search_db],
                force_update=True,
                query_database_dir=(
                    data_dir / "dbs/weekly_delta" if source_label == "overlay" else None
                ),
            )
    missing_raw = [
        f"{search_db}_{alignment_type}/{pdb_id}"
        for alignment_type in ("foldseek", "mmseqs")
        for pdb_id in sorted(full_queries)
        if not (
            data_dir
            / "dbs"
            / "subdbs"
            / f"{search_db}_{alignment_type}"
            / "aln"
            / f"{pdb_id}.parquet"
        ).is_file()
    ]
    if missing_raw:
        raise FileNotFoundError(f"missing weekly alignment results: {missing_raw[:10]}")
    full_shards = {pdb_id[1:3] for pdb_id in full_queries | inactive}
    if full_shards:
        _refresh_alignment_shards(
            data_dir,
            search_db=search_db,
            shards=full_shards,
            replacement_query_ids=full_queries | inactive,
            scorer_cfg=selected_cfg,
            scratch_dir=scratch_dir / f"map-{search_db}",
            threads=threads,
        )
    _rebase_unchanged_alignment_manifests(
        data_dir,
        search_db=search_db,
        repaired_shards=full_shards | {pdb_id[1:3] for pdb_id in targeted_queries},
    )
    return full_queries


def _write_pdb_manifest(path: Path, pdb_ids: Iterable[str]) -> Path:
    path.parent.mkdir(exist_ok=True, parents=True)
    identifiers = sorted(set(map(str, pdb_ids)))
    if path.is_file() and pd.read_parquet(path)["pdb_id"].tolist() == identifiers:
        return path
    frame = pd.DataFrame({"pdb_id": identifiers})
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)
    return path


def _refresh_foldseek_source_manifest(
    data_dir: Path, *, nextgen_root: Path, affected: set[str]
) -> None:
    """Point reused Foldseek source paths at the current managed snapshot."""
    manifest = data_dir / "manifests/foldseek_createdb_inputs.tsv"
    if not manifest.is_file() or not affected:
        return
    lines = manifest.read_text().splitlines()
    updated: list[str] = []
    for line in lines:
        old = Path(line)
        pdb_id = old.parent.name.removeprefix("pdb_0000")
        if pdb_id in affected:
            current = (
                nextgen_root
                / "data/entries/divided"
                / pdb_id[1:3]
                / f"pdb_0000{pdb_id}"
                / old.name
            )
            if current.is_file():
                line = str(current)
        updated.append(line)
    if updated != lines:
        temporary = manifest.with_suffix(".tmp")
        temporary.write_text("\n".join(updated) + "\n")
        temporary.replace(manifest)


def _changed_ligand_archive_entries(
    base_data_dir: Path, data_dir: Path, affected: set[str]
) -> set[str]:
    """Identify entries whose ligand poses or features actually changed."""
    by_shard: dict[str, list[str]] = {}
    for pdb_id in sorted(affected):
        by_shard.setdefault(pdb_id[1:3], []).append(pdb_id)

    changed: set[str] = set()
    for shard, pdb_ids in by_shard.items():

        def archive_rows(root: Path) -> dict[str, list[dict[str, Any]]]:
            path = root / "ligand_archives" / f"{shard}.parquet"
            if not path.is_file():
                return {}
            rows = pq.read_table(
                path,
                columns=["pdb_id", "ligand_asym_id", "sdf", "pharmacophore_features"],
                filters=[("pdb_id", "in", pdb_ids)],
            ).to_pylist()
            grouped: dict[str, list[dict[str, Any]]] = {}
            for row in rows:
                grouped.setdefault(str(row["pdb_id"]), []).append(row)
            return grouped

        before = archive_rows(base_data_dir)
        after = archive_rows(data_dir)
        for pdb_id in pdb_ids:
            ordered_before = sorted(
                before.get(pdb_id, []), key=lambda row: row["ligand_asym_id"]
            )
            ordered_after = sorted(
                after.get(pdb_id, []), key=lambda row: row["ligand_asym_id"]
            )
            if ordered_before != ordered_after:
                changed.add(pdb_id)
    return changed


def _unchanged_scoring_entries(
    base: Path, data_dir: Path, affected: set[str]
) -> set[str]:
    """Find revised entries whose coordinate and scoring inputs are unchanged."""
    if not affected:
        return set()

    def source_paths(root: Path) -> dict[str, Path]:
        manifest = root / "manifests/foldseek_createdb_inputs.tsv"
        if not manifest.is_file():
            return {}
        paths: dict[str, Path] = {}
        with manifest.open() as lines:
            for line in lines:
                path = Path(line.strip())
                pdb_id = path.parent.name.removeprefix("pdb_0000")
                if pdb_id in affected:
                    paths[pdb_id] = path
        return paths

    old_sources = source_paths(base)
    new_sources = source_paths(data_dir)

    def source_digest(path: Path) -> str:
        digest = hashlib.sha256()
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    unchanged = {
        pdb_id
        for pdb_id in affected & old_sources.keys() & new_sources.keys()
        if old_sources[pdb_id].is_file()
        and new_sources[pdb_id].is_file()
        and (
            old_sources[pdb_id] == new_sources[pdb_id]
            or source_digest(old_sources[pdb_id]) == source_digest(new_sources[pdb_id])
        )
    }
    if not unchanged:
        return set()

    columns_by_table = {
        "entry_chains": [
            "chain_asym_id",
            "chain_auth_id",
            "chain_type",
            "chain_receptor_type",
            "chain_sequence",
            "chain_sequence_noncanonical",
            "chain_is_holo",
            "chain_is_ligand_like",
        ],
        "entry_biounit_chains": [
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
            "chain_num_contacting_ions",
            "chain_num_contacting_artifacts",
            "chain_num_contacting_other_ligands",
            "chain_num_contacting_proteins",
        ],
        "annotation_table": [
            "ligand_id",
            "system_id",
            "ligand_smiles",
            "ligand_is_proper",
            "ligand_neighboring_residues",
            "ligand_interactions",
        ],
        "interface_annotation_table": [
            "system_id",
            "interface_chain_1_residue_numbers",
            "interface_chain_2_residue_numbers",
        ],
    }
    for table_name, columns in columns_by_table.items():

        def rows(root: Path) -> dict[str, list[dict[str, Any]]] | None:
            path = root / "index" / f"{table_name}.parquet"
            if not path.is_file():
                return None
            selected = ["entry_pdb_id", *columns]
            if not set(selected).issubset(pq.read_schema(path).names):
                return None
            records = pq.read_table(
                path,
                columns=selected,
                filters=[("entry_pdb_id", "in", sorted(unchanged))],
            ).to_pylist()
            grouped: dict[str, list[dict[str, Any]]] = {}
            for record in records:
                grouped.setdefault(str(record.pop("entry_pdb_id")), []).append(record)
            return grouped

        before, after = rows(base), rows(data_dir)
        if before is None or after is None:
            return set()
        unchanged = {
            pdb_id
            for pdb_id in unchanged
            if sorted(before.get(pdb_id, []), key=str)
            == sorted(after.get(pdb_id, []), key=str)
        }
        if not unchanged:
            return set()
    return unchanged - _changed_ligand_archive_entries(base, data_dir, unchanged)


def _entries_with_interfaces(
    base: Path, data_dir: Path, affected: set[str]
) -> set[str]:
    """Return changed entries with interfaces in either release."""
    if not affected:
        return set()
    present: set[str] = set()
    for root in (base, data_dir):
        path = root / "index/interface_annotation_table.parquet"
        if path.is_file():
            present.update(
                map(
                    str,
                    pq.read_table(
                        path,
                        columns=["entry_pdb_id"],
                        filters=[("entry_pdb_id", "in", sorted(affected))],
                    )
                    .column("entry_pdb_id")
                    .to_pylist(),
                )
            )
    return present


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
    ordered = sorted(set(map(str, shards)))
    workers = min(4, threads, len(ordered))
    if not workers:
        return 0

    def remove(shard: str) -> int:
        cached = data_dir / "scores/ligand_3d_by_query" / f"{shard}.parquet"
        if not cached.is_file():
            return 0
        shard_scratch = scratch_dir / shard
        shard_scratch.mkdir(exist_ok=True)
        temporary = shard_scratch / f"{shard}.parquet"
        previous_rows = pq.ParquetFile(cached).metadata.num_rows
        with duckdb.connect() as connection:
            connection.sql(f"SET threads={max(1, threads // workers)}")
            connection.sql(f"SET temp_directory='{shard_scratch.as_posix()}'")
            connection.register("affected_entries", affected_entries)
            temporary.unlink(missing_ok=True)
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
        return int(previous_rows - current_rows)

    if workers == 1:
        return sum(remove(shard) for shard in ordered)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return sum(pool.map(remove, ordered))


def _restore_score_query_caches(
    *, base_data_dir: Path, data_dir: Path, query_ids: set[str]
) -> None:
    """Restore only the old queries needed for target-only score repair."""
    by_shard: dict[str, list[str]] = {}
    for pdb_id in sorted(query_ids):
        source = base_data_dir / "dbs/subdbs/search_db=holo" / f"{pdb_id}.parquet"
        if not source.is_file():
            continue
        destination = data_dir / "dbs/subdbs/search_db=holo" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.is_file():
            link_or_copy_file(source, destination)
        by_shard.setdefault(pdb_id[1:3], []).append(pdb_id)

    sidecars = (
        (
            "ligand_3d_candidate_shards",
            "ligand_3d_candidates/search_db=holo",
            schemas.LIGAND_3D_CANDIDATE_SCHEMA,
        ),
        (
            "ligand_pair_score_shards",
            "ligand_pair_scores/search_db=holo",
            schemas.LIGAND_PAIR_SCORE_SCHEMA,
        ),
    )
    for shard, pdb_ids in by_shard.items():
        for packed_name, query_name, schema in sidecars:
            destinations = {
                pdb_id: data_dir
                / "scores"
                / query_name
                / f"shard={shard}"
                / f"{pdb_id}.parquet"
                for pdb_id in pdb_ids
            }
            missing = {
                pdb_id: path
                for pdb_id, path in destinations.items()
                if not path.is_file()
            }
            if not missing:
                continue
            packed = base_data_dir / "scores" / packed_name / f"shard={shard}.parquet"
            rows = pq.read_table(
                packed,
                columns=schema.names,
                filters=[("query_entry", "in", list(missing))],
            ).cast(schema)
            for pdb_id, destination in missing.items():
                destination.parent.mkdir(parents=True, exist_ok=True)
                selected = rows.filter(pc.equal(rows["query_entry"], pdb_id))
                temporary = destination.with_suffix(".parquet.tmp")
                pq.write_table(selected, temporary, compression="zstd")
                temporary.replace(destination)


def _repair_score_batches(
    data_dir: Path,
    repairs: pd.DataFrame,
    scorer_cfg: DictConfig,
    scratch_dir: Path,
    threads: int,
) -> None:
    batches = [
        (str(index), batch.to_dict("records"))
        for index, batch in repairs.groupby("repair_batch_index", sort=True)
    ]
    workers = min(4, max(1, threads // 4), len(batches))
    if workers < 2:
        for index, batch in batches:
            tasks.repair_batch_scores(
                data_dir=data_dir,
                repairs=batch,
                scorer_cfg=scorer_cfg,
                scratch_dir=scratch_dir / index,
                threads=threads,
            )
        return

    worker_threads = max(1, threads // workers)
    thread_vars = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    previous = {name: os.environ.get(name) for name in thread_vars}
    try:
        for name in thread_vars:
            os.environ[name] = str(worker_threads)
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=get_context("spawn")
        ) as pool:
            futures = {
                pool.submit(
                    tasks.repair_batch_scores,
                    data_dir=data_dir,
                    repairs=batch,
                    scorer_cfg=scorer_cfg,
                    scratch_dir=scratch_dir / index,
                    threads=worker_threads,
                ): index
                for index, batch in batches
            }
            for future in as_completed(futures):
                future.result()
                LOG.info("weekly score repair batch %s complete", futures[future])
    finally:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def _collate_repaired_candidate_shards(
    data_dir: Path,
    shards: list[str],
    scratch_dir: Path,
    threads: int,
    replacement_query_ids: set[str] | None = None,
    source_query_ids: set[str] | None = None,
) -> None:
    workers = min(4, threads, len(shards))
    if not workers:
        return

    def collate(shard: str) -> None:
        tasks.collate_ligand_3d_candidates(
            data_dir=data_dir,
            shards=[shard],
            scratch_dir=scratch_dir / shard,
            threads=max(1, threads // workers),
            replacement_query_ids=replacement_query_ids,
            source_query_ids=source_query_ids,
        )

    if workers == 1:
        for shard in shards:
            collate(shard)
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(collate, shards))


def _completed_holo_query_repair(
    data_dir: Path,
    *,
    affected_manifest: Path,
    full_manifest: Path,
    target_manifest: Path,
) -> dict[str, Any] | None:
    repair_manifest = data_dir / score.SCORE_REPAIR_RELATIVE
    plan_path = repair_manifest.with_suffix(".json")
    completion_path = data_dir / score.SCORE_REPAIR_DROP_RELATIVE.with_suffix(".json")
    if not all(
        path.is_file() for path in (repair_manifest, plan_path, completion_path)
    ):
        return None
    plan = json.loads(plan_path.read_text())
    completion = json.loads(completion_path.read_text())
    expected = {
        "affected_manifest": affected_manifest,
        "additional_full_query_manifest": full_manifest,
        "additional_target_query_manifest": target_manifest,
        "output": repair_manifest,
    }
    if completion.get("status") != "complete" or any(
        plan.get(key) != score._source_signature(path) for key, path in expected.items()
    ):
        return None
    if completion_path.stat().st_mtime_ns < repair_manifest.stat().st_mtime_ns:
        return None
    return cast(dict[str, Any], plan)


def repair_holo_scores(
    data_dir: Path,
    *,
    base_data_dir: Path,
    affected: set[str],
    full_alignment_queries: set[str],
    targeted_alignment_queries: set[str],
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
    _restore_score_query_caches(
        base_data_dir=base_data_dir,
        data_dir=data_dir,
        query_ids=targeted_alignment_queries,
    )
    affected_manifest = _write_pdb_manifest(
        data_dir / "manifests" / "weekly_affected_entries.parquet", affected
    )
    full_manifest = _write_pdb_manifest(
        data_dir / "manifests" / "weekly_full_score_queries.parquet",
        full_alignment_queries,
    )
    target_manifest = _write_pdb_manifest(
        data_dir / "manifests" / "weekly_target_score_queries.parquet",
        targeted_alignment_queries,
    )
    repair_manifest = data_dir / score.SCORE_REPAIR_RELATIVE
    plan = _completed_holo_query_repair(
        data_dir,
        affected_manifest=affected_manifest,
        full_manifest=full_manifest,
        target_manifest=target_manifest,
    )
    if plan is None:
        score.plan_score_batches(
            data_dir,
            batch_size=score_batch_size,
            threads=threads,
            scratch_dir=scratch_dir / "plan-scores",
            max_query_protein_chains=int(scorer_cfg.max_query_protein_chains),
            max_query_proper_ligand_chains=int(
                scorer_cfg.max_query_proper_ligand_chains
            ),
            reuse_mapped_alignments=True,
        )
        try:
            plan = score.plan_score_repair(
                data_dir,
                affected_manifest=affected_manifest,
                additional_full_query_manifest=full_manifest,
                additional_target_query_manifest=target_manifest,
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
        repairs = pd.read_parquet(repair_manifest)
        _repair_score_batches(
            data_dir,
            repairs,
            holo_cfg,
            scratch_dir / "queries",
            threads,
        )
        score.finalize_score_repair_queries(
            data_dir,
            repair_manifest=repair_manifest,
            max_new_drops=0,
        )
    else:
        repairs = pd.read_parquet(repair_manifest)

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
        _collate_repaired_candidate_shards(
            data_dir,
            existing_shards,
            scratch_dir / "candidate-shards",
            threads,
            replacement_query_ids=replacement_queries,
            source_query_ids=active_replacements,
        )
    if new_shards:
        _collate_repaired_candidate_shards(
            data_dir,
            new_shards,
            scratch_dir / "candidate-shards",
            threads,
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
        affected=_changed_ligand_archive_entries(base_data_dir, data_dir, affected),
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
    if int(ligand_plan["pair_count"]):
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
            data_dir=data_dir,
            shards=shards,
            scratch_dir=scratch_dir / "merge-ligand-3d",
            threads=threads,
        )
    repair_report = score.finalize_score_repair_artifacts(
        data_dir,
        repair_manifest=repair_manifest,
        require_packed_scores=False,
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
    output = data_dir / "scores" / f"search_db={search_db}" / f"{search_db}.parquet"
    if search_db == "apo":
        _merge_apo_scores(
            data_dir,
            output=output,
            new_queries=queries,
            replaced_queries=set(queries) | affected.difference(active),
            scratch_dir=scratch_dir / "apo-partition",
            threads=threads,
            memory_limit=memory_limit,
        )
    else:
        source_dir = data_dir / "dbs" / "subdbs" / f"search_db={search_db}"
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


def _merge_apo_scores(
    data_dir: Path,
    *,
    output: Path,
    new_queries: list[str],
    replaced_queries: set[str],
    scratch_dir: Path,
    threads: int,
    memory_limit: str,
) -> None:
    """Replace changed apo queries in the compact release score table."""
    if not replaced_queries:
        return
    source_dir = data_dir / "dbs/subdbs/search_db=apo"
    new_paths = [
        source_dir / f"{pdb_id}.parquet"
        for pdb_id in new_queries
        if (source_dir / f"{pdb_id}.parquet").is_file()
    ]
    if not output.is_file() and not new_paths:
        return
    scratch_dir.mkdir(parents=True, exist_ok=True)
    local_output = scratch_dir / "apo.parquet"
    local_output.unlink(missing_ok=True)

    def sql_path(path: Path) -> str:
        return path.as_posix().replace("'", "''")

    selects = []
    if output.is_file():
        selects.append(
            f"""
            SELECT old.* FROM read_parquet('{sql_path(output)}',
                                           hive_partitioning=false) AS old
            ANTI JOIN replaced_entries
              ON split_part(old.query_system, '__', 1) = replaced_entries.pdb_id
            """
        )
    if new_paths:
        paths = ", ".join(f"'{sql_path(path)}'" for path in new_paths)
        selects.append(
            f"SELECT scores.*, 'apo'::VARCHAR AS search_db "
            f"FROM read_parquet([{paths}], union_by_name=true, "
            "hive_partitioning=false) AS scores"
        )
    with duckdb.connect() as connection:
        connection.execute(f"SET threads={threads}")
        connection.execute(f"SET memory_limit='{memory_limit}'")
        connection.execute(f"SET temp_directory='{sql_path(scratch_dir)}'")
        connection.execute("SET preserve_insertion_order=false")
        connection.register(
            "replaced_entries", pd.DataFrame({"pdb_id": sorted(replaced_queries)})
        )
        connection.execute(
            f"""
            COPY (
                SELECT
                    query_system::VARCHAR AS query_system,
                    query_ligand_id::VARCHAR AS query_ligand_id,
                    target_system::VARCHAR AS target_system,
                    target_ligand_id::VARCHAR AS target_ligand_id,
                    protein_mapping::VARCHAR AS protein_mapping,
                    mapping::VARCHAR AS mapping,
                    protein_mapper::VARCHAR AS protein_mapper,
                    source::VARCHAR AS source,
                    metric::VARCHAR AS metric,
                    similarity::TINYINT AS similarity,
                    search_db::VARCHAR AS search_db
                FROM ({" UNION ALL BY NAME ".join(selects)})
            ) TO '{sql_path(local_output)}' (
                FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500000
            )
            """
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    install = output.with_suffix(".parquet.tmp")
    copyfile(local_output, install)
    install.replace(output)
    local_output.unlink(missing_ok=True)


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
        "interface_apo_structures": str(
            data_dir / RELEASE_PATHS["interface_apo_structures"]
        ),
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
    replacements = affected | full_alignment_queries
    representatives = (
        pq.read_table(
            data_dir / tasks.INTERFACE_REPRESENTATIVES_RELATIVE,
            columns=["entry_pdb_id"],
            filters=[("entry_pdb_id", "in", sorted(replacements))],
        )
        .column("entry_pdb_id")
        .to_pylist()
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
    current_queries = set(map(str, representatives)) & planned_queries - skipped_queries
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
        if not removed and repair_path is None:
            continue
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


def _ligand_chemistry_unchanged(base: Path, data_dir: Path, affected: set[str]) -> bool:
    """Check whether changed entries retain the same proper ligands and SMILES."""
    columns = [
        "entry_pdb_id",
        "ligand_id",
        "ligand_smiles",
        "ligand_is_proper",
        "system_type",
    ]

    def ligands(path: Path) -> set[tuple[str, str, str]]:
        rows = pd.read_parquet(
            path,
            columns=columns,
            filters=[("entry_pdb_id", "in", sorted(affected))],
        )
        rows = rows[
            rows["ligand_is_proper"].fillna(False) & rows["system_type"].eq("holo")
        ]
        return set(
            zip(
                rows["entry_pdb_id"].astype(str),
                rows["ligand_id"].astype(str),
                rows["ligand_smiles"].astype(str),
                strict=True,
            )
        )

    return ligands(base / RELEASE_PATHS["annotation_table"]) == ligands(
        data_dir / RELEASE_PATHS["annotation_table"]
    )


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
    shard_dir = data_dir / "exports" / "ligand_similarity_scores"
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


def _isolated_ligand_score_export(
    data_dir: Path, *, scratch_dir: Path, threads: int, memory_limit: str
) -> dict[str, Any]:
    """Release the large export's memory before continuing the update."""
    with ProcessPoolExecutor(max_workers=1, mp_context=get_context("spawn")) as pool:
        return pool.submit(
            refresh_ligand_score_export,
            data_dir,
            scratch_dir=scratch_dir,
            threads=threads,
            memory_limit=memory_limit,
        ).result()


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
        ("interface", list(clusters.INTERFACE_CLUSTER_METRICS)),
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


def _rebase_reused_scoring_artifacts(data_dir: Path, report: dict[str, Any]) -> None:
    """Keep unchanged scoring artifacts valid after entry tables are rewritten."""
    if not report.get("search_inputs", {}).get("reused_base"):
        return
    lookup = data_dir / tasks.ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    manifest = read_json_cache(
        data_dir / tasks.ALIGNMENT_CHAIN_LOOKUP_MANIFEST_RELATIVE
    )
    if manifest is None or manifest.get("output") != (
        tasks._interface_representative_output_signature(lookup)
    ):
        raise ValueError("reused alignment lookup changed during the weekly update")
    tasks._refresh_representative_source_manifests(data_dir)
    tasks._write_alignment_chain_lookup_manifest(data_dir)


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
        if report.get("search_inputs", {}).get("reused_base"):
            planned = json.loads((plan_dir / "plan.json").read_text())
            _refresh_foldseek_source_manifest(
                workspace,
                nextgen_root=Path(planned["nextgen_root"]),
                affected=set(report["affected_pdb_ids"]),
            )
        marker_path = workspace / "index" / "collation.json"
        marker = json.loads(marker_path.read_text())
        if marker.get("status") != "complete":
            if report["status"] == "complete":
                raise ValueError("completed weekly update has an incomplete index")
        else:
            _rebase_reused_scoring_artifacts(workspace, report)
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
    alignment_databases = list(search_databases)
    if "apo" in alignment_databases:
        alignment_databases.insert(
            alignment_databases.index("apo") + 1, "interface_apo"
        )
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

    _refresh_foldseek_source_manifest(
        workspace,
        nextgen_root=Path(plan["nextgen_root"]),
        affected=affected,
    )
    unchanged_scoring = _unchanged_scoring_entries(base, workspace, affected)
    scoring_affected = affected - unchanged_scoring
    report["unchanged_scoring_entries"] = sorted(unchanged_scoring)
    write_json_atomic(state_path, report)
    if not scoring_affected and set(report["completed_stages"]) == {
        "entries_and_archives"
    }:
        _finish_stage(state_path, report, "search_inputs", {"reused_base": True})
        _finish_stage(
            state_path,
            report,
            "alignments",
            {"full_queries": {search_db: [] for search_db in alignment_databases}},
        )
        _finish_stage(state_path, report, "scores", {"reused_base": True})
        _finish_stage(state_path, report, "ligand_chemistry", {"reused_base": True})
        _finish_stage(
            state_path,
            report,
            "cluster_assignments",
            {
                "ligand": str(workspace / "index/ligand_clusters.parquet"),
                "interface": str(workspace / "index/interface_clusters.parquet"),
            },
        )

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
        if pdb_search_databases:
            score.make_foldseek_input_manifest(
                workspace, Path(plan["nextgen_root"]) / "data/entries/divided"
            )
            score.make_mmseqs_input_fasta(workspace)
        shadowed = _prepare_weekly_search_overlay(
            base=base,
            data_dir=workspace,
            affected=scoring_affected,
            nextgen_root=Path(plan["nextgen_root"]),
            search_databases=[db for db in alignment_databases if db != "pred"],
            scratch_dir=scratch_dir / "weekly-search-overlay",
            threads=threads,
        )
        sequence_clusters = extend_protein_clusters(
            base=base,
            data_dir=workspace,
            affected=scoring_affected,
            backend="mmseqs",
            scratch_dir=scratch_dir / "protein-sequence-clusters",
            threads=threads,
            threshold=float(cfg.flow.protein_sequence_cluster_identity),
            coverage=float(cfg.flow.protein_cluster_coverage),
        )
        structure_clusters = extend_protein_clusters(
            base=base,
            data_dir=workspace,
            affected=scoring_affected,
            backend="foldseek",
            scratch_dir=scratch_dir / "protein-structure-clusters",
            threads=threads,
            threshold=float(cfg.flow.protein_structure_cluster_lddt),
            coverage=float(cfg.flow.protein_cluster_coverage),
        )
        search_details: dict[str, Any] = {
            "alignment_chain_lookup": str(lookup),
            "protein_plan": protein_plan,
            "protein_sequence_clusters": str(sequence_clusters),
            "protein_structure_clusters": str(structure_clusters),
            "search_databases": search_databases,
            "shadowed_entries": len(shadowed),
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
            search_databases=alignment_databases,
            affected=scoring_affected,
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
                    affected=scoring_affected,
                    full_queries=planned_queries[search_db] & scoring_affected,
                    targeted_queries=planned_queries[search_db] - scoring_affected,
                    scorer_cfg=cfg.scorer,
                    foldseek_cfg=cfg.foldseek,
                    mmseqs_cfg=cfg.mmseqs,
                    scratch_dir=scratch_dir / "alignment-repair" / search_db,
                    threads=threads,
                    batch_size=search_batch_size,
                )
            )
            for search_db in alignment_databases
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
    holo_repair_queries = set(
        map(str, report["alignment_repair_queries"].get("holo", []))
    )

    completed = set(report["completed_stages"])
    if "scores" not in completed:
        score_details: dict[str, Any] = {}
        if "holo" in search_databases:
            score_details["ligand"] = repair_holo_scores(
                workspace,
                base_data_dir=base,
                affected=scoring_affected,
                full_alignment_queries=full_query_sets["holo"],
                targeted_alignment_queries=holo_repair_queries.difference(
                    scoring_affected
                ),
                scorer_cfg=cfg.scorer,
                scratch_dir=scratch_dir / "holo-score-repair",
                threads=threads,
                memory_limit=memory_limit,
                score_batch_size=int(cfg.flow.make_batch_scores_batch_size),
                ligand_batch_size=int(cfg.flow.make_ligand_3d_scores_batch_size),
            )
            score_details["interface"] = (
                repair_interface_scores(
                    workspace,
                    affected=scoring_affected,
                    full_alignment_queries=holo_repair_queries,
                    scratch_dir=scratch_dir / "interface-score-repair",
                    threads=threads,
                    memory_limit=memory_limit,
                )
                if _entries_with_interfaces(
                    base, workspace, scoring_affected | holo_repair_queries
                )
                else {"reused_base": True}
            )
        if "apo" in search_databases:
            score_details["apo"] = repair_apo_scores(
                workspace,
                affected=scoring_affected,
                full_alignment_queries=set(
                    map(str, report["alignment_repair_queries"]["apo"])
                ),
                scorer_cfg=cfg.scorer,
                scratch_dir=scratch_dir / "apo-score-repair",
                threads=threads,
                memory_limit=memory_limit,
            )
        if "pred" in search_databases:
            score_details["pred"] = repair_non_holo_scores(
                workspace,
                search_db="pred",
                affected=scoring_affected,
                full_alignment_queries=full_query_sets["pred"],
                scorer_cfg=cfg.scorer,
                scratch_dir=scratch_dir / "pred-score-repair",
                threads=threads,
                memory_limit=memory_limit,
            )
        if "holo" in search_databases:
            score_details["ligand_export"] = _isolated_ligand_score_export(
                workspace,
                scratch_dir=scratch_dir / "ligand-score-export",
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
        if _ligand_chemistry_unchanged(base, workspace, scoring_affected):
            chemistry = {"reused_base": True}
        else:
            chemistry = refresh_ligand_chemistry(
                workspace,
                cfg=cfg,
                scratch_dir=scratch_dir / "ligand-chemistry",
                threads=threads,
            )
        _finish_stage(state_path, report, "ligand_chemistry", chemistry)
    if report["ligand_chemistry"].get("reused_base"):
        for name in (
            "ligands_per_smiles.parquet",
            get_similarity_scores.MHFP6_FINGERPRINT_FILE,
            "ligand_similarity_annotations.parquet",
        ):
            source = base / "fingerprints" / name
            destination = workspace / "fingerprints" / name
            if source.is_file() and not destination.is_file():
                link_or_copy_file(source, destination)

    completed = set(report["completed_stages"])
    if "cluster_assignments" not in completed:
        staged = workspace / "index/.staging/weekly_clusters"
        staged.mkdir(parents=True, exist_ok=True)
        ligand_path = staged / "ligand.parquet"
        interface_path = staged / "interface.parquet"
        extend_ligand_clusters(
            base=base, data_dir=workspace, affected=scoring_affected
        ).to_parquet(ligand_path, index=False)
        if (workspace / "index/interface_annotation_table.parquet").is_file():
            extend_interface_clusters(
                base=base, data_dir=workspace, affected=scoring_affected
            ).to_parquet(interface_path, index=False)
        _finish_stage(
            state_path,
            report,
            "cluster_assignments",
            {"ligand": str(ligand_path), "interface": str(interface_path)},
        )

    if "final_tables" not in set(report["completed_stages"]):
        report["status"] = "finalizing"
        write_json_atomic(state_path, report)
        staged = report["cluster_assignments"]
        interface_path = Path(staged["interface"])
        has_interfaces = (
            workspace / "index/interface_annotation_table.parquet"
        ).is_file()
        if has_interfaces and not interface_path.is_file():
            raise FileNotFoundError(interface_path)
        tasks.finalize_index(
            data_dir=workspace,
            weekly_ligand_clusters=Path(staged["ligand"]),
            weekly_interface_clusters=interface_path if has_interfaces else None,
        )
        _rebase_reused_scoring_artifacts(workspace, report)
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
