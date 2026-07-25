# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Reconstruct bounded similarity subsets from mapped search alignments."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path

import pandas as pd

from plinder.core.scores.entries import EntryView, LigandView, load_entry_views
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import PROTEIN_SIMILARITY_SCHEMA

ALIGNMENT_TYPES = ("foldseek", "mmseqs")
LOG = setup_logger(__name__)


def _pdb_id(identifier: str) -> str:
    return identifier.split("__", maxsplit=1)[0]


def _pdb_shard(pdb_id: str) -> str:
    return pdb_id[-3:-1]


def _require_file(path: Path, *, description: str) -> Path:
    if path.is_file():
        return path
    mode = "offline cache" if cpl.is_offline() else "release cache"
    raise FileNotFoundError(f"missing {description} in {mode}: {path}")


def _release_file(*, relative: str, data_dir: Path | None) -> Path:
    if data_dir is not None:
        return Path(data_dir) / relative
    return cpl.get_plinder_path(rel=relative)


def _alignment_relative_path(
    *, search_db: str, alignment_type: str, pdb_id: str
) -> str:
    cfg = get_config()
    shard = _pdb_shard(pdb_id)
    return (
        f"{cfg.data.alignments}/search_db={search_db}/"
        f"alignment_type={alignment_type}/shard={shard}.parquet"
    )


def prefetch_similarity_alignments(
    query_system_ids: Iterable[str],
    *,
    search_db: str = "holo",
    data_dir: Path | None = None,
) -> dict[str, dict[str, Path]]:
    """Download or verify only the mapped alignments needed by query systems.

    Calling this on an online node populates the normal PLINDER cache. Calling
    it later with offline mode enabled only checks those same local files.
    """
    if search_db != "holo":
        raise NotImplementedError(
            "subset reconstruction currently supports search_db='holo' only"
        )
    query_pdb_ids = sorted({_pdb_id(value) for value in query_system_ids})
    if not query_pdb_ids:
        raise ValueError("query_system_ids must not be empty")

    paths: dict[str, dict[str, Path]] = {}
    resolved: dict[str, Path | None] = {}
    resolution_errors: dict[str, FileNotFoundError] = {}
    for pdb_id in query_pdb_ids:
        paths[pdb_id] = {}
        for alignment_type in ALIGNMENT_TYPES:
            relative = _alignment_relative_path(
                search_db=search_db,
                alignment_type=alignment_type,
                pdb_id=pdb_id,
            )
            if relative not in resolved:
                try:
                    path = _release_file(relative=relative, data_dir=data_dir)
                    resolved[relative] = _require_file(
                        path,
                        description=(
                            f"{search_db} {alignment_type} mapped alignment "
                            f"shard {_pdb_shard(pdb_id)}"
                        ),
                    )
                except FileNotFoundError as exc:
                    resolved[relative] = None
                    resolution_errors[relative] = exc
            resolved_path = resolved[relative]
            if resolved_path is not None:
                paths[pdb_id][alignment_type] = resolved_path

        if not paths[pdb_id]:
            failures = []
            for alignment_type in ALIGNMENT_TYPES:
                relative = _alignment_relative_path(
                    search_db=search_db,
                    alignment_type=alignment_type,
                    pdb_id=pdb_id,
                )
                if relative in resolution_errors:
                    failures.append(str(resolution_errors[relative]))
            raise FileNotFoundError(
                f"no mapped alignment backend is available for {pdb_id}: "
                + "; ".join(failures)
            )
        if len(paths[pdb_id]) != len(ALIGNMENT_TYPES):
            missing = sorted(set(ALIGNMENT_TYPES) - paths[pdb_id].keys())
            LOG.warning(
                f"reconstructing {pdb_id} without unavailable alignment "
                f"backend(s): {', '.join(missing)}"
            )
    return paths


def _load_entry_subset(
    *, pdb_ids: set[str], data_dir: Path | None
) -> tuple[dict[str, EntryView], Path]:
    cfg = get_config()
    index_relative = f"{cfg.data.index}/{cfg.data.index_file}"
    index_path = _require_file(
        _release_file(relative=index_relative, data_dir=data_dir),
        description="annotation index",
    )
    entries = load_entry_views(pdb_ids=pdb_ids, data_dir=data_dir)
    if not entries:
        raise KeyError(f"no annotation rows found for PDB IDs {sorted(pdb_ids)}")
    return entries, index_path.parent.parent


def _selected_ligand_ids(
    entries: dict[str, EntryView], system_ids: set[str]
) -> set[str]:
    ligand_ids: set[str] = set()
    for entry in entries.values():
        for system_id, system in entry.systems.items():
            if system_id in system_ids:
                ligand_ids.update(ligand.id for ligand in system.ligands.values())
    return ligand_ids


def _validate_selection(
    *,
    entries: dict[str, EntryView],
    query_system_ids: set[str],
    target_system_ids: set[str],
    query_ligand_ids: set[str] | None,
    target_ligand_ids: set[str] | None,
) -> None:
    available_system_ids = {
        system_id for entry in entries.values() for system_id in entry.systems
    }
    missing_system_ids = (query_system_ids | target_system_ids) - available_system_ids
    if missing_system_ids:
        raise KeyError(f"unknown system IDs: {sorted(missing_system_ids)}")

    for name, selected, systems in [
        ("query_ligand_ids", query_ligand_ids, query_system_ids),
        ("target_ligand_ids", target_ligand_ids, target_system_ids),
    ]:
        if selected is None:
            continue
        missing = selected - _selected_ligand_ids(entries, systems)
        if missing:
            raise KeyError(
                f"{name} do not belong to selected systems: {sorted(missing)}"
            )


def _canonical_ligand_resolver(
    *, release_root: Path, data_dir: Path | None
) -> Callable[[LigandView], Path | None]:
    def resolve(ligand: LigandView) -> Path | None:
        candidates = [
            release_root
            / "raw_entries"
            / ligand.pdb_id[-3:-1]
            / ligand.pdb_id
            / "ligand_files"
            / f"{ligand.asym_id}.sdf",
            release_root
            / "ligand_archives"
            / ligand.pdb_id
            / "ligand_files"
            / f"{ligand.asym_id}.sdf",
        ]
        existing = next((path for path in candidates if path.is_file()), None)
        if existing is not None:
            return existing

        cfg = get_config()
        code = ligand.pdb_id[-3:-1]
        archive_relative = f"{cfg.data.ligand_archives}/{code}.parquet"
        archive = _require_file(
            _release_file(relative=archive_relative, data_dir=data_dir),
            description=f"canonical ligand archive for {code}",
        )
        packed = pd.read_parquet(
            archive,
            columns=["sdf"],
            filters=[
                ("pdb_id", "==", ligand.pdb_id),
                ("ligand_asym_id", "==", ligand.asym_id),
            ],
        )
        if len(packed) != 1:
            return None
        extracted = (
            archive.parent / ligand.pdb_id / "ligand_files" / f"{ligand.asym_id}.sdf"
        )
        extracted.parent.mkdir(exist_ok=True, parents=True)
        temporary = extracted.with_suffix(".tmp.sdf")
        temporary.write_bytes(packed.iloc[0]["sdf"])
        temporary.replace(extracted)
        return extracted if extracted.is_file() else None

    return resolve


def reconstruct_similarity_scores(
    query_system_ids: Iterable[str],
    target_system_ids: Iterable[str],
    *,
    query_ligand_ids: Iterable[str] | None = None,
    target_ligand_ids: Iterable[str] | None = None,
    search_db: str = "holo",
    include_shape: bool = True,
    data_dir: Path | None = None,
) -> pd.DataFrame:
    """Reconstruct directed ligand-level scores for a bounded system subset.

    Only mapped Foldseek/MMseqs shards and annotation rows for the requested
    systems are loaded. Canonical ligand SDF archives are accessed lazily, and
    only after a ligand pair has positive pocket coverage.
    """
    query_systems = set(query_system_ids)
    target_systems = set(target_system_ids)
    if not query_systems or not target_systems:
        raise ValueError("query_system_ids and target_system_ids must not be empty")
    query_ligands = set(query_ligand_ids) if query_ligand_ids is not None else None
    target_ligands = set(target_ligand_ids) if target_ligand_ids is not None else None

    alignment_paths = prefetch_similarity_alignments(
        query_systems,
        search_db=search_db,
        data_dir=data_dir,
    )
    selected_pdb_ids = {
        _pdb_id(system_id) for system_id in query_systems | target_systems
    }
    entries, release_root = _load_entry_subset(
        pdb_ids=selected_pdb_ids,
        data_dir=data_dir,
    )
    _validate_selection(
        entries=entries,
        query_system_ids=query_systems,
        target_system_ids=target_systems,
        query_ligand_ids=query_ligands,
        target_ligand_ids=target_ligands,
    )

    # Imported lazily so ordinary score parquet queries do not import ingest
    # machinery or RDKit shape-alignment code.
    from plinder.data.annotations.get_similarity_scores import Scorer

    scorer = Scorer(
        entries=entries,
        source_to_full_db_file={},
        db_dir=release_root,
        scores_dir=release_root,
        ligand_sdf_resolver=_canonical_ligand_resolver(
            release_root=release_root,
            data_dir=data_dir,
        ),
    )
    frames: list[pd.DataFrame] = []
    for pdb_id, paths in alignment_paths.items():
        query_ids = {
            system_id for system_id in query_systems if _pdb_id(system_id) == pdb_id
        }
        frame = scorer.aggregate_scores(
            pdb_id,
            search_db=search_db,
            data_dir=release_root if include_shape else None,
            source_to_aln_file={
                f"{search_db}_{alignment_type}": path
                for alignment_type, path in paths.items()
            },
            query_system_ids=query_ids,
            query_ligand_ids=query_ligands,
            target_system_ids=target_systems,
            target_ligand_ids=target_ligands,
        )
        if frame is not None and not frame.empty:
            frames.append(frame)

    if not frames:
        return pd.DataFrame(columns=PROTEIN_SIMILARITY_SCHEMA.names)
    return (
        pd.concat(frames, ignore_index=True)
        .sort_values(
            [
                "similarity",
                "query_system",
                "query_ligand_id",
                "target_system",
                "target_ligand_id",
            ],
            ascending=[False, True, True, True, True],
        )
        .reset_index(drop=True)
    )
