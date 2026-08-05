# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Resolve the minimal release assets needed to score custom structures."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from plinder.core.utils import cpl
from plinder.core.utils.config import get_config

SEARCH_BACKENDS = ("foldseek", "mmseqs")
EXACT_CLUSTER_CONTRACT: Mapping[str, object] = {
    "identity": 1.0,
    "coverage": 1.0,
    "coverage_mode": 0,
    "compressed_search_target": False,
}


@dataclass(frozen=True)
class SearchDatabaseBundle:
    """Portable exact-clustered target database for one search backend."""

    backend: str
    root: Path
    search_target: Path
    conversion_target: Path
    cluster_alignments: Path | None
    manifest: Mapping[str, object]


@dataclass(frozen=True)
class CustomScoringAssets:
    """Release files required by custom-structure scoring."""

    annotation_table: Path
    entry_chains: Path
    interface_annotations: Path
    alignment_chain_lookup: Path
    search_databases: Mapping[str, SearchDatabaseBundle]
    ligand_archives: Mapping[str, Path]


def _require_file(path: Path, *, description: str) -> Path:
    if path.is_file():
        return path
    mode = "offline cache" if cpl.is_offline() else "release cache"
    raise FileNotFoundError(f"missing {description} in {mode}: {path}")


def _release_path(*, relative: str, data_dir: Path | None) -> Path:
    if data_dir is not None:
        return Path(data_dir) / relative
    return cpl.get_plinder_path(rel=relative)


def _manifest_member(root: Path, value: object, *, field: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"invalid {field} in {root / 'exact_cluster.json'}")
    relative = Path(value)
    if relative.is_absolute() or relative.name != value:
        raise ValueError(f"unsafe {field} in {root / 'exact_cluster.json'}: {value}")
    return root / relative


def _require_database_prefix(
    prefix: Path,
    *,
    description: str,
    indexed: bool = False,
) -> None:
    _require_file(prefix.with_suffix(".dbtype"), description=description)
    if indexed:
        _require_file(Path(f"{prefix}.idx.dbtype"), description=f"indexed {description}")


def _validate_portable_links(root: Path) -> None:
    """Reject database links that would break after moving the bundle."""
    resolved_root = root.resolve()
    for path in root.rglob("*"):
        if not path.is_symlink():
            continue
        link = path.readlink()
        try:
            target = path.resolve(strict=True)
        except FileNotFoundError as exc:
            raise ValueError(f"broken database link: {path} -> {link}") from exc
        if link.is_absolute() or not target.is_relative_to(resolved_root):
            raise ValueError(f"non-portable database link: {path} -> {path.readlink()}")


def _local_search_database_root(data_dir: Path, backend: str) -> Path:
    """Resolve either a published bundle or an unmodified ingest output."""
    cfg = get_config()
    published = data_dir / str(cfg.data.search_databases) / f"holo_{backend}"
    if published.is_dir():
        return published
    ingest = data_dir / "dbs" / "subdbs" / f"holo_{backend}"
    return ingest if ingest.is_dir() else published


def resolve_search_database(
    backend: str,
    *,
    data_dir: Path | None = None,
) -> SearchDatabaseBundle:
    """Download or validate one portable PLINDER search database bundle."""
    if backend not in SEARCH_BACKENDS:
        raise ValueError(
            f"unsupported search backend {backend!r}; expected one of {SEARCH_BACKENDS}"
        )
    cfg = get_config()
    if data_dir is None:
        root = cpl.get_plinder_path(
            rel=f"{cfg.data.search_databases}/holo_{backend}"
        )
    else:
        root = _local_search_database_root(Path(data_dir), backend)
    root = Path(root)
    manifest_path = _require_file(
        root / "exact_cluster.json",
        description=f"{backend} search database manifest",
    )
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid search database manifest: {manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise ValueError(f"invalid search database manifest: {manifest_path}")
    if manifest.get("alignment_type") != backend:
        raise ValueError(
            f"search database backend mismatch in {manifest_path}: "
            f"{manifest.get('alignment_type')!r}"
        )
    if manifest.get("portable") is not True:
        raise ValueError(f"search database is not portable: {manifest_path}")
    mismatches = {
        key: manifest.get(key)
        for key, expected in EXACT_CLUSTER_CONTRACT.items()
        if manifest.get(key) != expected
    }
    if mismatches:
        raise ValueError(
            f"search database does not satisfy the exact-cluster contract in "
            f"{manifest_path}: {mismatches}"
        )

    search_target = _manifest_member(
        root, manifest.get("search_target"), field="search_target"
    )
    conversion_target = _manifest_member(
        root, manifest.get("conversion_target"), field="conversion_target"
    )
    cluster_alignments = None
    if backend == "mmseqs":
        cluster_alignments = _manifest_member(
            root,
            manifest.get("cluster_alignments"),
            field="cluster_alignments",
        )

    _require_database_prefix(
        search_target,
        description=f"{backend} representative search target",
        indexed=True,
    )
    _require_database_prefix(
        conversion_target,
        description=f"{backend} conversion target",
    )
    if cluster_alignments is not None:
        _require_database_prefix(
            cluster_alignments,
            description="MMseqs exact-cluster expansion database",
        )
    _validate_portable_links(root)
    return SearchDatabaseBundle(
        backend=backend,
        root=root,
        search_target=search_target,
        conversion_target=conversion_target,
        cluster_alignments=cluster_alignments,
        manifest=manifest,
    )


def _pdb_shard(identifier: str) -> str:
    pdb_id = str(identifier).split("__", maxsplit=1)[0].strip().lower()
    if len(pdb_id) < 4 or not pdb_id.replace("_", "").isalnum():
        raise ValueError(f"invalid PDB or system identifier: {identifier!r}")
    return pdb_id[-3:-1]


def resolve_ligand_archives(
    pdb_or_system_ids: Iterable[str],
    *,
    data_dir: Path | None = None,
) -> dict[str, Path]:
    """Download or validate only ligand-coordinate shards needed by targets."""
    cfg = get_config()
    codes = sorted({_pdb_shard(value) for value in pdb_or_system_ids})
    archives: dict[str, Path] = {}
    for code in codes:
        relative = f"{cfg.data.ligand_archives}/{code}.parquet"
        archives[code] = _require_file(
            _release_path(relative=relative, data_dir=data_dir),
            description=f"canonical ligand archive for shard {code}",
        )
    return archives


def resolve_custom_scoring_assets(
    *,
    data_dir: Path | None = None,
    backends: Iterable[str] = SEARCH_BACKENDS,
    ligand_pdb_ids: Iterable[str] = (),
) -> CustomScoringAssets:
    """Download or validate the bounded asset set for custom scoring.

    Ligand archives are intentionally omitted unless ``ligand_pdb_ids`` are
    supplied. A caller can therefore run protein search first and fetch only
    the coordinate shards containing target ligands with positive pocket
    coverage.
    """
    cfg = get_config()
    index_files = {
        "annotation_table": cfg.data.index_file,
        "entry_chains": cfg.data.entry_chain_file,
        "interface_annotations": cfg.data.interface_file,
        "alignment_chain_lookup": cfg.data.alignment_chain_lookup_file,
    }
    resolved_index: dict[str, Path] = {}
    for name, filename in index_files.items():
        relative = f"{cfg.data.index}/{filename}"
        resolved_index[name] = _require_file(
            _release_path(relative=relative, data_dir=data_dir),
            description=name.replace("_", " "),
        )

    selected_backends = tuple(dict.fromkeys(backends))
    databases = {
        backend: resolve_search_database(backend, data_dir=data_dir)
        for backend in selected_backends
    }
    ligand_archives = resolve_ligand_archives(
        ligand_pdb_ids,
        data_dir=data_dir,
    )
    return CustomScoringAssets(
        annotation_table=resolved_index["annotation_table"],
        entry_chains=resolved_index["entry_chains"],
        interface_annotations=resolved_index["interface_annotations"],
        alignment_chain_lookup=resolved_index["alignment_chain_lookup"],
        search_databases=databases,
        ligand_archives=ligand_archives,
    )
