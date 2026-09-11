# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Reconstruct bounded similarity subsets from mapped search alignments."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.scores.entries import (
    EntryView,
    InterfaceView,
    LigandView,
    load_entry_views,
)
from plinder.core.utils import cpl
from plinder.core.utils.log import setup_logger
from plinder.core.utils.schemas import (
    INTERFACE_SIMILARITY_SCHEMA,
    PROTEIN_SIMILARITY_SCHEMA,
)

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


def _release_file(
    name: str,
    *,
    data_dir: Path | None,
    description: str,
    **parameters: str,
) -> Path:
    release = PlinderRelease(data_dir)
    if data_dir is not None:
        return _require_file(
            release.path(name, **parameters),
            description=description,
        )
    try:
        return release.fetch(name, **parameters)
    except FileNotFoundError as exc:
        mode = "offline cache" if cpl.is_offline() else "release cache"
        raise FileNotFoundError(
            f"missing {description} in {mode}: {release.path(name, **parameters)}"
        ) from exc


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
    resolved: dict[tuple[str, str], Path | None] = {}
    resolution_errors: dict[tuple[str, str], FileNotFoundError] = {}
    for pdb_id in query_pdb_ids:
        paths[pdb_id] = {}
        shard = _pdb_shard(pdb_id)
        for alignment_type in ALIGNMENT_TYPES:
            key = (alignment_type, shard)
            if key not in resolved:
                try:
                    resolved[key] = _release_file(
                        "alignment_shard",
                        data_dir=data_dir,
                        description=(
                            f"{search_db} {alignment_type} mapped alignment "
                            f"shard {shard}"
                        ),
                        search_db=search_db,
                        alignment_type=alignment_type,
                        shard=shard,
                    )
                except FileNotFoundError as exc:
                    resolved[key] = None
                    resolution_errors[key] = exc
            resolved_path = resolved[key]
            if resolved_path is not None:
                paths[pdb_id][alignment_type] = resolved_path

        if not paths[pdb_id]:
            failures = []
            for alignment_type in ALIGNMENT_TYPES:
                key = (alignment_type, shard)
                if key in resolution_errors:
                    failures.append(str(resolution_errors[key]))
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
    index_path = _release_file(
        "annotation_table",
        data_dir=data_dir,
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


def _selected_interfaces(
    entries: Mapping[str, EntryView], interface_ids: set[str]
) -> dict[str, InterfaceView]:
    available = {
        interface_id: interface
        for entry in entries.values()
        for interface_id, interface in entry.interfaces.items()
    }
    missing = interface_ids - available.keys()
    if missing:
        raise KeyError(f"unknown interface system IDs: {sorted(missing)}")
    return {interface_id: available[interface_id] for interface_id in interface_ids}


def _load_interface_alignments(
    alignment_paths: Mapping[str, Mapping[str, Path]],
    *,
    target_pdb_ids: set[str],
) -> pd.DataFrame:
    """Load only compact residue maps needed by an interface subset."""
    columns = [
        "query_entry",
        "target_entry",
        "query_chain_mapped",
        "target_chain_mapped",
        "source",
        "query_selected_residue_numbers",
        "target_selected_residue_numbers",
    ]
    frames: list[pd.DataFrame] = []
    targets = sorted(target_pdb_ids)
    for query_pdb_id, paths in alignment_paths.items():
        for alignment_type, path in paths.items():
            frame = pd.read_parquet(
                path,
                columns=columns,
                filters=[
                    ("query_entry", "==", query_pdb_id),
                    ("target_entry", "in", targets),
                ],
            )
            if frame.empty:
                continue
            frame["source"] = alignment_type
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True)


def _interface_side_index(
    interfaces: Mapping[str, InterfaceView],
) -> dict[tuple[str, str], list[tuple[str, int, frozenset[int]]]]:
    sides: dict[tuple[str, str], list[tuple[str, int, frozenset[int]]]] = defaultdict(
        list
    )
    for interface in interfaces.values():
        for side, (chain, residues) in enumerate(
            (
                (interface.chain_1, interface.chain_1_residue_number_to_index),
                (interface.chain_2, interface.chain_2_residue_number_to_index),
            ),
            start=1,
        ):
            asym_id = chain.split(".", maxsplit=1)[-1]
            sides[(interface.pdb_id, asym_id)].append((
                interface.id,
                side,
                frozenset(residues),
            ))
    return sides


def _as_int_list(value: object) -> list[int]:
    if (
        value is None
        or isinstance(value, (str, bytes))
        or not isinstance(value, Iterable)
    ):
        return []
    return [int(item) for item in value]


def calculate_interface_similarity_scores(
    alignments: pd.DataFrame,
    *,
    query_interfaces: Mapping[str, InterfaceView],
    target_interfaces: Mapping[str, InterfaceView],
) -> pd.DataFrame:
    """Calculate directional two-chain interface coverage from compact maps.

    For each backend, the direct and swapped chain assignments are evaluated
    as the product of the two directed residue coverages. The better complete
    assignment is retained, followed by the better Foldseek/MMseqs backend.
    """
    query_sides = _interface_side_index(query_interfaces)
    target_sides = _interface_side_index(target_interfaces)
    coverage: dict[tuple[str, str, str, int, int], float] = {}
    for row in alignments.itertuples(index=False):
        source = str(row.source)
        query_key = (str(row.query_entry), str(row.query_chain_mapped))
        target_key = (str(row.target_entry), str(row.target_chain_mapped))
        matching_query_sides = query_sides.get(query_key, [])
        matching_target_sides = target_sides.get(target_key, [])
        if not matching_query_sides or not matching_target_sides:
            continue
        residue_pairs = list(
            zip(
                _as_int_list(row.query_selected_residue_numbers),
                _as_int_list(row.target_selected_residue_numbers),
            )
        )
        if not residue_pairs:
            continue
        for query_id, query_side, query_residues in matching_query_sides:
            if not query_residues:
                continue
            for target_id, target_side, target_residues in matching_target_sides:
                matched_query_residues = {
                    query_number
                    for query_number, target_number in residue_pairs
                    if query_number in query_residues
                    and target_number in target_residues
                }
                if not matched_query_residues:
                    continue
                value = len(matched_query_residues) / len(query_residues)
                key = (source, query_id, target_id, query_side, target_side)
                coverage[key] = max(coverage.get(key, 0.0), value)

    def assignment(
        source: str,
        query: InterfaceView,
        target: InterfaceView,
        *,
        swapped: bool,
    ) -> tuple[float, float, float, str] | None:
        target_sides_order = (2, 1) if swapped else (1, 2)
        values = [
            coverage.get((source, query.id, target.id, query_side, target_side))
            for query_side, target_side in zip((1, 2), target_sides_order)
        ]
        first, second = values
        if first is None or second is None:
            return None
        target_chains = (target.chain_2, target.chain_1) if swapped else target.chains
        mapping = ";".join(
            f"{query_chain}:{target_chain}"
            for query_chain, target_chain in zip(query.chains, target_chains)
        )
        return first * second, first, second, mapping

    records: list[dict[str, object]] = []
    for query in query_interfaces.values():
        for target in target_interfaces.values():
            if query.id == target.id:
                continue
            backend_scores: dict[str, tuple[float, float, float, str]] = {}
            for source in ALIGNMENT_TYPES:
                candidates = [
                    result
                    for swapped in (False, True)
                    if (result := assignment(source, query, target, swapped=swapped))
                    is not None
                ]
                if candidates:
                    backend_scores[source] = max(
                        candidates,
                        key=lambda result: result[0],
                    )
            if not backend_scores:
                continue
            best_score = max(value[0] for value in backend_scores.values())
            best_sources = [
                source
                for source, value in backend_scores.items()
                if abs(value[0] - best_score) < 1e-12
            ]
            source = (
                "both" if len(best_sources) == len(ALIGNMENT_TYPES) else best_sources[0]
            )
            mapping_source = (
                "foldseek" if "foldseek" in best_sources else best_sources[0]
            )
            records.append({
                "query_system": query.id,
                "target_system": target.id,
                "mapping": backend_scores[mapping_source][3],
                "source": source,
                "metric": "interface_qcov",
                "iface1_qcov": backend_scores[mapping_source][1],
                "iface2_qcov": backend_scores[mapping_source][2],
                "similarity": max(0, min(100, round(best_score * 100))),
            })
    if not records:
        return pd.DataFrame(columns=INTERFACE_SIMILARITY_SCHEMA.names)
    return (
        pd.DataFrame
        .from_records(records, columns=INTERFACE_SIMILARITY_SCHEMA.names)
        .sort_values(
            ["similarity", "query_system", "target_system"],
            ascending=[False, True, True],
        )
        .reset_index(drop=True)
    )


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

        code = ligand.pdb_id[-3:-1]
        archive = _release_file(
            "ligand_archive",
            data_dir=data_dir,
            description=f"canonical ligand archive for {code}",
            shard=code,
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
        pd
        .concat(frames, ignore_index=True)
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


def reconstruct_interface_similarity_scores(
    query_interface_ids: Iterable[str],
    target_interface_ids: Iterable[str],
    *,
    search_db: str = "holo",
    data_dir: Path | None = None,
) -> pd.DataFrame:
    """Reconstruct directed interface coverage for a bounded interface subset.

    Only the mapped alignment shards for the requested query PDB entries and
    chain and interface rows are loaded. The result is therefore
    reconstructable from the compact public artifacts without the private
    all-vs-all score table used for release clustering.
    """
    query_ids = set(query_interface_ids)
    target_ids = set(target_interface_ids)
    if not query_ids or not target_ids:
        raise ValueError(
            "query_interface_ids and target_interface_ids must not be empty"
        )
    alignment_paths = prefetch_similarity_alignments(
        query_ids,
        search_db=search_db,
        data_dir=data_dir,
    )
    selected_pdb_ids = {_pdb_id(value) for value in query_ids | target_ids}
    entries, _ = _load_entry_subset(
        pdb_ids=selected_pdb_ids,
        data_dir=data_dir,
    )
    query_interfaces = _selected_interfaces(entries, query_ids)
    target_interfaces = _selected_interfaces(entries, target_ids)
    alignments = _load_interface_alignments(
        alignment_paths,
        target_pdb_ids={interface.pdb_id for interface in target_interfaces.values()},
    )
    return calculate_interface_similarity_scores(
        alignments,
        query_interfaces=query_interfaces,
        target_interfaces=target_interfaces,
    )
