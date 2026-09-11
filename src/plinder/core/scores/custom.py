# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Annotate and score custom mmCIF structures against PLINDER."""

from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
from bisect import bisect_left
from collections import Counter
from collections.abc import Iterable as IterableABC
from collections.abc import Iterator
from dataclasses import dataclass, replace
from gzip import open as gzip_open
from pathlib import Path
from threading import Lock
from typing import Any, Iterable, Literal, Mapping

import numpy as np
import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.utils import cpl

LOG = logging.getLogger(__name__)
CIF_SEARCH_BACKENDS = ("foldseek", "mmseqs")
SEQUENCE_SEARCH_BACKENDS = ("mmseqs",)
SEARCH_BACKENDS = CIF_SEARCH_BACKENDS
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


@dataclass(frozen=True)
class CustomQueryInputs:
    """Written per-chain coordinates, sequences, and stable metadata."""

    root: Path
    chain_cif_dir: Path
    sequence_fasta: Path
    chain_manifest: Path


@dataclass(frozen=True)
class CustomQueryDatabases:
    """Foldseek/MMseqs query DBs and their backend identifier mapping."""

    root: Path
    inputs: CustomQueryInputs
    databases: Mapping[str, Path]
    identifier_map: Path


@dataclass(frozen=True)
class CustomProteinSearchConfig:
    """Shared filters for custom queries against the PLINDER protein DBs."""

    evalue: float = 0.01
    sensitivity: float = 11.0
    max_seqs: int = 10_000
    coverage: float = 0.0
    min_seq_id: float = 0.0

    def __post_init__(self) -> None:
        if self.max_seqs < 1:
            raise ValueError("max_seqs must be positive")
        if not 0 <= self.coverage <= 1:
            raise ValueError("coverage must be in [0, 1]")
        if not 0 <= self.min_seq_id <= 1:
            raise ValueError("min_seq_id must be in [0, 1]")


@dataclass(frozen=True)
class CustomStructureAnnotations:
    """Custom ligand/interface annotations and their generated ligand SDFs."""

    root: Path
    ligand_sdf_root: Path
    annotation_table: Path
    entry_chains: Path
    interface_annotations: Path
    entries_by_structure: Mapping[str, Any]


@dataclass(frozen=True)
class CustomScoringResult:
    """Files produced by the complete custom-CIF scoring workflow."""

    annotations: CustomStructureAnnotations
    query_inputs: CustomQueryInputs
    query_databases: CustomQueryDatabases
    protein_hits: Mapping[str, Path]
    score_alignments: Mapping[str, Path]
    protein_score_alignments: Mapping[str, Path]
    protein_scores: Path
    aligned_pocket_residues: Path | None
    ligand_scores: Path | None
    interface_scores: Path | None


@dataclass(frozen=True)
class CustomSequenceScoringResult:
    """Files produced by protein-sequence scoring against PLINDER pockets."""

    query_inputs: CustomQueryInputs
    query_databases: CustomQueryDatabases
    protein_hits: Mapping[str, Path]
    protein_score_alignments: Mapping[str, Path]
    protein_scores: Path
    aligned_pocket_residues: Path | None
    sequence_links: Path
    best_sequence_links: Path


@dataclass(frozen=True)
class _SequenceEntry:
    """Small entry stand-in needed while reversing custom alignments."""

    pdb_id: str


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
        _require_file(
            Path(f"{prefix}.idx.dbtype"), description=f"indexed {description}"
        )


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
    published = PlinderRelease(data_dir).path("search_database", backend=backend)
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
    if data_dir is None:
        root = _release_file(
            "search_database",
            data_dir=None,
            description=f"{backend} search database",
            backend=backend,
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
    codes = sorted({_pdb_shard(value) for value in pdb_or_system_ids})
    archives: dict[str, Path] = {}
    for code in codes:
        archives[code] = _release_file(
            "ligand_archive",
            data_dir=data_dir,
            description=f"canonical ligand archive for shard {code}",
            shard=code,
        )
    return archives


def resolve_custom_scoring_assets(
    *,
    data_dir: Path | None = None,
    backends: Iterable[str] = CIF_SEARCH_BACKENDS,
    ligand_pdb_ids: Iterable[str] = (),
) -> CustomScoringAssets:
    """Download or validate the bounded asset set for custom scoring.

    Ligand archives are intentionally omitted unless ``ligand_pdb_ids`` are
    supplied. A caller can therefore run protein search first and fetch only
    the coordinate shards containing target ligands with positive pocket
    coverage.
    """
    index_files = (
        "annotation_table",
        "entry_chains",
        "interface_annotations",
        "alignment_chain_lookup",
    )
    resolved_index: dict[str, Path] = {}
    for name in index_files:
        resolved_index[name] = _release_file(
            name,
            data_dir=data_dir,
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


def _mmcif_stem(path: Path) -> str:
    """Return a stable structure ID while accepting only mmCIF suffixes."""
    name = path.name
    lowered = name.lower()
    for suffix in (".mmcif.gz", ".cif.gz", ".mmcif", ".cif"):
        if lowered.endswith(suffix):
            structure_id = name[: -len(suffix)]
            if not structure_id:
                break
            return structure_id
    raise ValueError(f"custom structure must be an mmCIF file: {path}")


def _clean_protein_sequence(sequence: str, *, context: str) -> str:
    cleaned = "".join(str(sequence).replace(";", "").split()).upper()
    cleaned = cleaned.replace("?", "X")
    if not cleaned or re.fullmatch(r"[A-Z]+", cleaned) is None:
        raise ValueError(f"invalid protein sequence for {context}: {sequence!r}")
    return cleaned


def _fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    """Yield identifier and sequence pairs from a plain or gzipped FASTA."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"missing protein FASTA: {source}")
    handle = (
        gzip_open(source, mode="rt")
        if source.name.lower().endswith(".gz")
        else source.open()
    )
    with handle:
        identifier: str | None = None
        parts: list[str] = []
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(parts)
                header = text[1:].strip()
                if not header:
                    raise ValueError(
                        f"empty FASTA header in {source} at line {line_number}"
                    )
                identifier = header.split(maxsplit=1)[0]
                parts = []
            else:
                if identifier is None:
                    raise ValueError(
                        f"sequence precedes the first FASTA header in {source} "
                        f"at line {line_number}"
                    )
                parts.append(text)
        if identifier is not None:
            yield identifier, "".join(parts)


def write_custom_sequence_query_files(
    sequence_fasta: Path,
    *,
    work_dir: Path,
    min_chain_length: int = 12,
) -> CustomQueryInputs:
    """Write stable internal query IDs for a protein FASTA.

    Each FASTA record represents one ligand-free protein chain. Original
    identifiers remain in the manifest and generated result tables; internal
    IDs keep search-database identifiers unambiguous.
    """
    if min_chain_length < 1:
        raise ValueError("min_chain_length must be positive")
    source = Path(sequence_fasta)
    root = Path(work_dir) / "query_inputs"
    if source.resolve().is_relative_to(root.resolve()):
        raise ValueError("the protein FASTA must be outside the generated work tree")
    records = list(_fasta_records(source))
    if not records:
        raise ValueError(f"protein FASTA contains no records: {source}")
    identifiers = [identifier for identifier, _ in records]
    duplicates = sorted(
        identifier for identifier, count in Counter(identifiers).items() if count > 1
    )
    if duplicates:
        raise ValueError(f"protein FASTA repeats identifiers: {duplicates[:10]}")

    if root.exists():
        shutil.rmtree(root)
    chain_cif_dir = root / "chains"
    chain_cif_dir.mkdir(parents=True)
    written_fasta = root / "query_sequences.fasta"
    chain_manifest = root / "query_chains.parquet"
    rows: list[dict[str, Any]] = []
    fasta_rows: list[str] = []
    skipped: dict[str, str] = {}
    for sequence_id, raw_sequence in records:
        sequence = _clean_protein_sequence(
            raw_sequence,
            context=f"FASTA record {sequence_id}",
        )
        if len(sequence) < min_chain_length:
            skipped[sequence_id] = (
                f"protein length {len(sequence)} is below {min_chain_length}"
            )
            continue
        query_id = f"cq{len(rows):08d}"
        rows.append({
            "query_id": query_id,
            "query_chain_id": f"{query_id}__A",
            "structure_id": query_id,
            "sequence_id": sequence_id,
            "source_fasta": str(source.resolve()),
            "chain_asym_id": "A",
            "sequence": sequence,
            "sequence_length": len(sequence),
            "sequence_source": "polymer",
            "resolved_residue_numbers": [],
        })
        fasta_rows.append(f">{query_id}\n{sequence}\n")
    if not rows:
        raise ValueError(
            f"no protein sequences of at least {min_chain_length} residues were "
            f"found; record diagnostics={skipped}"
        )
    written_fasta.write_text("".join(fasta_rows))
    pd.DataFrame(rows).to_parquet(chain_manifest, index=False)
    return CustomQueryInputs(
        root=root,
        chain_cif_dir=chain_cif_dir,
        sequence_fasta=written_fasta,
        chain_manifest=chain_manifest,
    )


def _protein_asym_sequences(block: Any) -> dict[str, str]:
    """Map label asym IDs to canonical sequences for protein entities only."""
    from plinder.data.annotations.cif_utils import get_label_asym_sequences

    if "struct_asym" not in block or "entity_poly" not in block:
        return {}
    struct_asym = block["struct_asym"]
    entity_poly = block["entity_poly"]
    if not {"id", "entity_id"}.issubset(struct_asym):
        return {}
    if not {"entity_id", "type"}.issubset(entity_poly):
        return {}
    protein_entities = {
        str(entity_id)
        for entity_id, polymer_type in zip(
            entity_poly["entity_id"].as_array(),
            entity_poly["type"].as_array(),
        )
        if "polypeptide" in str(polymer_type).lower()
    }
    sequences = get_label_asym_sequences(block)
    return {
        str(asym_id): _clean_protein_sequence(
            sequences[str(asym_id)], context=f"chain {asym_id}"
        )
        for asym_id, entity_id in zip(
            struct_asym["id"].as_array(),
            struct_asym["entity_id"].as_array(),
        )
        if str(entity_id) in protein_entities and str(asym_id) in sequences
    }


def _select_assembly_ids(
    available: Iterable[str], selected: Iterable[str] | None
) -> list[str]:
    available_ids = list(dict.fromkeys(str(value) for value in available))
    if selected is None:
        return available_ids
    selected_values = [selected] if isinstance(selected, str) else selected
    selected_ids = list(dict.fromkeys(str(value) for value in selected_values))
    if not selected_ids:
        raise ValueError("assembly_ids must not be empty when provided")
    missing = sorted(set(selected_ids).difference(available_ids))
    if missing:
        raise ValueError(
            f"requested assembly IDs are absent from the mmCIF: {missing}; "
            f"available={available_ids}"
        )
    return selected_ids


def _query_chain_atoms(
    cif_file: Any,
    *,
    structure_mode: str,
    assembly_ids: Iterable[str] | None,
) -> dict[str, Any]:
    """Return one coordinate representative per label asym protein chain."""
    import biotite.structure as struc
    from biotite.structure import filter_heavy

    from plinder.data.annotations.cif_utils import (
        build_biounit,
        get_structure_with_altloc,
    )

    if structure_mode == "as_is":
        if assembly_ids is not None:
            raise ValueError("assembly_ids are only valid in pdb mode")
        atoms = get_structure_with_altloc(
            cif_file,
            model=1,
            use_author_fields=False,
            include_bonds=False,
        )
        atoms = atoms[filter_heavy(atoms)]
        return {
            str(chain_id): atoms[atoms.chain_id == chain_id]
            for chain_id in np.unique(atoms.chain_id)
        }
    if structure_mode != "pdb":
        raise ValueError("structure_mode must be 'as_is' or 'pdb'")

    block = list(cif_file.values())[0]
    available_assemblies = block["pdbx_struct_assembly_gen"]["assembly_id"].as_array(
        str
    )
    selected = _select_assembly_ids(available_assemblies, assembly_ids)
    if not selected:
        raise ValueError("pdb mode requires at least one deposited assembly")
    representatives: dict[str, Any] = {}
    for assembly_id in selected:
        assembly = build_biounit(cif_file, assembly_id)
        for instance_chain in np.unique(assembly.chain_id):
            chain = str(instance_chain)
            if "." not in chain:
                continue
            asym_id = chain.split(".", maxsplit=1)[1]
            representatives.setdefault(
                asym_id,
                assembly[assembly.chain_id == instance_chain],
            )
    # Ensure every representative contains at least one residue. This also
    # narrows the otherwise broad Any returned by Biotite's dynamic CIF API.
    return {
        asym_id: atoms
        for asym_id, atoms in representatives.items()
        if struc.get_residue_count(atoms) > 0
    }


def _resolved_protein_sequence(atoms: Any, *, context: str) -> str | None:
    """Derive a sequence for coordinate-only mmCIFs without entity tables."""
    import biotite.structure as struc
    from biotite.sequence import ProteinSequence

    try:
        sequences, _ = struc.to_sequence(atoms, allow_hetero=True)
    except (IndexError, TypeError, ValueError, struc.BadStructureError):
        return None
    if len(sequences) != 1:
        return None
    sequence = sequences[0]
    if not isinstance(sequence, ProteinSequence):
        return None
    return _clean_protein_sequence(str(sequence), context=context)


def write_custom_query_files(
    cif_files: Iterable[Path],
    *,
    work_dir: Path,
    structure_mode: str = "as_is",
    assembly_ids: Iterable[str] | None = None,
    min_chain_length: int = 12,
) -> CustomQueryInputs:
    """Write one Foldseek-ready CIF and one FASTA record per protein chain.

    ``as_is`` treats each supplied coordinate file as an already assembled
    structure. ``pdb`` applies the selected deposited assembly definitions,
    but searches only one rigid-transform-equivalent copy per label asym ID.
    Nucleic-acid entities and protein chains shorter than ``min_chain_length``
    are excluded from both search backends.
    """
    import biotite.structure as struc
    from biotite.file import DeserializationError, InvalidFileError

    from plinder.data.annotations.cif_utils import (
        check_custom_mmcif_fields,
        read_mmcif_file,
    )
    from plinder.data.annotations.save_utils import save_cif_file

    if min_chain_length < 1:
        raise ValueError("min_chain_length must be positive")
    sources = tuple(sorted((Path(path) for path in cif_files), key=str))
    if not sources:
        raise ValueError("at least one custom mmCIF is required")
    requested_assemblies = (
        tuple(str(value) for value in assembly_ids)
        if assembly_ids is not None and not isinstance(assembly_ids, str)
        else assembly_ids
    )
    root = Path(work_dir) / "query_inputs"
    resolved_root = root.resolve()
    for source in sources:
        _mmcif_stem(source)
        if not source.is_file():
            raise FileNotFoundError(f"missing custom mmCIF: {source}")
        if source.resolve().is_relative_to(resolved_root):
            raise ValueError(
                "custom input files must be outside the generated work tree"
            )
    if root.exists():
        shutil.rmtree(root)
    chain_cif_dir = root / "chains"
    chain_cif_dir.mkdir(parents=True)
    sequence_fasta = root / "query_sequences.fasta"
    chain_manifest = root / "query_chains.parquet"

    structure_ids: set[str] = set()
    rows: list[dict[str, Any]] = []
    fasta_records: list[str] = []
    skipped_chains: dict[str, str] = {}
    for source in sources:
        structure_id = _mmcif_stem(source)
        if structure_id in structure_ids:
            raise ValueError(f"duplicate custom structure ID: {structure_id}")
        structure_ids.add(structure_id)
        try:
            cif_file = read_mmcif_file(source)
            block = list(cif_file.values())[0]
        except (DeserializationError, IndexError, InvalidFileError, OSError) as exc:
            raise ValueError(f"cannot parse custom mmCIF {source}: {exc}") from exc
        check_custom_mmcif_fields(
            block,
            source=source,
            structure_mode=structure_mode,
        )
        sequences = _protein_asym_sequences(block)
        try:
            atoms_by_asym = _query_chain_atoms(
                cif_file,
                structure_mode=structure_mode,
                assembly_ids=requested_assemblies,
            )
        except (
            IndexError,
            InvalidFileError,
            KeyError,
            TypeError,
            ValueError,
        ) as exc:
            raise ValueError(
                f"failed to read {structure_mode} coordinates from custom "
                f"mmCIF {source}: {exc}"
            ) from exc
        for asym_id in sorted(atoms_by_asym):
            sequence = sequences.get(asym_id)
            sequence_source = "polymer"
            if sequence is None:
                sequence = _resolved_protein_sequence(
                    atoms_by_asym[asym_id],
                    context=f"{structure_id} chain {asym_id}",
                )
                sequence_source = "coordinates"
            if sequence is None:
                skipped_chains[f"{structure_id}__{asym_id}"] = (
                    "not identifiable as a protein"
                )
                continue
            if len(sequence) < min_chain_length:
                skipped_chains[f"{structure_id}__{asym_id}"] = (
                    f"protein length {len(sequence)} is below {min_chain_length}"
                )
                continue
            atoms = atoms_by_asym[asym_id].copy()
            residue_starts = struc.get_residue_starts(atoms, add_exclusive_stop=False)
            if len(residue_starts) < 1:
                continue
            query_id = f"cq{len(rows):08d}"
            query_chain_id = f"{structure_id}__{asym_id}"
            resolved_residue_numbers = [
                int(atoms.res_id[index]) for index in residue_starts
            ]
            if sequence_source == "coordinates":
                if len(resolved_residue_numbers) != len(sequence):
                    raise ValueError(
                        "resolved residue count does not match the emitted FASTA "
                        f"sequence for {structure_id} chain {asym_id}"
                    )
                if len(set(resolved_residue_numbers)) != len(resolved_residue_numbers):
                    raise ValueError(
                        f"resolved residue numbers are ambiguous for {structure_id} "
                        f"chain {asym_id}; coordinate-derived FASTA positions cannot "
                        "be mapped uniquely"
                    )
            atoms.chain_id = np.full(len(atoms), "A")
            chain_path = chain_cif_dir / f"{query_id}.cif"
            save_cif_file(
                atoms,
                query_id,
                chain_path,
                source_block=block,
                source_asym_ids={"A": asym_id},
                protein_sequences={"A": sequence},
            )
            rows.append({
                "query_id": query_id,
                "query_chain_id": query_chain_id,
                "structure_id": structure_id,
                "source_mmcif": str(source.resolve()),
                "chain_asym_id": asym_id,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "sequence_source": sequence_source,
                "resolved_residue_numbers": resolved_residue_numbers,
            })
            fasta_records.append(f">{query_id}\n{sequence}\n")
    if not rows:
        raise ValueError(
            f"no protein chains of at least {min_chain_length} residues were found; "
            f"chain diagnostics={skipped_chains}"
        )
    sequence_fasta.write_text("".join(fasta_records))
    pd.DataFrame(rows).to_parquet(chain_manifest, index=False)
    return CustomQueryInputs(
        root=root,
        chain_cif_dir=chain_cif_dir,
        sequence_fasta=sequence_fasta,
        chain_manifest=chain_manifest,
    )


def annotate_custom_cif_files(
    cif_files: Iterable[Path],
    *,
    work_dir: Path,
    structure_mode: Literal["as_is", "pdb"] = "as_is",
    assembly_ids: Iterable[str] | None = None,
    ligand_smiles_dict: Mapping[str, str] | None = None,
    ligand_ccd_code_dict: Mapping[str, str] | None = None,
    include_ligands: bool = True,
    include_interfaces: bool = True,
    interface_annotate_prodigy: bool = True,
    data_dir: Path | None = None,
) -> CustomStructureAnnotations:
    """Annotate custom structures and write the scorer's three index tables.

    Ligand chemistry overrides are passed directly to
    :meth:`Entry.from_custom_cif_file`. The same component mapping is applied
    to every input file; callers with different meanings for a repeated
    component name such as ``LIG`` should score those files in separate calls.
    Generated SDFs live under ``ligand_sdf_root/<entry>/ligand_files`` and
    contain the bond orders assigned during custom-CIF ingest.
    """
    from plinder.core.scores.entries import entry_views_from_df
    from plinder.data.annotations.aggregate_annotations import Entry

    sources = tuple(sorted((Path(path) for path in cif_files), key=str))
    if not sources:
        raise ValueError("at least one custom mmCIF is required")
    structure_ids = [_mmcif_stem(source) for source in sources]
    duplicates = sorted(
        structure_id
        for structure_id in set(structure_ids)
        if structure_ids.count(structure_id) > 1
    )
    if duplicates:
        raise ValueError(f"duplicate custom structure IDs: {duplicates}")

    root = Path(work_dir) / "custom_annotations"
    root.mkdir(exist_ok=True, parents=True)
    ligand_sdf_root = root / "ligands"
    ligand_sdf_root.mkdir(exist_ok=True, parents=True)
    requested_assemblies = (
        tuple(str(value) for value in assembly_ids)
        if assembly_ids is not None and not isinstance(assembly_ids, str)
        else assembly_ids
    )

    annotated: dict[str, Any] = {}
    entry_ids: set[str] = set()
    for structure_id, source in zip(structure_ids, sources, strict=True):
        entry = Entry.from_custom_cif_file(
            None if structure_mode == "pdb" else structure_id,
            source,
            ligand_smiles_dict=(
                dict(ligand_smiles_dict) if ligand_smiles_dict is not None else None
            ),
            ligand_ccd_code_dict=(
                dict(ligand_ccd_code_dict) if ligand_ccd_code_dict is not None else None
            ),
            save_folder=ligand_sdf_root if include_ligands else None,
            structure_mode=structure_mode,
            assembly_ids=requested_assemblies,
            include_ligands=include_ligands,
            include_interfaces=include_interfaces,
            interface_annotate_prodigy=interface_annotate_prodigy,
            data_dir=data_dir,
        )
        entry_id = str(entry.pdb_id)
        if entry_id in entry_ids:
            raise ValueError(
                f"custom inputs produce duplicate entry ID {entry_id!r}; "
                "score them in separate calls"
            )
        entry_ids.add(entry_id)
        protein_chains = [
            chain
            for chain in entry.chains.values()
            if "polypeptide" in str(chain.chain_type_str).lower()
        ]
        if not protein_chains:
            raise ValueError(
                f"custom mmCIF {source} produced no protein chains to score"
            )
        annotated[structure_id] = entry

    ligand_frames = [entry.to_df() for entry in annotated.values()]
    ligand_frames = [frame for frame in ligand_frames if not frame.empty]
    annotation = (
        pd.concat(ligand_frames, ignore_index=True)
        if ligand_frames
        else pd.DataFrame({
            "entry_pdb_id": pd.Series(dtype="string"),
            "system_id": pd.Series(dtype="string"),
        })
    )
    entry_chains = pd.concat(
        [entry.chains_to_df() for entry in annotated.values()],
        ignore_index=True,
    )
    interface_rows = [
        interface.to_row()
        for entry in annotated.values()
        for interface in entry.interfaces
    ]
    interface_annotations = pd.DataFrame.from_records(interface_rows)
    if interface_annotations.empty:
        interface_annotations = pd.DataFrame({
            "entry_pdb_id": pd.Series(dtype="string")
        })

    annotation_table = root / "annotation_table.parquet"
    entry_chain_table = root / "entry_chains.parquet"
    interface_table = root / "interface_annotation_table.parquet"
    annotation.to_parquet(annotation_table, index=False)
    entry_chains.to_parquet(entry_chain_table, index=False)
    interface_annotations.to_parquet(interface_table, index=False)
    entry_views = entry_views_from_df(
        annotation,
        entry_chains=entry_chains,
        interface_annotations=interface_annotations,
    )
    entries_by_structure = {
        structure_id: entry_views[str(entry.pdb_id)]
        for structure_id, entry in annotated.items()
    }
    return CustomStructureAnnotations(
        root=root,
        ligand_sdf_root=ligand_sdf_root,
        annotation_table=annotation_table,
        entry_chains=entry_chain_table,
        interface_annotations=interface_table,
        entries_by_structure=entries_by_structure,
    )


def _run_command(command: list[str]) -> None:
    LOG.info("running custom search command: %s", " ".join(command))
    subprocess.check_call(command)


def _database_lookup_identifiers(database: Path) -> list[str]:
    lookup = database.with_suffix(".lookup")
    if not lookup.is_file():
        raise FileNotFoundError(f"search query database has no lookup: {lookup}")
    identifiers: list[str] = []
    for line in lookup.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            identifiers.append(fields[1])
    return identifiers


def _map_backend_query_identifiers(
    identifiers: Iterable[str], query_ids: Iterable[str], *, backend: str
) -> dict[str, str]:
    expected = tuple(query_ids)
    output: dict[str, str] = {}
    matched: dict[str, list[str]] = {query_id: [] for query_id in expected}
    for identifier in identifiers:
        candidates = [
            query_id
            for query_id in expected
            if identifier == query_id
            or identifier.startswith(f"{query_id}_")
            or re.search(rf"(?:^|/){re.escape(query_id)}(?:[_.]|$)", identifier)
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"cannot map {backend} query identifier {identifier!r}: {candidates}"
            )
        query_id = candidates[0]
        output[identifier] = query_id
        matched[query_id].append(identifier)
    invalid = {
        query_id: values for query_id, values in matched.items() if len(values) != 1
    }
    if invalid:
        raise ValueError(
            f"{backend} query database must contain one record per chain: {invalid}"
        )
    return output


def create_custom_query_databases(
    inputs: CustomQueryInputs,
    *,
    work_dir: Path,
    backends: Iterable[str] = CIF_SEARCH_BACKENDS,
    threads: int = 1,
) -> CustomQueryDatabases:
    """Create unindexed, scratch-local query DBs for selected search backends."""
    if threads < 1:
        raise ValueError("threads must be positive")
    selected = tuple(dict.fromkeys(backends))
    if not selected:
        raise ValueError("at least one custom search backend is required")
    unsupported = sorted(set(selected).difference(SEARCH_BACKENDS))
    if unsupported:
        raise ValueError(f"unsupported search backends: {unsupported}")
    for backend in selected:
        if shutil.which(backend) is None:
            raise FileNotFoundError(
                f"{backend} executable is required for custom protein searches"
            )
    root = Path(work_dir) / "query_databases"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    databases = {backend: root / backend for backend in selected}
    if "foldseek" in selected:
        _run_command([
            "foldseek",
            "createdb",
            str(inputs.chain_cif_dir),
            str(databases["foldseek"]),
            "--threads",
            str(threads),
            "--chain-name-mode",
            "1",
            "--coord-store-mode",
            "2",
        ])
    if "mmseqs" in selected:
        _run_command([
            "mmseqs",
            "createdb",
            str(inputs.sequence_fasta),
            str(databases["mmseqs"]),
            "--threads",
            str(threads),
        ])
    manifest = pd.read_parquet(inputs.chain_manifest)
    query_ids = manifest["query_id"].astype(str).tolist()
    mapping_rows: list[dict[str, str]] = []
    for backend, database in databases.items():
        mapping = _map_backend_query_identifiers(
            _database_lookup_identifiers(database),
            query_ids,
            backend=backend,
        )
        mapping_rows.extend(
            {
                "backend": backend,
                "backend_query_id": identifier,
                "query_id": query_id,
            }
            for identifier, query_id in mapping.items()
        )
    identifier_map = root / "query_identifier_map.parquet"
    pd.DataFrame(mapping_rows).sort_values(["backend", "query_id"]).to_parquet(
        identifier_map, index=False
    )
    return CustomQueryDatabases(
        root=root,
        inputs=inputs,
        databases=databases,
        identifier_map=identifier_map,
    )


def _plinder_entry_subset(values: Iterable[str] | None) -> tuple[str, ...] | None:
    if values is None:
        return None
    if isinstance(values, str):
        values = (values,)
    selected = tuple(
        dict.fromkeys(
            str(value).strip().lower() for value in values if str(value).strip()
        )
    )
    if not selected:
        raise ValueError("plinder_entry_ids must contain at least one PDB ID")
    invalid = [
        value for value in selected if re.fullmatch(r"[a-z0-9]{4}", value) is None
    ]
    if invalid:
        raise ValueError(f"invalid PLINDER PDB IDs: {invalid[:10]}")
    return selected


def _build_mmseqs_target_subset(
    entry_chains: Path,
    interface_annotations: Path,
    *,
    entry_ids: Iterable[str],
    output_dir: Path,
    threads: int = 1,
) -> SearchDatabaseBundle:
    """Build an MMseqs target from scoreable chains in selected entries."""
    selected = _plinder_entry_subset(entry_ids)
    assert selected is not None
    if threads < 1:
        raise ValueError("threads must be positive")
    if shutil.which("mmseqs") is None:
        raise FileNotFoundError("mmseqs executable is required for custom searches")

    chains = pd.read_parquet(
        entry_chains,
        columns=[
            "entry_pdb_id",
            "chain_asym_id",
            "chain_auth_id",
            "chain_receptor_type",
            "chain_is_holo",
            "chain_sequence",
        ],
        filters=[("entry_pdb_id", "in", list(selected))],
    )
    interfaces = pd.read_parquet(
        interface_annotations,
        columns=["entry_pdb_id", "interface_chain_1", "interface_chain_2"],
        filters=[("entry_pdb_id", "in", list(selected))],
    )
    interface_keys = pd.concat(
        [
            interfaces[["entry_pdb_id", column]].rename(
                columns={column: "chain_instance"}
            )
            for column in ["interface_chain_1", "interface_chain_2"]
        ],
        ignore_index=True,
    )
    interface_keys["chain_asym_id"] = (
        interface_keys.pop("chain_instance").astype(str).str.split(".", n=1).str[-1]
    )
    interface_keys = interface_keys.drop_duplicates()
    interface_keys["chain_is_interface"] = True
    chains = chains.merge(
        interface_keys,
        on=["entry_pdb_id", "chain_asym_id"],
        how="left",
        validate="one_to_one",
    )
    chains = chains.loc[
        chains["chain_receptor_type"].fillna("").astype(str).eq("protein")
        & (
            chains["chain_is_holo"].fillna(False).astype(bool)
            | chains["chain_is_interface"].eq(True)
        )
        & chains["chain_auth_id"].notna()
        & chains["chain_sequence"].notna()
    ].copy()
    found = set(chains["entry_pdb_id"].astype(str))
    missing = sorted(set(selected).difference(found))
    if missing:
        raise ValueError(
            "selected PLINDER entries have no scoreable receptor or interface "
            f"chains: {missing[:10]}"
        )

    chains["target_id"] = (
        chains["entry_pdb_id"].astype(str) + "_" + chains["chain_auth_id"].astype(str)
    )
    duplicate_ids = chains.loc[
        chains["target_id"].duplicated(keep=False), "target_id"
    ].drop_duplicates()
    if not duplicate_ids.empty:
        raise ValueError(
            "selected PLINDER chains have ambiguous author IDs: "
            f"{duplicate_ids.tolist()[:10]}"
        )

    records: list[str] = []
    for row in chains.sort_values("target_id").itertuples(index=False):
        if any(character.isspace() for character in str(row.target_id)):
            raise ValueError(f"invalid PLINDER target identifier: {row.target_id!r}")
        sequence = _clean_protein_sequence(
            str(row.chain_sequence),
            context=f"{row.entry_pdb_id} chain {row.chain_auth_id}",
        )
        records.append(f">{row.target_id}\n{sequence}\n")

    root = Path(output_dir) / "mmseqs"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(exist_ok=True, parents=True)
    fasta = root / "targets.fasta"
    fasta.write_text("".join(records))
    database = root / "targets"
    _run_command([
        "mmseqs",
        "createdb",
        str(fasta),
        str(database),
        "--threads",
        str(threads),
    ])
    return SearchDatabaseBundle(
        backend="mmseqs",
        root=root,
        search_target=database,
        conversion_target=database,
        cluster_alignments=None,
        manifest={"plinder_entry_ids": list(selected)},
    )


def _resolve_workflow_assets(
    *,
    data_dir: Path | None,
    backends: tuple[str, ...],
    plinder_entry_ids: tuple[str, ...] | None,
    work_dir: Path,
    threads: int,
) -> CustomScoringAssets:
    if plinder_entry_ids is not None and backends != ("mmseqs",):
        raise ValueError("plinder_entry_ids currently requires backends=('mmseqs',)")
    assets = resolve_custom_scoring_assets(
        data_dir=data_dir,
        backends=() if plinder_entry_ids is not None else backends,
    )
    if plinder_entry_ids is None:
        return assets
    subset = _build_mmseqs_target_subset(
        assets.entry_chains,
        assets.interface_annotations,
        entry_ids=plinder_entry_ids,
        output_dir=work_dir / "plinder_target_subset",
        threads=threads,
    )
    return replace(assets, search_databases={"mmseqs": subset})


def _parse_plinder_target_identifier(
    identifier: str, *, backend: str
) -> tuple[str, str]:
    cleaned = str(identifier)
    for value in (
        "_xyz-enrich.cif.gz",
        "_xyz-enrich.cif",
        "_xyz-enrich",
        "pdb_0000",
        ".cif.gz",
        ".cif",
    ):
        cleaned = cleaned.replace(value, "")
    entry_id, separator, author_chain = cleaned.partition("_")
    if not separator or not entry_id or not author_chain:
        raise ValueError(f"invalid {backend} target identifier: {identifier!r}")
    if backend == "foldseek":
        author_chain = re.sub(r"^MODEL_\d+_", "", author_chain)
    return entry_id.lower(), author_chain


def _load_target_chain_mapping(
    lookup_path: Path, target_identifiers: Iterable[str], *, backend: str
) -> pd.DataFrame:
    parsed = pd.DataFrame([
        {
            "target_backend_id": identifier,
            "target_entry": entry_id,
            "target_chain_auth_id": auth_id,
        }
        for identifier in sorted(set(target_identifiers))
        for entry_id, auth_id in [
            _parse_plinder_target_identifier(identifier, backend=backend)
        ]
    ])
    if parsed.empty:
        return parsed.assign(
            target_chain_asym_id=pd.Series(dtype="string"),
            target_selected_residue_numbers=pd.Series(dtype="object"),
            target_selected_residue_indices=pd.Series(dtype="object"),
        )
    entry_ids = sorted(set(parsed["target_entry"]))
    lookup = pd.read_parquet(
        lookup_path,
        columns=[
            "entry_pdb_id",
            "chain_asym_id",
            "chain_auth_id",
            "selected_residue_numbers",
            "selected_residue_indices",
        ],
        filters=[("entry_pdb_id", "in", entry_ids)],
    ).rename(
        columns={
            "entry_pdb_id": "target_entry",
            "chain_asym_id": "target_chain_asym_id",
            "chain_auth_id": "target_chain_auth_id",
            "selected_residue_numbers": "target_selected_residue_numbers",
            "selected_residue_indices": "target_selected_residue_indices",
        }
    )
    key_columns = ["target_entry", "target_chain_auth_id"]
    ambiguous = (
        lookup
        .groupby(key_columns, dropna=False)["target_chain_asym_id"]
        .nunique()
        .loc[lambda values: values > 1]
    )
    if not ambiguous.empty:
        raise ValueError(
            "alignment chain lookup has ambiguous author-chain mappings: "
            f"{list(ambiguous.index[:10])}"
        )
    lookup = lookup.drop_duplicates(key_columns)
    mapped = parsed.merge(lookup, on=key_columns, how="left", validate="many_to_one")
    missing = mapped.loc[
        mapped["target_chain_asym_id"].isna(), "target_backend_id"
    ].tolist()
    if missing:
        raise ValueError(
            f"PLINDER alignment lookup cannot map target chains: {missing[:10]}"
        )
    return mapped


def map_custom_alignment_hits(
    *,
    raw_alignment: Path,
    backend: str,
    query_databases: CustomQueryDatabases,
    alignment_chain_lookup: Path,
    output_path: Path,
) -> Path:
    """Replace backend identifiers with custom and PLINDER label-asym IDs."""
    if backend not in SEARCH_BACKENDS:
        raise ValueError(f"unsupported search backend: {backend}")
    raw = pd.read_parquet(raw_alignment)
    identifier_map = pd.read_parquet(
        query_databases.identifier_map,
        filters=[("backend", "=", backend)],
    )
    chains = pd.read_parquet(query_databases.inputs.chain_manifest)
    query_mapping = identifier_map.merge(
        chains,
        on="query_id",
        how="left",
        validate="many_to_one",
    )
    output_path.parent.mkdir(exist_ok=True, parents=True)
    if raw.empty:
        columns = [
            "query_id",
            "query_chain_id",
            "structure_id",
            "query_chain_asym_id",
            "target_entry",
            "target_chain_asym_id",
            "source",
        ]
        pd.DataFrame({
            column: pd.Series(dtype="string") for column in columns
        }).to_parquet(output_path, index=False)
        return output_path
    query_columns = query_mapping[
        [
            "backend_query_id",
            "query_id",
            "query_chain_id",
            "structure_id",
            "chain_asym_id",
            "sequence_source",
            "resolved_residue_numbers",
        ]
    ].rename(
        columns={
            "backend_query_id": "query_backend_id",
            "chain_asym_id": "query_chain_asym_id",
            "sequence_source": "query_sequence_source",
            "resolved_residue_numbers": "query_resolved_residue_numbers",
        }
    )
    raw = raw.rename(
        columns={"query": "query_backend_id", "target": "target_backend_id"}
    )
    mapped = raw.merge(
        query_columns,
        on="query_backend_id",
        how="left",
        validate="many_to_one",
    )
    missing_queries = mapped.loc[
        mapped["query_id"].isna(), "query_backend_id"
    ].drop_duplicates()
    if not missing_queries.empty:
        raise ValueError(
            f"cannot map {backend} query chains: {missing_queries.tolist()[:10]}"
        )
    targets = _load_target_chain_mapping(
        alignment_chain_lookup,
        mapped["target_backend_id"],
        backend=backend,
    )
    mapped = mapped.merge(
        targets,
        on="target_backend_id",
        how="left",
        validate="many_to_one",
    )
    mapped["source"] = backend
    mapped = mapped.drop(columns=["query_pdb_id", "target_pdb_id"], errors="ignore")
    leading = [
        "query_id",
        "query_chain_id",
        "structure_id",
        "query_chain_asym_id",
        "target_entry",
        "target_chain_asym_id",
        "source",
    ]
    remainder = [column for column in mapped.columns if column not in leading]
    mapped = mapped[leading + remainder].sort_values(
        ["query_id", "target_entry", "target_chain_asym_id"],
        ignore_index=True,
    )
    mapped.to_parquet(output_path, index=False, compression="zstd")
    return output_path


def _int_values(value: object, *, column: str) -> list[int]:
    if value is None or isinstance(value, (str, bytes)):
        return []
    if not isinstance(value, IterableABC):
        if pd.isna(value):
            return []
        raise ValueError(f"custom alignment column {column} is not a list")
    return [int(item) for item in value]


def _selected_positions_for_custom_hit(
    row: Any,
    *,
    backend: str,
    query_entry: Any,
) -> tuple[list[int], list[int], bytes]:
    """Map selected custom residues through one raw pairwise alignment."""
    qaln = str(row.qaln).upper()
    taln = str(row.taln).upper()
    alignment_length = min(len(qaln), len(taln))
    qaln = qaln[:alignment_length]
    taln = taln[:alignment_length]
    query_chain = str(row.query_chain_asym_id)
    query_index_to_number = query_entry.selected_index_to_number_per_chain.get(
        query_chain, {}
    )
    if not query_index_to_number or alignment_length == 0:
        return [], [], b""

    target_numbers = _int_values(
        row.target_selected_residue_numbers,
        column="target_selected_residue_numbers",
    )
    target_indices = _int_values(
        row.target_selected_residue_indices,
        column="target_selected_residue_indices",
    )
    if len(target_numbers) != len(target_indices):
        raise ValueError(
            "target selected-residue numbers and indices have different lengths "
            f"for {row.target_entry} chain {row.target_chain_asym_id}"
        )
    target_index_to_number = dict(zip(target_indices, target_numbers, strict=True))

    query_start = int(row.qstart) - 1
    target_start = int(row.tstart) - 1
    query_alignment_positions = [
        position for position, residue in enumerate(qaln) if residue != "-"
    ]
    target_alignment_positions = [
        position for position, residue in enumerate(taln) if residue != "-"
    ]
    if backend in SEQUENCE_SEARCH_BACKENDS:
        query_numbers = set(query_index_to_number.values())
        if str(row.query_sequence_source) == "coordinates":
            resolved_numbers = _int_values(
                row.query_resolved_residue_numbers,
                column="query_resolved_residue_numbers",
            )
            if len(resolved_numbers) != len(set(resolved_numbers)):
                raise ValueError(
                    "coordinate-derived FASTA has ambiguous residue numbers for "
                    f"{row.structure_id} chain {query_chain}"
                )
            residue_number_to_fasta_position = {
                number: position
                for position, number in enumerate(resolved_numbers, start=1)
            }
            missing_numbers = sorted(
                query_numbers.difference(residue_number_to_fasta_position)
            )
            if missing_numbers:
                raise ValueError(
                    "coordinate-derived FASTA is missing selected residues for "
                    f"{row.structure_id} chain {query_chain}: {missing_numbers}"
                )
            query_candidates = [
                (
                    residue_number_to_fasta_position[number] - (query_start + 1),
                    number,
                )
                for number in query_numbers
            ]
        elif str(row.query_sequence_source) == "polymer":
            query_candidates = [
                (number - (query_start + 1), number) for number in query_numbers
            ]
        else:
            raise ValueError(
                "unknown query FASTA sequence source for "
                f"{row.structure_id} chain {query_chain}: "
                f"{row.query_sequence_source!r}"
            )
        target_selected_numbers = set(target_index_to_number.values())
    else:
        query_candidates = [
            (index - query_start, number)
            for index, number in query_index_to_number.items()
        ]
        target_selected_numbers = set()

    aligned_query_numbers: list[int] = []
    aligned_target_numbers: list[int] = []
    residue_identity = bytearray()
    for query_offset, query_number in sorted(query_candidates):
        if not 0 <= query_offset < len(query_alignment_positions):
            continue
        alignment_position = query_alignment_positions[query_offset]
        if taln[alignment_position] == "-":
            continue
        target_index = target_start + bisect_left(
            target_alignment_positions, alignment_position
        )
        if backend in SEQUENCE_SEARCH_BACKENDS:
            target_number = (
                target_index + 1
                if target_index + 1 in target_selected_numbers
                else None
            )
        else:
            target_number = target_index_to_number.get(target_index)
        aligned_query_numbers.append(int(query_number))
        aligned_target_numbers.append(
            int(target_number) if target_number is not None else -1
        )
        residue_identity.append(qaln[alignment_position] == taln[alignment_position])
    return aligned_query_numbers, aligned_target_numbers, bytes(residue_identity)


def prepare_custom_score_alignments(
    protein_hits: Mapping[str, Path],
    *,
    entries_by_structure: Mapping[str, Any],
    output_dir: Path,
) -> dict[str, Path]:
    """Write hit tables in the compact alignment format used by PLINDER scores."""
    from plinder.data.annotations.get_similarity_scores import (
        get_sequence_similarity_helper,
    )

    unsupported = sorted(set(protein_hits).difference(SEARCH_BACKENDS))
    if unsupported:
        raise ValueError(f"unsupported search backends: {unsupported}")
    query_entry_ids = [str(entry.pdb_id) for entry in entries_by_structure.values()]
    if len(query_entry_ids) != len(set(query_entry_ids)):
        raise ValueError("custom scoring entries must have distinct entry IDs")
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    outputs: dict[str, Path] = {}
    required = {
        "structure_id",
        "query_chain_asym_id",
        "target_entry",
        "target_chain_asym_id",
        "target_selected_residue_numbers",
        "target_selected_residue_indices",
        "qstart",
        "tstart",
        "qcov",
        "fident",
        "qaln",
        "taln",
    }
    for backend, hit_path in protein_hits.items():
        hits = pd.read_parquet(hit_path)
        output = output_dir / f"{backend}.parquet"
        if hits.empty:
            pd.DataFrame({
                column: pd.Series(dtype="string")
                for column in (
                    "query_entry",
                    "target_entry",
                    "query_chain_mapped",
                    "target_chain_mapped",
                    "source",
                )
            }).to_parquet(output, index=False)
            outputs[backend] = output
            continue
        backend_required = required | (
            {"query_sequence_source", "query_resolved_residue_numbers"}
            if backend in SEQUENCE_SEARCH_BACKENDS
            else set()
        )
        missing = sorted(backend_required.difference(hits.columns))
        if missing:
            raise ValueError(f"{backend} custom hit table is missing columns {missing}")
        unknown_structures = sorted(
            set(hits["structure_id"].astype(str)).difference(entries_by_structure)
        )
        if unknown_structures:
            raise ValueError(
                f"{backend} hits reference unknown custom structures: "
                f"{unknown_structures}"
            )
        hits = hits.loc[~hits["target_entry"].astype(str).isin(query_entry_ids)].copy()
        if hits.empty:
            pd.DataFrame({
                column: pd.Series(dtype="string")
                for column in (
                    "query_entry",
                    "target_entry",
                    "query_chain_mapped",
                    "target_chain_mapped",
                    "source",
                )
            }).to_parquet(output, index=False)
            outputs[backend] = output
            continue
        selected = [
            _selected_positions_for_custom_hit(
                row,
                backend=backend,
                query_entry=entries_by_structure[str(row.structure_id)],
            )
            for row in hits.itertuples(index=False)
        ]
        hits["query_selected_residue_numbers"] = [value[0] for value in selected]
        hits["target_selected_residue_numbers"] = [value[1] for value in selected]
        hits["selected_residue_identity"] = [value[2] for value in selected]
        hits["query_entry"] = [
            str(entries_by_structure[str(value)].pdb_id)
            for value in hits["structure_id"]
        ]
        hits["query_chain_mapped"] = hits["query_chain_asym_id"].astype(str)
        hits["target_chain_mapped"] = hits["target_chain_asym_id"].astype(str)
        hits["source"] = backend
        hits["seqsim"] = [
            get_sequence_similarity_helper(str(qaln).upper(), str(taln).upper())
            for qaln, taln in zip(hits["qaln"], hits["taln"], strict=True)
        ]
        hits["fident_qcov"] = hits["fident"] * hits["qcov"]
        hits["seqsim_qcov"] = hits["seqsim"] * hits["qcov"]
        if backend == "foldseek":
            if "lddt" not in hits:
                raise ValueError("foldseek custom hit table is missing column lddt")
            hits["lddt_qcov"] = hits["lddt"] * hits["qcov"]
        hits = hits.drop(
            columns=[
                "qaln",
                "taln",
                "target_selected_residue_indices",
            ],
            errors="ignore",
        )
        leading = [
            "query_entry",
            "target_entry",
            "query_chain_mapped",
            "target_chain_mapped",
            "source",
        ]
        hits = hits[leading + [column for column in hits if column not in leading]]
        hits = hits.sort_values(
            [
                "query_entry",
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
            ],
            ignore_index=True,
        )
        hits.to_parquet(output, index=False, compression="zstd")
        outputs[backend] = output
    return outputs


def _release_positions_for_custom_hit(
    row: Any,
    *,
    backend: str,
) -> tuple[list[int], list[int], bytes]:
    """Project PLINDER selected residues through a custom-chain alignment."""
    qaln = str(row.qaln).upper()
    taln = str(row.taln).upper()
    alignment_length = min(len(qaln), len(taln))
    qaln = qaln[:alignment_length]
    taln = taln[:alignment_length]
    if alignment_length == 0:
        return [], [], b""

    selected_numbers = _int_values(
        row.target_selected_residue_numbers,
        column="target_selected_residue_numbers",
    )
    selected_indices = _int_values(
        row.target_selected_residue_indices,
        column="target_selected_residue_indices",
    )
    if len(selected_numbers) != len(selected_indices):
        raise ValueError(
            "target selected-residue numbers and indices have different lengths "
            f"for {row.target_entry} chain {row.target_chain_asym_id}"
        )
    if not selected_numbers:
        return [], [], b""

    target_start = int(row.tstart) - 1
    query_start = int(row.qstart) - 1
    query_alignment_positions = [
        position for position, residue in enumerate(qaln) if residue != "-"
    ]
    resolved_numbers = _int_values(
        row.query_resolved_residue_numbers,
        column="query_resolved_residue_numbers",
    )
    target_alignment_positions = [
        position for position, residue in enumerate(taln) if residue != "-"
    ]
    if backend in SEQUENCE_SEARCH_BACKENDS:
        selected = [
            (number - (target_start + 1), number) for number in selected_numbers
        ]
    else:
        selected = [
            (index - target_start, number)
            for index, number in zip(
                selected_indices,
                selected_numbers,
                strict=True,
            )
        ]

    query_numbers: list[int] = []
    custom_numbers: list[int] = []
    residue_identity = bytearray()
    for target_offset, target_number in sorted(selected):
        if not 0 <= target_offset < len(target_alignment_positions):
            continue
        alignment_position = target_alignment_positions[target_offset]
        if qaln[alignment_position] == "-":
            continue
        query_index = query_start + bisect_left(
            query_alignment_positions, alignment_position
        )
        if backend == "foldseek" or str(row.query_sequence_source) == "coordinates":
            custom_number = (
                resolved_numbers[query_index]
                if 0 <= query_index < len(resolved_numbers)
                else None
            )
        elif str(row.query_sequence_source) == "polymer":
            custom_number = query_index + 1
        else:
            raise ValueError(
                "unknown query sequence source for "
                f"{row.structure_id} chain {row.query_chain_asym_id}: "
                f"{row.query_sequence_source!r}"
            )
        query_numbers.append(int(target_number))
        custom_numbers.append(int(custom_number) if custom_number is not None else -1)
        residue_identity.append(qaln[alignment_position] == taln[alignment_position])
    return query_numbers, custom_numbers, bytes(residue_identity)


def prepare_custom_protein_score_alignments(
    protein_hits: Mapping[str, Path],
    *,
    entries_by_structure: Mapping[str, Any],
    output_dir: Path,
) -> dict[str, Path]:
    """Write reverse compact alignments for PLINDER-to-custom protein scores."""
    from plinder.data.annotations.get_similarity_scores import (
        get_sequence_similarity_helper,
    )

    unsupported = sorted(set(protein_hits).difference(SEARCH_BACKENDS))
    if unsupported:
        raise ValueError(f"unsupported search backends: {unsupported}")
    entry_id_by_structure = {
        str(structure_id): str(entry.pdb_id)
        for structure_id, entry in entries_by_structure.items()
    }
    if len(entry_id_by_structure) != len(set(entry_id_by_structure.values())):
        raise ValueError("custom scoring entries must have distinct entry IDs")
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    outputs: dict[str, Path] = {}
    required = {
        "structure_id",
        "query_chain_asym_id",
        "query_sequence_source",
        "query_resolved_residue_numbers",
        "target_entry",
        "target_chain_asym_id",
        "target_selected_residue_numbers",
        "target_selected_residue_indices",
        "qstart",
        "tstart",
        "qcov",
        "tcov",
        "fident",
        "qaln",
        "taln",
    }
    empty_columns = (
        "query_entry",
        "target_entry",
        "query_chain_mapped",
        "target_chain_mapped",
        "source",
    )
    for backend, hit_path in protein_hits.items():
        hits = pd.read_parquet(hit_path)
        output = output_dir / f"{backend}.parquet"
        if hits.empty:
            pd.DataFrame({
                column: pd.Series(dtype="string") for column in empty_columns
            }).to_parquet(output, index=False)
            outputs[backend] = output
            continue
        missing = sorted(required.difference(hits.columns))
        if missing:
            raise ValueError(f"{backend} custom hit table is missing columns {missing}")
        unknown_structures = sorted(
            set(hits["structure_id"].astype(str)).difference(entry_id_by_structure)
        )
        if unknown_structures:
            raise ValueError(
                f"{backend} hits reference unknown custom structures: "
                f"{unknown_structures}"
            )
        custom_entry_ids = hits["structure_id"].astype(str).map(entry_id_by_structure)
        hits = hits.loc[
            hits["target_entry"].astype(str).to_numpy() != custom_entry_ids.to_numpy()
        ].copy()
        if hits.empty:
            pd.DataFrame({
                column: pd.Series(dtype="string") for column in empty_columns
            }).to_parquet(output, index=False)
            outputs[backend] = output
            continue

        selected = [
            _release_positions_for_custom_hit(row, backend=backend)
            for row in hits.itertuples(index=False)
        ]
        custom_entries = hits["structure_id"].astype(str).map(entry_id_by_structure)
        release_entries = hits["target_entry"].astype(str).copy()
        custom_chains = hits["query_chain_asym_id"].astype(str).copy()
        release_chains = hits["target_chain_asym_id"].astype(str).copy()
        custom_coverage = hits["qcov"].copy()
        hits["query_entry"] = release_entries
        hits["target_entry"] = custom_entries
        hits["query_chain_mapped"] = release_chains
        hits["target_chain_mapped"] = custom_chains
        hits["source"] = backend
        hits["qcov"] = hits["tcov"]
        hits["tcov"] = custom_coverage
        hits["query_selected_residue_numbers"] = [value[0] for value in selected]
        hits["target_selected_residue_numbers"] = [value[1] for value in selected]
        hits["selected_residue_identity"] = [value[2] for value in selected]
        hits["seqsim"] = [
            get_sequence_similarity_helper(str(qaln).upper(), str(taln).upper())
            for qaln, taln in zip(hits["qaln"], hits["taln"], strict=True)
        ]
        hits["fident_qcov"] = hits["fident"] * hits["qcov"]
        hits["seqsim_qcov"] = hits["seqsim"] * hits["qcov"]
        if backend == "foldseek":
            if "lddt" not in hits:
                raise ValueError("foldseek custom hit table is missing column lddt")
            hits["lddt_qcov"] = hits["lddt"] * hits["qcov"]
        hits = hits.drop(
            columns=[
                "qaln",
                "taln",
                "target_selected_residue_indices",
            ],
            errors="ignore",
        )
        leading = list(empty_columns)
        hits = hits[leading + [column for column in hits if column not in leading]]
        hits = hits.sort_values(
            [
                "query_entry",
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
            ],
            ignore_index=True,
        )
        hits.to_parquet(output, index=False, compression="zstd")
        outputs[backend] = output
    return outputs


def run_custom_protein_searches(
    *,
    query_databases: CustomQueryDatabases,
    assets: CustomScoringAssets,
    output_dir: Path,
    scratch_dir: Path,
    config: CustomProteinSearchConfig | None = None,
    backends: Iterable[str] = CIF_SEARCH_BACKENDS,
    threads: int = 1,
) -> dict[str, Path]:
    """Search custom protein chains and write backend-independent hit tables."""
    from plinder.data.annotations.get_similarity_scores import run_alignment
    from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig

    if threads < 1:
        raise ValueError("threads must be positive")
    config = config or CustomProteinSearchConfig()
    selected = tuple(dict.fromkeys(backends))
    unsupported = sorted(set(selected).difference(SEARCH_BACKENDS))
    if unsupported:
        raise ValueError(f"unsupported search backends: {unsupported}")
    output_dir = Path(output_dir)
    scratch_dir = Path(scratch_dir)
    raw_root = output_dir / "raw"
    raw_root.mkdir(exist_ok=True, parents=True)
    outputs: dict[str, Path] = {}
    for backend in selected:
        if backend not in assets.search_databases:
            raise ValueError(f"custom scoring assets omit {backend}")
        bundle = assets.search_databases[backend]
        backend_scratch = scratch_dir / backend
        backend_scratch.mkdir(exist_ok=True, parents=True)
        alignment_config: FoldseekConfig | MMSeqsConfig
        if backend == "foldseek":
            alignment_config = FoldseekConfig(
                evalue=config.evalue,
                sensitivity=config.sensitivity,
                max_seqs=config.max_seqs,
                coverage=config.coverage,
                min_seq_id=config.min_seq_id,
            )
        else:
            alignment_config = MMSeqsConfig(
                evalue=config.evalue,
                sensitivity=config.sensitivity,
                max_seqs=config.max_seqs,
                coverage=config.coverage,
                min_seq_id=config.min_seq_id,
            )
        raw_prefix = raw_root / backend
        run_alignment(
            aln_type=backend,
            query_db=query_databases.databases[backend],
            target_db=bundle.conversion_target,
            search_target_db=bundle.search_target,
            cluster_alignment_db=bundle.cluster_alignments,
            search_db=backend_scratch / "result",
            aln_file=raw_prefix,
            alignment_config=alignment_config,
            tmp_dir=backend_scratch / "tmp",
            remove_tmp=True,
            threads=threads,
        )
        output = output_dir / f"{backend}.parquet"
        map_custom_alignment_hits(
            raw_alignment=raw_prefix.with_suffix(".parquet"),
            backend=backend,
            query_databases=query_databases,
            alignment_chain_lookup=assets.alignment_chain_lookup,
            output_path=output,
        )
        outputs[backend] = output
    return outputs


def _target_entry_ids(score_alignments: Mapping[str, Path]) -> set[str]:
    target_ids: set[str] = set()
    for path in score_alignments.values():
        frame = pd.read_parquet(path, columns=["target_entry"])
        target_ids.update(frame["target_entry"].dropna().astype(str))
    return target_ids


def _query_entry_ids(score_alignments: Mapping[str, Path]) -> set[str]:
    query_ids: set[str] = set()
    for path in score_alignments.values():
        frame = pd.read_parquet(path, columns=["query_entry"])
        query_ids.update(frame["query_entry"].dropna().astype(str))
    return query_ids


def _load_release_entry_views(
    assets: CustomScoringAssets,
    *,
    pdb_ids: Iterable[str],
) -> dict[str, Any]:
    from plinder.core.scores.entries import entry_views_from_df

    selected = sorted(set(map(str, pdb_ids)))
    if not selected:
        return {}
    annotation = pd.read_parquet(
        assets.annotation_table,
        filters=[("entry_pdb_id", "in", selected)],
    )
    entry_chains = pd.read_parquet(
        assets.entry_chains,
        filters=[("entry_pdb_id", "in", selected)],
    )
    interfaces = pd.read_parquet(
        assets.interface_annotations,
        filters=[("entry_pdb_id", "in", selected)],
    )
    return entry_views_from_df(
        annotation,
        entry_chains=entry_chains,
        interface_annotations=interfaces,
    )


def _custom_ligand_sdf_resolver(
    *,
    query_entry_ids: set[str],
    query_ligand_root: Path,
    assets: CustomScoringAssets,
    extracted_target_root: Path,
    data_dir: Path | None,
) -> Any:
    """Return an SDF resolver that opens target archives only when called."""
    archive_cache = dict(assets.ligand_archives)
    packed_sdf_cache: dict[str, dict[tuple[str, str], bytes]] = {}
    cache_lock = Lock()

    def resolve(ligand: Any) -> Path | None:
        pdb_id = str(ligand.pdb_id)
        asym_id = str(ligand.asym_id)
        if pdb_id in query_entry_ids:
            custom_path = query_ligand_root / pdb_id / "ligand_files" / f"{asym_id}.sdf"
            return custom_path if custom_path.is_file() else None

        if data_dir is not None:
            raw_path = (
                Path(data_dir)
                / "raw_entries"
                / _pdb_shard(pdb_id)
                / pdb_id
                / "ligand_files"
                / f"{asym_id}.sdf"
            )
            if raw_path.is_file():
                return raw_path
        extracted = extracted_target_root / pdb_id / "ligand_files" / f"{asym_id}.sdf"
        if extracted.is_file():
            return extracted

        shard = _pdb_shard(pdb_id)
        with cache_lock:
            if extracted.is_file():
                return extracted
            if shard not in packed_sdf_cache:
                archive = archive_cache.get(shard)
                if archive is None:
                    archive = _release_file(
                        "ligand_archive",
                        data_dir=data_dir,
                        description=f"canonical ligand archive for shard {shard}",
                        shard=shard,
                    )
                    archive_cache[shard] = archive
                packed = pd.read_parquet(
                    archive,
                    columns=["pdb_id", "ligand_asym_id", "sdf"],
                )
                packed_sdf_cache[shard] = {
                    (str(row.pdb_id), str(row.ligand_asym_id)): bytes(row.sdf)
                    for row in packed.itertuples(index=False)
                }
            sdf = packed_sdf_cache[shard].get((pdb_id, asym_id))
        if sdf is None:
            LOG.warning(
                "canonical ligand archive has no row for %s/%s",
                pdb_id,
                asym_id,
            )
            return None
        extracted.parent.mkdir(exist_ok=True, parents=True)
        temporary = extracted.with_suffix(".tmp.sdf")
        temporary.write_bytes(sdf)
        temporary.replace(extracted)
        return extracted

    return resolve


def calculate_custom_protein_similarity_scores(
    protein_score_alignments: Mapping[str, Path],
    *,
    assets: CustomScoringAssets,
    work_dir: Path,
    plinder_system_ids: Iterable[str] | None = None,
    custom_chain_ids: Iterable[str] | None = None,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Score PLINDER receptors and ligand pockets against custom protein chains."""
    from plinder.core.utils.schemas import PROTEIN_SIMILARITY_SCHEMA
    from plinder.data.annotations.get_similarity_scores import Scorer

    entries = _load_release_entry_views(
        assets,
        pdb_ids=_query_entry_ids(protein_score_alignments),
    )
    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True, parents=True)
    scorer = Scorer(
        entries=entries,
        source_to_full_db_file={},
        db_dir=work_dir,
        scores_dir=work_dir,
        include_pli_fident=True,
    )
    selected_plinder_systems = (
        set(map(str, plinder_system_ids)) if plinder_system_ids is not None else None
    )
    selected_custom_chains = (
        set(map(str, custom_chain_ids)) if custom_chain_ids is not None else None
    )
    source_to_alignment = {
        f"apo_{backend}": path for backend, path in protein_score_alignments.items()
    }
    alignments = scorer.load_alignments(
        source_to_aln_file=source_to_alignment,
        search_db="apo",
    )
    alignments_by_query_entry = (
        {
            str(entry_id): group.droplevel("query_entry")
            for entry_id, group in alignments.groupby(level="query_entry", sort=False)
        }
        if not alignments.empty
        else {}
    )
    frames: list[pd.DataFrame] = []
    for entry_id in sorted(entries):
        query_entry_alignments = alignments_by_query_entry.get(entry_id)
        if query_entry_alignments is None:
            continue
        frame = scorer.aggregate_scores(
            entry_id,
            search_db="apo",
            source_to_aln_file=source_to_alignment,
            query_entry_alignments=query_entry_alignments,
            query_system_ids=selected_plinder_systems,
            target_system_ids=selected_custom_chains,
        )
        if frame is not None and not frame.empty:
            frames.append(frame)
    result = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=PROTEIN_SIMILARITY_SCHEMA.names)
    )
    if not result.empty:
        result = result.sort_values(
            [
                "similarity",
                "query_system",
                "query_ligand_id",
                "target_system",
            ],
            ascending=[False, True, True, True],
            ignore_index=True,
        )
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(exist_ok=True, parents=True)
        result.to_parquet(output_path, index=False, compression="zstd")
    return result


def write_custom_aligned_pocket_residues(
    protein_score_alignments: Mapping[str, Path],
    *,
    protein_scores: Path,
    assets: CustomScoringAssets,
    output_path: Path,
) -> Path:
    """Write residue pairs supporting custom-chain pocket identity scores."""
    columns = [
        "plinder_system_id",
        "plinder_ligand_id",
        "plinder_entry_id",
        "plinder_chain_instance",
        "plinder_chain_asym_id",
        "plinder_residue_number",
        "custom_structure_id",
        "custom_chain_asym_id",
        "custom_residue_number",
        "residue_identical",
        "source",
    ]
    scores = pd.read_parquet(
        protein_scores,
        columns=[
            "query_system",
            "query_ligand_id",
            "target_system",
            "protein_mapping",
            "source",
            "metric",
        ],
    )
    pocket_scores = scores.loc[
        scores["metric"].astype(str).eq("pocket_fident"),
        [
            "query_system",
            "query_ligand_id",
            "target_system",
            "protein_mapping",
            "source",
        ],
    ].drop_duplicates()
    accepted: set[tuple[str, str, str, str, str, str]] = set()
    for score in pocket_scores.itertuples(index=False):
        if not isinstance(score.protein_mapping, str):
            continue
        score_sources = (
            tuple(protein_score_alignments)
            if str(score.source) == "both"
            else (str(score.source),)
        )
        for pair in score.protein_mapping.split(";"):
            if ":" not in pair:
                continue
            release_instance, custom_instance = pair.split(":", maxsplit=1)
            release_chain = release_instance.split(".", maxsplit=1)[-1]
            custom_chain = custom_instance.split(".", maxsplit=1)[-1]
            for source in score_sources:
                accepted.add((
                    str(score.query_system),
                    str(score.query_ligand_id),
                    str(score.target_system),
                    source,
                    release_chain,
                    custom_chain,
                ))
    if not accepted:
        result = pd.DataFrame(columns=columns)
    else:
        entries = _load_release_entry_views(
            assets,
            pdb_ids=_query_entry_ids(protein_score_alignments),
        )
        pocket_membership: dict[tuple[str, str, int], list[tuple[str, str, str]]] = {}
        for entry_id, entry in entries.items():
            for system in entry.systems.values():
                for ligand in system.ligands.values():
                    if not ligand.is_proper:
                        continue
                    for (
                        instance_chain,
                        number_to_index,
                    ) in ligand.pocket_residue_number_to_index.items():
                        asym_id = instance_chain.split(".", maxsplit=1)[-1]
                        for residue_number in number_to_index:
                            pocket_membership.setdefault(
                                (str(entry_id), asym_id, int(residue_number)), []
                            ).append((
                                str(system.id),
                                str(ligand.id),
                                str(instance_chain),
                            ))

        rows: list[dict[str, Any]] = []
        for backend, alignment_path in protein_score_alignments.items():
            alignments = pd.read_parquet(alignment_path)
            for row in alignments.itertuples(index=False):
                release_entry = str(row.query_entry)
                release_chain = str(row.query_chain_mapped)
                custom_entry = str(row.target_entry)
                custom_chain = str(row.target_chain_mapped)
                target_system = f"{custom_entry}_{custom_chain}"
                release_numbers = _int_values(
                    row.query_selected_residue_numbers,
                    column="query_selected_residue_numbers",
                )
                custom_numbers = _int_values(
                    row.target_selected_residue_numbers,
                    column="target_selected_residue_numbers",
                )
                identities = bytes(row.selected_residue_identity)
                if not (len(release_numbers) == len(custom_numbers) == len(identities)):
                    raise ValueError(
                        "custom pocket residue alignment columns have different "
                        f"lengths for {release_entry} chain {release_chain}"
                    )
                for release_number, custom_number, identical in zip(
                    release_numbers,
                    custom_numbers,
                    identities,
                    strict=True,
                ):
                    for system_id, ligand_id, instance_chain in pocket_membership.get(
                        (release_entry, release_chain, release_number), []
                    ):
                        if (
                            system_id,
                            ligand_id,
                            target_system,
                            backend,
                            release_chain,
                            custom_chain,
                        ) not in accepted:
                            continue
                        rows.append({
                            "plinder_system_id": system_id,
                            "plinder_ligand_id": ligand_id,
                            "plinder_entry_id": release_entry,
                            "plinder_chain_instance": instance_chain,
                            "plinder_chain_asym_id": release_chain,
                            "plinder_residue_number": release_number,
                            "custom_structure_id": custom_entry,
                            "custom_chain_asym_id": custom_chain,
                            "custom_residue_number": (
                                custom_number if custom_number >= 0 else pd.NA
                            ),
                            "residue_identical": bool(identical),
                            "source": backend,
                        })
        result = pd.DataFrame.from_records(rows, columns=columns)
        if not result.empty:
            result = result.drop_duplicates().sort_values(
                [
                    "custom_structure_id",
                    "custom_chain_asym_id",
                    "plinder_system_id",
                    "plinder_ligand_id",
                    "source",
                    "plinder_chain_instance",
                    "plinder_residue_number",
                ],
                ignore_index=True,
            )
    result["plinder_residue_number"] = result["plinder_residue_number"].astype("Int64")
    result["custom_residue_number"] = result["custom_residue_number"].astype("Int64")
    result["residue_identical"] = result["residue_identical"].astype("boolean")
    output_path = Path(output_path)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    result.to_parquet(output_path, index=False, compression="zstd")
    return output_path


def write_custom_sequence_link_tables(
    *,
    protein_scores: Path,
    chain_manifest: Path,
    annotation_table: Path,
    output_path: Path,
    best_output_path: Path,
) -> tuple[Path, Path]:
    """Attach sequence IDs and ligand chemistry to pocket-identity scores."""
    manifest = pd.read_parquet(
        chain_manifest,
        columns=[
            "structure_id",
            "chain_asym_id",
            "sequence_id",
            "sequence_length",
        ],
    )
    manifest["target_system"] = (
        manifest["structure_id"].astype(str)
        + "_"
        + manifest["chain_asym_id"].astype(str)
    )
    target_map = manifest[
        ["target_system", "sequence_id", "sequence_length"]
    ].drop_duplicates()
    if target_map.duplicated("target_system").any():
        raise ValueError("sequence manifest has ambiguous target-system IDs")

    scores = pd.read_parquet(protein_scores)
    score_identity_columns = [
        "query_system",
        "query_ligand_id",
        "target_system",
        "protein_mapping",
        "protein_mapper",
        "source",
    ]
    score_pair_columns = [
        "query_system",
        "query_ligand_id",
        "target_system",
    ]
    pocket = scores.loc[
        scores["metric"].astype(str).eq("pocket_fident"),
        [
            *score_identity_columns,
            "similarity",
        ],
    ].copy()
    pli = scores.loc[
        scores["metric"].astype(str).eq("pli_fident"),
        [*score_pair_columns, "similarity"],
    ].rename(columns={"similarity": "pli_fident"})
    if pli.duplicated(score_pair_columns).any():
        raise ValueError("protein scores contain duplicate PLI identity rows")
    pocket = pocket.merge(
        pli,
        on=score_pair_columns,
        how="left",
        validate="one_to_one",
    )
    pocket = pocket.merge(
        target_map,
        on="target_system",
        how="left",
        validate="many_to_one",
    )
    if pocket["sequence_id"].isna().any():
        missing = pocket.loc[pocket["sequence_id"].isna(), "target_system"].unique()
        raise ValueError(f"protein scores reference unknown sequences: {missing[:10]}")

    chemistry_columns = [
        "system_id",
        "ligand_id",
        "ligand_ccd_code",
        "ligand_unique_ccd_code",
        "ligand_smiles",
    ]
    scored_pairs = (
        pocket[["query_system", "query_ligand_id"]]
        .drop_duplicates()
        .rename(
            columns={
                "query_system": "system_id",
                "query_ligand_id": "ligand_id",
            }
        )
    )
    scored_system_ids = scored_pairs["system_id"].astype(str).drop_duplicates().tolist()
    chemistry = (
        pd.read_parquet(
            annotation_table,
            columns=chemistry_columns,
            filters=[("system_id", "in", scored_system_ids)],
        ).drop_duplicates()
        if scored_system_ids
        else pd.DataFrame(columns=chemistry_columns)
    )
    chemistry = chemistry.merge(
        scored_pairs,
        on=["system_id", "ligand_id"],
        how="inner",
        validate="many_to_one",
    )
    if chemistry.duplicated(["system_id", "ligand_id"]).any():
        raise ValueError("release annotation has conflicting ligand chemistry rows")
    pocket = pocket.merge(
        chemistry,
        left_on=["query_system", "query_ligand_id"],
        right_on=["system_id", "ligand_id"],
        how="left",
        validate="many_to_one",
        indicator="_chemistry_merge",
    )
    missing_chemistry = pocket["_chemistry_merge"].ne("both")
    if missing_chemistry.any():
        missing = pocket.loc[
            missing_chemistry, ["query_system", "query_ligand_id"]
        ].drop_duplicates()
        raise ValueError(
            "release annotation has no row for scored PLINDER ligands: "
            f"{missing.head(10).to_dict('records')}"
        )
    pocket["plinder_pdb_id"] = (
        pocket["query_system"].astype(str).str.split("__", n=1).str[0]
    )
    pocket = pocket.rename(
        columns={
            "query_system": "plinder_system_id",
            "query_ligand_id": "plinder_ligand_id",
            "similarity": "pocket_fident",
            "ligand_ccd_code": "plinder_ligand_ccd_code",
            "ligand_unique_ccd_code": "plinder_ligand_unique_ccd_code",
            "ligand_smiles": "plinder_ligand_smiles",
        }
    ).drop(columns=["system_id", "ligand_id", "target_system", "_chemistry_merge"])
    link_columns = [
        "sequence_id",
        "sequence_length",
        "plinder_pdb_id",
        "plinder_system_id",
        "plinder_ligand_id",
        "pocket_fident",
        "pli_fident",
        "protein_mapping",
        "protein_mapper",
        "source",
        "plinder_ligand_ccd_code",
        "plinder_ligand_unique_ccd_code",
        "plinder_ligand_smiles",
    ]
    pocket = pocket[link_columns].sort_values(
        ["sequence_id", "pocket_fident", "plinder_pdb_id", "plinder_ligand_id"],
        ascending=[True, False, True, True],
        ignore_index=True,
    )
    output_path = Path(output_path)
    best_output_path = Path(best_output_path)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    pocket.to_parquet(output_path, index=False, compression="zstd")

    best_hits = pocket.drop_duplicates("sequence_id", keep="first")
    best = manifest[["sequence_id", "sequence_length"]].merge(
        best_hits.drop(columns="sequence_length"),
        on="sequence_id",
        how="left",
        validate="one_to_one",
    )
    best = best[link_columns]
    best.to_parquet(best_output_path, index=False, compression="zstd")
    return output_path, best_output_path


def add_sequence_ids_to_aligned_pocket_residues(
    aligned_pocket_residues: Path,
    *,
    chain_manifest: Path,
) -> Path:
    """Add original FASTA identifiers to sequence-query residue mappings."""
    path = Path(aligned_pocket_residues)
    residues = pd.read_parquet(path)
    identifiers = pd.read_parquet(
        chain_manifest,
        columns=["structure_id", "sequence_id"],
    ).drop_duplicates()
    if identifiers.duplicated("structure_id").any():
        raise ValueError("sequence manifest has ambiguous structure IDs")
    residues = residues.merge(
        identifiers,
        left_on="custom_structure_id",
        right_on="structure_id",
        how="left",
        validate="many_to_one",
    ).drop(columns="structure_id")
    if not residues.empty and residues["sequence_id"].isna().any():
        missing = residues.loc[
            residues["sequence_id"].isna(), "custom_structure_id"
        ].unique()
        raise ValueError(
            f"residue mappings reference unknown sequences: {missing[:10]}"
        )
    leading = ["sequence_id"]
    residues = residues[
        leading + [column for column in residues if column not in leading]
    ]
    temporary = path.with_suffix(".tmp.parquet")
    residues.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(path)
    return path


def score_custom_sequence_file(
    sequence_fasta: Path,
    *,
    work_dir: Path,
    scratch_dir: Path | None = None,
    data_dir: Path | None = None,
    search_config: CustomProteinSearchConfig | None = None,
    backends: Iterable[str] = SEQUENCE_SEARCH_BACKENDS,
    plinder_entry_ids: Iterable[str] | None = None,
    threads: int = 1,
    store_aligned_pocket_residues: bool = False,
) -> CustomSequenceScoringResult:
    """Search protein sequences and score PLINDER ligand-pocket identity.

    ``plinder_entry_ids`` builds a small MMseqs target from the scoreable
    receptor and interface chains in the selected release entries instead of
    fetching the complete search database. This bounded mode currently requires
    ``backends=("mmseqs",)``.
    """
    selected_backends = tuple(dict.fromkeys(backends))
    if not selected_backends:
        raise ValueError("at least one sequence search backend is required")
    unsupported = sorted(set(selected_backends).difference(SEQUENCE_SEARCH_BACKENDS))
    if unsupported:
        raise ValueError(
            f"unsupported sequence search backends: {unsupported}; "
            f"expected a subset of {SEQUENCE_SEARCH_BACKENDS}"
        )
    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True, parents=True)
    selected_plinder_entries = _plinder_entry_subset(plinder_entry_ids)
    query_inputs = write_custom_sequence_query_files(
        Path(sequence_fasta),
        work_dir=work_dir,
    )
    query_databases = create_custom_query_databases(
        query_inputs,
        work_dir=work_dir,
        backends=selected_backends,
        threads=threads,
    )
    assets = _resolve_workflow_assets(
        data_dir=data_dir,
        backends=selected_backends,
        plinder_entry_ids=selected_plinder_entries,
        work_dir=work_dir,
        threads=threads,
    )
    protein_hits = run_custom_protein_searches(
        query_databases=query_databases,
        assets=assets,
        output_dir=work_dir / "protein_hits",
        scratch_dir=(
            Path(scratch_dir) if scratch_dir is not None else work_dir / "scratch"
        ),
        config=search_config,
        backends=selected_backends,
        threads=threads,
    )
    manifest = pd.read_parquet(query_inputs.chain_manifest)
    entries_by_structure = {
        str(structure_id): _SequenceEntry(pdb_id=str(structure_id))
        for structure_id in manifest["structure_id"]
    }
    protein_score_alignments = prepare_custom_protein_score_alignments(
        protein_hits,
        entries_by_structure=entries_by_structure,
        output_dir=work_dir / "protein_score_alignments",
    )
    protein_score_path = work_dir / "protein_scores.parquet"
    calculate_custom_protein_similarity_scores(
        protein_score_alignments,
        assets=assets,
        work_dir=work_dir / "protein_score_work",
        output_path=protein_score_path,
    )
    aligned_pocket_residue_path = (
        work_dir / "aligned_pocket_residues.parquet"
        if store_aligned_pocket_residues
        else None
    )
    if aligned_pocket_residue_path is not None:
        write_custom_aligned_pocket_residues(
            protein_score_alignments,
            protein_scores=protein_score_path,
            assets=assets,
            output_path=aligned_pocket_residue_path,
        )
        add_sequence_ids_to_aligned_pocket_residues(
            aligned_pocket_residue_path,
            chain_manifest=query_inputs.chain_manifest,
        )
    sequence_links, best_sequence_links = write_custom_sequence_link_tables(
        protein_scores=protein_score_path,
        chain_manifest=query_inputs.chain_manifest,
        annotation_table=assets.annotation_table,
        output_path=work_dir / "sequence_links.parquet",
        best_output_path=work_dir / "best_sequence_links.parquet",
    )
    return CustomSequenceScoringResult(
        query_inputs=query_inputs,
        query_databases=query_databases,
        protein_hits=protein_hits,
        protein_score_alignments=protein_score_alignments,
        protein_scores=protein_score_path,
        aligned_pocket_residues=aligned_pocket_residue_path,
        sequence_links=sequence_links,
        best_sequence_links=best_sequence_links,
    )


def calculate_custom_similarity_scores(
    score_alignments: Mapping[str, Path],
    *,
    annotations: CustomStructureAnnotations,
    assets: CustomScoringAssets,
    work_dir: Path,
    include_shape: bool = True,
    shape_score_threads: int = 1,
    target_system_ids: Iterable[str] | None = None,
    target_ligand_ids: Iterable[str] | None = None,
    data_dir: Path | None = None,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Calculate directed ligand-pocket scores from custom structures."""
    from plinder.core.utils.schemas import PROTEIN_SIMILARITY_SCHEMA
    from plinder.data.annotations.get_similarity_scores import Scorer

    query_entries = {
        str(entry.pdb_id): entry for entry in annotations.entries_by_structure.values()
    }
    target_entries = _load_release_entry_views(
        assets,
        pdb_ids=_target_entry_ids(score_alignments).difference(query_entries),
    )
    entries = {**target_entries, **query_entries}
    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True, parents=True)
    scorer = Scorer(
        entries=entries,
        source_to_full_db_file={},
        db_dir=work_dir,
        scores_dir=work_dir,
        ligand_sdf_resolver=_custom_ligand_sdf_resolver(
            query_entry_ids=set(query_entries),
            query_ligand_root=annotations.ligand_sdf_root,
            assets=assets,
            extracted_target_root=work_dir / "target_ligands",
            data_dir=data_dir,
        ),
        shape_score_threads=shape_score_threads,
        include_pli_fident=True,
    )
    selected_systems = (
        set(map(str, target_system_ids)) if target_system_ids is not None else None
    )
    selected_ligands = (
        set(map(str, target_ligand_ids)) if target_ligand_ids is not None else None
    )
    source_to_alignment = {
        f"holo_{backend}": path for backend, path in score_alignments.items()
    }
    frames: list[pd.DataFrame] = []
    for query_entry_id, query_entry in query_entries.items():
        if not query_entry.systems:
            continue
        frame = scorer.aggregate_scores(
            query_entry_id,
            search_db="holo",
            data_dir=work_dir if include_shape else None,
            source_to_aln_file=source_to_alignment,
            target_system_ids=selected_systems,
            target_ligand_ids=selected_ligands,
            include_holo_protein_scores=True,
        )
        if frame is not None and not frame.empty:
            frames.append(frame)
    result = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=PROTEIN_SIMILARITY_SCHEMA.names)
    )
    if not result.empty:
        result = result.sort_values(
            [
                "similarity",
                "query_system",
                "query_ligand_id",
                "target_system",
                "target_ligand_id",
            ],
            ascending=[False, True, True, True, True],
            ignore_index=True,
        )
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(exist_ok=True, parents=True)
        result.to_parquet(output_path, index=False, compression="zstd")
    return result


def calculate_custom_interface_similarity_scores(
    score_alignments: Mapping[str, Path],
    *,
    annotations: CustomStructureAnnotations,
    assets: CustomScoringAssets,
    target_interface_ids: Iterable[str] | None = None,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Calculate directed protein-interface coverage from custom structures."""
    from plinder.core.scores.reconstruct import (
        calculate_interface_similarity_scores,
    )

    query_entry_ids = {
        str(entry.pdb_id) for entry in annotations.entries_by_structure.values()
    }
    target_entries = _load_release_entry_views(
        assets,
        pdb_ids=_target_entry_ids(score_alignments).difference(query_entry_ids),
    )
    query_interfaces = {
        interface_id: interface
        for entry in annotations.entries_by_structure.values()
        for interface_id, interface in entry.interfaces.items()
    }
    target_interfaces = {
        interface_id: interface
        for entry in target_entries.values()
        for interface_id, interface in entry.interfaces.items()
    }
    if target_interface_ids is not None:
        selected = set(map(str, target_interface_ids))
        missing = selected.difference(target_interfaces)
        if missing:
            raise KeyError(f"unknown target interface system IDs: {sorted(missing)}")
        target_interfaces = {
            interface_id: target_interfaces[interface_id] for interface_id in selected
        }
    alignment_frames = [pd.read_parquet(path) for path in score_alignments.values()]
    alignments = (
        pd.concat(alignment_frames, ignore_index=True)
        if alignment_frames
        else pd.DataFrame()
    )
    result = calculate_interface_similarity_scores(
        alignments,
        query_interfaces=query_interfaces,
        target_interfaces=target_interfaces,
    )
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(exist_ok=True, parents=True)
        result.to_parquet(output_path, index=False, compression="zstd")
    return result


def score_custom_cif_files(
    cif_files: Iterable[Path],
    *,
    work_dir: Path,
    scratch_dir: Path | None = None,
    data_dir: Path | None = None,
    structure_mode: Literal["as_is", "pdb"] = "as_is",
    assembly_ids: Iterable[str] | None = None,
    ligand_smiles_dict: Mapping[str, str] | None = None,
    ligand_ccd_code_dict: Mapping[str, str] | None = None,
    include_ligands: bool | None = True,
    include_interfaces: bool | None = True,
    include_shape: bool = True,
    interface_annotate_prodigy: bool = True,
    search_config: CustomProteinSearchConfig | None = None,
    backends: Iterable[str] = CIF_SEARCH_BACKENDS,
    plinder_entry_ids: Iterable[str] | None = None,
    threads: int = 1,
    shape_score_threads: int = 1,
    store_aligned_pocket_residues: bool = False,
) -> CustomScoringResult:
    """Run custom-CIF annotation, protein search, and PLINDER scoring.

    Unknown component IDs such as ``LIG`` can be assigned chemistry with
    either ``ligand_smiles_dict`` or ``ligand_ccd_code_dict``. The resulting
    bond-aware SDFs are used as the query conformers for shape, color, and
    SuCOS scores. Target ligand archives are first opened after a positive
    pocket match, and each required archive shard is read at most once.

    The returned object points to the written custom annotation tables,
    backend hit tables, compact score alignments, and the final protein,
    ligand, and/or interface score parquet files under ``work_dir``. Protein
    scores use each PLINDER ligand pocket as the directed query and each
    custom protein chain as a ligand-free target, so ``include_ligands=False``
    still yields protein metrics and ``pocket_fident``. Passing ``None`` for
    ligand or interface inclusion annotates that feature and emits its score
    table only when the custom structures contain a proper ligand or protein
    interface, respectively. ``plinder_entry_ids`` builds a small MMseqs target
    from scoreable receptor and interface chains in selected release entries
    rather than fetching the complete search database; this bounded mode
    currently requires ``backends=("mmseqs",)``.
    """
    sources = tuple(Path(path) for path in cif_files)
    if not sources:
        raise ValueError("at least one custom mmCIF is required")
    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True, parents=True)
    selected_plinder_entries = _plinder_entry_subset(plinder_entry_ids)
    selected_backends = tuple(dict.fromkeys(backends))
    if not selected_backends:
        raise ValueError("at least one custom search backend is required")
    requested_assemblies = (
        tuple(str(value) for value in assembly_ids)
        if assembly_ids is not None and not isinstance(assembly_ids, str)
        else assembly_ids
    )
    annotate_ligands = include_ligands is not False
    annotate_interfaces = include_interfaces is not False
    annotations = annotate_custom_cif_files(
        sources,
        work_dir=work_dir,
        structure_mode=structure_mode,
        assembly_ids=requested_assemblies,
        ligand_smiles_dict=ligand_smiles_dict,
        ligand_ccd_code_dict=ligand_ccd_code_dict,
        include_ligands=annotate_ligands,
        include_interfaces=annotate_interfaces,
        interface_annotate_prodigy=interface_annotate_prodigy,
        data_dir=data_dir,
    )
    query_inputs = write_custom_query_files(
        sources,
        work_dir=work_dir,
        structure_mode=structure_mode,
        assembly_ids=requested_assemblies,
    )
    query_databases = create_custom_query_databases(
        query_inputs,
        work_dir=work_dir,
        backends=selected_backends,
        threads=threads,
    )
    assets = _resolve_workflow_assets(
        data_dir=data_dir,
        backends=selected_backends,
        plinder_entry_ids=selected_plinder_entries,
        work_dir=work_dir,
        threads=threads,
    )
    protein_hits = run_custom_protein_searches(
        query_databases=query_databases,
        assets=assets,
        output_dir=work_dir / "protein_hits",
        scratch_dir=(
            Path(scratch_dir) if scratch_dir is not None else work_dir / "scratch"
        ),
        config=search_config,
        backends=selected_backends,
        threads=threads,
    )
    score_alignments = prepare_custom_score_alignments(
        protein_hits,
        entries_by_structure=annotations.entries_by_structure,
        output_dir=work_dir / "score_alignments",
    )
    protein_score_alignments = prepare_custom_protein_score_alignments(
        protein_hits,
        entries_by_structure=annotations.entries_by_structure,
        output_dir=work_dir / "protein_score_alignments",
    )
    protein_score_path = work_dir / "protein_scores.parquet"
    calculate_custom_protein_similarity_scores(
        protein_score_alignments,
        assets=assets,
        work_dir=work_dir / "protein_score_work",
        output_path=protein_score_path,
    )
    aligned_pocket_residue_path = (
        work_dir / "aligned_pocket_residues.parquet"
        if store_aligned_pocket_residues
        else None
    )
    if aligned_pocket_residue_path is not None:
        write_custom_aligned_pocket_residues(
            protein_score_alignments,
            protein_scores=protein_score_path,
            assets=assets,
            output_path=aligned_pocket_residue_path,
        )
    has_proper_ligands = any(
        ligand.is_proper
        for entry in annotations.entries_by_structure.values()
        for system in entry.systems.values()
        for ligand in system.ligands.values()
    )
    score_ligands = include_ligands is True or (
        include_ligands is None and has_proper_ligands
    )
    ligand_score_path = work_dir / "ligand_scores.parquet" if score_ligands else None
    if ligand_score_path is not None:
        calculate_custom_similarity_scores(
            score_alignments,
            annotations=annotations,
            assets=assets,
            work_dir=work_dir / "score_work",
            include_shape=include_shape,
            shape_score_threads=shape_score_threads,
            data_dir=data_dir,
            output_path=ligand_score_path,
        )
    has_interfaces = any(
        entry.interfaces for entry in annotations.entries_by_structure.values()
    )
    score_interfaces = include_interfaces is True or (
        include_interfaces is None and has_interfaces
    )
    interface_score_path = (
        work_dir / "interface_scores.parquet" if score_interfaces else None
    )
    if interface_score_path is not None:
        calculate_custom_interface_similarity_scores(
            score_alignments,
            annotations=annotations,
            assets=assets,
            output_path=interface_score_path,
        )
    return CustomScoringResult(
        annotations=annotations,
        query_inputs=query_inputs,
        query_databases=query_databases,
        protein_hits=protein_hits,
        score_alignments=score_alignments,
        protein_score_alignments=protein_score_alignments,
        protein_scores=protein_score_path,
        aligned_pocket_residues=aligned_pocket_residue_path,
        ligand_scores=ligand_score_path,
        interface_scores=interface_score_path,
    )
