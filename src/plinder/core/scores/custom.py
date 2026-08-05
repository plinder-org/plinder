# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Resolve the minimal release assets needed to score custom structures."""

from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from plinder.core.utils import cpl
from plinder.core.utils.config import get_config

LOG = logging.getLogger(__name__)
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


def _check_custom_mmcif_fields(
    block: Any,
    *,
    source: Path,
    structure_mode: str,
) -> None:
    """Check only the mmCIF fields required by the selected ingest mode."""
    if "atom_site" not in block:
        raise ValueError(f"custom mmCIF {source} has no _atom_site category")
    atom_site = block["atom_site"]
    required_columns = {
        "group_PDB",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    }
    missing = sorted(required_columns.difference(atom_site))
    identifier_pairs = {
        "atom name": ("label_atom_id", "auth_atom_id"),
        "residue name": ("label_comp_id", "auth_comp_id"),
        "chain ID": ("label_asym_id", "auth_asym_id"),
        "residue number": ("label_seq_id", "auth_seq_id"),
    }
    missing_identifiers = [
        f"{description} ({first} or {second})"
        for description, (first, second) in identifier_pairs.items()
        if first not in atom_site and second not in atom_site
    ]
    if missing or missing_identifiers:
        details = [f"_atom_site.{column}" for column in missing]
        details.extend(missing_identifiers)
        raise ValueError(
            f"custom mmCIF {source} is missing required coordinate fields: "
            + ", ".join(details)
        )

    # These fields have unambiguous defaults for single-model custom input.
    atom_count = atom_site.row_count
    if "pdbx_PDB_model_num" not in atom_site:
        atom_site["pdbx_PDB_model_num"] = np.ones(atom_count, dtype=np.int32)
    if "pdbx_PDB_ins_code" not in atom_site:
        atom_site["pdbx_PDB_ins_code"] = ["."] * atom_count

    if structure_mode != "pdb":
        return
    if "label_asym_id" not in atom_site:
        raise ValueError(
            f"custom mmCIF {source} needs _atom_site.label_asym_id in pdb mode "
            "because assembly definitions reference label asym IDs"
        )
    assembly_fields = {
        "pdbx_struct_assembly_gen": {
            "assembly_id",
            "oper_expression",
            "asym_id_list",
        },
        "pdbx_struct_oper_list": {
            "id",
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[1]",
            "vector[2]",
            "vector[3]",
        },
    }
    missing_assembly_fields: list[str] = []
    for category_name, columns in assembly_fields.items():
        if category_name not in block:
            missing_assembly_fields.append(f"_{category_name}")
            continue
        missing_assembly_fields.extend(
            f"_{category_name}.{column}"
            for column in sorted(columns.difference(block[category_name]))
        )
    if missing_assembly_fields:
        raise ValueError(
            f"custom mmCIF {source} is missing fields required by pdb mode: "
            + ", ".join(missing_assembly_fields)
            + "; use structure_mode='as_is' for an already assembled model"
        )


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

    from plinder.core.structure.atoms import is_hydrogen_isotope
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
        atoms = atoms[~is_hydrogen_isotope(atoms.element)]
        return {
            str(chain_id): atoms[atoms.chain_id == chain_id]
            for chain_id in np.unique(atoms.chain_id)
        }
    if structure_mode != "pdb":
        raise ValueError("structure_mode must be 'as_is' or 'pdb'")

    block = list(cif_file.values())[0]
    available_assemblies = block["pdbx_struct_assembly_gen"][
        "assembly_id"
    ].as_array(str)
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

    from plinder.data.annotations.cif_utils import read_mmcif_file
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
            raise ValueError("custom input files must be outside the generated work tree")
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
        _check_custom_mmcif_fields(
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
            sequence = sequences.get(asym_id) or _resolved_protein_sequence(
                atoms_by_asym[asym_id],
                context=f"{structure_id} chain {asym_id}",
            )
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
            residue_starts = struc.get_residue_starts(
                atoms, add_exclusive_stop=False
            )
            if len(residue_starts) < 1:
                continue
            query_id = f"cq{len(rows):08d}"
            query_chain_id = f"{structure_id}__{asym_id}"
            resolved_residue_numbers = [
                int(atoms.res_id[index]) for index in residue_starts
            ]
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
            rows.append(
                {
                    "query_id": query_id,
                    "query_chain_id": query_chain_id,
                    "structure_id": structure_id,
                    "source_mmcif": str(source.resolve()),
                    "chain_asym_id": asym_id,
                    "sequence": sequence,
                    "sequence_length": len(sequence),
                    "resolved_residue_numbers": resolved_residue_numbers,
                }
            )
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
        query_id: values
        for query_id, values in matched.items()
        if len(values) != 1
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
    threads: int = 1,
) -> CustomQueryDatabases:
    """Create unindexed, scratch-local query DBs for both search backends."""
    if threads < 1:
        raise ValueError("threads must be positive")
    for backend in SEARCH_BACKENDS:
        if shutil.which(backend) is None:
            raise FileNotFoundError(
                f"{backend} executable is required for custom protein searches"
            )
    root = Path(work_dir) / "query_databases"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    databases = {backend: root / backend for backend in SEARCH_BACKENDS}
    _run_command(
        [
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
        ]
    )
    _run_command(
        [
            "mmseqs",
            "createdb",
            str(inputs.sequence_fasta),
            str(databases["mmseqs"]),
            "--threads",
            str(threads),
        ]
    )
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
    pd.DataFrame(mapping_rows).sort_values(
        ["backend", "query_id"]
    ).to_parquet(identifier_map, index=False)
    return CustomQueryDatabases(
        root=root,
        inputs=inputs,
        databases=databases,
        identifier_map=identifier_map,
    )


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
    parsed = pd.DataFrame(
        [
            {
                "target_backend_id": identifier,
                "target_entry": entry_id,
                "target_chain_auth_id": auth_id,
            }
            for identifier in sorted(set(target_identifiers))
            for entry_id, auth_id in [
                _parse_plinder_target_identifier(identifier, backend=backend)
            ]
        ]
    )
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
        lookup.groupby(key_columns, dropna=False)["target_chain_asym_id"]
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
            "PLINDER alignment lookup cannot map target chains: "
            f"{missing[:10]}"
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
        pd.DataFrame({column: pd.Series(dtype="string") for column in columns}).to_parquet(
            output_path, index=False
        )
        return output_path
    query_columns = query_mapping[
        [
            "backend_query_id",
            "query_id",
            "query_chain_id",
            "structure_id",
            "chain_asym_id",
        ]
    ].rename(
        columns={
            "backend_query_id": "query_backend_id",
            "chain_asym_id": "query_chain_asym_id",
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
    mapped = mapped.drop(
        columns=["query_pdb_id", "target_pdb_id"], errors="ignore"
    )
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


def run_custom_protein_searches(
    *,
    query_databases: CustomQueryDatabases,
    assets: CustomScoringAssets,
    output_dir: Path,
    scratch_dir: Path,
    config: CustomProteinSearchConfig | None = None,
    backends: Iterable[str] = SEARCH_BACKENDS,
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
