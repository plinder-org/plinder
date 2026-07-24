# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Protocol

import biotite.structure as struc
import biotite.structure.io.mol as mol
import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from rdkit import Chem

WaterSelection = Literal["none", "interacting", "all"]


class AnnotationRow(Protocol):
    """Minimal interface shared by dicts and pandas Series."""

    def get(self, key: str, default: Any = None) -> Any:
        ...


@dataclass(frozen=True)
class SystemReconstructionOptions:
    """Configure which atoms are present in reconstructed mmCIF views."""

    system_waters: WaterSelection = "interacting"
    receptor_waters: WaterSelection = "interacting"
    system_include_ligands: bool = True
    system_include_other_protein_chains: bool = False
    system_include_other_ligand_chains: bool = False
    receptor_include_other_protein_chains: bool = False
    validate_biounit_chain_ids: bool = True


@dataclass(frozen=True)
class SystemReconstructionOutputs:
    """Explicit output paths; ``None`` means that file is not written."""

    system_cif: Path | None = None
    receptor_cif: Path | None = None
    sequences_fasta: Path | None = None


@dataclass(frozen=True)
class ReconstructedSystem:
    """In-memory views rebuilt from a source PDB mmCIF."""

    biounit: struc.AtomArray
    system: struc.AtomArray
    receptor: struc.AtomArray


@dataclass(frozen=True)
class _ChainSelections:
    """Chain instances selected for reconstruction views."""

    expected_biounit: set[str]
    system: set[str]
    receptor: set[str]


def save_ligands(
    atoms: struc.AtomArray,
    ligand_chain_ids: list[str] | dict[str, list[str]],
    output_folder: str | Path,
) -> None:
    """Save one canonical ASU SDF per retained ligand.

    RDKit supplies chemically normalized bond/aromaticity information when it
    can sanitize the ligand.  If it cannot, Biotite serializes the source bond
    graph directly so the canonical coordinate archive remains complete.

    Parameters
    ----------
    atoms : AtomArray
        Full system atoms with bonds.
    ligand_chain_ids : list[str] or dict[str, list[str]]
        Either a list of chain IDs (one SDF per chain), or a mapping of
        output name -> member chain IDs. The mapping form writes a single
        SDF spanning all member chains, so covalently-linked ligand chains
        (a macrocycle deposited as several chains) are saved as one molecule.
    output_folder : str or Path
        Directory to write SDF files.
    """
    import logging

    from plinder.data.utils.annotations.cif_utils import atoms_to_rdkit_mol

    log = logging.getLogger(__name__)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # Normalize to {output_name: [member_chain_ids]}.
    groups = (
        ligand_chain_ids
        if isinstance(ligand_chain_ids, dict)
        else {chain_id: [chain_id] for chain_id in ligand_chain_ids}
    )

    for chain_id, member_chain_ids in groups.items():
        lig_mask = np.isin(atoms.chain_id, member_chain_ids)
        if not np.any(lig_mask):
            raise ValueError(f"No atoms found for ligand chain(s) {member_chain_ids!r}")
        lig_atoms = atoms[lig_mask]
        output_file = output_folder / f"{chain_id}.sdf"
        try:
            rdkit_mol = atoms_to_rdkit_mol(lig_atoms)
        except Exception as exc:
            log.warning(
                "save_ligands: RDKit conversion failed for chain %s; "
                "writing the unsanitized Biotite structure: %s",
                chain_id,
                exc,
            )
            sdf_file = mol.SDFile()
            mol.set_structure(sdf_file, lig_atoms, record_name=chain_id)
            sdf_file.write(output_file)
        else:
            rdkit_mol.SetProp("_Name", chain_id)
            with Chem.SDWriter(str(output_file)) as writer:
                writer.write(rdkit_mol)


def save_cif_file(
    atoms: struc.AtomArray,
    name: str,
    output_cif_file: str | Path,
) -> None:
    """Save structure as mmCIF.

    Parameters
    ----------
    atoms : AtomArray
        Atoms to save.
    name : str
        Data block name.
    output_cif_file : str or Path
        Output path.
    """
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, atoms, data_block=name, include_bonds=True)
    cif_file.write(str(output_cif_file))


def _string_list(value: Any) -> list[str]:
    """Normalize Arrow/Pandas/list values from an annotation row."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value]
    try:
        if bool(np.isnan(value)):
            return []
    except (TypeError, ValueError):
        pass
    raise TypeError(f"Expected a list-like annotation value, got {type(value)!r}")


def _annotation_chains(annotation: AnnotationRow, column: str) -> list[str]:
    return sorted(set(_string_list(annotation.get(column))))


def _water_mask(
    biounit: struc.AtomArray,
    annotation: AnnotationRow,
    selection: WaterSelection,
) -> NDArray[np.bool_]:
    if selection == "none":
        empty_mask: NDArray[np.bool_] = np.zeros(biounit.array_length(), dtype=np.bool_)
        return empty_mask
    if selection == "all":
        solvent_mask: NDArray[np.bool_] = struc.filter_solvent(biounit)
        return solvent_mask
    if selection != "interacting":
        raise ValueError(f"Unknown water selection: {selection!r}")

    mask: NDArray[np.bool_] = np.zeros(biounit.array_length(), dtype=np.bool_)
    for encoded_residue in _string_list(annotation.get("system_water_residues")):
        try:
            chain_id, residue_number = encoded_residue.rsplit("_", maxsplit=1)
            residue_id = int(residue_number)
        except ValueError as exc:
            raise ValueError(
                "system_water_residues values must have the form "
                f"'<instance>.<asym>_<residue>', got {encoded_residue!r}"
            ) from exc
        mask |= (biounit.chain_id == chain_id) & (biounit.res_id == residue_id)
    return mask


def _require_chains(
    biounit: struc.AtomArray,
    chain_ids: set[str],
    *,
    view_name: str,
) -> None:
    available = set(str(chain_id) for chain_id in np.unique(biounit.chain_id))
    missing = chain_ids - available
    if missing:
        raise ValueError(
            f"Cannot reconstruct {view_name}: assembly is missing chains "
            f"{sorted(missing)}"
        )


def _biounit_chain_roles(
    annotation: AnnotationRow,
    biounit_chains: pd.DataFrame | None,
) -> dict[str, str]:
    """Return chain-instance roles for the annotation's biological assembly."""
    if biounit_chains is None:
        # Existing V2 rows carry repeated assembly membership. New V3 ingest
        # never writes these columns, but retaining a read-only fallback keeps
        # source reconstruction usable for already-published V2 data.
        all_chains = set(
            _annotation_chains(annotation, "system_biounit_chains_asym_id")
        )
        if not all_chains:
            return {}
        waters = set(
            _annotation_chains(annotation, "system_biounit_water_chains_asym_id")
        )
        ligands = set(_annotation_chains(annotation, "system_ligand_chains"))
        ligands.update(
            _annotation_chains(annotation, "system_other_ligand_chains_asym_id")
        )
        return {
            chain: (
                "water"
                if chain in waters
                else "ligand"
                if chain in ligands
                else "receptor"
            )
            for chain in all_chains
        }
    required = {
        "entry_pdb_id",
        "biounit_id",
        "chain_instance",
        "chain_asym_id",
        "chain_role",
    }
    missing = required.difference(biounit_chains.columns)
    if missing:
        raise ValueError(
            "biounit-chain metadata is missing columns " f"{sorted(missing)}"
        )
    pdb_id = str(annotation.get("entry_pdb_id", ""))
    biounit_id = str(annotation.get("system_biounit_id", ""))
    selected = biounit_chains[
        biounit_chains["entry_pdb_id"].astype(str).eq(pdb_id)
        & biounit_chains["biounit_id"].astype(str).eq(biounit_id)
    ].copy()
    if selected.empty:
        raise ValueError(
            "No biological-assembly chain metadata for "
            f"entry={pdb_id!r}, biounit={biounit_id!r}"
        )
    if selected["chain_instance"].duplicated().any():
        duplicates = sorted(
            selected.loc[
                selected["chain_instance"].duplicated(keep=False),
                "chain_instance",
            ]
            .astype(str)
            .unique()
        )
        raise ValueError(f"Duplicate biological-assembly chains: {duplicates}")
    roles = dict(
        zip(
            selected["chain_instance"].astype(str),
            selected["chain_role"].astype(str),
        )
    )
    invalid_roles = sorted(set(roles.values()) - {"receptor", "ligand", "water"})
    if invalid_roles:
        raise ValueError(f"Unknown biological-assembly chain roles: {invalid_roles}")
    mismatched_asym_ids = selected[
        selected["chain_instance"]
        .astype(str)
        .str.split(".", n=1)
        .str[-1]
        .ne(selected["chain_asym_id"].astype(str))
    ]
    if not mismatched_asym_ids.empty:
        raise ValueError(
            "biounit-chain metadata has inconsistent chain instance/asym IDs"
        )
    return roles


def _select_chains(
    annotation: AnnotationRow,
    biounit_chains: pd.DataFrame | None,
    options: SystemReconstructionOptions,
) -> _ChainSelections:
    protein_chains = set(
        _annotation_chains(annotation, "system_protein_chains_asym_id")
    )
    ligand_chains = set(
        _annotation_chains(annotation, "system_ligand_chains_asym_id")
        or _annotation_chains(annotation, "system_ligand_chains")
    )
    if not protein_chains:
        raise ValueError("annotation contains no system protein chains")

    roles = _biounit_chain_roles(annotation, biounit_chains)
    include_other = (
        options.system_include_other_protein_chains
        or options.system_include_other_ligand_chains
        or options.receptor_include_other_protein_chains
    )
    if include_other and not roles:
        raise ValueError(
            "biounit-chain metadata is required when including other chains"
        )
    system_members = protein_chains | ligand_chains
    other_receptor_chains = {
        chain
        for chain, role in roles.items()
        if role == "receptor" and chain not in system_members
    }
    other_ligand_chains = {
        chain
        for chain, role in roles.items()
        if role == "ligand" and chain not in system_members
    }

    system_chains = set(protein_chains)
    if options.system_include_ligands:
        system_chains.update(ligand_chains)
    if options.system_include_other_protein_chains:
        system_chains.update(other_receptor_chains)
    if options.system_include_other_ligand_chains:
        system_chains.update(other_ligand_chains)

    receptor_chains = set(protein_chains)
    if options.receptor_include_other_protein_chains:
        receptor_chains.update(other_receptor_chains)
    return _ChainSelections(
        expected_biounit=set(roles),
        system=system_chains,
        receptor=receptor_chains,
    )


def reconstruct_system(
    source_mmcif: Path | str,
    annotation: AnnotationRow,
    *,
    biounit_chains: pd.DataFrame | None = None,
    options: SystemReconstructionOptions = SystemReconstructionOptions(),
) -> ReconstructedSystem:
    """Rebuild system and receptor views from a PDB mmCIF and parquet row.

    The annotation row must be one ligand-level row from the system to
    reconstruct.  System-level columns are repeated for every ligand row, so
    no grouping or coordinate data from the parquet is required.

    Parameters
    ----------
    source_mmcif : Path or str
        Original PDB mmCIF used for annotation (plain or gzip-compressed).
    annotation : mapping-like
        A dictionary or pandas Series containing the system selection columns.
    biounit_chains : pandas.DataFrame, optional
        Normalized biological-assembly membership rows. Required when any
        ``include_other_*`` option is enabled and used for assembly validation
        when supplied.
    options : SystemReconstructionOptions
        Atom-content choices for the two returned views.
    """
    from plinder.data.utils.annotations.cif_utils import (
        build_biounit,
        read_mmcif_file,
    )

    assembly_id = str(annotation.get("system_biounit_id", ""))
    if not assembly_id:
        raise ValueError("annotation is missing system_biounit_id")
    biounit = build_biounit(read_mmcif_file(source_mmcif), assembly_id)

    selections = _select_chains(annotation, biounit_chains, options)
    actual_biounit_chains = set(
        str(chain_id) for chain_id in np.unique(biounit.chain_id)
    )
    if (
        options.validate_biounit_chain_ids
        and selections.expected_biounit
        and selections.expected_biounit != actual_biounit_chains
    ):
        raise ValueError(
            "Source mmCIF assembly chains differ from the annotation: "
            f"expected={sorted(selections.expected_biounit)}, "
            f"actual={sorted(actual_biounit_chains)}"
        )

    _require_chains(biounit, selections.system, view_name="system")
    _require_chains(biounit, selections.receptor, view_name="receptor")
    system_mask = np.isin(biounit.chain_id, list(selections.system))
    system_mask |= _water_mask(biounit, annotation, options.system_waters)
    receptor_mask = np.isin(biounit.chain_id, list(selections.receptor))
    receptor_mask |= _water_mask(biounit, annotation, options.receptor_waters)
    return ReconstructedSystem(
        biounit=biounit,
        system=biounit[system_mask],
        receptor=biounit[receptor_mask],
    )


def save_reconstructed_system(
    source_mmcif: Path | str,
    annotation: AnnotationRow,
    *,
    outputs: SystemReconstructionOutputs,
    biounit_chains: pd.DataFrame | None = None,
    options: SystemReconstructionOptions = SystemReconstructionOptions(),
    overwrite: bool = False,
    reconstructed: ReconstructedSystem | None = None,
) -> dict[str, Path]:
    """Reconstruct a system and write exactly the requested mmCIF/FASTA files.

    No PDB files or assembly-rotated ligand SDFs are produced.  Output parents
    are created automatically.  Existing files are rejected unless
    ``overwrite=True``.  A previously reconstructed in-memory view may be
    supplied to avoid parsing and expanding the source assembly again.
    """
    requested: dict[str, Path] = {
        name: Path(path)
        for name, path in {
            "system_cif": outputs.system_cif,
            "receptor_cif": outputs.receptor_cif,
            "sequences_fasta": outputs.sequences_fasta,
        }.items()
        if path is not None
    }
    if not requested:
        raise ValueError("At least one reconstruction output path is required")
    if len(set(requested.values())) != len(requested):
        raise ValueError("Reconstruction output paths must be unique")
    existing = [path for path in requested.values() if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Refusing to overwrite reconstruction outputs: "
            + ", ".join(str(path) for path in existing)
        )

    if reconstructed is None:
        reconstructed = reconstruct_system(
            source_mmcif,
            annotation,
            biounit_chains=biounit_chains,
            options=options,
        )
    system_id = str(annotation.get("system_id", "plinder_system"))
    for path in requested.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    if "system_cif" in requested:
        save_cif_file(reconstructed.system, system_id, requested["system_cif"])
    if "receptor_cif" in requested:
        save_cif_file(
            reconstructed.receptor,
            system_id,
            requested["receptor_cif"],
        )
    if "sequences_fasta" in requested:
        from plinder.data.utils.annotations.cif_utils import (
            get_label_asym_sequences,
            read_mmcif_container,
        )

        asym_to_sequence = get_label_asym_sequences(
            read_mmcif_container(Path(source_mmcif))
        )
        receptor_chain_ids = _select_chains(
            annotation,
            biounit_chains,
            options,
        ).receptor
        missing_sequences = sorted(
            chain_id
            for chain_id in receptor_chain_ids
            if chain_id.split(".", maxsplit=1)[-1] not in asym_to_sequence
        )
        if missing_sequences:
            raise ValueError(
                "Source mmCIF has no entity_poly sequence for reconstructed "
                f"receptor chains {missing_sequences}"
            )
        with requested["sequences_fasta"].open("w") as fasta:
            for chain_id in sorted(receptor_chain_ids):
                asym_id = chain_id.split(".", maxsplit=1)[-1]
                fasta.write(f">{chain_id}\n{asym_to_sequence[asym_id]}\n")
    return requested
