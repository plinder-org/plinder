# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import re
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


def _output_asym_ids(chain_ids: list[str]) -> dict[str, str]:
    """Return unique mmCIF asym IDs while preserving assembly chain IDs."""
    output: dict[str, str] = {}
    used: set[str] = set()
    for chain_id in chain_ids:
        candidate = chain_id
        supported = re.fullmatch(r"(?:[A-Za-z0-9]+|[0-9]+\.[A-Za-z0-9]+)", candidate)
        if supported is None or candidate in used:
            candidate = "C" + chain_id.encode("utf-8").hex()
        suffix = 2
        base = candidate
        while candidate in used:
            candidate = f"{base}X{suffix}"
            suffix += 1
        output[chain_id] = candidate
        used.add(candidate)
    return output


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

    from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

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
    *,
    source_block: pdbx.CIFBlock | None = None,
    source_asym_ids: dict[str, str] | None = None,
    protein_sequences: dict[str, str] | None = None,
) -> None:
    """Save a self-contained structure as PDBx/mmCIF.

    Parameters
    ----------
    atoms : AtomArray
        Atoms to save.
    name : str
        Data block name.
    output_cif_file : str or Path
        Output path.
    source_block : CIFBlock, optional
        Source metadata. Relevant entity, polymer, and chemical-component
        rows are retained for the atoms in the output view.
    source_asym_ids : dict, optional
        Map output chain IDs to label asym IDs in ``source_block``. Assembly
        instance IDs default to their suffix after the first ``.``.
    protein_sequences : dict, optional
        Canonical protein sequences keyed by output chain ID. These are used
        when the source has no polymer metadata.
    """
    if atoms.array_length() == 0:
        raise ValueError("cannot write an empty mmCIF structure")
    atoms = atoms.copy()
    input_chain_ids = list(dict.fromkeys(atoms.chain_id.astype(str)))
    chain_id_map = _output_asym_ids(input_chain_ids)
    source_asym_ids = source_asym_ids or {}
    output_to_source_asym = {
        chain_id_map[chain_id]: source_asym_ids.get(
            chain_id, chain_id.split(".", maxsplit=1)[-1]
        )
        for chain_id in input_chain_ids
    }
    protein_sequences = {
        chain_id_map.get(chain_id, chain_id): sequence
        for chain_id, sequence in (protein_sequences or {}).items()
    }
    atoms.chain_id = np.asarray(
        [chain_id_map[str(chain_id)] for chain_id in atoms.chain_id]
    )
    chain_ids = list(dict.fromkeys(atoms.chain_id.astype(str)))

    source_asym_to_entity: dict[str, str] = {}
    if source_block is not None and "struct_asym" in source_block:
        struct_asym = source_block["struct_asym"]
        if {"id", "entity_id"}.issubset(struct_asym):
            source_asym_to_entity = dict(
                zip(
                    struct_asym["id"].as_array(str),
                    struct_asym["entity_id"].as_array(str),
                )
            )

    chain_to_entity: dict[str, str] = {}
    used_entity_ids: set[str] = set()
    for chain_id in chain_ids:
        source_asym_id = output_to_source_asym[chain_id]
        entity_id = source_asym_to_entity.get(source_asym_id)
        if entity_id is not None:
            chain_to_entity[chain_id] = entity_id
            used_entity_ids.add(entity_id)
    next_entity_id = 1
    for chain_id in chain_ids:
        if chain_id in chain_to_entity:
            continue
        while str(next_entity_id) in used_entity_ids:
            next_entity_id += 1
        entity_id = str(next_entity_id)
        chain_to_entity[chain_id] = entity_id
        used_entity_ids.add(entity_id)
        next_entity_id += 1
    atoms.set_annotation(
        "label_entity_id",
        np.asarray([chain_to_entity[str(chain_id)] for chain_id in atoms.chain_id]),
    )

    source_polymer_metadata_ids: set[str] = set()
    if (
        source_block is not None
        and "entity_poly" in source_block
        and "entity_poly_seq" in source_block
        and "entity_id" in source_block["entity_poly"]
        and "entity_id" in source_block["entity_poly_seq"]
    ):
        source_polymer_metadata_ids = set(
            source_block["entity_poly"]["entity_id"].as_array(str)
        ).intersection(source_block["entity_poly_seq"]["entity_id"].as_array(str))
    for chain_id, sequence in protein_sequences.items():
        if chain_id not in chain_to_entity:
            raise ValueError(f"protein sequence refers to absent chain {chain_id!r}")
        if chain_to_entity[chain_id] in source_polymer_metadata_ids:
            continue
        chain_mask = atoms.chain_id == chain_id
        chain = atoms[chain_mask]
        residue_starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
        residue_count = len(residue_starts) - 1
        if residue_count > len(sequence):
            raise ValueError(
                f"chain {chain_id!r} has {residue_count} resolved residues but "
                f"its supplied sequence has length {len(sequence)}"
            )
        new_residue_ids = np.empty(chain.array_length(), dtype=int)
        for index, (start, stop) in enumerate(
            zip(residue_starts[:-1], residue_starts[1:]), start=1
        ):
            new_residue_ids[start:stop] = index
        atoms.res_id[chain_mask] = new_residue_ids
        atoms.ins_code[chain_mask] = ""

    annotation_categories = set(atoms.get_annotation_categories())
    if "atom_id" not in annotation_categories:
        atoms.set_annotation("atom_id", np.arange(1, atoms.array_length() + 1))
    if "occupancy" not in annotation_categories:
        atoms.set_annotation("occupancy", np.ones(atoms.array_length()))
    if "b_factor" not in annotation_categories:
        atoms.set_annotation("b_factor", np.zeros(atoms.array_length()))
    if "charge" not in annotation_categories:
        atoms.set_annotation("charge", np.zeros(atoms.array_length(), dtype=int))

    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, atoms, data_block=name)
    block = cif_file[name]
    atom_site = block["atom_site"]
    # Keep biological-assembly copies distinct for readers that default to
    # author fields (including Biotite itself).
    atom_site["auth_asym_id"] = atom_site["label_asym_id"].as_array(str)
    if "cell" in block:
        block["cell"]["entry_id"] = [name]
    block["entry"] = pdbx.CIFCategory({"id": [name]})
    block["struct_asym"] = pdbx.CIFCategory(
        {
            "id": chain_ids,
            "entity_id": [chain_to_entity[chain_id] for chain_id in chain_ids],
        }
    )

    source_entity_types: dict[str, str] = {}
    if source_block is not None and "entity" in source_block:
        source_entity = source_block["entity"]
        if {"id", "type"}.issubset(source_entity):
            source_entity_types = dict(
                zip(
                    source_entity["id"].as_array(str),
                    source_entity["type"].as_array(str),
                )
            )
    entity_types: dict[str, str] = {}
    for chain_id in chain_ids:
        entity_id = chain_to_entity[chain_id]
        chain = atoms[atoms.chain_id == chain_id]
        entity_types.setdefault(
            entity_id,
            source_entity_types.get(
                entity_id,
                "polymer"
                if chain_id in protein_sequences
                or np.any(struc.filter_amino_acids(chain))
                or np.any(struc.filter_nucleotides(chain))
                else "water"
                if np.all(struc.filter_solvent(chain))
                else "non-polymer",
            ),
        )
    entity_ids = list(dict.fromkeys(chain_to_entity.values()))
    block["entity"] = pdbx.CIFCategory(
        {
            "id": entity_ids,
            "type": [entity_types[entity_id] for entity_id in entity_ids],
        }
    )

    def copy_source_rows(
        category_name: str,
        key_name: str,
        selected_values: set[str],
    ) -> pdbx.CIFCategory | None:
        if source_block is None or category_name not in source_block:
            return None
        source_category = source_block[category_name]
        if key_name not in source_category:
            return None
        row_mask = np.isin(
            source_category[key_name].as_array(str), list(selected_values)
        )
        if not np.any(row_mask):
            return None
        columns: dict[str, pdbx.CIFColumn] = {}
        for column_name, column in source_category.items():
            mask = column.mask.array[row_mask] if column.mask is not None else None
            columns[column_name] = pdbx.CIFColumn(
                column.data.array[row_mask], mask=mask
            )
        return pdbx.CIFCategory(columns)

    polymer_entity_ids = {
        entity_id
        for entity_id, entity_type in entity_types.items()
        if entity_type == "polymer"
    }
    nonpoly_entity_ids = set(entity_ids).difference(polymer_entity_ids)
    if nonpoly_entity_ids:
        label_seq_id = atom_site["label_seq_id"]
        nonpoly_atom_mask = np.isin(
            atom_site["label_entity_id"].as_array(str),
            list(nonpoly_entity_ids),
        )
        label_seq_mask = (
            label_seq_id.mask.array.copy()
            if label_seq_id.mask is not None
            else np.full(atom_site.row_count, pdbx.MaskValue.PRESENT)
        )
        label_seq_mask[nonpoly_atom_mask] = pdbx.MaskValue.INAPPLICABLE
        atom_site["label_seq_id"] = pdbx.CIFColumn(
            label_seq_id.data.array, mask=label_seq_mask
        )
        if "struct_conn" in block:
            struct_conn = block["struct_conn"]
            nonpoly_chain_ids = {
                chain_id
                for chain_id in chain_ids
                if chain_to_entity[chain_id] in nonpoly_entity_ids
            }
            for partner in (1, 2):
                asym_column = f"ptnr{partner}_label_asym_id"
                seq_column = f"ptnr{partner}_label_seq_id"
                if asym_column not in struct_conn or seq_column not in struct_conn:
                    continue
                seq_id = struct_conn[seq_column]
                seq_mask = (
                    seq_id.mask.array.copy()
                    if seq_id.mask is not None
                    else np.full(struct_conn.row_count, pdbx.MaskValue.PRESENT)
                )
                seq_mask[
                    np.isin(
                        struct_conn[asym_column].as_array(str),
                        list(nonpoly_chain_ids),
                    )
                ] = pdbx.MaskValue.INAPPLICABLE
                struct_conn[seq_column] = pdbx.CIFColumn(
                    seq_id.data.array, mask=seq_mask
                )
    entity_poly = copy_source_rows("entity_poly", "entity_id", polymer_entity_ids)
    entity_poly_seq = copy_source_rows(
        "entity_poly_seq", "entity_id", polymer_entity_ids
    )
    if entity_poly is not None and entity_poly_seq is not None:
        if "pdbx_strand_id" in entity_poly:
            entity_poly["pdbx_strand_id"] = [
                ",".join(
                    chain_id
                    for chain_id in chain_ids
                    if chain_to_entity[chain_id] == entity_id
                )
                for entity_id in entity_poly["entity_id"].as_array(str)
            ]
        block["entity_poly"] = entity_poly
        block["entity_poly_seq"] = entity_poly_seq
    elif polymer_entity_ids:
        from plinder.core.utils.constants import ONE_TO_THREE

        poly_rows: list[tuple[str, str, str]] = []
        sequence_rows: list[tuple[str, str, int, str]] = []
        for entity_id in sorted(polymer_entity_ids):
            entity_chains = [
                chain_id
                for chain_id in chain_ids
                if chain_to_entity[chain_id] == entity_id
            ]
            chain_id = entity_chains[0]
            chain = atoms[atoms.chain_id == chain_id]
            protein_sequence = protein_sequences.get(chain_id)
            residue_starts = struc.get_residue_starts(chain, add_exclusive_stop=False)
            residue_names = chain.res_name[residue_starts].astype(str).tolist()
            if protein_sequence is None:
                if not np.any(struc.filter_amino_acids(chain)):
                    raise ValueError(
                        "source mmCIF metadata is required to write a "
                        f"non-protein polymer chain {chain_id!r}"
                    )
                from plinder.core.utils.constants import THREE_TO_ONE

                protein_sequence = "".join(
                    THREE_TO_ONE.get(residue_name, "X")
                    for residue_name in residue_names
                )
            monomers = [ONE_TO_THREE.get(symbol, "UNK") for symbol in protein_sequence]
            for residue_id, residue_name in zip(
                chain.res_id[residue_starts], residue_names
            ):
                if 1 <= residue_id <= len(monomers):
                    monomers[int(residue_id) - 1] = residue_name
            poly_rows.append((entity_id, protein_sequence, ",".join(entity_chains)))
            sequence_rows.extend(
                (entity_id, monomer, index, "n")
                for index, monomer in enumerate(monomers, start=1)
            )
        block["entity_poly"] = pdbx.CIFCategory(
            {
                "entity_id": [row[0] for row in poly_rows],
                "type": ["polypeptide(L)"] * len(poly_rows),
                "nstd_linkage": ["no"] * len(poly_rows),
                "nstd_monomer": ["yes" if "X" in row[1] else "no" for row in poly_rows],
                "pdbx_seq_one_letter_code": [row[1] for row in poly_rows],
                "pdbx_seq_one_letter_code_can": [row[1] for row in poly_rows],
                "pdbx_strand_id": [row[2] for row in poly_rows],
            }
        )
        block["entity_poly_seq"] = pdbx.CIFCategory(
            {
                "entity_id": [row[0] for row in sequence_rows],
                "mon_id": [row[1] for row in sequence_rows],
                "num": [row[2] for row in sequence_rows],
                "hetero": [row[3] for row in sequence_rows],
            }
        )

    pdbx_entity_nonpoly = copy_source_rows(
        "pdbx_entity_nonpoly", "entity_id", nonpoly_entity_ids
    )
    if pdbx_entity_nonpoly is not None:
        block["pdbx_entity_nonpoly"] = pdbx_entity_nonpoly
    elif nonpoly_entity_ids:
        nonpoly_rows: list[tuple[str, str]] = []
        for entity_id in sorted(nonpoly_entity_ids):
            chain_id = next(
                chain_id
                for chain_id in chain_ids
                if chain_to_entity[chain_id] == entity_id
            )
            chain = atoms[atoms.chain_id == chain_id]
            nonpoly_rows.append((entity_id, str(chain.res_name[0])))
        block["pdbx_entity_nonpoly"] = pdbx.CIFCategory(
            {
                "entity_id": [row[0] for row in nonpoly_rows],
                "name": [row[1] for row in nonpoly_rows],
                "comp_id": [row[1] for row in nonpoly_rows],
            }
        )

    component_ids = set(atoms.res_name.astype(str))
    polymer_component_ids: set[str] = set()
    if "entity_poly_seq" in block:
        polymer_component_ids.update(block["entity_poly_seq"]["mon_id"].as_array(str))
        component_ids.update(polymer_component_ids)
    source_component_types: dict[str, str] = {}
    if source_block is not None and "chem_comp" in source_block:
        source_chem_comp = source_block["chem_comp"]
        if {"id", "type"}.issubset(source_chem_comp):
            source_component_types = dict(
                zip(
                    source_chem_comp["id"].as_array(str),
                    source_chem_comp["type"].as_array(str),
                )
            )
    sorted_component_ids = sorted(component_ids)
    block["chem_comp"] = pdbx.CIFCategory(
        {
            "id": sorted_component_ids,
            "type": [
                source_component_types.get(
                    component_id,
                    "PEPTIDE LINKING"
                    if component_id == "GLY"
                    else "L-PEPTIDE LINKING"
                    if component_id in polymer_component_ids
                    else "NON-POLYMER",
                )
                for component_id in sorted_component_ids
            ],
        }
    )
    block["atom_type"] = pdbx.CIFCategory(
        {"symbol": sorted(set(atoms.element.astype(str)))}
    )
    if "chem_comp_bond" in block:
        bond_order = block["chem_comp_bond"]["value_order"]
        block["chem_comp_bond"]["value_order"] = pdbx.CIFColumn(
            np.char.lower(bond_order.data.array.astype(str)),
            mask=bond_order.mask,
        )
    if "struct_conn" in block:
        connection_types = sorted(
            set(block["struct_conn"]["conn_type_id"].as_array(str))
        )
        block["struct_conn_type"] = pdbx.CIFCategory({"id": connection_types})
    cif_file.write(str(output_cif_file))


def save_reconstructed_chain(
    source_mmcif: Path | str,
    *,
    assembly_id: str,
    chain_instance: str,
    source_asym_id: str,
    output_cif: Path | str,
    structure_id: str,
    superpose_to: struc.AtomArray | None = None,
    overwrite: bool = False,
) -> Path:
    """Write one biological-assembly chain as a self-contained mmCIF.

    If ``superpose_to`` is provided, the reconstructed chain is fitted to that
    protein chain using sequence-matched C-alpha atoms before it is written.
    """
    from plinder.data.annotations.cif_utils import (
        build_biounit,
        get_label_asym_sequences,
        read_mmcif_container,
        read_mmcif_file,
    )

    if not assembly_id:
        raise ValueError("assembly_id must not be empty")
    if not chain_instance:
        raise ValueError("chain_instance must not be empty")
    if not source_asym_id:
        raise ValueError("source_asym_id must not be empty")
    output_cif = Path(output_cif)
    if output_cif.exists() and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite reconstruction output: {output_cif}"
        )

    biounit = build_biounit(read_mmcif_file(source_mmcif), assembly_id)
    chain_mask = biounit.chain_id.astype(str) == chain_instance
    if not np.any(chain_mask):
        available = sorted(set(biounit.chain_id.astype(str)))
        raise ValueError(
            f"Assembly {assembly_id!r} has no chain {chain_instance!r}; "
            f"available chains are {available}"
        )

    source_block = read_mmcif_container(Path(source_mmcif))
    source_sequences = get_label_asym_sequences(source_block)
    protein_sequences = {}
    if source_asym_id in source_sequences:
        protein_sequences[chain_instance] = source_sequences[source_asym_id]
    chain_atoms = biounit[chain_mask]
    if superpose_to is not None:
        from plinder.core.structure.superimpose import superimpose_chain

        chain_atoms, _, _, _ = superimpose_chain(superpose_to, chain_atoms)
    output_cif.parent.mkdir(parents=True, exist_ok=True)
    save_cif_file(
        chain_atoms,
        structure_id,
        output_cif,
        source_block=source_block,
        source_asym_ids={chain_instance: source_asym_id},
        protein_sequences=protein_sequences,
    )
    return output_cif


def reconstruct_interface(
    source_mmcif: Path | str,
    annotation: AnnotationRow,
) -> struc.AtomArray:
    """Rebuild one annotated protein-interface chain pair from a source mmCIF."""
    from plinder.data.annotations.cif_utils import build_biounit, read_mmcif_file

    assembly_id = str(annotation.get("system_biounit_id", ""))
    chain_1 = str(annotation.get("interface_chain_1", ""))
    chain_2 = str(annotation.get("interface_chain_2", ""))
    if not assembly_id:
        raise ValueError("interface annotation is missing system_biounit_id")
    if not chain_1 or not chain_2:
        raise ValueError("interface annotation is missing one or both chains")
    if chain_1 == chain_2:
        raise ValueError("an interface requires two distinct chain instances")

    biounit = build_biounit(read_mmcif_file(source_mmcif), assembly_id)
    selected_chains = {chain_1, chain_2}
    _require_chains(biounit, selected_chains, view_name="interface")
    return biounit[np.isin(biounit.chain_id, list(selected_chains))]


def save_reconstructed_interface(
    source_mmcif: Path | str,
    annotation: AnnotationRow,
    *,
    output_cif: Path | str,
    overwrite: bool = False,
    reconstructed: struc.AtomArray | None = None,
) -> Path:
    """Write one annotated biological-assembly interface as a valid mmCIF."""
    from plinder.data.annotations.cif_utils import (
        get_label_asym_sequences,
        read_mmcif_container,
    )

    output_cif = Path(output_cif)
    if output_cif.exists() and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite reconstruction output: {output_cif}"
        )
    atoms = (
        reconstruct_interface(source_mmcif, annotation)
        if reconstructed is None
        else reconstructed
    )
    chain_ids = [
        str(annotation.get("interface_chain_1", "")),
        str(annotation.get("interface_chain_2", "")),
    ]
    source_asym_ids = {
        chain_id: chain_id.split(".", maxsplit=1)[-1] for chain_id in chain_ids
    }
    source_block = read_mmcif_container(Path(source_mmcif))
    source_sequences = get_label_asym_sequences(source_block)
    protein_sequences = {
        chain_id: source_sequences[source_asym_id]
        for chain_id, source_asym_id in source_asym_ids.items()
        if source_asym_id in source_sequences
    }
    structure_id = str(annotation.get("system_id", "interface"))
    output_cif.parent.mkdir(parents=True, exist_ok=True)
    save_cif_file(
        atoms,
        structure_id,
        output_cif,
        source_block=source_block,
        source_asym_ids=source_asym_ids,
        protein_sequences=protein_sequences,
    )
    return output_cif


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
    from plinder.data.annotations.cif_utils import (
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
    from plinder.data.annotations.cif_utils import (
        get_label_asym_sequences,
        read_mmcif_container,
    )

    source_block = read_mmcif_container(Path(source_mmcif))
    source_sequences = get_label_asym_sequences(source_block)

    def output_protein_sequences(atoms: struc.AtomArray) -> dict[str, str]:
        return {
            str(chain_id): source_sequences[source_asym_id]
            for chain_id in np.unique(atoms.chain_id)
            for source_asym_id in [str(chain_id).split(".", maxsplit=1)[-1]]
            if source_asym_id in source_sequences
        }

    if "system_cif" in requested:
        save_cif_file(
            reconstructed.system,
            system_id,
            requested["system_cif"],
            source_block=source_block,
            protein_sequences=output_protein_sequences(reconstructed.system),
        )
    if "receptor_cif" in requested:
        save_cif_file(
            reconstructed.receptor,
            system_id,
            requested["receptor_cif"],
            source_block=source_block,
            protein_sequences=output_protein_sequences(reconstructed.receptor),
        )
    if "sequences_fasta" in requested:
        asym_to_sequence = source_sequences
        receptor_chain_ids = _select_chains(
            annotation,
            biounit_chains,
            options,
        ).receptor
        output_chain_ids = _output_asym_ids(
            list(dict.fromkeys(reconstructed.receptor.chain_id.astype(str)))
        )
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
                fasta.write(
                    f">{output_chain_ids[chain_id]}\n" f"{asym_to_sequence[asym_id]}\n"
                )
    return requested
