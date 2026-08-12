# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""mmCIF I/O utilities using biotite.

Generic helpers for reading CIF blocks, extracting scalar values and
category rows, plus ligand bond-order detection and assignment from
SMILES templates.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from math import prod
from pathlib import Path
from typing import TypedDict

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from rdkit import Chem

from plinder.core.structure.smallmols_utils import (
    mol_assigned_bond_orders_by_template,
)

LOG = logging.getLogger(__name__)


class EntryTaxonomy(TypedDict):
    """Taxonomy fields extracted from deposited entry metadata."""

    source_taxonomy_ids: list[int]
    source_organism_names: list[str]
    host_taxonomy_ids: list[int]
    host_organism_names: list[str]


# Single source of truth lives in ``plinder.core.structure.atoms`` so
# both ``plinder.core`` and ``plinder.data`` filter H/D/T isotopes
# consistently.
from plinder.core.structure.atoms import is_hydrogen_isotope  # noqa: E402

# ---------------------------------------------------------------------------
# Generic CIF I/O helpers
# ---------------------------------------------------------------------------


def read_mmcif_file(mmcif_filename: Path | str) -> pdbx.CIFFile:
    """Read an mmCIF file and add unambiguous optional atom-site defaults."""
    import gzip

    path = str(mmcif_filename)
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as f:
            cif_file = pdbx.CIFFile.read(f)
    else:
        cif_file = pdbx.CIFFile.read(path)
    for block in cif_file.values():
        if "atom_site" not in block:
            continue
        atom_site = block["atom_site"]
        atom_count = atom_site.row_count
        if "pdbx_PDB_model_num" not in atom_site:
            atom_site["pdbx_PDB_model_num"] = np.ones(atom_count, dtype=np.int32)
        if "pdbx_PDB_ins_code" not in atom_site:
            atom_site["pdbx_PDB_ins_code"] = ["."] * atom_count
    return cif_file


def check_custom_mmcif_fields(
    block: pdbx.CIFBlock,
    *,
    source: Path,
    structure_mode: str,
    require_label_ids: bool = False,
) -> None:
    """Check fields required to read custom coordinates and assemblies.

    Model number and insertion code are optional for a single-model custom
    structure and receive unambiguous defaults. Deposited-PDB mode additionally
    requires the assembly operation tables. Other deposition metadata (authors,
    citations, experimental details, and validation categories) is not required.
    """
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
    if require_label_ids:
        label_columns = {
            "label_atom_id",
            "label_comp_id",
            "label_asym_id",
            "label_seq_id",
        }
        missing.extend(sorted(label_columns.difference(atom_site)))
        missing_identifiers: list[str] = []
    else:
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

    atom_count = atom_site.row_count
    if "pdbx_PDB_model_num" not in atom_site:
        atom_site["pdbx_PDB_model_num"] = np.ones(atom_count, dtype=np.int32)
    if "pdbx_PDB_ins_code" not in atom_site:
        atom_site["pdbx_PDB_ins_code"] = ["."] * atom_count

    if structure_mode != "pdb":
        if structure_mode != "as_is":
            raise ValueError("structure_mode must be 'as_is' or 'pdb'")
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


def read_mmcif_container(mmcif_filename: Path) -> pdbx.CIFBlock:
    """Parse mmcif file and return the first data block."""
    cif_file = read_mmcif_file(mmcif_filename)
    return list(cif_file.values())[0]


def _alphabetic_id(index: int) -> str:
    """Return a deterministic spreadsheet-style alphabetic identifier."""
    value = index + 1
    result = ""
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


@contextmanager
def _alphabetic_altloc_ids(
    cif_file: pdbx.CIFFile | pdbx.CIFBlock,
) -> Iterator[dict[str, str]]:
    """Temporarily map non-alphabetic altloc IDs for Biotite filtering."""
    block = (
        cif_file if isinstance(cif_file, pdbx.CIFBlock) else list(cif_file.values())[0]
    )
    specifications = [
        ("atom_site", "label_alt_id"),
        ("struct_conn", "pdbx_ptnr1_label_alt_id"),
        ("struct_conn", "pdbx_ptnr2_label_alt_id"),
    ]
    columns: list[tuple[pdbx.CIFCategory, str, pdbx.CIFColumn]] = []
    all_ids: set[str] = set()
    for category_name, column_name in specifications:
        if category_name not in block or column_name not in block[category_name]:
            continue
        category = block[category_name]
        column = category[column_name]
        columns.append((category, column_name, column))
        all_ids.update(column.as_array(str))

    no_altloc = {".", "?", " ", ""}
    invalid_ids = sorted(
        altloc_id
        for altloc_id in all_ids
        if altloc_id not in no_altloc and not altloc_id.isalpha()
    )
    if not invalid_ids:
        yield {}
        return

    used_ids = {altloc_id for altloc_id in all_ids if altloc_id.isalpha()}
    mapping: dict[str, str] = {}
    candidate_index = 0
    for invalid_id in invalid_ids:
        while (candidate := _alphabetic_id(candidate_index)) in used_ids:
            candidate_index += 1
        mapping[invalid_id] = candidate
        used_ids.add(candidate)
        candidate_index += 1

    try:
        for category, column_name, original in columns:
            values = original.as_array(str)
            mapped = [mapping.get(value, value) for value in values]
            category[column_name] = pdbx.CIFColumn(mapped, mask=original.mask)
        yield mapping
    finally:
        for category, column_name, original in columns:
            category[column_name] = original


def get_structure_with_altloc(
    cif_file: pdbx.CIFFile | pdbx.CIFBlock,
    *,
    model: int = 1,
    use_author_fields: bool = False,
    include_bonds: bool = False,
    extra_fields: list[str] | None = None,
) -> struc.AtomArray:
    """Load one model using its deposited-first alternate conformers.

    Biotite only recognizes alphabetic alternate-location IDs while filtering.
    Non-alphabetic source IDs are therefore mapped temporarily, then restored
    on the returned ``selected_altloc_id`` annotation so validation can select
    the exact same deposited conformer.
    """
    requested_extra_fields = list(extra_fields or [])
    include_label_alt_id = "label_alt_id" in requested_extra_fields
    if not include_label_alt_id:
        requested_extra_fields.append("label_alt_id")
    with _alphabetic_altloc_ids(cif_file) as mapping:
        atoms = pdbx.get_structure(
            cif_file,
            model=model,
            altloc="first",
            use_author_fields=use_author_fields,
            include_bonds=include_bonds,
            extra_fields=requested_extra_fields,
        )
    if not isinstance(atoms, struc.AtomArray):
        raise TypeError("loading one mmCIF model must return an AtomArray")
    inverse_mapping = {mapped: source for source, mapped in mapping.items()}
    source_altloc_ids = np.asarray(
        [inverse_mapping.get(value, value) for value in atoms.label_alt_id]
    )
    # ``label_alt_id`` is atom-level: atoms shared by all conformers retain
    # ``.``.  Record the selected source conformer on every atom in the
    # residue so the choice survives subsequent hydrogen removal and can be
    # used to select the matching validation record.
    selected_altloc_ids = np.full(len(atoms), ".", dtype=source_altloc_ids.dtype)
    residue_starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    for start, stop in zip(residue_starts[:-1], residue_starts[1:]):
        selected = next(
            (
                altloc_id
                for altloc_id in source_altloc_ids[start:stop]
                if altloc_id not in {".", "?", " ", ""}
            ),
            ".",
        )
        selected_altloc_ids[start:stop] = selected
    atoms.set_annotation("selected_altloc_id", selected_altloc_ids)
    if not include_label_alt_id:
        atoms.del_annotation("label_alt_id")
    return atoms


def get_unit_cell_with_altloc(
    cif_file: pdbx.CIFFile,
    *,
    model: int = 1,
    use_author_fields: bool = False,
) -> struc.AtomArray:
    """Build a unit cell using normalized deposited-first altlocs."""
    with _alphabetic_altloc_ids(cif_file):
        atoms = pdbx.get_unit_cell(
            cif_file,
            model=model,
            altloc="first",
            use_author_fields=use_author_fields,
        )
    if not isinstance(atoms, struc.AtomArray):
        raise TypeError("loading one mmCIF unit-cell model must return an AtomArray")
    return atoms


def get_label_asym_sequences(block: pdbx.CIFBlock) -> dict[str, str]:
    """Extract polymer sequences keyed by label asym ID.

    ``entity_poly.pdbx_strand_id`` contains author chain IDs and therefore
    cannot be used for reconstructed biological assemblies, whose chain IDs
    are ``<instance>.<label_asym_id>``.  The stable mapping is
    ``struct_asym.id -> struct_asym.entity_id -> entity_poly.entity_id``.
    """
    if "struct_asym" not in block or "entity_poly" not in block:
        return {}
    struct_asym = block["struct_asym"]
    entity_poly = block["entity_poly"]
    if not {"id", "entity_id"}.issubset(struct_asym):
        return {}
    if not {"entity_id", "pdbx_seq_one_letter_code_can"}.issubset(entity_poly):
        return {}

    entity_sequences = {
        str(entity_id): "".join(str(sequence).replace(";", "").split())
        for entity_id, sequence in zip(
            entity_poly["entity_id"].as_array(),
            entity_poly["pdbx_seq_one_letter_code_can"].as_array(),
        )
    }
    return {
        str(asym_id): entity_sequences[str(entity_id)]
        for asym_id, entity_id in zip(
            struct_asym["id"].as_array(),
            struct_asym["entity_id"].as_array(),
        )
        if str(entity_id) in entity_sequences
    }


def get_entry_taxonomy(
    block: pdbx.CIFBlock,
) -> EntryTaxonomy:
    """Extract distinct source and expression-host organisms for an entry."""
    source_taxonomy_ids: set[int] = set()
    source_organism_names: set[str] = set()
    host_taxonomy_ids: set[int] = set()
    host_organism_names: set[str] = set()

    def add_values(
        category_name: str,
        taxonomy_column: str,
        organism_column: str,
        taxonomy_ids: set[int],
        organism_names: set[str],
    ) -> None:
        if category_name not in block:
            return
        category = block[category_name]
        if taxonomy_column in category:
            for value in category[taxonomy_column].as_array(str):
                text = str(value).strip()
                if text in {"", ".", "?"}:
                    continue
                try:
                    taxonomy_ids.add(int(text))
                except ValueError:
                    LOG.warning(
                        "ignoring non-integer %s.%s value %r",
                        category_name,
                        taxonomy_column,
                        text,
                    )
        if organism_column in category:
            organism_names.update(
                text
                for value in category[organism_column].as_array(str)
                if (text := str(value).strip()) not in {"", ".", "?"}
            )

    add_values(
        "entity_src_gen",
        "pdbx_gene_src_ncbi_taxonomy_id",
        "pdbx_gene_src_scientific_name",
        source_taxonomy_ids,
        source_organism_names,
    )
    add_values(
        "entity_src_nat",
        "pdbx_ncbi_taxonomy_id",
        "pdbx_organism_scientific",
        source_taxonomy_ids,
        source_organism_names,
    )
    add_values(
        "pdbx_entity_src_syn",
        "ncbi_taxonomy_id",
        "organism_scientific",
        source_taxonomy_ids,
        source_organism_names,
    )
    add_values(
        "entity_src_gen",
        "pdbx_host_org_ncbi_taxonomy_id",
        "pdbx_host_org_scientific_name",
        host_taxonomy_ids,
        host_organism_names,
    )
    return {
        "source_taxonomy_ids": sorted(source_taxonomy_ids),
        "source_organism_names": sorted(source_organism_names),
        "host_taxonomy_ids": sorted(host_taxonomy_ids),
        "host_organism_names": sorted(host_organism_names),
    }


def get_mmcif_revision(block: pdbx.CIFBlock) -> tuple[int, int]:
    """Return the latest structure-model major/minor revision in an mmCIF."""
    category_name = "pdbx_audit_revision_history"
    if category_name not in block:
        raise ValueError(f"mmCIF has no {category_name} category")
    category = block[category_name]
    required = {"major_revision", "minor_revision"}
    if not required.issubset(category):
        raise ValueError(f"mmCIF {category_name} has no revision columns")

    major = category["major_revision"].as_array()
    minor = category["minor_revision"].as_array()
    rows = list(range(len(major)))
    if "data_content_type" in category:
        content_type = category["data_content_type"].as_array()
        rows = [
            index
            for index, value in enumerate(content_type)
            if str(value).lower() == "structure model"
        ]
    if not rows:
        raise ValueError(f"mmCIF {category_name} has no structure-model revisions")
    return max((int(major[index]), int(minor[index])) for index in rows)


def get_model_count(cif_file: pdbx.CIFFile) -> int:
    """Return the number of models in a CIF (1 if no model column present)."""
    block = list(cif_file.values())[0]
    if "atom_site" not in block:
        return 0
    atom_site = block["atom_site"]
    if "pdbx_PDB_model_num" not in atom_site:
        return 1
    return int(len(set(atom_site["pdbx_PDB_model_num"].as_array())))


def _operation_expression_count(expression: str) -> int:
    """Count transformations encoded by an assembly operation expression."""
    groups = re.findall(r"\(([^()]*)\)", expression)
    if not groups:
        groups = [expression]
    group_sizes = []
    for group in groups:
        size = 0
        for token in group.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                first, last = token.split("-", maxsplit=1)
                size += int(last) - int(first) + 1
            else:
                size += 1
        if size == 0:
            raise ValueError(f"empty assembly operation group in {expression!r}")
        group_sizes.append(size)
    return prod(group_sizes)


def get_legacy_chain_instance_mapping(
    cif_file: pdbx.CIFFile | pdbx.CIFBlock,
    assembly_id: str,
) -> dict[str, str]:
    """Map canonical per-asym instance IDs to historical global-copy IDs.

    The canonical Biotite representation numbers transformed copies separately
    for each label asym ID.  Earlier Plinder releases used OpenStructure, which
    numbered every generated transformation globally across assembly rows.
    This mapping preserves those historical identifiers as lookup aliases.
    """
    block = (
        cif_file if isinstance(cif_file, pdbx.CIFBlock) else list(cif_file.values())[0]
    )
    category_name = "pdbx_struct_assembly_gen"
    if category_name not in block:
        return {}
    category = block[category_name]
    required = {"assembly_id", "oper_expression", "asym_id_list"}
    if not required.issubset(category):
        return {}

    mapping: dict[str, str] = {}
    per_asym_offset: defaultdict[str, int] = defaultdict(int)
    global_offset = 0
    for row_assembly_id, expression, asym_id_list in zip(
        category["assembly_id"].as_array(str),
        category["oper_expression"].as_array(str),
        category["asym_id_list"].as_array(str),
    ):
        if str(row_assembly_id) != str(assembly_id):
            continue
        operation_count = _operation_expression_count(str(expression))
        asym_ids = [value.strip() for value in str(asym_id_list).split(",")]
        for asym_id in asym_ids:
            if not asym_id:
                continue
            local_offset = per_asym_offset[asym_id]
            for operation_index in range(operation_count):
                canonical = f"{local_offset + operation_index + 1}.{asym_id}"
                legacy = f"{global_offset + operation_index + 1}.{asym_id}"
                mapping[canonical] = legacy
            per_asym_offset[asym_id] += operation_count
        global_offset += operation_count
    return mapping


def build_biounit(
    cif_file: pdbx.CIFFile,
    assembly_id: str,
) -> struc.AtomArray:
    """Build a biological assembly with stable ``instance.asym`` chain IDs.

    Biotite's ``sym_id`` enumerates transformed copies independently for each
    source asym chain.  Using it avoids assuming that an assembly consists of
    complete, contiguous ASU-sized blocks, which is false when operators apply
    to only a subset of chains.
    """
    with _alphabetic_altloc_ids(cif_file):
        biounit = pdbx.get_assembly(
            cif_file,
            assembly_id=assembly_id,
            model=1,
            altloc="first",
            use_author_fields=False,
            include_bonds=True,
            extra_fields=["label_asym_id"],
        )
        biounit = biounit[~is_hydrogen_isotope(biounit.element)]
        if biounit.bonds is None:
            raise ValueError(
                f"assembly {assembly_id}: biotite returned no bonds despite "
                "include_bonds=True"
            )
        biounit.chain_id = np.asarray(
            [
                f"{int(sym_id) + 1}.{asym_id}"
                for sym_id, asym_id in zip(biounit.sym_id, biounit.label_asym_id)
            ]
        )
        legacy_mapping = get_legacy_chain_instance_mapping(cif_file, assembly_id)
        biounit.set_annotation(
            "legacy_chain_id",
            np.asarray(
                [
                    legacy_mapping.get(chain_id, chain_id)
                    for chain_id in biounit.chain_id
                ]
            ),
        )
        apply_struct_conn_bonds(biounit, list(cif_file.values())[0])
    return biounit


def _cif_scalar(block: pdbx.CIFBlock, category: str, column: str) -> str | None:
    """Read a single scalar value from a CIF category, or None."""
    if category not in block:
        return None
    cat = block[category]
    if column not in cat:
        return None
    val = cat[column].as_array()[0]
    if val in ("?", "."):
        return None
    return str(val)


def _iter_category_rows(
    block: pdbx.CIFBlock, category: str, columns: list[str]
) -> list[dict[str, str]]:
    """Iterate over rows of a CIF category as dicts."""
    if category not in block:
        return []
    cat = block[category]
    arrays = {}
    for col in columns:
        if col not in cat:
            return []
        arrays[col] = cat[col].as_array()
    n = len(next(iter(arrays.values())))
    return [{col: arrays[col][i] for col in columns} for i in range(n)]


def get_entry_info(data: pdbx.CIFBlock) -> dict[str, str | None]:
    """Get entry-level information from a CIF block.

    Parameters
    ----------
    data : pdbx.CIFBlock
    Returns
    -------
    dict[str, str | None]
    """
    entry_info = {}
    mappings = [
        ("entry_oligomeric_state", "pdbx_struct_assembly", "oligomeric_details"),
        ("entry_determination_method", "exptl", "method"),
        ("entry_keywords", "struct_keywords", "pdbx_keywords"),
        ("entry_pH", "exptl_crystal_grow", "pH"),
    ]
    for key, cat_name, col_name in mappings:
        entry_info[key] = _cif_scalar(data, cat_name, col_name)
    resolution_options = [
        ("refine", "ls_d_res_high"),
        # ("em_3d_reconstruction", "resolution"), # TODO: add this back for next annotation rerun
    ]
    resolution = None
    for cat_name, col_name in resolution_options:
        r = _cif_scalar(data, cat_name, col_name)
        if r is not None:
            resolution = r
            break
    entry_info["entry_resolution"] = resolution
    return entry_info


def get_chain_external_mappings(
    data: pdbx.CIFBlock,
) -> dict[str, dict[str, dict[str, list[tuple[str, str] | None]]]]:
    """Get additional metadata directory from nextgen mmcif."""
    per_chain: dict[str, dict[str, dict[str, set[tuple[str, str] | None]]]] = {}

    # SIFTS mapping
    for row in _iter_category_rows(
        data,
        "pdbx_sifts_xref_db_segments",
        ["asym_id", "xref_db", "xref_db_acc", "seq_id_start", "seq_id_end"],
    ):
        if row["asym_id"] not in per_chain:
            per_chain[row["asym_id"]] = defaultdict(lambda: defaultdict(set))
        per_chain[row["asym_id"]][row["xref_db"]][row["xref_db_acc"]].add(
            (
                row["seq_id_start"],
                row["seq_id_end"],
            )
        )

    # UniProt mapping
    for row in _iter_category_rows(
        data,
        "pdbx_sifts_unp_segments",
        ["asym_id", "unp_acc", "seq_id_start", "seq_id_end"],
    ):
        if row["asym_id"] not in per_chain:
            per_chain[row["asym_id"]] = defaultdict(lambda: defaultdict(set))
        per_chain[row["asym_id"]]["UniProt"][row["unp_acc"]].add(
            (
                row["seq_id_start"],
                row["seq_id_end"],
            )
        )

    # BIRD entries with PRD codes
    for row in _iter_category_rows(data, "pdbx_molecule", ["asym_id"]):
        if row["asym_id"] not in per_chain:
            per_chain[row["asym_id"]] = defaultdict(lambda: defaultdict(set))
        per_chain[row["asym_id"]]["BIRD"][row["asym_id"]].add(None)

    per_chain_list: dict[str, dict[str, dict[str, list[tuple[str, str] | None]]]] = {}
    for chain in per_chain:
        per_chain_list[chain] = {}
        for mapping in per_chain[chain]:
            per_chain_list[chain][mapping] = {
                k: list(v) for k, v in per_chain[chain][mapping].items()
            }
    return per_chain_list


# ---------------------------------------------------------------------------
# CIF -> RDKit conversion
# ---------------------------------------------------------------------------


def atoms_to_rdkit_mol(
    atoms: "struc.AtomArray",
    assign_stereo: bool = True,
) -> "Chem.Mol":
    """Convert a biotite AtomArray to a sanitized RDKit Mol.

    Stereochemistry is, optionally, assigned from 3D coordinates before
    the final ``RemoveAllHs`` so chiral tags are stamped on heavy atoms
    and survive hydrogen removal.

    Parameters
    ----------
    atoms : AtomArray
        Atoms with bonds (e.g. from ``include_bonds=True`` or set
        explicitly by the caller). Multi-atom inputs must carry
        bonds — missing / empty bonds raise ``ValueError``. Single
        atoms (ions) are allowed to have no bonds.
    assign_stereo : bool
        If True, call ``AssignAtomChiralTagsFromStructure`` on the result.

    Returns
    -------
    Chem.Mol
        Sanitized RDKit molecule with 3D coordinates and PDB atom info,
        heavy atoms only.

    Raises
    ------
    ValueError
        If the input has no bonds, or if RDKit conversion fails.

    Notes
    -----
    Hydrogen atoms *and isotopes* (D, T) are removed. biotite's
    ``element`` is a string, so a naive ``element != "H"`` filter would
    leak deuterium/tritium into the mol; we pre-filter the common
    mass-1 isotopes and additionally call ``RemoveAllHs`` as a
    belt-and-braces catch for anything RDKit still classifies as
    hydrogen via atomic number.

    Warnings
    --------
    The input **must carry bonds** (``atoms.bonds`` non-empty for
    multi-atom inputs). Callers are expected to have either loaded the
    CIF with ``include_bonds=True`` (which reads ``_chem_comp_bond``
    and ``_struct_conn``) or to have populated bonds themselves. The
    function will not re-derive bonds via
    ``connect_via_residue_names`` because that fallback silently drops
    inter-residue peptide bonds for non-standard residues in
    multi-residue ligands — better to fail loudly than hand back a
    structurally-wrong mol.
    """
    from biotite.interface import rdkit as rdkit_interface
    from biotite.structure import BondList

    # TODO(peppr): temporary local sanitize carrying boron/main-group
    # over-valence fixes not yet in a released peppr. Revert to
    # `from peppr import sanitize as peppr_sanitize` once upstream.
    # See plinder.core.utils.sanitize.
    from plinder.core.utils.sanitize import sanitize as peppr_sanitize

    heavy = atoms[~is_hydrogen_isotope(atoms.element)]

    # Multi-atom inputs must carry bonds; single atoms (ions) don't need any.
    if heavy.bonds is None or heavy.bonds.as_array().shape[0] == 0:
        if heavy.array_length() == 1:
            # add empty bondlist for single atoms to convert to RDKit mol
            heavy.bonds = BondList(1)
        else:
            raise ValueError(
                "atoms_to_rdkit_mol requires bonds on multi-atom inputs. "
                "Load the CIF with include_bonds=True (which parses "
                "_chem_comp_bond + _struct_conn) or populate atoms.bonds "
                "before calling. A connect_via_residue_names fallback was "
                "removed because it silently drops inter-residue peptide "
                "bonds for non-standard residues in multi-residue ligands."
            )

    mol = rdkit_interface.to_mol(heavy, kekulize=True, use_dative_bonds=True)
    if mol is None:
        raise ValueError("Failed to convert AtomArray to RDKit Mol")

    peppr_sanitize(mol)
    if assign_stereo:
        Chem.AssignAtomChiralTagsFromStructure(mol)
    # RDKit's RemoveAllHs keys on atomic number, so it strips any
    # hydrogen isotope atom that survived the element-string filter.
    # Safe after stereo assignment — chiral tags live on heavy atoms.
    #
    # sanitize=False: peppr.sanitize already sanitized the mol and deliberately
    # tolerates over-valent main-group centres (boron cages, Be, …) that RDKit's
    # valence check rejects. RemoveAllHs re-sanitizes by default, which would
    # re-raise AtomValenceException for those molecules and undo peppr's
    # tolerance — so strip Hs without re-validating valence.
    return Chem.RemoveAllHs(mol, sanitize=False)


# ---------------------------------------------------------------------------
# CIF ligand parsing
# ---------------------------------------------------------------------------


def parse_struct_conn(
    block: pdbx.CIFBlock,
) -> list[dict[str, str]]:
    """Parse ``_struct_conn`` into a list of connection dicts."""
    if "struct_conn" not in block:
        return []
    conn = block["struct_conn"]
    cols = {
        "conn_type_id": "conn_type",
        "ptnr1_label_asym_id": "chain1",
        "ptnr1_label_seq_id": "seq1",
        "ptnr1_label_atom_id": "atom1",
        "ptnr1_label_comp_id": "comp1",
        "ptnr2_label_asym_id": "chain2",
        "ptnr2_label_seq_id": "seq2",
        "ptnr2_label_atom_id": "atom2",
        "ptnr2_label_comp_id": "comp2",
        "ptnr1_auth_seq_id": "auth_seq1",
        "ptnr2_auth_seq_id": "auth_seq2",
    }
    arrays = {}
    for cif_col, key in cols.items():
        if cif_col not in conn:
            return []
        arrays[key] = conn[cif_col].as_array()
    n = len(arrays["conn_type"])
    return [{k: arrays[k][i] for k in arrays} for i in range(n)]


def apply_struct_conn_bonds(
    atoms: "struc.AtomArray",
    block: pdbx.CIFBlock,
) -> None:
    """Add inter-residue covalent bonds from ``_struct_conn`` in-place."""

    parsed_connections = parse_struct_conn(block)
    connections = []
    partner_keys: set[tuple[str, int, str]] = set()
    for connection in parsed_connections:
        if connection["conn_type"] != "covale":
            continue
        try:
            res_id1 = int(connection["seq1"]) if connection["seq1"] != "." else -1
            res_id2 = int(connection["seq2"]) if connection["seq2"] != "." else -1
        except ValueError:
            continue
        key1 = (connection["chain1"], res_id1, connection["atom1"])
        key2 = (connection["chain2"], res_id2, connection["atom2"])
        connections.append((key1, key2))
        partner_keys.update((key1, key2))
    if not connections:
        return
    if atoms.bonds is None:
        atoms.bonds = struc.BondList(atoms.array_length())

    categories = set(atoms.get_annotation_categories())
    if "label_asym_id" in categories:
        label_ids = atoms.get_annotation("label_asym_id")
    else:
        label_ids = np.asarray(
            [str(chain_id).split(".", maxsplit=1)[-1] for chain_id in atoms.chain_id]
        )
    if "sym_id" in categories:
        instance_ids = atoms.get_annotation("sym_id")
    else:
        chain_parts = [
            str(chain_id).split(".", maxsplit=1) for chain_id in atoms.chain_id
        ]
        instance_ids = np.asarray(
            [parts[0] if len(parts) == 2 else "" for parts in chain_parts]
        )

    candidate_mask = (
        np.isin(label_ids, [key[0] for key in partner_keys])
        & np.isin(atoms.res_id, [key[1] for key in partner_keys])
        & np.isin(atoms.atom_name, [key[2] for key in partner_keys])
    )
    partner_indices: dict[tuple[str, int, str], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for index in np.flatnonzero(candidate_mask):
        key = (
            str(label_ids[index]),
            int(atoms.res_id[index]),
            str(atoms.atom_name[index]),
        )
        if key in partner_keys:
            partner_indices[key][str(instance_ids[index])].append(int(index))

    for key1, key2 in connections:
        partner1 = partner_indices.get(key1, {})
        partner2 = partner_indices.get(key2, {})
        # A source row applies within each transformed copy, not to the
        # Cartesian product of all biological-assembly copies.
        for instance_id in set(partner1) & set(partner2):
            for index1 in partner1[instance_id]:
                bonded_indices, _ = atoms.bonds.get_bonds(index1)
                existing_neighbors = set(int(index) for index in bonded_indices)
                for index2 in partner2[instance_id]:
                    if index2 not in existing_neighbors:
                        atoms.bonds.add_bond(
                            index1,
                            index2,
                            struc.BondType.SINGLE,
                        )
                        existing_neighbors.add(index2)


# ---------------------------------------------------------------------------
# Ligand bond order detection and assignment
# ---------------------------------------------------------------------------


class MissingBondOrderError(ValueError):
    """Raised when a CIF file has ligands with unresolvable bond orders."""

    pass


# Minimum fraction of CCD heavy atoms that must be present in a CIF
# for connect_via_residue_names to produce reliable bonds.
_MIN_CCD_ATOM_OVERLAP = 0.5


def _get_hetatm_comp_ids(block: pdbx.CIFBlock) -> set[str]:
    """Extract non-polymer component IDs from atom_site."""
    if "atom_site" not in block:
        return set()
    atom_site = block["atom_site"]
    group_pdb = atom_site["group_PDB"].as_array()
    comp_ids = atom_site["label_comp_id"].as_array()
    return {comp_ids[i] for i in range(len(group_pdb)) if group_pdb[i] == "HETATM"}


def _get_cif_bond_comp_ids(block: pdbx.CIFBlock) -> set[str]:
    """Return the set of comp_ids that already have _chem_comp_bond entries."""
    if "chem_comp_bond" not in block:
        return set()
    return set(block["chem_comp_bond"]["comp_id"].as_array())


def _is_known_compound(comp_id: str, atom_names: set[str] | None = None) -> bool:
    """Check if a component ID is known to the CCD compound library.

    "Known" resolves via :func:`_get_ccd_atomarray`, so a code present only in
    the downloaded ``components.cif`` (not biotite's bundled CCD) counts as
    known — correct for reference SMILES/stereo. Note the biounit bond path
    (biotite ``connect_via_residue_names``) is still bundled-CCD-based, so such
    a code only gets its intra-residue bonds if the CIF carries
    ``_chem_comp_bond`` (deposited entries always do) or via the components.cif
    bond fallback in :func:`ligand_utils._fill_missing_ccd_bonds`.

    If *atom_names* is provided, also verify that the CIF atom names
    overlap with the CCD entry. Bond assignment via
    ``connect_via_residue_names`` relies on atom-name matching, so a
    compound whose names don't match CCD will get wrong bonds even if
    the comp_id exists in the dictionary (e.g. Boltz ``LIG`` =/= CCD
    ``LIG``).
    """
    try:
        from plinder.data.annotations.ligand_utils import _get_ccd_atomarray

        ref = _get_ccd_atomarray(comp_id)
        if ref is None:
            return False
        if atom_names is not None:
            ref_heavy = ref[~is_hydrogen_isotope(ref.element)]
            ref_names = set(ref_heavy.atom_name)
            if not ref_names or not atom_names:
                return False
            # All CIF atom names must exist in the CCD entry
            unknown_names = atom_names - ref_names
            if unknown_names:
                return False
            # Enough CCD atoms must be present for reliable bond assignment
            if len(atom_names & ref_names) < _MIN_CCD_ATOM_OVERLAP * len(ref_names):
                return False
        return True
    except Exception as e:
        LOG.warning(f"CCD lookup failed for {comp_id}: {e}")
        return False


def get_unknown_ligand_ids(cif_input: pdbx.CIFFile | Path | str) -> set[str]:
    """Return HETATM comp_ids not in CCD and missing ``_chem_comp_bond``.

    Parameters
    ----------
    cif_input : CIFFile, Path, or str
        Biotite CIFFile or path to an mmCIF file.

    Returns
    -------
    set[str]
        Component IDs requiring user-supplied SMILES.
    """
    if not isinstance(cif_input, pdbx.CIFFile):
        cif_input = pdbx.CIFFile.read(str(cif_input))
    block = list(cif_input.values())[0]

    hetatm_ids = _get_hetatm_comp_ids(block)
    if not hetatm_ids:
        return set()

    cif_bond_ids = _get_cif_bond_comp_ids(block)

    # Collect heavy-atom names per comp_id for validation
    atom_names_per_comp: dict[str, set[str]] = {}
    if "atom_site" in block:
        atom_site = block["atom_site"]
        comp_ids = atom_site["label_comp_id"].as_array()
        a_names = atom_site["label_atom_id"].as_array()
        elements = (
            atom_site["type_symbol"].as_array() if "type_symbol" in atom_site else None
        )
        for comp_id in hetatm_ids:
            if elements is not None:
                mask = (comp_ids == comp_id) & ~is_hydrogen_isotope(elements)
            else:
                mask = comp_ids == comp_id
            atom_names_per_comp[comp_id] = set(a_names[mask])

    unknown = set()
    for comp_id in hetatm_ids:
        if comp_id in cif_bond_ids:
            # if bonds defined in cif - consider chemistry as known
            continue
        if _is_known_compound(comp_id, atom_names=atom_names_per_comp.get(comp_id)):
            continue
        unknown.add(comp_id)
    return unknown


def _rdkit_bond_to_cif(bond: Chem.rdchem.Bond) -> tuple[str, str]:
    """Map an RDKit bond to ``(value_order, pdbx_aromatic_flag)``.

    Biotite's ``_parse_intra_residue_bonds`` needs both columns; the
    ``(order, flag)`` pair keys into
    :data:`biotite.structure.io.pdbx.convert.COMP_BOND_ORDER_TO_TYPE`
    — without the aromatic flag biotite silently falls back to the
    CCD library, which fails for custom residues.
    """
    aromatic_flag = "Y" if bond.GetIsAromatic() else "N"
    order_map = {
        Chem.rdchem.BondType.SINGLE: "sing",
        Chem.rdchem.BondType.DOUBLE: "doub",
        Chem.rdchem.BondType.TRIPLE: "trip",
        Chem.rdchem.BondType.AROMATIC: "arom",
    }
    order = order_map.get(bond.GetBondType(), "sing")
    return order, aromatic_flag


def check_cif_bond_orders(cif_input: pdbx.CIFFile | Path | str) -> None:
    """Raise if any ligand has unresolvable bond orders.

    Parameters
    ----------
    cif_input : CIFFile, Path, or str
        Biotite CIFFile or path to an mmCIF file.

    Raises
    ------
    MissingBondOrderError
        If any ligand is unknown to CCD and has no ``_chem_comp_bond``.
    """
    unknown = get_unknown_ligand_ids(cif_input)
    if unknown:
        raise MissingBondOrderError(
            f"CIF file contains unknown ligands {unknown} with "
            "no _chem_comp_bond category and no CCD library match. "
            "Provide ligand SMILES to assign bond orders."
        )


def _bonds_by_position(
    comp_id: str,
    template_heavy: Chem.Mol,
    lig_heavy: struc.AtomArray,
) -> list[tuple[str, str, str, str]]:
    """Assign bonds by trusting positional atom-order correspondence.

    Assumes CIF heavy atoms appear in the same order as heavy atoms in
    the SMILES template (the convention used by Boltz, AlphaFold3,
    Chai-1, etc.). Verifies by comparing elements at each position and
    raises ``ValueError`` on any mismatch, pointing at the offending
    position so the caller can diagnose it quickly.

    Returns a list of ``(atom_name_1, atom_name_2, value_order,
    pdbx_aromatic_flag)`` tuples ready to be written to
    ``_chem_comp_bond``.
    """
    n_template = template_heavy.GetNumAtoms()
    n_cif = lig_heavy.array_length()
    if n_template != n_cif:
        raise ValueError(
            f"Atom count mismatch for {comp_id}: CIF has {n_cif} heavy atoms, "
            f"SMILES has {n_template}."
        )

    cif_elements = [str(e).upper() for e in lig_heavy.element]
    for i, (cif_el, tmpl_atom) in enumerate(
        zip(cif_elements, template_heavy.GetAtoms())
    ):
        tmpl_el = tmpl_atom.GetSymbol().upper()
        if cif_el != tmpl_el:
            raise ValueError(
                f"Element mismatch for {comp_id} at position {i}: "
                f"CIF has {cif_el}, SMILES has {tmpl_el}. Set "
                "force_substructure_match=True if the CIF doesn't "
                "preserve SMILES atom order."
            )

    atom_names = list(lig_heavy.atom_name)
    out: list[tuple[str, str, str, str]] = []
    for bond in template_heavy.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        order, aromatic_flag = _rdkit_bond_to_cif(bond)
        out.append((atom_names[idx1], atom_names[idx2], order, aromatic_flag))
    return out


def _bonds_by_substructure_match(
    comp_id: str,
    template: Chem.Mol,
    lig_heavy: struc.AtomArray,
) -> list[tuple[str, str, str, str]]:
    """Assign bonds via RDKit substructure matching.

    Opt-in alternative to :func:`_bonds_by_position` for CIFs whose
    atom order does not match SMILES parse order. Invoked only when
    ``force_substructure_match=True`` is passed to
    :func:`assign_bond_orders_from_smiles` — there is no automatic
    fallback between the two paths.

    Substructure matching needs CIF connectivity (RDKit can't search
    a graph that has no edges). If ``lig_heavy.bonds`` is empty, bonds
    are inferred from interatomic distances
    (``connect_via_distances``); bond orders are then reassigned from
    the SMILES template via ``AssignBondOrdersFromTemplate``. The
    positional path doesn't need this fallback because it never reads
    the CIF's bond list — it copies bonds straight from the SMILES
    template using positional atom-name lookup.

    Raises
    ------
    ValueError
        If the ligand cannot be parsed as an RDKit mol, or if
        :func:`peppr.sanitize` fails — a half-sanitized mol has
        undefined aromaticity perception, and feeding it to
        ``AssignBondOrdersFromTemplate`` can silently match the wrong
        substructure. Better to fail loudly than to emit chemically
        wrong bond orders.
    """
    from biotite.interface import rdkit as rdkit_interface

    # TODO(peppr): temporary local sanitize carrying boron/main-group
    # over-valence fixes not yet in a released peppr. Revert to
    # `from peppr import sanitize as peppr_sanitize` once upstream.
    # See plinder.core.utils.sanitize.
    from plinder.core.utils.sanitize import sanitize as peppr_sanitize

    if lig_heavy.bonds is None or lig_heavy.bonds.as_array().shape[0] == 0:
        # Unknown residue — infer bonds from distances
        lig_heavy.bonds = struc.connect_via_distances(lig_heavy)
    # Ensure bond types are SINGLE (1), not ANY/UNSPECIFIED (0),
    # so RDKit template matching can reassign proper orders
    bond_arr = lig_heavy.bonds.as_array()
    bond_arr[:, 2] = np.where(bond_arr[:, 2] == 0, 1, bond_arr[:, 2])
    lig_heavy.bonds = struc.BondList(lig_heavy.array_length(), bond_arr)
    rdkit_mol = rdkit_interface.to_mol(lig_heavy)
    if rdkit_mol is None:
        raise ValueError(f"Could not parse ligand {comp_id} as RDKit mol")
    try:
        peppr_sanitize(rdkit_mol)
    except Exception as e:
        raise ValueError(
            f"peppr_sanitize failed for {comp_id}: {e}. "
            "Proceeding to substructure matching with a half-sanitized "
            "mol can silently produce wrong bond orders, so aborting."
        ) from e
    fixed_mol = mol_assigned_bond_orders_by_template(template, rdkit_mol)

    atom_names = [
        a.GetPDBResidueInfo().GetName().strip()
        if a.GetPDBResidueInfo()
        else lig_heavy.atom_name[a.GetIdx()]
        for a in fixed_mol.GetAtoms()
    ]

    out: list[tuple[str, str, str, str]] = []
    for bond in fixed_mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        if idx1 < len(atom_names) and idx2 < len(atom_names):
            order, aromatic_flag = _rdkit_bond_to_cif(bond)
            out.append((atom_names[idx1], atom_names[idx2], order, aromatic_flag))
    return out


def _component_reference_instance(
    atoms: struc.AtomArray,
    comp_id: str,
) -> struc.AtomArray:
    """Return one component instance after checking all copies agree."""
    component_atoms = atoms[atoms.res_name == comp_id]
    if component_atoms.array_length() == 0:
        raise ValueError(f"No atoms found for component {comp_id} in CIF")

    instance_keys = list(
        dict.fromkeys(
            (str(chain), int(res_id))
            for chain, res_id in zip(component_atoms.chain_id, component_atoms.res_id)
        )
    )
    instances = [
        (
            key,
            component_atoms[
                (component_atoms.chain_id == key[0])
                & (component_atoms.res_id == key[1])
            ],
        )
        for key in instance_keys
    ]
    ref_key, ref_atoms = instances[0]
    ref_names = tuple(ref_atoms.atom_name)
    ref_elements = tuple(ref_atoms.element)
    for key, instance_atoms in instances[1:]:
        names = tuple(instance_atoms.atom_name)
        elements = tuple(instance_atoms.element)
        if names != ref_names or elements != ref_elements:
            raise ValueError(
                f"{comp_id}: instances disagree on heavy-atom naming/order. "
                f"Instance {ref_key} has {len(ref_names)} atoms starting with "
                f"{ref_names[:5]}; instance {key} has {len(names)} atoms "
                f"starting with {names[:5]}. mmCIF _chem_comp_bond is keyed "
                "by component ID, so every copy must use the same chemistry "
                "and atom names. Use distinct component IDs when copies differ."
            )
    if len(instances) > 1:
        LOG.info(
            "%s: %d instances share one atom naming scheme; writing one "
            "_chem_comp_bond definition",
            comp_id,
            len(instances),
        )
    return ref_atoms


def _existing_chem_comp_bonds(
    block: pdbx.CIFBlock,
    *,
    replace_components: set[str] | None = None,
) -> list[tuple[str, str, str, str, str]]:
    """Read existing component bonds except explicitly replaced components."""
    if "chem_comp_bond" not in block:
        return []
    replace_components = replace_components or set()
    category = block["chem_comp_bond"]
    required = {"comp_id", "atom_id_1", "atom_id_2", "value_order"}
    missing = sorted(required.difference(category))
    if missing:
        raise ValueError(
            "existing _chem_comp_bond is missing required fields: "
            + ", ".join(f"_chem_comp_bond.{name}" for name in missing)
        )
    comp_ids = category["comp_id"].as_array(str)
    atom_ids_1 = category["atom_id_1"].as_array(str)
    atom_ids_2 = category["atom_id_2"].as_array(str)
    orders = category["value_order"].as_array(str)
    aromatic_flags = (
        category["pdbx_aromatic_flag"].as_array(str)
        if "pdbx_aromatic_flag" in category
        else np.full(category.row_count, "N")
    )
    return [
        (comp_id, atom_id_1, atom_id_2, order.lower(), aromatic_flag)
        for comp_id, atom_id_1, atom_id_2, order, aromatic_flag in zip(
            comp_ids,
            atom_ids_1,
            atom_ids_2,
            orders,
            aromatic_flags,
        )
        if comp_id not in replace_components
    ]


def _set_chem_comp_bonds(
    block: pdbx.CIFBlock,
    bonds: list[tuple[str, str, str, str, str]],
) -> None:
    """Write dictionary-valid component bond rows."""
    block["chem_comp_bond"] = pdbx.CIFCategory(
        {
            "comp_id": [bond[0] for bond in bonds],
            "atom_id_1": [bond[1] for bond in bonds],
            "atom_id_2": [bond[2] for bond in bonds],
            "value_order": [bond[3].lower() for bond in bonds],
            "pdbx_aromatic_flag": [bond[4] for bond in bonds],
        }
    )


def enrich_cif_with_ccd_bonds(
    cif_file: pdbx.CIFFile,
    ligand_ccd_codes: dict[str, str],
) -> dict[str, str]:
    """Assign custom component bonds from Chemical Component Dictionary entries.

    The input mapping uses the custom component ID in ``_atom_site`` as its
    key and a reference CCD code as its value, for example ``{"LIG": "ATP"}``.
    Atom names are matched first. If they differ, provisional connectivity is
    inferred from the supplied coordinates and graph-matched to the CCD
    template. The function refuses partial or ambiguous mappings instead of
    assigning bonds from element order alone.

    Returns
    -------
    dict[str, str]
        Canonical CCD SMILES keyed by the custom component ID.
    """
    from plinder.data.annotations.ligand_utils import (
        _get_ccd_atomarray,
        _get_ccd_mol,
        _get_ccd_smiles,
    )

    if not ligand_ccd_codes:
        return {}
    block = list(cif_file.values())[0]
    atoms = pdbx.get_structure(
        cif_file, model=1, use_author_fields=False, include_bonds=False
    )
    atoms = atoms[~is_hydrogen_isotope(atoms.element)]
    bonds = _existing_chem_comp_bonds(
        block,
        replace_components=set(ligand_ccd_codes),
    )
    smiles_by_component: dict[str, str] = {}
    bond_type_to_cif = {
        struc.BondType.SINGLE: ("sing", "N"),
        struc.BondType.DOUBLE: ("doub", "N"),
        struc.BondType.TRIPLE: ("trip", "N"),
        struc.BondType.QUADRUPLE: ("quad", "N"),
        struc.BondType.AROMATIC_SINGLE: ("sing", "Y"),
        struc.BondType.AROMATIC_DOUBLE: ("doub", "Y"),
        struc.BondType.AROMATIC_TRIPLE: ("trip", "Y"),
        struc.BondType.AROMATIC: ("arom", "Y"),
    }
    for custom_comp_id, reference_code in ligand_ccd_codes.items():
        custom_comp_id = str(custom_comp_id).strip()
        reference_code = str(reference_code).strip().upper()
        if not custom_comp_id or not reference_code:
            raise ValueError("ligand component IDs and CCD codes must not be empty")
        custom_atoms = _component_reference_instance(atoms, custom_comp_id)
        ccd_atoms = _get_ccd_atomarray(reference_code)
        if ccd_atoms is None:
            raise ValueError(
                f"CCD code {reference_code!r} supplied for component "
                f"{custom_comp_id!r} was not found in the Chemical Component "
                "Dictionary"
            )
        ccd_atoms = ccd_atoms[~is_hydrogen_isotope(ccd_atoms.element)]
        if ccd_atoms.bonds is None or len(ccd_atoms.bonds.as_array()) == 0:
            raise ValueError(
                f"CCD code {reference_code!r} has no heavy-atom bonds to assign "
                f"to component {custom_comp_id!r}"
            )
        if custom_atoms.array_length() != ccd_atoms.array_length():
            raise ValueError(
                f"cannot map component {custom_comp_id!r} to CCD code "
                f"{reference_code!r}: custom CIF has "
                f"{custom_atoms.array_length()} heavy atoms but CCD has "
                f"{ccd_atoms.array_length()}"
            )

        custom_names = [str(name) for name in custom_atoms.atom_name]
        ccd_names = [str(name) for name in ccd_atoms.atom_name]
        if len(set(custom_names)) != len(custom_names):
            raise ValueError(
                f"component {custom_comp_id!r} has duplicate heavy-atom names; "
                "CCD bond assignment requires unique atom names"
            )
        bonds_to_emit: list[tuple[str, str, str, str]] | None = None
        if set(custom_names) == set(ccd_names):
            ccd_element_by_name = dict(zip(ccd_names, ccd_atoms.element))
            element_mismatches = [
                name
                for name, element in zip(custom_names, custom_atoms.element)
                if element != ccd_element_by_name[name]
            ]
            if element_mismatches:
                raise ValueError(
                    f"cannot map component {custom_comp_id!r} to CCD code "
                    f"{reference_code!r}: atom names match but elements differ "
                    f"for {element_mismatches[:5]}"
                )
            ccd_index_to_custom_name = {
                index: ccd_name for index, ccd_name in enumerate(ccd_names)
            }
        else:
            ccd_template = _get_ccd_mol(reference_code)
            if ccd_template is None:
                raise ValueError(
                    f"CCD code {reference_code!r} could not provide an RDKit "
                    f"template for component {custom_comp_id!r}"
                )
            try:
                bonds_to_emit = _bonds_by_substructure_match(
                    custom_comp_id,
                    ccd_template,
                    custom_atoms,
                )
            except ValueError as exc:
                raise ValueError(
                    f"cannot safely map component {custom_comp_id!r} to CCD "
                    f"code {reference_code!r} from atom names or inferred "
                    f"connectivity: {exc}. Provide a SMILES whose parse order "
                    "matches the custom CIF if its geometry does not support "
                    "distance-based bond inference."
                ) from exc

        if bonds_to_emit is None:
            bonds_to_emit = []
            for atom_index_1, atom_index_2, bond_type in ccd_atoms.bonds.as_array():
                bond_type = struc.BondType(bond_type)
                if bond_type not in bond_type_to_cif:
                    raise ValueError(
                        f"CCD code {reference_code!r} contains unsupported bond "
                        f"type {bond_type.name} for component {custom_comp_id!r}"
                    )
                order, aromatic_flag = bond_type_to_cif[bond_type]
                bonds_to_emit.append(
                    (
                        ccd_index_to_custom_name[int(atom_index_1)],
                        ccd_index_to_custom_name[int(atom_index_2)],
                        order,
                        aromatic_flag,
                    )
                )
        bonds.extend(
            (custom_comp_id, atom_1, atom_2, order, aromatic_flag)
            for atom_1, atom_2, order, aromatic_flag in bonds_to_emit
        )
        smiles = _get_ccd_smiles(reference_code)
        if smiles is None:
            raise ValueError(
                f"CCD code {reference_code!r} could not be converted to RDKit "
                f"chemistry for component {custom_comp_id!r}"
            )
        smiles_by_component[custom_comp_id] = smiles

    _set_chem_comp_bonds(block, bonds)
    return smiles_by_component


def enrich_cif_with_smiles_bonds(
    cif_file: pdbx.CIFFile,
    ligand_smiles: dict[str, str],
    force_substructure_match: bool = False,
) -> None:
    """Add ``_chem_comp_bond`` rows to a CIFFile in-memory.

    Mutates ``cif_file`` by appending bond entries for unknown ligands
    using the provided SMILES templates. Known CCD compounds are
    skipped. Existing ``_chem_comp_bond`` rows are preserved.

    See :func:`assign_bond_orders_from_smiles` for the full description
    of the atom-order assumption and the ``force_substructure_match``
    opt-in.

    Parameters
    ----------
    cif_file : pdbx.CIFFile
        CIF object to mutate in place.
    ligand_smiles : dict[str, str]
        Mapping of component ID (e.g. ``LIG``) to SMILES.
    force_substructure_match : bool, default=False
        If ``True``, skip the positional element check entirely and
        assign bonds via RDKit substructure matching instead.

    Raises
    ------
    MissingBondOrderError
        If unknown ligands remain without SMILES.
    ValueError
        If SMILES is invalid, atom counts differ, element order does
        not match (default path), sanitize / template matching fails
        (substructure path — this path fails loudly rather than risk
        emitting chemically wrong bond orders), or multiple instances
        of the same comp_id disagree on heavy-atom naming/order
        (mmCIF ``_chem_comp_bond`` is keyed by comp_id so all
        instances must share atom naming for biotite to apply the
        single bond definition correctly).
    """
    block = list(cif_file.values())[0]

    unknown_ids = get_unknown_ligand_ids(cif_file)
    if not unknown_ids:
        LOG.info("All ligands are known or already have bond orders, nothing to do")
        return

    missing_smiles = unknown_ids - set(ligand_smiles.keys())
    if missing_smiles:
        raise MissingBondOrderError(
            f"Unknown ligands {missing_smiles} need SMILES but none were provided"
        )

    to_process = {k: v for k, v in ligand_smiles.items() if k in unknown_ids}
    skipped = set(ligand_smiles.keys()) - unknown_ids
    if skipped:
        LOG.info(f"Skipping known compounds: {skipped}")

    atoms = pdbx.get_structure(
        cif_file, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[~is_hydrogen_isotope(atoms.element)]
    bonds = _existing_chem_comp_bonds(block)

    for comp_id, smiles in to_process.items():
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError(f"Invalid SMILES for {comp_id}: {smiles}")
        template_heavy = Chem.RemoveHs(template, sanitize=False)

        lig_heavy = _component_reference_instance(atoms, comp_id)

        if force_substructure_match:
            bonds_to_emit = _bonds_by_substructure_match(comp_id, template, lig_heavy)
        else:
            bonds_to_emit = _bonds_by_position(comp_id, template_heavy, lig_heavy)

        for atom_name_1, atom_name_2, value_order, aromatic_flag in bonds_to_emit:
            bonds.append(
                (
                    comp_id,
                    atom_name_1,
                    atom_name_2,
                    value_order,
                    aromatic_flag,
                )
            )

    _set_chem_comp_bonds(block, bonds)


def assign_bond_orders_from_smiles(
    cif_path: Path,
    ligand_smiles: dict[str, str],
    output_path: Path | None = None,
    force_substructure_match: bool = False,
) -> Path:
    """Disk-based wrapper around :func:`enrich_cif_with_smiles_bonds`.

    Reads ``cif_path``, enriches the CIF in memory, and writes the
    result to ``output_path`` (or overwrites ``cif_path`` when
    ``output_path`` is ``None``). Callers that already have a
    ``pdbx.CIFFile`` object in memory should use
    :func:`enrich_cif_with_smiles_bonds` directly to avoid the read /
    write round-trip.

    Atom-order assumption
    ---------------------
    By default this function assumes that the heavy-atom order in the
    CIF exactly matches the heavy-atom parse order of the SMILES. This
    is the convention produced by structure-prediction tools that
    accept SMILES input (e.g. Boltz, AlphaFold3, Chai-1): their output
    CIF writes ligand atoms in the same order that the SMILES was
    parsed. Under this assumption the mapping from CIF atom -> SMILES
    atom is the identity, and bond orders can be copied directly from
    the SMILES template with zero ambiguity.

    The function verifies the assumption by comparing the element at
    each position. If counts or elements don't match, ``ValueError``
    is raised pointing at the first mismatch.

    Set ``force_substructure_match=True`` to fully replace the default
    path with RDKit substructure matching. This does NOT fall back on
    failure — it is the only method used when the flag is set. Slower,
    can be ambiguous for symmetric molecules, and should only be used
    for CIFs from tools that don't preserve SMILES atom order.

    Parameters
    ----------
    cif_path : Path
        Input mmCIF file.
    ligand_smiles : dict[str, str]
        Mapping of component ID (e.g. ``LIG``) to SMILES.
    output_path : Path | None
        Output path. Defaults to overwriting *cif_path*.
    force_substructure_match : bool, default=False
        If ``True``, skip the positional element check entirely and
        assign bonds via RDKit substructure matching instead. Use only
        when CIF atom order is not guaranteed to match SMILES order.

    Returns
    -------
    Path
        Path to the written CIF file.

    Raises
    ------
    MissingBondOrderError, ValueError
        Propagated from :func:`enrich_cif_with_smiles_bonds`. See
        that function's docstring for the full list of failure modes.
    """
    if output_path is None:
        output_path = cif_path
    cif_file = pdbx.CIFFile.read(str(cif_path))
    enrich_cif_with_smiles_bonds(
        cif_file,
        ligand_smiles=ligand_smiles,
        force_substructure_match=force_substructure_match,
    )
    cif_file.write(str(output_path))
    return output_path
