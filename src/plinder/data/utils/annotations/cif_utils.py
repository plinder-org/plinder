# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""mmCIF I/O utilities using biotite.

Generic helpers for reading CIF blocks, extracting scalar values and
category rows, plus ligand bond-order detection and assignment from
SMILES templates.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import biotite.structure as struc
import biotite.structure.info as bt_info
import biotite.structure.io.pdbx as pdbx
import numpy as np
from rdkit import Chem

from plinder.core.structure.smallmols_utils import (
    mol_assigned_bond_orders_by_template,
)

LOG = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Generic CIF I/O helpers
# ---------------------------------------------------------------------------


def read_mmcif_file(mmcif_filename: Path | str) -> pdbx.CIFFile:
    """Read an mmCIF file, handling .gz transparently."""
    import gzip

    path = str(mmcif_filename)
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as f:
            return pdbx.CIFFile.read(f)
    return pdbx.CIFFile.read(path)


def read_mmcif_container(mmcif_filename: Path) -> pdbx.CIFBlock:
    """Parse mmcif file and return the first data block."""
    cif_file = read_mmcif_file(mmcif_filename)
    return list(cif_file.values())[0]


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
# CIF → RDKit conversion
# ---------------------------------------------------------------------------


def atoms_to_rdkit_mol(
    atoms: "struc.AtomArray",
    assign_stereo: bool = True,
) -> "Chem.Mol":
    """Convert a biotite AtomArray to a sanitized RDKit Mol.

    Hydrogen atoms are removed.  Bonds are assigned from CCD residue
    names if not already present.  Stereochemistry is optionally
    assigned from 3D coordinates.

    Parameters
    ----------
    atoms : AtomArray
        Heavy atoms with optional bonds (e.g. from ``include_bonds=True``).
        If bonds are missing, ``connect_via_residue_names`` is used.
    assign_stereo : bool
        If True, call ``AssignStereochemistryFrom3D`` on the result.

    Returns
    -------
    Chem.Mol
        Sanitized RDKit molecule with 3D coordinates and PDB atom info.

    Raises
    ------
    ValueError
        If the conversion fails.
    """
    from biotite.interface import rdkit as rdkit_interface
    from peppr import sanitize as peppr_sanitize

    heavy = atoms[atoms.element != "H"]
    if heavy.bonds is None or heavy.bonds.as_array().shape[0] == 0:
        heavy.bonds = struc.connect_via_residue_names(heavy)
    mol = rdkit_interface.to_mol(heavy)
    if mol is None:
        raise ValueError("Failed to convert AtomArray to RDKit Mol")
    peppr_sanitize(mol)
    if assign_stereo:
        Chem.AssignStereochemistryFrom3D(mol)
    return mol


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

    connections = parse_struct_conn(block)
    if not connections:
        return
    if atoms.bonds is None:
        atoms.bonds = struc.BondList(atoms.array_length())

    existing = set(
        (min(b[0], b[1]), max(b[0], b[1])) for b in atoms.bonds.as_array()[:, :2]
    )
    label_ids = np.array([c.split(".")[-1] if "." in c else c for c in atoms.chain_id])

    for c in connections:
        if c["conn_type"] != "covale":
            continue
        try:
            r1 = int(c["seq1"]) if c["seq1"] != "." else -1
            r2 = int(c["seq2"]) if c["seq2"] != "." else -1
        except ValueError:
            continue

        mask1 = (
            (label_ids == c["chain1"])
            & (atoms.res_id == r1)
            & (atoms.atom_name == c["atom1"])
        )
        mask2 = (
            (label_ids == c["chain2"])
            & (atoms.res_id == r2)
            & (atoms.atom_name == c["atom2"])
        )

        for i1 in np.where(mask1)[0]:
            for i2 in np.where(mask2)[0]:
                pair = (min(int(i1), int(i2)), max(int(i1), int(i2)))
                if pair not in existing:
                    atoms.bonds.add_bond(int(i1), int(i2), struc.BondType.SINGLE)
                    existing.add(pair)


# ---------------------------------------------------------------------------
# Bridged interaction detection (synced with peppr-internal)
# TODO: remove once peppr >= 0.14 is released with these methods.
# ---------------------------------------------------------------------------

# Water bridge lower bound: 0.75 * VdW_sum (~2.28 A for O-O)
# avoids clashes but allows short water-mediated H-bonds.
# Upper bound: 1.15 * VdW_sum (~3.50 A for O-O), standard H-bond max.
_WATER_BRIDGE_DISTANCE_SCALING = (0.75, 1.15)

# Metals that form coordination bonds (not spectator ions like Na/Cl/K)
_COORDINATION_METALS = frozenset(
    {
        "MG",
        "CA",
        "ZN",
        "FE",
        "FE2",  # Fe(II)
        "MN",
        "CO",
        "CU",
        "CU1",  # Cu(I)
        "NI",
        "CD",
        "MO",
        "4MO",  # Mo(IV)
        "6MO",  # Mo(VI)
        "W",
        "V",
    }
)
_METAL_ACCEPTOR_PATTERN = (
    "["
    "$([O]),"
    "$([#7;!$([nX3]);!$([NX3]-*=[!#6]);!$([NX3]-[a]);!$([NX4])]),"
    "$([#16]),"
    "$([*;-{1-};!+{1-}])"
    "]"
)


def _find_bridged_interactions(
    receptor: "struc.AtomArray",
    ligand: "struc.AtomArray",
    bridge_atoms: "struc.AtomArray",
    receptor_pattern: str,
    ligand_pattern: str,
    distance_scaling: tuple[float, float],
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Find interactions bridged by intermediary atoms (water or metal).

    TODO: remove once peppr has ContactMeasurement.find_bridged_interactions.
    """
    import biotite.structure.info as info
    from peppr.contacts import ContactMeasurement, find_atoms_by_pattern

    if bridge_atoms.array_length() == 0:
        return []

    try:
        cm = ContactMeasurement(receptor, ligand)
    except Exception as e:
        LOG.warning(f"ContactMeasurement setup failed: {e}")
        return []

    receptor_matched = find_atoms_by_pattern(cm._binding_site_mol, receptor_pattern)
    ligand_matched = find_atoms_by_pattern(cm._ligand_mol, ligand_pattern)
    if len(receptor_matched) == 0 or len(ligand_matched) == 0:
        return []

    receptor_coords = cm._binding_site.coord[receptor_matched]
    ligand_coords = cm._ligand.coord[ligand_matched]
    lo, hi = sorted(distance_scaling)

    r_vdw = np.array(
        [info.vdw_radius_single(e) for e in cm._binding_site.element[receptor_matched]]
    )
    l_vdw = np.array(
        [info.vdw_radius_single(e) for e in cm._ligand.element[ligand_matched]]
    )

    bridges: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for bi in range(bridge_atoms.array_length()):
        b_coord = bridge_atoms.coord[bi]
        b_vdw = info.vdw_radius_single(bridge_atoms.element[bi])

        r_dists = np.linalg.norm(receptor_coords - b_coord, axis=1)
        r_thresholds = r_vdw + b_vdw
        r_contacts = receptor_matched[
            (r_dists >= lo * r_thresholds) & (r_dists <= hi * r_thresholds)
        ]
        if len(r_contacts) == 0:
            continue

        l_dists = np.linalg.norm(ligand_coords - b_coord, axis=1)
        l_thresholds = l_vdw + b_vdw
        l_contacts = ligand_matched[
            (l_dists >= lo * l_thresholds) & (l_dists <= hi * l_thresholds)
        ]
        if len(l_contacts) == 0:
            continue

        for ri in r_contacts:
            for li in l_contacts:
                bridges.append(
                    (
                        cm._binding_site_indices[ri : ri + 1],
                        np.array([li], dtype=int),
                        np.array([bi], dtype=int),
                    )
                )

    return bridges


def find_water_bridges(
    receptor: "struc.AtomArray",
    ligand: "struc.AtomArray",
    waters: "struc.AtomArray",
    distance_scaling: tuple[float, float] = _WATER_BRIDGE_DISTANCE_SCALING,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Find water-mediated hydrogen bonds between receptor and ligand."""
    from peppr.common import ACCEPTOR_PATTERN, DONOR_PATTERN

    water_oxygens = waters[waters.element == "O"]
    hbond_pattern = "[" + DONOR_PATTERN[1:-1] + "," + ACCEPTOR_PATTERN[1:-1] + "]"
    return _find_bridged_interactions(
        receptor,
        ligand,
        water_oxygens,
        hbond_pattern,
        hbond_pattern,
        distance_scaling,
    )


def find_metal_bridges(
    receptor: "struc.AtomArray",
    ligand: "struc.AtomArray",
    metals: "struc.AtomArray",
    cutoff: float = 3.0,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Find metal-mediated coordination between receptor and ligand."""
    from peppr.contacts import ContactMeasurement, find_atoms_by_pattern

    coord_mask = np.isin(metals.res_name, list(_COORDINATION_METALS))
    if not np.any(coord_mask):
        return []
    coord_metals = metals[coord_mask]

    try:
        cm = ContactMeasurement(receptor, ligand)
    except Exception as e:
        LOG.warning(f"ContactMeasurement setup failed for metal bridges: {e}")
        return []

    receptor_matched = find_atoms_by_pattern(
        cm._binding_site_mol, _METAL_ACCEPTOR_PATTERN
    )
    ligand_matched = find_atoms_by_pattern(cm._ligand_mol, _METAL_ACCEPTOR_PATTERN)
    if len(receptor_matched) == 0 or len(ligand_matched) == 0:
        return []

    bridges: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for bi in range(coord_metals.array_length()):
        b_coord = coord_metals.coord[bi]
        r_dists = np.linalg.norm(
            cm._binding_site.coord[receptor_matched] - b_coord, axis=1
        )
        r_contacts = receptor_matched[r_dists < cutoff]
        if len(r_contacts) == 0:
            continue
        l_dists = np.linalg.norm(cm._ligand.coord[ligand_matched] - b_coord, axis=1)
        l_contacts = ligand_matched[l_dists < cutoff]
        if len(l_contacts) == 0:
            continue
        for ri in r_contacts:
            for li in l_contacts:
                bridges.append(
                    (
                        cm._binding_site_indices[ri : ri + 1],
                        np.array([li], dtype=int),
                        np.array([bi], dtype=int),
                    )
                )
    return bridges


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

    If *atom_names* is provided, also verify that the CIF atom names
    overlap with the CCD entry. Bond assignment via
    ``connect_via_residue_names`` relies on atom-name matching, so a
    compound whose names don't match CCD will get wrong bonds even if
    the comp_id exists in the dictionary (e.g. Boltz ``LIG`` =/= CCD
    ``LIG``).
    """
    try:
        ref = bt_info.residue(comp_id)
        if atom_names is not None:
            ref_heavy = ref[ref.element != "H"]
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
                mask = (comp_ids == comp_id) & (elements != "H")
            else:
                mask = comp_ids == comp_id
            atom_names_per_comp[comp_id] = set(a_names[mask])

    unknown = set()
    for comp_id in hetatm_ids:
        if comp_id in cif_bond_ids:
            continue
        if _is_known_compound(comp_id, atom_names=atom_names_per_comp.get(comp_id)):
            continue
        unknown.add(comp_id)
    return unknown


def _rdkit_bond_order_to_cif(bond_type: Chem.rdchem.BondType) -> str:
    mapping = {
        Chem.rdchem.BondType.SINGLE: "SING",
        Chem.rdchem.BondType.DOUBLE: "DOUB",
        Chem.rdchem.BondType.TRIPLE: "TRIP",
        Chem.rdchem.BondType.AROMATIC: "AROM",
    }
    return mapping.get(bond_type, "SING")


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


def assign_bond_orders_from_smiles(
    cif_path: Path,
    ligand_smiles: dict[str, str],
    output_path: Path | None = None,
) -> Path:
    """Assign bond orders to unknown ligands using SMILES templates.

    Known CCD compounds are skipped. Writes ``_chem_comp_bond`` into
    the CIF, preserving any existing entries.

    Parameters
    ----------
    cif_path : Path
        Input mmCIF file.
    ligand_smiles : dict[str, str]
        Mapping of component ID (e.g. ``LIG``) to SMILES.
    output_path : Path | None
        Output path. Defaults to overwriting *cif_path*.

    Returns
    -------
    Path
        Path to the written CIF file.

    Raises
    ------
    MissingBondOrderError
        If unknown ligands remain without SMILES.
    ValueError
        If SMILES is invalid or template matching fails.
    """
    if output_path is None:
        output_path = cif_path

    cif_file = pdbx.CIFFile.read(str(cif_path))
    block = list(cif_file.values())[0]

    # Determine which ligands actually need bond order assignment
    unknown_ids = get_unknown_ligand_ids(cif_file)
    if not unknown_ids:
        LOG.info("All ligands are known or already have bond orders, nothing to do")
        cif_file.write(str(output_path))
        return output_path

    # Check that user provided SMILES for all unknown ligands
    missing_smiles = unknown_ids - set(ligand_smiles.keys())
    if missing_smiles:
        raise MissingBondOrderError(
            f"Unknown ligands {missing_smiles} need SMILES but none were provided"
        )

    # Filter to only process unknown ligands
    to_process = {k: v for k, v in ligand_smiles.items() if k in unknown_ids}
    skipped = set(ligand_smiles.keys()) - unknown_ids
    if skipped:
        LOG.info(f"Skipping known compounds: {skipped}")

    atoms = pdbx.get_structure(
        cif_file, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[atoms.element != "H"]

    # Preserve existing _chem_comp_bond rows
    comp_id_list: list[str] = []
    atom_id_1_list: list[str] = []
    atom_id_2_list: list[str] = []
    value_order_list: list[str] = []

    if "chem_comp_bond" in block:
        existing = block["chem_comp_bond"]
        for i in range(existing.row_count):
            comp_id_list.append(existing["comp_id"].as_array()[i])
            atom_id_1_list.append(existing["atom_id_1"].as_array()[i])
            atom_id_2_list.append(existing["atom_id_2"].as_array()[i])
            value_order_list.append(existing["value_order"].as_array()[i])

    from biotite.interface import rdkit as rdkit_interface
    from peppr import sanitize as peppr_sanitize

    for comp_id, smiles in to_process.items():
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError(f"Invalid SMILES for {comp_id}: {smiles}")

        lig_mask = atoms.res_name == comp_id
        if not np.any(lig_mask):
            raise ValueError(f"No atoms found for component {comp_id} in CIF")

        lig_atoms = atoms[lig_mask]
        lig_heavy = lig_atoms[lig_atoms.element != "H"]

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
            LOG.warning(f"peppr_sanitize failed for {comp_id}: {e}")
        fixed_mol = mol_assigned_bond_orders_by_template(template, rdkit_mol)

        atom_names = [
            a.GetPDBResidueInfo().GetName().strip()
            if a.GetPDBResidueInfo()
            else lig_heavy.atom_name[a.GetIdx()]
            for a in fixed_mol.GetAtoms()
        ]

        for bond in fixed_mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            if idx1 < len(atom_names) and idx2 < len(atom_names):
                comp_id_list.append(comp_id)
                atom_id_1_list.append(atom_names[idx1])
                atom_id_2_list.append(atom_names[idx2])
                value_order_list.append(_rdkit_bond_order_to_cif(bond.GetBondType()))

    bond_cat = pdbx.CIFCategory(
        {
            "comp_id": comp_id_list,
            "atom_id_1": atom_id_1_list,
            "atom_id_2": atom_id_2_list,
            "value_order": value_order_list,
        }
    )
    block["chem_comp_bond"] = bond_cat

    cif_file.write(str(output_path))
    return output_path
