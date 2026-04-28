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

# Single source of truth lives in ``plinder.core.structure.atoms`` so
# both ``plinder.core`` and ``plinder.data`` filter H/D/T isotopes
# consistently.
from plinder.core.structure.atoms import is_hydrogen_isotope  # noqa: E402

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


def get_model_count(cif_file: pdbx.CIFFile) -> int:
    """Return the number of models in a CIF (1 if no model column present)."""
    block = list(cif_file.values())[0]
    if "atom_site" not in block:
        return 0
    atom_site = block["atom_site"]
    if "pdbx_PDB_model_num" not in atom_site:
        return 1
    return int(len(set(atom_site["pdbx_PDB_model_num"].as_array())))


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
        If True, call ``AssignStereochemistryFrom3D`` on the result.

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
    from peppr import sanitize as peppr_sanitize

    heavy = atoms[~is_hydrogen_isotope(atoms.element)]
    # Multi-atom inputs must carry bonds; single atoms (ions) don't need any.
    if heavy.array_length() > 1 and (
        heavy.bonds is None or heavy.bonds.as_array().shape[0] == 0
    ):
        raise ValueError(
            "atoms_to_rdkit_mol requires bonds on multi-atom inputs. "
            "Load the CIF with include_bonds=True (which parses "
            "_chem_comp_bond + _struct_conn) or populate atoms.bonds "
            "before calling. A connect_via_residue_names fallback was "
            "removed because it silently drops inter-residue peptide "
            "bonds for non-standard residues in multi-residue ligands."
        )
    mol = rdkit_interface.to_mol(heavy)
    if mol is None:
        raise ValueError("Failed to convert AtomArray to RDKit Mol")
    peppr_sanitize(mol)
    if assign_stereo:
        Chem.AssignStereochemistryFrom3D(mol)
    # RDKit's RemoveAllHs keys on atomic number, so it strips any
    # hydrogen isotope atom that survived the element-string filter.
    # Safe after stereo assignment — chiral tags live on heavy atoms.
    return Chem.RemoveAllHs(mol)


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
        Chem.rdchem.BondType.SINGLE: "SING",
        Chem.rdchem.BondType.DOUBLE: "DOUB",
        Chem.rdchem.BondType.TRIPLE: "TRIP",
        Chem.rdchem.BondType.AROMATIC: "AROM",
    }
    order = order_map.get(bond.GetBondType(), "SING")
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
    from peppr import sanitize as peppr_sanitize

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

    # Preserve existing _chem_comp_bond rows. biotite's parser requires
    # pdbx_aromatic_flag to consume the category — default to "N" when
    # absent so pre-existing rows remain parseable.
    comp_id_list: list[str] = []
    atom_id_1_list: list[str] = []
    atom_id_2_list: list[str] = []
    value_order_list: list[str] = []
    aromatic_flag_list: list[str] = []

    if "chem_comp_bond" in block:
        existing = block["chem_comp_bond"]
        existing_flag = (
            existing["pdbx_aromatic_flag"].as_array()
            if "pdbx_aromatic_flag" in existing
            else None
        )
        for i in range(existing.row_count):
            comp_id_list.append(existing["comp_id"].as_array()[i])
            atom_id_1_list.append(existing["atom_id_1"].as_array()[i])
            atom_id_2_list.append(existing["atom_id_2"].as_array()[i])
            value_order_list.append(existing["value_order"].as_array()[i])
            aromatic_flag_list.append(
                existing_flag[i] if existing_flag is not None else "N"
            )

    for comp_id, smiles in to_process.items():
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError(f"Invalid SMILES for {comp_id}: {smiles}")
        template_heavy = Chem.RemoveHs(template, sanitize=False)

        lig_mask = atoms.res_name == comp_id
        if not np.any(lig_mask):
            raise ValueError(f"No atoms found for component {comp_id} in CIF")

        # mmCIF schema keys ``_chem_comp_bond`` by ``comp_id``, not by
        # instance — biotite applies a single bond definition to every
        # copy via atom-name lookup. So multi-instance custom residues
        # (docking ensembles, multi-copy systems) require that all
        # instances share the same heavy-atom naming, otherwise the
        # bonds we emit from instance 1 won't be findable in the others.
        # We validate that explicitly and emit bonds once from the
        # reference instance — refuse to silently produce wrong bonds.
        all_lig_atoms = atoms[lig_mask]
        instances: list[tuple[tuple[str, int], struc.AtomArray]] = []
        seen_keys: dict[tuple[str, int], None] = {}
        for chain, res_id in zip(all_lig_atoms.chain_id, all_lig_atoms.res_id):
            seen_keys.setdefault((str(chain), int(res_id)), None)
        for chain, res_id in seen_keys:
            inst_mask = (all_lig_atoms.chain_id == chain) & (
                all_lig_atoms.res_id == res_id
            )
            inst = all_lig_atoms[inst_mask]
            inst_heavy = inst[~is_hydrogen_isotope(inst.element)]
            instances.append(((chain, res_id), inst_heavy))

        ref_key, ref_heavy = instances[0]
        ref_names = tuple(ref_heavy.atom_name)
        for key, inst_heavy in instances[1:]:
            inst_names = tuple(inst_heavy.atom_name)
            if inst_names != ref_names:
                raise ValueError(
                    f"{comp_id}: instances disagree on heavy-atom naming/order. "
                    f"Instance {ref_key} has {len(ref_names)} atoms "
                    f"starting with {ref_names[:5]}; instance {key} "
                    f"has {len(inst_names)} atoms starting with "
                    f"{inst_names[:5]}. mmCIF ``_chem_comp_bond`` is "
                    "keyed by comp_id and biotite applies bonds to all "
                    "copies via atom-name match — every instance must "
                    "share identical heavy-atom naming. Use distinct "
                    "comp_ids if instances differ chemically."
                )
        if len(instances) > 1:
            LOG.info(
                f"{comp_id}: {len(instances)} instances with consistent "
                "atom naming, defining _chem_comp_bond once "
                "(biotite applies to all copies via atom-name match)."
            )
        lig_heavy = ref_heavy

        if force_substructure_match:
            bonds_to_emit = _bonds_by_substructure_match(comp_id, template, lig_heavy)
        else:
            bonds_to_emit = _bonds_by_position(comp_id, template_heavy, lig_heavy)

        for atom_name_1, atom_name_2, value_order, aromatic_flag in bonds_to_emit:
            comp_id_list.append(comp_id)
            atom_id_1_list.append(atom_name_1)
            atom_id_2_list.append(atom_name_2)
            value_order_list.append(value_order)
            aromatic_flag_list.append(aromatic_flag)

    block["chem_comp_bond"] = pdbx.CIFCategory(
        {
            "comp_id": comp_id_list,
            "atom_id_1": atom_id_1_list,
            "atom_id_2": atom_id_2_list,
            "value_order": value_order_list,
            "pdbx_aromatic_flag": aromatic_flag_list,
        }
    )


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
