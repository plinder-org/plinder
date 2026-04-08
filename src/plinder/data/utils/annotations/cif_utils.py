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

import biotite.structure.io.pdbx as pdbx
import numpy as np
from ost import conop, io, mol
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol

from plinder.core.structure.smallmols_utils import (
    mol_assigned_bond_orders_by_template,
)

LOG = logging.getLogger(__name__)
_COMPOUND_LIB = conop.GetDefaultLib()


# ---------------------------------------------------------------------------
# Generic CIF I/O helpers
# ---------------------------------------------------------------------------


def read_mmcif_container(mmcif_filename: Path) -> pdbx.CIFBlock:
    """Parse mmcif file and return the first data block.

    Parameters
    ----------
    mmcif_filename : Path
    Returns
    -------
    pdbx.CIFBlock
    """
    path = str(mmcif_filename)
    if path.endswith(".gz"):
        import gzip

        with gzip.open(path, "rt", encoding="utf-8") as f:
            cif_file = pdbx.CIFFile.read(f)
    else:
        cif_file = pdbx.CIFFile.read(path)
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
            (row["seq_id_start"], row["seq_id_end"])
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
            (row["seq_id_start"], row["seq_id_end"])
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
# CIF ligand parsing
# ---------------------------------------------------------------------------


def get_ligand_chainid_comp_id_map(data: pdbx.CIFBlock) -> dict[str, set[str]]:
    """Map chain IDs to their non-polymer component IDs."""
    if "atom_site" not in data:
        return {}
    atom_site = data["atom_site"]
    group_pdb = atom_site["group_PDB"].as_array()
    comp_ids = atom_site["label_comp_id"].as_array()
    asym_ids = atom_site["label_asym_id"].as_array()

    chain_comp_id_map: dict[str, set[str]] = defaultdict(set)
    for i in range(len(group_pdb)):
        if group_pdb[i] == "HETATM":
            chain_comp_id_map[asym_ids[i]].add(comp_ids[i])
    return chain_comp_id_map


def get_bond_info(
    data: pdbx.CIFBlock, comp_ids: set[str]
) -> dict[str, list[tuple[str, str, str]]]:
    """Extract _chem_comp_bond info for given component IDs."""
    if "chem_comp_bond" not in data:
        return {}
    bond_cat = data["chem_comp_bond"]
    cids = bond_cat["comp_id"].as_array()
    a1s = bond_cat["atom_id_1"].as_array()
    a2s = bond_cat["atom_id_2"].as_array()
    orders = bond_cat["value_order"].as_array()

    bonds_dict: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for i in range(len(cids)):
        if cids[i] not in comp_ids or cids[i] == "HOH":
            continue
        bonds_dict[cids[i]].append((a1s[i], a2s[i], orders[i]))
    return bonds_dict


def bond_pdb_order(value_order: str) -> Chem.rdchem.BondType | None:
    """Convert PDB bond order string to RDKit BondType."""
    if value_order.casefold() == "sing":
        return Chem.rdchem.BondType(1)
    if value_order.casefold() == "doub":
        return Chem.rdchem.BondType(2)
    if value_order.casefold() == "trip":
        return Chem.rdchem.BondType(3)
    return None


def get_rdkit_mol_from_pdb_block(
    pdb_block: str, bonds_dict: dict[str, list[tuple[str, str, str]]]
) -> str:
    """Build SMILES from PDB block using _chem_comp_bond info."""
    rdmol = Chem.MolFromPDBBlock(pdb_block)
    atoms_ids = [
        f"{atm.GetPDBResidueInfo().GetResidueName().strip()}"
        + f":{atm.GetPDBResidueInfo().GetName().strip()}"
        for atm in rdmol.GetAtoms()
    ]

    rw_mol = RWMol(rdmol)
    for comp_id, bonds in bonds_dict.items():
        for row in bonds:
            atom_1, atom_2 = row[0], row[1]
            if (f"{comp_id}:{atom_1}" not in atoms_ids) | (
                f"{comp_id}:{atom_2}" not in atoms_ids
            ):
                pass
            if atom_1.startswith("H") | atom_2.startswith("H"):
                pass
            else:
                try:
                    atom_1_ids = _get_all_indices(atoms_ids, f"{comp_id}:{atom_1}")
                    atom_2_ids = _get_all_indices(atoms_ids, f"{comp_id}:{atom_2}")
                    for a1, a2, order in zip(
                        atom_1_ids, atom_2_ids, np.repeat(row[2], len(atom_1_ids))
                    ):
                        bo = bond_pdb_order(order)
                        rw_mol.RemoveBond(int(a1), int(a2))
                        rw_mol.AddBond(int(a1), int(a2), bo)
                except ValueError:
                    LOG.warning(f"Error perceiving {atom_1}-{atom_2} bond")
                except RuntimeError:
                    LOG.warning(f"Duplicate bond {atom_1}-{atom_2}")

    return str(Chem.MolToSmiles(rw_mol.GetMol()))


def _get_all_indices(lst: list[str], item: str) -> list[int]:
    arr = np.array(lst)
    return [int(i) for i in np.where(arr == item)[0]]


def get_smiles_from_cif(
    data: pdbx.CIFBlock, ent: io.EntityHandle, polymer_cutoff: int = 20
) -> dict[str, str]:
    """Extract SMILES for each ligand chain using _chem_comp_bond."""
    from plinder.data.utils.annotations.interaction_utils import pdbize

    rdk_mols = {}
    chain_id_comp_id_map = get_ligand_chainid_comp_id_map(data)
    for chain_id, list_of_comp_ids in chain_id_comp_id_map.items():
        bonds_dict = get_bond_info(data, list_of_comp_ids)
        mol_ent = mol.CreateEntityFromView(
            ent.Select(f"chain='{chain_id}'"),
            True,
        )
        if len(mol_ent.residues) < polymer_cutoff:
            pdb_block = io.EntityToPDBStr(pdbize(ent, mol_ent)[0])
            rdk_mols[chain_id] = get_rdkit_mol_from_pdb_block(pdb_block, bonds_dict)
        elif sum([res.name == "HOH" for res in mol_ent.residues]) > 0:
            continue
    return rdk_mols


def get_rdkit_mol_with_bond_order_from_cif(
    rdk_smiles_dict: dict[str, str], chain_id: str
) -> str:
    return rdk_smiles_dict.get(chain_id, "")


def ost_ent_to_rdkit_mol(ent: mol.EntityHandle) -> Chem.Mol | None:
    """Convert an OST entity to an RDKit Mol via PDB block, with SDF fallback."""
    pdbstring = io.EntityToPDBStr(ent).strip()
    rdkit_mol = Chem.MolFromPDBBlock(pdbstring, sanitize=False, removeHs=False)
    if rdkit_mol is None:
        sdfstring = io.EntityToSDFStr(ent).strip()
        rdkit_mol = Chem.MolFromMolBlock(sdfstring, sanitize=False)
    if rdkit_mol is not None:
        rdkit_mol = Chem.RemoveAllHs(rdkit_mol, sanitize=False)
    return rdkit_mol


# ---------------------------------------------------------------------------
# Ligand bond order detection and assignment
# ---------------------------------------------------------------------------


class MissingBondOrderError(ValueError):
    """Raised when a CIF file has ligands with unresolvable bond orders."""

    pass


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


def _is_known_compound(comp_id: str) -> bool:
    """Check if a component ID is known to the CCD compound library."""
    return _COMPOUND_LIB.FindCompound(comp_id) is not None


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

    unknown = set()
    for comp_id in hetatm_ids:
        if comp_id in cif_bond_ids:
            continue
        if _is_known_compound(comp_id):
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

    ent = io.LoadMMCIF(str(cif_path)).Select("")

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

    for comp_id, smiles in to_process.items():
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError(f"Invalid SMILES for {comp_id}: {smiles}")

        ligand_view = ent.Select(f"rname={comp_id}")
        if not ligand_view.IsValid() or ligand_view.GetAtomCount() == 0:
            raise ValueError(f"No atoms found for component {comp_id} in CIF")

        ligand_ent = mol.CreateEntityFromView(ligand_view, True)

        rdkit_mol = ost_ent_to_rdkit_mol(ligand_ent)
        if rdkit_mol is None:
            raise ValueError(f"Could not parse ligand {comp_id} as RDKit mol")
        fixed_mol = mol_assigned_bond_orders_by_template(template, rdkit_mol)

        atom_names = [a.name.strip() for a in ligand_ent.atoms]

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
