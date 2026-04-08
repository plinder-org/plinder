# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Check and assign ligand bond orders in mmCIF files.

Cofolding tools (AlphaFold3, Boltz, Chai-1) output mmCIF files without
``_chem_comp_bond``. This module detects missing bond orders and assigns
them from user-supplied SMILES templates.

Only ligands unknown to the CCD library *and* missing from
``_chem_comp_bond`` are processed. Known compounds (ATP, NAD, HEM, etc.)
are skipped automatically.
"""

from __future__ import annotations

import logging
from pathlib import Path

import biotite.structure.io.pdbx as pdbx
from openbabel import pybel
from ost import conop, io, mol
from rdkit import Chem

from plinder.core.structure.smallmols_utils import (
    mol_assigned_bond_orders_by_template,
    params_removeHs,
)

LOG = logging.getLogger(__name__)
_COMPOUND_LIB = conop.GetDefaultLib()


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

        # Convert to RDKit mol via OpenBabel bond perception
        pdbstring = io.EntityToPDBStr(ligand_ent).strip()
        sdfstring = pybel.readstring("pdb", pdbstring).write("sdf")
        rdkit_mol = Chem.MolFromMolBlock(sdfstring, sanitize=False)
        if rdkit_mol is None:
            raise ValueError(f"Could not parse ligand {comp_id} as RDKit mol")

        rdkit_mol = params_removeHs(rdkit_mol)
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
