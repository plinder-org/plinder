# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb_io
import biotite.structure.io.pdbx as pdbx
import numpy as np
from rdkit import Chem

# Define available names for receptor (protein/NA) and ligand chains in PDB format
PDB_RECEPTOR_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
PDB_LIGAND_CHAINS = PDB_RECEPTOR_CHAINS.lower() + "0123456789"
WATER_CHAIN_NAME = "_"


def save_ligands(
    atoms: struc.AtomArray,
    ligand_chain_ids: list[str],
    ligand_smiles: list[str],
    ligand_num_unresolved_heavy_atoms: list[int | None],
    output_folder: str | Path,
) -> None:
    """Save ligand SDF files from AtomArray.

    Parameters
    ----------
    atoms : AtomArray
        Full system atoms with bonds.
    ligand_chain_ids : list[str]
        Chain IDs identifying each ligand.
    ligand_smiles : list[str]
        Reference SMILES for each ligand.
    ligand_num_unresolved_heavy_atoms : list[int | None]
        Number of unresolved heavy atoms per ligand.
    output_folder : str or Path
        Directory to write SDF files.
    """
    from plinder.data.utils.annotations.cif_utils import atoms_to_rdkit_mol

    for chain_id, smiles, num_unresolved in zip(
        ligand_chain_ids,
        ligand_smiles,
        ligand_num_unresolved_heavy_atoms,
    ):
        lig_mask = atoms.chain_id == chain_id
        if not np.any(lig_mask):
            continue
        lig_atoms = atoms[lig_mask]
        try:
            rdkit_mol = atoms_to_rdkit_mol(lig_atoms)
        except Exception:
            continue
        if rdkit_mol is None:
            continue
        rdkit_mol.SetProp("_Name", chain_id)
        with Chem.SDWriter(str(Path(output_folder) / f"{chain_id}.sdf")) as w:
            w.write(rdkit_mol)


def save_pdb_file(
    full_system: struc.AtomArray,
    receptor_chains: list[str],
    ligand_chains: list[str],
    output_pdb_file: str | Path,
    output_mapping_file: str | Path,
    waters: dict[str, list[int]],
    water_mapping_file: str | Path,
) -> None:
    """Rename chains to PDB single-letter convention and save.

    Parameters
    ----------
    full_system : AtomArray
        System atoms (receptor + ligand, no waters yet).
    receptor_chains : list[str]
        Original receptor chain IDs (protein and/or nucleic acid).
    ligand_chains : list[str]
        Original ligand chain IDs.
    output_pdb_file : str or Path
        Path to output PDB file.
    output_mapping_file : str or Path
        Path to output chain mapping JSON.
    waters : dict[str, list[int]]
        Water chain IDs mapped to residue numbers.
    water_mapping_file : str or Path
        Path to output water mapping JSON.
    """
    # Remove waters from the main structure (added back separately)
    atoms = full_system[~struc.filter_solvent(full_system)]

    # Build chain renaming
    receptor_chain_index = 0
    ligand_chain_index = 0
    name_mapping: dict[str, str] = {}

    for original_name in np.unique(atoms.chain_id):
        if original_name in receptor_chains:
            final_name = PDB_RECEPTOR_CHAINS[receptor_chain_index]
            receptor_chain_index += 1
        elif original_name in ligand_chains:
            final_name = PDB_LIGAND_CHAINS[ligand_chain_index]
            ligand_chain_index += 1
        else:
            continue
        name_mapping[original_name] = final_name

    # Apply renaming
    new_chain_ids = atoms.chain_id.copy()
    for old, new in name_mapping.items():
        new_chain_ids[atoms.chain_id == old] = new
    atoms.chain_id = new_chain_ids

    # Add water residues
    water_mapping: dict[str, dict[int, int]] = {}
    if waters:
        water_atoms_list = []
        index = 1
        for chain_name, resnums in waters.items():
            water_mapping[chain_name] = {}
            chain_mask = full_system.chain_id == chain_name
            for resnum in resnums:
                res_mask = chain_mask & (full_system.res_id == resnum)
                if not np.any(res_mask):
                    continue
                water_res = full_system[res_mask].copy()
                water_res.chain_id[:] = WATER_CHAIN_NAME
                water_res.res_id[:] = index
                water_atoms_list.append(water_res)
                water_mapping[chain_name][int(resnum)] = index
                index += 1
        if water_atoms_list:
            water_arr = water_atoms_list[0]
            for wa in water_atoms_list[1:]:
                water_arr = water_arr + wa
            atoms = atoms + water_arr

    # Write PDB
    pdb_file = pdb_io.PDBFile()
    pdb_file.set_structure(atoms)
    pdb_file.write(str(output_pdb_file))

    with open(output_mapping_file, "w") as f:
        json.dump(name_mapping, f)
    if waters:
        with open(water_mapping_file, "w") as f:
            json.dump(water_mapping, f)


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
