# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Generator

import numpy as np
import pandas as pd
from biotite.structure.io import load_structure
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from plinder.core.structure import diffdock_utils
from plinder.core.structure.contacts import get_atom_neighbors
from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.rdkit_utils import fix_valency_issues

if TYPE_CHECKING:
    from plinder.data.utils.annotations.aggregate_annotations import Entry

LOG = setup_logger(__name__)


def ligand_is_rdkit_loadable(sdf_path: Path) -> bool:
    """Check if structure is loadable by rdkit

    Parameters
    ----------
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    bool
        True if loadable, False otherwise.
    """
    try:
        mol = next(Chem.SDMolSupplier(str(sdf_path)))
        if mol is not None:
            return True
        else:
            return False
    except RuntimeError:
        return False


def ligand_is_rdkit_loadable_with_fix(sdf_path: Path) -> bool:
    """Check if structure is loadable by rdkit after fixing

    Parameters
    ----------
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    bool
        True if loadable, False otherwise.
    """
    mol = next(Chem.SDMolSupplier(str(sdf_path), sanitize=False))
    try:
        mol = fix_valency_issues(mol)
        if mol is not None:
            return True
        else:
            return False

    except RuntimeError:
        return False


def ligand_is_obabel_loadable(sdf_path: Path) -> bool:
    """Check if structure is loadable by openbabel

    Parameters
    ----------
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    bool
        True if loadable, False otherwise.
    """

    obconversion = ob.OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = ob.OBMol()
    return bool(obconversion.ReadFile(obmol, str(sdf_path)))


def ligand_is_obabel_loadable_with_rdkit_fix(sdf_path: Path) -> bool:
    """Check if structure is loadable by openbabel after fixing

    Parameters
    ----------
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    bool
        True if loadable, False otherwise.
    """
    obconversion = ob.OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = ob.OBMol()
    if obconversion.ReadFile(obmol, str(sdf_path)):
        return True
    else:
        mol = next(Chem.SDMolSupplier(str(sdf_path), sanitize=False))
        try:
            mol = fix_valency_issues(mol)
            if mol is not None:
                fixed_sdf_str = Chem.MolToMolBlock(mol)
                return bool(obconversion.ReadString(obmol, fixed_sdf_str))
            else:
                return False
        except Exception:
            return False


def ligand_matches_smiles_atom_num(smiles: str, sdf_path: Path) -> bool:
    """Check if atom number in ligand sdf matches the one from annotation smiles

    Parameters
    ----------
    smiles : str
        smiles from annotation foile
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    bool
        True if matching, False otherwise.
    """
    mol = next(Chem.SDMolSupplier(str(sdf_path), sanitize=False))
    try:
        mol = fix_valency_issues(mol)
    except Exception:
        return False
    if mol is None:
        return False
    try:
        target_mol = Chem.MolFromSmiles(smiles, sanitize=False)
        target_mol = fix_valency_issues(target_mol)
    except Exception:
        return False
    if target_mol is None:
        return False
    else:
        # atms = [i.GetSymbol() for i in mol.GetAtoms() if i.GetSymbol() != "H"]
        # target_atms = [j.GetSymbol() for j in target_mol.GetAtoms()]
        return bool(
            len(Chem.RemoveHs(mol).GetAtoms())
            == len(Chem.RemoveHs(target_mol).GetAtoms())
        )


def get_molvs_ligand_validation(sdf_path: Path) -> list[str]:
    """Check for errors that could be detected by MolVS

    Parameters
    ----------
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    list[str]
        [] if not validation error.
    """
    validations = [
        rdMolStandardize.NoAtomValidation(),
        rdMolStandardize.FragmentValidation(),
        rdMolStandardize.NeutralValidation(),
    ]
    mol = fix_valency_issues(next(Chem.SDMolSupplier(str(sdf_path), sanitize=False)))
    vm = rdMolStandardize.MolVSValidation(validations)
    return list(vm.validate(mol))


def get_rdkit_ligand_validation(sdf_path: Path) -> list[str]:
    """Check for errors that could be detected by RDKit

    Parameters
    ----------
    sdf_path : Path
        Path to ligand sdf

    Returns
    -------
    list[str]
        [] if not validation error.
    """
    mol = fix_valency_issues(next(Chem.SDMolSupplier(str(sdf_path), sanitize=False)))
    vm = rdMolStandardize.RDKitValidation()
    return list(vm.validate(mol))


def ligand_positions_correct(
    center_of_mass: np.ndarray[Any, Any], sdf_path: Path
) -> bool:
    """Check if ligand position is maintained

    Parameters
    ----------
    center_of_mass : np.ndarray[Any, Any]
        Center of mass from annotation
    sdf_path : Path
        Path to sdf ligand

    Returns
    -------
    bool
        True if position is maintained, otherwise False
    """
    mol = next(Chem.SDMolSupplier(str(sdf_path), sanitize=False))
    mol = fix_valency_issues(mol)
    conf = mol.GetConformer()
    return bool(
        np.allclose(
            np.mean(conf.GetPositions(), axis=0), center_of_mass, atol=0.1, rtol=0.1
        )
    )


def file_loadbable_via_biotite(structure_path: Path) -> bool:
    """Ligand is loadable by biotite

    Parameters
    ----------
    structure_path : Path
        Path to stucture (could be protein/ligand/complex)

    Returns
    -------
    bool
       True if loadable, otherwise False
    """

    try:
        if structure_path.suffix == ".sdf":
            arr = load_structure(structure_path)
        else:
            arr = load_structure(structure_path, use_author_fields=False)
        return bool(arr.coord.shape[0] > 0)
    except Exception:
        return False


def ligand_protein_neighbor_still_preserved(
    neighbor_chains: set[str], protein_file: Path, ligand_file: Path, radius: int = 6
) -> bool:
    """
    Check if all ligand protein neighbors are saved
    in protein structure

    Parameters
    ----------
    neighbor_chains : set[str]
        Set of protein chain ids neighboring to ligand
    protein_file : Path
        Path to protein structure
    ligand_file : Path
        Path to ligand sdf
    radius : int
        Radius around ligand to extract neighbors

    Returns
    -------
    bool
        True if all neighbors are present, otherwise False
    """

    protein_arr = load_structure(protein_file, use_author_fields=False)
    ligand_arr = load_structure(ligand_file)
    neighbor_arr = get_atom_neighbors(protein_arr, ligand_arr, radius=radius)
    return {neighbor.chain_id for neighbor in neighbor_arr} == neighbor_chains


def ligand_protein_neighbor_still_preserved_complex(
    neighbor_chains: set[str],
    complex_file: Path,
    protein_chains: set[str],
    ligand_chains: set[str],
    radius: int = 10,
) -> bool:
    """
    Check if all ligand protein neighbors are saved
    in complex structure

    Parameters
    ----------
    neighbor_chains : set[str]
        Set of protein chain ids neighboring to ligand
    complex_file : Path
        Path to complex structure
    radius : int
        Radius around ligand to extract neighbors

    Returns
    -------
    bool
        True if all neighbors are present, otherwise False
    """
    complex_arr = load_structure(complex_file, use_author_fields=False)
    # protein_arr = complex_arr[~complex_arr.hetero]
    # ligand_arr = complex_arr[complex_arr.hetero]
    protein_arr = complex_arr[np.isin(complex_arr.chain_id, list(protein_chains))]
    ligand_arr = complex_arr[np.isin(complex_arr.chain_id, list(ligand_chains))]
    neighbor_arr = get_atom_neighbors(protein_arr, ligand_arr, radius=radius)
    return {str(neighbor.chain_id) for neighbor in neighbor_arr} == neighbor_chains


def all_ligand_chains_present(ligand_chains: set[str], complex_file: Path) -> bool:
    """
    Check if all ligand chains are correctly saved

    Parameters
    ----------
    ligand_chains : set[str]
        Set of ligand chain ids from annotation
    complex_file : Path
        Path to complex structure

    Returns
    -------
    bool
        True if all ligand chains are present, otherwise False
    """
    complex_arr = load_structure(complex_file, use_author_fields=False)
    # ligand_arr = complex_arr[complex_arr.hetero]
    # chains = {ch for ch in ligand_arr.chain_id}
    chains = {ch for ch in complex_arr.chain_id}
    return bool(ligand_chains.issubset(chains))


def all_protein_chains_present(protein_chains: set[str], complex_file: Path) -> bool:
    """
    Check if all protein chains are correctly saved

    Parameters
    ----------
    protein_chains : set[str]
        Set of protein chain ids from annotation
    complex_file : Path
        Path to complex structure

    Returns
    -------
    bool
        True if all protein chains are present, otherwise False
    """
    complex_arr = load_structure(complex_file, use_author_fields=False)
    # protein_arr = complex_arr[~complex_arr.hetero]
    # chains = {ch for ch in protein_arr.chain_id}
    chains = {ch for ch in complex_arr.chain_id}
    return bool(protein_chains.issubset(chains))


def ligand_is_diffdock_loadable(ligand_file: Path) -> bool:
    try:
        lig = diffdock_utils.read_molecule(str(ligand_file))
        diffdock_utils.get_lig_graph_with_matching(lig)
        return True
    except Exception:
        return False


def run_structure_checks(system_dict: dict[str, Any]) -> list[dict[str, Any]]:
    """Aggregate all the checks for a single system

    Parameters
    ----------
    system_dict : dict[str, Any]
        System dictionary

    Returns
    -------
    list[dict[str, Any]]
    """
    system_checks = []
    complex_path: Any = system_dict["complex_path"]
    for ligand_chain, ligand in system_dict["ligands"].items():
        ligand_path: Path = ligand["ligand_path"]
        protein_chains: Any = system_dict["protein_chains"]
        ligand_chains: Any = system_dict["ligand_chains"]
        system_checks.append(
            {
                "system_id": complex_path.parent.name,
                "ligand_instance": ligand["ligand_instance"],
                "ligand_asym_id": ligand["ligand_asym_id"],
                "ligand_is_rdkit_loadable": ligand_is_rdkit_loadable(ligand_path),
                "ligand_is_rdkit_loadable_with_fix": ligand_is_rdkit_loadable_with_fix(
                    ligand_path
                ),
                "ligand_is_obabel_loadable": ligand_is_obabel_loadable(ligand_path),
                "ligand_is_obabel_loadable_with_rdkit_fix": ligand_is_obabel_loadable_with_rdkit_fix(
                    ligand_path
                ),
                # "ligand_is_diffdock_loadable": ligand_is_diffdock_loadable(ligand_path),
                "ligand_matches_smiles_atom_num": ligand_matches_smiles_atom_num(
                    ligand["ligand_resolved_smiles"], ligand_path
                ),
                "ligand_positions_correct": ligand_positions_correct(
                    ligand["center_of_mass"], ligand_path
                ),
                "ligand_loadbable_via_biotite": file_loadbable_via_biotite(ligand_path),
                "ligand_molvs_validation": get_molvs_ligand_validation(ligand_path),
                "ligand_rdkit_validation": get_rdkit_ligand_validation(ligand_path),
                "complex_loadbable_via_biotite": file_loadbable_via_biotite(
                    complex_path
                ),
                "ligand_protein_neighbor_still_preserved_complex": ligand_protein_neighbor_still_preserved_complex(
                    protein_chains, complex_path, protein_chains, ligand_chains
                ),
                "all_ligand_chains_present": all_ligand_chains_present(
                    ligand_chains, complex_path
                ),
                "all_protein_chains_present": all_protein_chains_present(
                    protein_chains, complex_path
                ),
            }
        )
    return system_checks


def prepare_system_dict(
    system_structure_path: Path,
    entry: "Entry",
) -> Generator[dict[str, Any], None, None]:
    for sys_id, annotation in entry.systems.items():
        system_folder = sys_id
        ligands = {}
        if (system_structure_path / system_folder).exists():
            for lig in annotation.ligands:
                ligand_chain = f"{lig.instance}.{lig.asym_id}"
                ligands[ligand_chain] = {
                    "ligand_path": system_structure_path
                    / system_folder
                    / "ligand_files"
                    / f"{ligand_chain}.sdf",
                    "ligand_instance": lig.instance,
                    "ligand_asym_id": lig.asym_id,
                    "ligand_smiles": lig.smiles,
                    "ligand_resolved_smiles": lig.resolved_smiles,
                    "center_of_mass": lig.centroid,
                }
            complex_file = system_structure_path / system_folder / "system.cif"
            receptor_path = system_structure_path / system_folder / "receptor.cif"
            system_dict = {
                "complex_path": complex_file,
                "receptor_path": receptor_path,
                "ligands": ligands,
                "protein_chains": set(annotation.protein_chains_asym_id),
                "ligand_chains": set(annotation.ligand_chains),
            }
        else:
            continue
        yield system_dict


def prepare_system_dicts(
    system_structure_path: Path, all_entries_json: Path
) -> Generator[dict[str, Any], None, None]:
    """
    Prepare dictionary input

    Parameters
    ----------
    system_structure_path : Path
        Folder for all system structure
    all_entries_json : Path
        All entries jsons

    Returns
    -------
    list[dict[str, Any]]
    """
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    with open(all_entries_json, "rb") as f:
        data = json.load(f)
    for entry in data:
        yield from prepare_system_dict(
            system_structure_path, Entry.model_validate(entry)
        )


def run_all_checks(system_structure_path: Path, all_entries_json: Path) -> pd.DataFrame:
    """
    Run all checks

    Parameters
    ----------
    system_structure_path : Path
        Folder for all system structure
    all_entries_json : Path
        All entries jsons

    Returns
    -------
    list[dict[str, Any]]
    """
    system_list = prepare_system_dicts(system_structure_path, all_entries_json)
    # TODO: Parallelize this, currently having issues with multiprocessing
    results = []
    for sys_dict in system_list:
        results.append(run_structure_checks(sys_dict))
    return pd.DataFrame([j for i in results for j in i])
