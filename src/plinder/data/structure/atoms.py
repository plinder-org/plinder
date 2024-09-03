# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
import os
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

import biotite.structure as struc
from biotite import TextFile
from biotite.structure import connect_via_residue_names
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from rdkit import Chem
from rdkit.Chem import AllChem, Mol

from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)

DOCKQ_BACKBONE_ATOMS = ["C", "CA", "N", "O"]

THREE_LETTER_AA = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "SER",
    "THR",
    "VAL",
    "TYR",
    "TRP",
    "PRO",
    "HIS",
    "PHE",
]

BINDING_SITE_METALS = [
    "MG",
    "K",
    "MN",
    "NA",
    "ZN",
    "MG",
    "CA",
    "CD",
    "FE",
    "CU",
    "4MO",
    "CO",
]

_AtomArrayOrStack = Union[AtomArray, AtomArrayStack]


def biotite_pdbfile() -> TextFile:
    from biotite.structure.io.pdb import PDBFile

    return PDBFile


def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import PDBxFile

    return PDBxFile


def atom_array_from_pdb_file(structure: Path | AtomArray) -> _AtomArrayOrStack | None:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_pdbfile()
        try:
            model = reader.read(str(structure))
            arr = model.get_structure(model=1)  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


def atom_array_from_cif_file(
    structure: Path | AtomArray, use_author_fields: bool = True
) -> AtomArray | None:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_ciffile()
        try:
            if structure.suffix == ".gz":
                with gzip.open(str(structure), "rt", encoding="utf-8") as f:
                    mod = reader.read(f)
            else:
                mod = reader.read(structure)
            arr = get_structure(mod, model=1, use_author_fields=use_author_fields)  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


def atom_array_to_rdkit_mol(
    lig_atom_array: AtomArray,
    smiles: str,
    lig_resn: str = "LIG",
    chain_mapping: Optional[Dict[str, str]] = None,
    add_h: bool = True,
) -> Mol:
    # Change chain from two letters to one
    if chain_mapping is None:
        lig_atom_array.chain_id[:] = "X"
    else:
        lig_ch_id = list(set(lig_atom_array.chain_id))[0]
        lig_atom_array.chain_id[:] = chain_mapping.get(lig_ch_id, "X")

    temp1 = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    file = struc.io.pdbx.PDBxFile()
    struc.io.pdbx.set_structure(file, lig_atom_array, data_block=lig_resn)
    struc.io.save_structure(temp1.name, lig_atom_array)
    pdb_block = temp1.read()
    os.unlink(temp1.name)
    ligand_pose = Chem.MolFromPDBBlock(pdb_block)
    try:
        if smiles != "":
            # If we have smiles, use it to assign bond orders
            # Load original molecule from smiles
            template = Chem.MolFromSmiles(smiles)
            new_mol = AllChem.AssignBondOrdersFromTemplate(template, ligand_pose)
        else:
            new_mol = Chem.Mol(ligand_pose)
    except (ValueError, TypeError) as e:
        log.warning(f"Bond assignment to ligand failed...: {e}")
        # If bond assignment or addition of H fails
        new_mol = Chem.Mol(ligand_pose)
    try:
        if add_h:
            log.info("Adding Hs ...")
            new_mol_h = Chem.AddHs(new_mol, addCoords=True)
            return new_mol_h
        else:
            new_mol_h = Chem.Mol(new_mol)
            return new_mol_h
    except (ValueError, TypeError) as e:
        log.warning(f"Addition of H to ligand failed...: {e}")
        new_mol_h = Chem.Mol(new_mol)
        return new_mol_h


def add_bond_order_to_protein_atom_array(atom_array: AtomArray) -> AtomArray:
    bonds = connect_via_residue_names(atom_array)
    atom_array.bonds = bonds
    return atom_array
