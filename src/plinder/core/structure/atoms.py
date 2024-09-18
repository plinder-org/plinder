# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
import gzip
import os
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

import biotite.structure as struc
import numpy as np
from biotite import TextFile
from biotite.structure import get_residues
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from numpy.typing import NDArray
from pinder.core.structure.atoms import (
    _convert_resn_to_sequence_and_numbering,
    _get_structure_and_res_info,
    align_sequences,
    apply_mask,
)
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdMolDescriptors, rdRascalMCES

from plinder.core.utils import constants as pc
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


# TODO: WHY DO WE NEED THIS?
# def biotite_pdbfile() -> TextFile:
#     from biotite.structure.io.pdb import PDBFile

#     return PDBFile
# def atom_array_from_pdb_file(structure: Path | AtomArray) -> _AtomArrayOrStack | None:
#     if isinstance(structure, str):
#         structure = Path(structure)

#     if isinstance(structure, Path):
#         reader = biotite_pdbfile()
#         try:
#             model = reader.read(str(structure))
#             arr = model.get_structure(model=1)  # noqa
#             return arr
#         except Exception as e:
#             log.error(f"Unable to parse {structure}! {str(e)}")
#     return structure


def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import PDBxFile

    return PDBxFile


def atom_array_from_cif_file(
    structure: Path | _AtomArrayOrStack, use_author_fields: bool = True
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


# TODO: WHY DO WE NEED THIS?
def atom_array_to_rdkit_mol(
    lig_atom_array: _AtomArrayOrStack,
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


def generate_input_conformer(template_mol: Mol, addHs: bool = True) -> Mol:
    _mol = copy.deepcopy(template_mol)
    if addHs:
        _mol = Chem.AddHs(_mol, addCoords=True)
    conformer_atom_matches = get_template_to_mol_matches(template_mol, _mol)
    ps = AllChem.ETKDGv2()
    ps.useRandomCoords = True
    AllChem.EmbedMolecule(_mol, ps)
    AllChem.MMFFOptimizeMolecule(_mol, confId=0)
    return _mol, conformer_atom_matches


# TODO: check - if this make sense!!!
# def match_ligands(
#     resolved_smiles: str,
#     unresolved_sdf: Path,
#     add_hydrogen: bool = False
# ) -> tuple[tuple[int, ...], Chem.rdchem.Mol, Chem.rdchem.Mol]:
#     try:
#         resolved_mol = Chem.MolFromSmiles(resolved_smiles)
#         unresolved_mol = next(Chem.SDMolSupplier(unresolved_sdf))
#         if add_hydrogen:
#             resolved_mol = Chem.AddHs(resolved_mol, addCoords=True)
#             unresolved_mol = Chem.AddHs(unresolved_mol, addCoords=True)
#         AllChem.ConstrainedEmbed(resolved_mol, unresolved_mol)
#     except Exception:
#         # TODO: Need to figure out how to handle failure, for now,
#         # set it to unresolved
#         resolved_mol = unresolved_mol = next(Chem.SDMolSupplier(unresolved_sdf))
#     # Returns all possible set of matches in case of symmetric molecules
#     matches = resolved_mol.GetSubstructMatches(unresolved_mol)
#     return matches[0], resolved_mol, unresolved_mol


def match_ligands(
    input_smiles: str,
    resolved_sdf: Path,
) -> tuple[Chem.Mol, Chem.Mol, list[dict[int, int]]]:
    template_mol = Chem.MolFromSmiles(input_smiles)
    resolved_mol = Chem.MolFromMolFile(resolved_sdf.__str__())
    atom_matches = get_template_to_mol_matches(template_mol, resolved_mol)
    return template_mol, resolved_mol, atom_matches


def get_template_to_mol_matches(
    template: Chem.Mol, mol: Chem.Mol
) -> list[dict[int, int]]:
    """
    Function that works a lot like get_matched_template but can better deal with fragmented molecules
    """
    rascal_opts = rdRascalMCES.RascalOptions()
    rascal_opts.similarityThreshold = 0.1
    rascal_opts.allBestMCESs = True
    rascal_opts.returnEmptyMCES = True
    rascal_opts.completeAromaticRings = False
    rascal_opts.ringMatchesRingOnly = False
    rascal_opts.timeout = 20

    results = rdRascalMCES.FindMCES(mol, template, rascal_opts)
    atom_matches = np.array(results[0].atomMatches())
    bond_matches = np.array(results[0].bondMatches())

    numHA_template = rdMolDescriptors.CalcNumHeavyAtoms(template)
    numHA_mol = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    print("AAA", mol, template)
    if len(atom_matches) < min(numHA_template, numHA_mol):
        # if not complete molecule is matched
        match_mol = copy.deepcopy(mol)
        ref_mol = copy.deepcopy(template)

        log.warning(
            "get_mol_to_template_matches: could not match template fully - retry with unmatched bonds set as UNSPECIFIED"
        )
        # set all unmatched bonds to UNSPECIFIED to help with the match
        if len(bond_matches):
            [
                b.SetBondType(Chem.BondType.UNSPECIFIED)
                for b in match_mol.GetBonds()
                if not b in bond_matches[:, 0]
            ]
            [
                b.SetBondType(Chem.BondType.UNSPECIFIED)
                for b in ref_mol.GetBonds()
                if not b in bond_matches[:, 1]
            ]
            # run again
            results2 = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)

        # if still not fully matched - attempt one more!
        if len(results2[0].atomMatches()) < min(numHA_template, numHA_mol):
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in match_mol.GetBonds()]
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in ref_mol.GetBonds()]
            # run again
            results3 = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)
            if len(results3[0].atomMatches()) > len(results[0].atomMatches()):
                if len(results3[0].atomMatches()) > len(results2[0].atomMatches()):
                    results = results3
                else:
                    results = results2
    return [
        {ix_smi: ix_mol for ix_smi, ix_mol in match.atomMatches()} for match in results
    ]


def get_residue_index_mapping_mask(
    ref_seqs: dict[str, str], subject_arr: _AtomArrayOrStack
) -> dict[str, list[int]]:
    mask_map = {}
    for reference_chain, ref_seq in ref_seqs.items():
        # refernece and subject chain are the same
        subject_chain = reference_chain
        subject_mask = subject_arr.chain_id == subject_chain
        subject_ch_arr = apply_mask(subject_arr, subject_mask)

        subj_info = _get_structure_and_res_info(subject_ch_arr)
        subj_seq, subj_numbering = _convert_resn_to_sequence_and_numbering(subj_info)
        alignments = align_sequences(
            ref_seq, subj_seq, ref_numbering=None, subject_numbering=subj_numbering
        )
        (
            _,
            _,
            ref_numbering_mapped,
            _,
        ) = alignments
        mask = np.zeros(len(ref_seq))
        for i in range(len(mask)):
            if (i + 1) in ref_numbering_mapped:
                mask[i] = 1
        mask_map[reference_chain] = mask
    return mask_map


def get_ligand_atom_index_mapping_mask(
    ref_mol: Mol, matching_indices: tuple[int]
) -> NDArray[np._int]:
    mask = np.zeros(len(ref_mol.GetAtoms()))
    for atm_idx in range(len(mask)):
        if atm_idx in matching_indices:
            mask[atm_idx] = 1
    return mask


def make_one_hot_atom_features(atom_name: str) -> list[int]:
    allowed_atom_names = ["C", "N", "O", "S", "P"]
    striped_atom_name = "".join(filter(lambda x: not x.isdigit(), atom_name))[0]
    return [1 if striped_atom_name == atm else 0 for atm in allowed_atom_names]


def get_per_residue_nask(
    residue_reference_atom_list: list[str], atom_list: list[int]
) -> list[int]:
    return [1 if i in atom_list else 0 for i in residue_reference_atom_list]


def get_per_residue_atoms(
    atom_array: _AtomArrayOrStack, resi: int, resn: str
) -> NDArray[np.str_]:
    return atom_array[
        (atom_array.res_id == resi) & (atom_array.res_name == resn)
    ].atom_name


def make_atom_mask(
    atom_array: _AtomArrayOrStack, seq_res: list[str], seq_mask: list[int]
) -> list[int]:
    seq_res = [pc.ONE_TO_THREE[aa] for aa in seq_res]
    resi, resn = get_residues(atom_array)
    residue_tuple = list(zip(tuple(resi), tuple(resn)))

    atom_mask = []
    resolved_residue_start = 0
    for seq_res, mask in zip(seq_res, seq_mask):
        if mask == 0:
            atom_mask.append([0 for i in range(len(pc.ORDERED_AA_FULL_ATOM[seq_res]))])
        else:
            resi, resn = residue_tuple[resolved_residue_start]
            atom_mask.append(
                get_per_residue_nask(
                    pc.ORDERED_AA_FULL_ATOM[resn],
                    get_per_residue_atoms(atom_array, resi, resn),
                )
            )
    return [atm for res in atom_mask for atm in res]
