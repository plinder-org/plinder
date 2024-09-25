# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
import gzip
from pathlib import Path
from typing import Any, Union

import numpy as np
from biotite import TextFile
from biotite.structure import get_residues
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdDepictor, rdMolDescriptors, rdRascalMCES

from plinder.core.structure.vendored import (
    _convert_resn_to_sequence_and_numbering,
    _get_structure_and_res_info,
    align_sequences,
    apply_mask,
)
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


def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import CIFFile

    return CIFFile


def atom_array_from_cif_file(
    structure: Path | _AtomArrayOrStack, use_author_fields: bool = True
) -> AtomArray | None:
    # TODO: this conversion may be unnecessary as the annotation is
    # Path | _AtomArrayOrStack
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
            arr = get_structure(
                mod, model=1, use_author_fields=use_author_fields, include_bonds=True
            )  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


def generate_input_conformer(
    template_mol: Chem.Mol, addHs: bool = False, minimize_maxIters: int = -1
) -> Chem.Mol:
    _mol = copy.deepcopy(template_mol)
    # need to add Hs to generate sensible conformers
    _mol = Chem.AddHs(_mol)
    # ps = AllChem.ETKDGv2()
    # try embedding molecule using ETKDGv2 (default)
    confid = AllChem.EmbedMolecule(
        _mol,
        # ps,
        useRandomCoords=True,
        useBasicKnowledge=True,
        maxAttempts=100,
        randomSeed=42,
    )
    if confid != -1:
        if minimize_maxIters > 0:
            # molecule successfully embedded - minimize
            success = AllChem.MMFFOptimizeMolecule(_mol, maxIters=minimize_maxIters)
            # 0 if the optimization converged,
            # -1 if the forcefield could not be set up,
            # 1 if more iterations are required.
            if success == 1:
                log.info(
                    f"generate_conformer: MMFFOptimizeMolecule - more iterations are required, doubling the steps (2x {minimize_maxIters})"
                )
                # extend optimization to double the steps (extends by the same amount)
                AllChem.MMFFOptimizeMolecule(_mol, maxIters=minimize_maxIters)
            elif success == -1:
                log.warning(
                    "generate_conformer: MMFFOptimizeMolecule - the forcefield could not be set up"
                )

    else:
        # this means EmbedMolecule failed
        log.warning(
            "generate_conformer: default EmbedMolecule - failed, try using useBasicKnowledge=False"
        )
        # try less optimal approach
        confid = AllChem.EmbedMolecule(
            _mol,
            useRandomCoords=True,
            useBasicKnowledge=False,
            maxAttempts=100,
            randomSeed=42,
        )
        if confid == -1:
            # if that still fails - try generating just 2D conformer
            log.warning(
                "generate_conformer: EmbedMolecule - failed, trying rdDepictor.Compute2DCoords instead"
            )
            confid = rdDepictor.Compute2DCoords(_mol)

    # verify that mol has conformers
    if _mol.GetNumConformers() == 0:
        raise ValueError("Could not generate conformer")

    if not addHs:
        # remove Hs if they should not be kept
        _mol = params_removeHs(_mol)

    return _mol


def params_removeHs(mol: Chem.Mol) -> Chem.Mol:
    params = Chem.rdmolops.RemoveHsParameters()
    params.removeIsotopes = True
    params.removeDegreeZero = True
    params.removeHigherDegrees = True
    params.removeOnlyHNeighbors = True
    params.removeNontetrahedralNeighbors = True
    return Chem.rdmolops.RemoveHs(mol, params, sanitize=False)


def match_ligands(
    input_smiles: str,
    resolved_sdf: str | Path,
) -> tuple[Chem.Mol, Chem.Mol, tuple[_AtomArrayOrStack, _AtomArrayOrStack]]:
    template_mol = Chem.MolFromSmiles(input_smiles)
    resolved_mol = Chem.MolFromMolFile(resolved_sdf.__str__())
    atom_order_stacks = get_template_to_mol_matches(template_mol, resolved_mol)
    return template_mol, resolved_mol, atom_order_stacks


def get_template_to_mol_matches(
    template: Chem.Mol, mol: Chem.Mol
) -> tuple[_AtomArrayOrStack, _AtomArrayOrStack]:
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
    if len(atom_matches) < min(numHA_template, numHA_mol):
        # if not complete molecule is matched
        match_mol = copy.deepcopy(mol)
        ref_mol = copy.deepcopy(template)

        log.warning(
            "get_template_to_mol_matches: could not match template fully - retry with unmatched bonds set as UNSPECIFIED"
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
            if len(results2[0].atomMatches()) > len(results[0].atomMatches()):
                results = results2

        # if still not fully matched - attempt one more!
        if len(results[0].atomMatches()) < min(numHA_template, numHA_mol):
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in match_mol.GetBonds()]
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in ref_mol.GetBonds()]
            # run again
            results3 = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)
            if len(results3[0].atomMatches()) > len(results[0].atomMatches()):
                results = results3

    # convert to atom order array stacks
    template_atom_order_stack1 = np.array(
        [[ix_templ for _, ix_templ in match.atomMatches()] for match in results]
    )
    mol_atom_order_stack2 = np.array(
        [[ix_mol for ix_mol, _ in match.atomMatches()] for match in results]
    )
    return template_atom_order_stack1, mol_atom_order_stack2


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
    ref_mol: Mol, matching_indices: tuple[int, ...]
) -> NDArray[np._int]:
    mask = np.zeros(len(ref_mol.GetAtoms()))
    for atm_idx in range(len(mask)):
        if atm_idx in matching_indices:
            mask[atm_idx] = 1
    return mask


def make_one_hot_atom_features(atom_name: list[str]) -> list[int]:
    allowed_atom_names = ["C", "N", "O", "S", "P"]
    striped_atom_name = "".join(filter(lambda x: not x.isdigit(), atom_name))[0]
    return [1 if striped_atom_name == atm else 0 for atm in allowed_atom_names]


def _convert_pdb_atom_name_to_elem_symbol(atom_name: str) -> str:
    return "".join(filter(lambda x: not x.isdigit(), atom_name))[0]


def get_per_residue_mask(
    residue_reference_atom_list: list[str], atom_list: list[str]
) -> list[int]:
    res_mask = [1 if i in atom_list else 0 for i in residue_reference_atom_list]
    return res_mask


def get_per_residue_atoms(
    atom_array: _AtomArrayOrStack, resi: int, resn: str
) -> NDArray[np.str_]:
    return atom_array[
        (atom_array.res_id == resi) & (atom_array.res_name == resn)
    ].atom_name


def make_atom_mask(
    atom_array: _AtomArrayOrStack, seq_res: str, seq_mask: list[int]
) -> list[int]:
    seq_res_three_aa = [pc.ONE_TO_THREE[aa] for aa in seq_res]
    resi, resn = get_residues(atom_array)
    residue_tuple = list(zip(tuple(resi), tuple(resn)))

    atom_mask = []
    resolved_residue_start = 0
    for seq_, mask in zip(seq_res_three_aa, seq_mask):
        if mask == 0:
            atom_mask.append([0 for i in range(len(pc.ORDERED_AA_FULL_ATOM[seq_]))])
        else:
            resi, resn = residue_tuple[resolved_residue_start]
            atom_mask.append(
                get_per_residue_mask(
                    pc.ORDERED_AA_FULL_ATOM[resn],
                    get_per_residue_atoms(atom_array, resi, resn),
                )
            )
    return [atm for res in atom_mask for atm in res]


def _stack_atom_array_features(
    atom_arr: _AtomArrayOrStack,
    atom_arr_feat: str,
    chain_order_list: list[str] | None,
) -> list[NDArray[np.int_ | np.str_ | np.float_]]:
    assert chain_order_list is not None
    return [
        getattr(atom_arr[atom_arr.chain_id == chain], atom_arr_feat)
        for chain in chain_order_list
    ]


def _stack_ligand_feat(
    feat_dict: dict[str, Any], chain_order_list: list[str] | None
) -> list[list[list[int]]]:
    assert chain_order_list is not None
    return [feat_dict[chain] for chain in chain_order_list]


def _one_hot_encode_stack(
    stack: list[NDArray],
    feature_dict: dict[str, int],
    unknown_name_filler: str,
) -> list[NDArray]:
    feat_array = []
    unknown_name_filler_value = feature_dict[unknown_name_filler]
    for per_chain_feat in stack:
        feat_array_by_chain = np.zeros(
            (len(per_chain_feat), len(set(list(feature_dict.values()))))
        )
        for index, value in enumerate(per_chain_feat):
            feat_array_by_chain[
                index, feature_dict.get(value, unknown_name_filler_value)
            ] = 1.0
        feat_array.append(feat_array_by_chain)
    return feat_array


def _sequence_full_atom_type_array(
    input_sequences: dict[str, str]
) -> dict[str, NDArray]:
    """Resolved sequence full atom features."""
    seq_atom_dict = {}
    for chain, sequence in input_sequences.items():
        feat = []
        for res in sequence:
            for atom in pc.ORDERED_AA_FULL_ATOM[pc.ONE_TO_THREE[res]]:
                feat.append(_convert_pdb_atom_name_to_elem_symbol(atom))
        seq_atom_dict[chain] = np.array(feat)
    return seq_atom_dict
