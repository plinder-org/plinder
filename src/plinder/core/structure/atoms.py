# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any, Union

import numpy as np
from biotite import TextFile
from biotite.structure import get_residues
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from numpy.typing import NDArray
from rdkit.Chem import Mol

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
