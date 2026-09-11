# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, Union

import biotite.sequence as seq
import biotite.sequence.align as align
import numpy as np
from biotite import TextFile
from biotite.structure import get_residues
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import CIFFile, get_structure, set_structure
from numpy.typing import NDArray
from rdkit.Chem import Mol

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


def write_cif(arr: AtomArray, filepath: Path) -> None:
    """Write an AtomArray to an mmCIF file whose data block is the file stem."""
    filepath = Path(filepath)
    if filepath.suffix != ".cif":
        raise ValueError(f"mmCIF output must end in .cif, got {filepath}")
    filepath.parent.mkdir(exist_ok=True, parents=True)
    cif_file = CIFFile()
    set_structure(cif_file, arr, data_block=filepath.stem)
    cif_file.write(str(filepath))


def apply_mask(atoms: _AtomArrayOrStack, mask: NDArray[np.bool_]) -> _AtomArrayOrStack:
    """Filter the atoms of an AtomArray or AtomArrayStack with a boolean mask."""
    if isinstance(atoms, AtomArray):
        return atoms[mask]
    if isinstance(atoms, AtomArrayStack):
        return atoms[..., mask]
    raise TypeError("atoms must be an AtomArray or AtomArrayStack")


def resn2seq(resn: Iterable[str]) -> str:
    """Convert residue names to one-letter code, with ``X`` for unknown names.

    Selenocysteine and pyrrolysine become cysteine and lysine because
    biotite's protein alphabet has no ``U`` or ``O``.
    """
    three_to_one = pc.three_to_one_noncanonical_mapping
    substitutions = {"U": "C", "O": "K"}
    return "".join(
        substitutions.get(one, one)
        for one in (three_to_one.get(three, "X") for three in resn)
    )


def align_sequences(
    ref_seq: str,
    subject_seq: str,
    ref_numbering: Sequence[int] | None = None,
    subject_numbering: Sequence[int] | None = None,
) -> tuple[str, str, list[int], list[int]]:
    """Globally align two sequences and return their matched residues.

    The alignment uses BLOSUM62 with affine gap penalties and no terminal
    penalty.  Only positions where both sequences carry a known residue are
    kept, and returned as the two mapped sequences plus their residue
    numbering (1-based unless numbering is given).
    """
    alignment = align.align_optimal(
        seq.ProteinSequence(subject_seq),
        seq.ProteinSequence(ref_seq),
        align.SubstitutionMatrix.std_protein_matrix(),
        gap_penalty=(-10.0, -1.0),
        terminal_penalty=False,
        local=False,
        max_number=1,
    )[0]
    subject_symbols, ref_symbols = align.get_symbols(alignment)
    if ref_numbering is None:
        ref_numbering = list(range(1, len(ref_seq) + 1))
    if subject_numbering is None:
        subject_numbering = list(range(1, len(subject_seq) + 1))
    ref_mapped = ""
    subject_mapped = ""
    ref_numbering_mapped: list[int] = []
    subject_numbering_mapped: list[int] = []
    ref_index = -1
    subject_index = -1
    for subject_symbol, ref_symbol in zip(subject_symbols, ref_symbols):
        if ref_symbol:
            ref_index += 1
        if subject_symbol:
            subject_index += 1
        if ref_symbol and subject_symbol and "X" not in (ref_symbol, subject_symbol):
            ref_numbering_mapped.append(int(ref_numbering[ref_index]))
            subject_numbering_mapped.append(int(subject_numbering[subject_index]))
            ref_mapped += ref_symbol
            subject_mapped += subject_symbol
    return ref_mapped, subject_mapped, ref_numbering_mapped, subject_numbering_mapped


def get_residue_index_mapping_mask(
    ref_seqs: dict[str, str], subject_arr: _AtomArrayOrStack
) -> dict[str, NDArray[np.float64]]:
    """Mark which reference-sequence positions are resolved in each chain."""
    mask_map = {}
    for chain, ref_seq in ref_seqs.items():
        # The reference and subject chain IDs are the same.
        subject_ch_arr = apply_mask(subject_arr, subject_arr.chain_id == chain)
        subj_numbering, subj_resn = get_residues(subject_ch_arr)
        _, _, ref_numbering_mapped, _ = align_sequences(
            ref_seq,
            resn2seq(subj_resn),
            subject_numbering=subj_numbering.tolist(),
        )
        resolved = set(ref_numbering_mapped)
        mask = np.zeros(len(ref_seq))
        for i in range(len(mask)):
            if (i + 1) in resolved:
                mask[i] = 1
        mask_map[chain] = mask
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
        feat_array_by_chain = np.zeros((
            len(per_chain_feat),
            len(set(list(feature_dict.values()))),
        ))
        for index, value in enumerate(per_chain_feat):
            feat_array_by_chain[
                index, feature_dict.get(value, unknown_name_filler_value)
            ] = 1.0
        feat_array.append(feat_array_by_chain)
    return feat_array


def _sequence_full_atom_type_array(
    input_sequences: dict[str, str],
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
