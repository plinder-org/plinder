# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
import gzip
import os
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
from biotite import TextFile
from biotite.sequence.align import Alignment
from biotite.structure import connect_via_residue_names
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Mol

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


def generate_conformer(mol):
    _mol = copy.deepcopy(mol)
    ps = AllChem.ETKDGv2()
    ps.useRandomCoords = True
    AllChem.EmbedMolecule(_mol, ps)
    AllChem.MMFFOptimizeMolecule(_mol, confId=0)
    return _mol


def match_ligands(
    resolved_smiles: str, unresolved_sdf: Path, add_hydrogen
) -> tuple[tuple[int, ...], Chem.rdchem.Mol, Chem.rdchem.Mol]:
    try:
        resolved_mol = Chem.MolFromSmiles(resolved_smiles)
        unresolved_mol = next(Chem.SDMolSupplier(unresolved_sdf))
        if add_hydrogen:
            resolved_mol = Chem.AddHs(resolved_mol, addCoords=True)
            unresolved_mol = Chem.AddHs(unresolved_mol, addCoords=True)
        AllChem.ConstrainedEmbed(resolved_mol, unresolved_mol)
    except Exception:
        # TODO: Need to figure out how to handle failure, for now,
        # set it to unresolved
        resolved_mol = unresolved_mol = next(Chem.SDMolSupplier(unresolved_sdf))
    # Returns all possible set of matches in case of symmetric molecules
    matches = resolved_mol.GetSubstructMatches(unresolved_mol)
    return matches[0], resolved_mol, unresolved_mol


def _get_seq_aligned_structures(
    ref: Path | _AtomArrayOrStack, subject: Path | _AtomArrayOrStack
) -> tuple[_AtomArrayOrStack, _AtomArrayOrStack]:
    """Find common sequence and set identical residue numbering for common seq.
    This method only safely works on a single chain.
    """

    ref_info = _get_structure_and_res_info(ref)
    subj_info = _get_structure_and_res_info(subject)

    subj_resid_map, alns = _align_and_map_sequences(ref_info, subj_info)

    # Need to remove ref residues in order to match subject
    ref_structure, _, _ = ref_info
    subj_structure, _, _ = subj_info

    assert isinstance(ref_structure, (AtomArray, AtomArrayStack))
    assert isinstance(subj_structure, (AtomArray, AtomArrayStack))

    ref_mask = mask_from_res_list(ref_structure, list(subj_resid_map.values()))
    subj_mask = mask_from_res_list(subj_structure, list(subj_resid_map.keys()))

    ref_aligned = apply_mask(ref_structure, ref_mask)
    subject_aligned = apply_mask(subj_structure, subj_mask)
    subject_aligned.res_id = np.array(
        [subj_resid_map[ri] for ri in subject_aligned.res_id]
    )
    return ref_aligned, subject_aligned


def _align_structure_to_ref_sequence(
    ref_seq: str,
    subject: Path | _AtomArrayOrStack,
) -> tuple[_AtomArrayOrStack, str, list[int]]:
    """ """
    subj_info = _get_structure_and_res_info(subject)
    subj_seq, subj_numbering = _convert_resn_to_sequence_and_numbering(subj_info)

    # Setting ref_numbering=None sets it to 1 -> len(ref_seq)
    alignments = align_sequences(
        ref_seq, subj_seq, ref_numbering=None, subject_numbering=subj_numbering
    )
    (
        ref_seq_mapped,
        subj_seq_mapped,
        ref_numbering_mapped,
        subj_numbering_mapped,
    ) = alignments

    subject_common = f"{len(subj_seq_mapped)}/{len(subj_seq)}"
    ref_common = f"{len(ref_seq_mapped)}/{len(ref_seq)}"
    log.debug(
        f"{subject_common} residues in subject matched to "
        f"{ref_common} residues in ref"
    )

    # Renumber subject residues to match aligned reference
    subj_resid_map = dict(zip(subj_numbering_mapped, ref_numbering_mapped))

    # Need to remove ref residues in order to match subject
    subj_structure, _, _ = subj_info

    assert isinstance(subj_structure, (AtomArray, AtomArrayStack))

    subj_mask = mask_from_res_list(subj_structure, list(subj_resid_map.keys()))

    subject_aligned = apply_mask(subj_structure, subj_mask)
    subject_aligned.res_id = np.array(
        [subj_resid_map[ri] for ri in subject_aligned.res_id]
    )
    return subject_aligned, subj_seq, subj_numbering


def get_per_chain_seq_to_structure_alignments(
    ref_seqs: dict[str, str],
    subject_arr: _AtomArrayOrStack,
) -> dict[str, dict[int, int]]:
    """ """
    subject_aligned_list: list[AtomArrayStack] = []
    subj_seq_dict: dict[str, str] = {}
    subj_numbering_dict: dict[str, list[int]] = {}
    for reference_chain, ref_seq in ref_seqs.items():
        # refernece and subject chain are the same
        subject_chain = reference_chain
        subject_mask = subject_arr.chain_id == subject_chain
        subject_ch_arr = apply_mask(subject_arr, subject_mask)

        subject_aligned, subj_seq, subj_numbering = _align_structure_to_ref_sequence(
            ref_seq, subject_ch_arr
        )
        subject_aligned_list.append(subject_aligned)
        subj_seq_dict[subject_chain] = subj_seq
        subj_numbering_dict[subject_chain] = subj_numbering
        combined_renumbered_arr = subject_aligned_list[0]
        for arr in subject_aligned_list[1:]:
            combined_renumbered_arr += arr
    return combined_renumbered_arr, subj_seq_dict, subj_numbering_dict


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
        for i in len(mask):
            if (i + 1) in ref_numbering_mapped:
                mask[i] = 1
        mask_map[reference_chain] = mask
        return mask_map


def _align_and_map_sequences(
    ref_info: tuple[_AtomArrayOrStack, list[int], list[str]],
    subj_info: tuple[_AtomArrayOrStack, list[int], list[str]],
) -> tuple[dict[int, int], tuple[str, str, list[int], list[int]]]:
    """Aligns two sequences and maps the numbering from the reference to the subject.

    This method only safely works on a single chain.
    """
    ref_seq, ref_numbering = _convert_resn_to_sequence_and_numbering(ref_info)
    subj_seq, subj_numbering = _convert_resn_to_sequence_and_numbering(subj_info)

    alignments = align_sequences(ref_seq, subj_seq, ref_numbering, subj_numbering)
    (
        ref_seq_mapped,
        subj_seq_mapped,
        ref_numbering_mapped,
        subj_numbering_mapped,
    ) = alignments

    subject_common = f"{len(subj_seq_mapped)}/{len(subj_seq)}"
    ref_common = f"{len(ref_seq_mapped)}/{len(ref_seq)}"
    log.debug(
        f"{subject_common} residues in subject matched to "
        f"{ref_common} residues in ref"
    )

    # Renumber subject residues to match aligned reference
    subj_resid_map = dict(zip(subj_numbering_mapped, ref_numbering_mapped))
    return subj_resid_map, alignments


def _get_paired_structures_and_chains(
    ref: Path | _AtomArrayOrStack,
    subject: Path | _AtomArrayOrStack,
) -> tuple[_AtomArrayOrStack, _AtomArrayOrStack, NDArray[np.str_], NDArray[np.str_]]:
    ref_arr: _AtomArrayOrStack = atom_array_from_pdb_file(ref)
    subject_arr: _AtomArrayOrStack = atom_array_from_pdb_file(subject)
    ref_chains = struc.get_chains(ref_arr)
    subject_chains = struc.get_chains(subject_arr)

    # Sometimes returns repeated chains
    _, sub_idx = np.unique(subject_chains, return_index=True)
    subject_chains = subject_chains[np.sort(sub_idx)]

    _, ref_idx = np.unique(ref_chains, return_index=True)
    ref_chains = ref_chains[np.sort(ref_idx)]
    if len(ref_chains) != len(subject_chains):
        raise ValueError("Ref and subject must have same number of chains!")

    if (ref_chains != subject_chains).any():
        log.error("Ref and subject have different chains!")
        log.warning(
            "Sequence alignment will assume this order: "
            f"{ref_chains}, {subject_chains} "
        )
    return ref_arr, subject_arr, ref_chains, subject_chains


def get_per_chain_seq_alignments(
    ref: Path | _AtomArrayOrStack,
    subject: Path | _AtomArrayOrStack,
) -> dict[str, dict[int, int]]:
    """Computes per-chain sequence alignments between reference and subject
    structures, and provides a mapping of residue numbers from the subject
    to the reference for each chain.

    This function processes each chain separately, aligns their sequences,
    and creates a mapping of residue numbering from the subject structure
    to the reference structure based on sequence alignment. It is designed
    to work with multi-chain structures, aligning each chain independently.

    Parameters
    ----------
    ref : Path | _AtomArrayOrStack
        The file path or structure array of the reference structure.
    subject : Path | _AtomArrayOrStack
        The file path or structure array of the subject structure to align with the reference.

    Returns
    -------
    dict[str, dict[int, int]]
        A dictionary where keys are chain identifiers from the subject structure and values are dictionaries mapping residue numbers in the subject structure to those in the reference structure for the aligned sequence segments.

    Raises
    ------
    ValueError
        If the reference and subject structures do not have the same number of chains or if the chains cannot be matched properly.

    """
    (
        ref_arr,
        subject_arr,
        ref_chains,
        subject_chains,
    ) = _get_paired_structures_and_chains(ref, subject)
    ch_res_maps: dict[str, dict[int, int]] = {}
    for ref_ch, subject_ch in zip(ref_chains, subject_chains):
        ref_mask = ref_arr.chain_id == ref_ch
        subject_mask = subject_arr.chain_id == subject_ch
        ref_ch_arr = apply_mask(ref_arr, ref_mask)
        subject_ch_arr = apply_mask(subject_arr, subject_mask)

        ref_info = _get_structure_and_res_info(ref_ch_arr)
        subj_info = _get_structure_and_res_info(subject_ch_arr)
        existing2new, alns = _align_and_map_sequences(ref_info, subj_info)
        ch_res_maps[subject_ch] = existing2new
    return ch_res_maps


def _get_paired_structures_and_chains(
    ref: Path | _AtomArrayOrStack,
    subject: Path | _AtomArrayOrStack,
) -> tuple[_AtomArrayOrStack, _AtomArrayOrStack, NDArray[np.str_], NDArray[np.str_]]:
    ref_arr: _AtomArrayOrStack = atom_array_from_pdb_file(ref)
    subject_arr: _AtomArrayOrStack = atom_array_from_pdb_file(subject)
    ref_chains = struc.get_chains(ref_arr)
    subject_chains = struc.get_chains(subject_arr)

    # Sometimes returns repeated chains
    _, sub_idx = np.unique(subject_chains, return_index=True)
    subject_chains = subject_chains[np.sort(sub_idx)]

    _, ref_idx = np.unique(ref_chains, return_index=True)
    ref_chains = ref_chains[np.sort(ref_idx)]
    if len(ref_chains) != len(subject_chains):
        raise ValueError("Ref and subject must have same number of chains!")

    if (ref_chains != subject_chains).any():
        log.error("Ref and subject have different chains!")
        log.warning(
            "Sequence alignment will assume this order: "
            f"{ref_chains}, {subject_chains} "
        )
    return ref_arr, subject_arr, ref_chains, subject_chains


def get_seq_alignments(
    ref_seq: str,
    subject_seq: str,
    gap_penalty: tuple[float, float] = (-10.0, -1.0),
) -> list[Alignment]:
    """Generate optimal global alignments between two sequences using BLOSUM62 matrix.

    Parameters
    ----------
    ref_seq : (str)
        The reference sequence.
    subject_seq : (str)
        The subject sequence.
    gap_penalty : (tuple[float, float], optional)
        A tuple consisting of the gap open penalty and the gap extension penalty.

    Returns
    -------
        list[Alignment]:
            A list of `Alignment` objects that represent the optimal global alignment(s).

    Note:
       The function uses the BLOSUM62 matrix by default for alignment scoring.
    """
    ref_protein_seq = seq.ProteinSequence(ref_seq)
    subject_protein_seq = seq.ProteinSequence(subject_seq)

    matrix = align.SubstitutionMatrix.std_protein_matrix()  # BLOSUM62 matrix
    # matrix = seq.align.SubstitutionMatrix(alph, alph, "IDENTITY")
    # alph = seq.ProteinSequence.alphabet
    alignments = align.align_optimal(
        subject_protein_seq,
        ref_protein_seq,
        matrix,
        gap_penalty=gap_penalty,
        terminal_penalty=False,
        local=False,
        max_number=1,
    )
    assert isinstance(alignments, list)
    return alignments


def get_seq_identity(
    ref_seq: str | None = None,
    subject_seq: str | None = None,
    alignments: list[Alignment] | None = None,
    gap_penalty: tuple[float, float] = (-10.0, -1.0),
) -> float:
    """Align an arbitrary sequence with the reference sequence
    Return sequence identity between the two.
    """

    all_none = all(arg is None for arg in [ref_seq, subject_seq, alignments])
    if all_none:
        raise ValueError("Must provide one of ref_seq and subject_seq or alignments!")

    if not alignments:
        assert isinstance(ref_seq, str)
        assert isinstance(subject_seq, str)
        alignments = get_seq_alignments(ref_seq, subject_seq, gap_penalty)

    aln = alignments[0]
    identity: float = align.get_sequence_identity(aln, "shortest")
    if identity < 0.90:
        log.debug(f"Alignment issue detected with identity of {identity}")
    return identity


def get_seq_aligned_structures(
    ref: Path | _AtomArrayOrStack,
    subject: Path | _AtomArrayOrStack,
) -> tuple[_AtomArrayOrStack, _AtomArrayOrStack]:
    """Computes per-chain sequence alignments between reference and subject
    structures, and creates a new set of structures where the residue numbers
    in the subject structure are renumbered to match those in the reference
    structure for the aligned sequence segments for each chain.

    This function processes each chain separately, aligns their sequences,
    and then constructs the multi-chain structure using the chain-aligned
    sequences.

    Parameters
    ----------
    ref : Path | _AtomArrayOrStack
        The file path or structure array of the reference structure.
    subject : Path | _AtomArrayOrStack
        The file path or structure array of the subject structure to align with the reference.


    Returns
    -------
    tuple
        A tuple containing:
        - Reference structure with aligned sequence and residue numbering (_AtomArrayOrStack).
        - Subject structure with aligned sequence and residue numbering (_AtomArrayOrStack).

    Raises
    ------
    ValueError
        If the reference and subject structures do not have the same number of chains or if the chains cannot be matched properly.
    """
    (
        ref_arr,
        subject_arr,
        ref_chains,
        subject_chains,
    ) = _get_paired_structures_and_chains(ref, subject)
    ref_arrays = []
    subject_arrays = []
    for ref_ch, subject_ch in zip(ref_chains, subject_chains):
        ref_mask = ref_arr.chain_id == ref_ch
        subject_mask = subject_arr.chain_id == subject_ch
        ref_ch_arr = apply_mask(ref_arr, ref_mask)
        subject_ch_arr = apply_mask(subject_arr, subject_mask)
        ref_aligned, subject_aligned = _get_seq_aligned_structures(
            ref_ch_arr, subject_ch_arr
        )
        ref_arrays.append(ref_aligned)
        subject_arrays.append(subject_aligned)

    # construct multi-chain objects again
    ref_aligned = ref_arrays.pop(0)
    subject_aligned = subject_arrays.pop(0)
    for ref_arr, sub_arr in zip(ref_arrays, subject_arrays):
        ref_aligned += ref_arr
        subject_aligned += sub_arr
    return ref_aligned, subject_aligned


def invert_chain_seq_map(
    seq_map: dict[str, dict[int, int]] | list[dict[str, dict[int, int]]],
) -> dict[str, dict[int, int]] | list[dict[str, dict[int, int]]]:
    def _invert(seq_map: dict[str, dict[int, int]]) -> dict[str, dict[int, int]]:
        rev_map: dict[str, dict[int, int]] = {ch: {} for ch in seq_map}
        for ch, rmap in seq_map.items():
            rev_map[ch] = {v: k for k, v in rmap.items()}
        return rev_map

    inverted: dict[str, dict[int, int]] | list[dict[str, dict[int, int]]]
    if isinstance(seq_map, list):
        inverted = []
        for sm in seq_map:
            inverted.append(_invert(sm))
    else:
        inverted = _invert(seq_map)
    return inverted


def resn2seq(resn: list[str]) -> str:
    three_to_one = pc.three_to_one_noncanonical_mapping
    seq_list = [three_to_one.get(three, "X") for three in resn]
    # Make sure we convert selenocysteine to cysteine
    # https://github.com/biotite-dev/biotite/blob/master/src/biotite/sequence/seqtypes.py#L364-L369
    # It appears that pyrrolysine (PYL->O) is also not supported - convert to lysine.
    seq_str = ""
    for s in seq_list:
        if s == "U":
            seq_str += "C"
        elif s == "O":
            seq_str += "K"
        else:
            seq_str += s
    return seq_str


def write_pdb(arr: AtomArray, filepath: Path) -> None:
    """Write AtomArray to a PDB file.

    Parameters
    ----------
    arr : AtomArray
        AtomArray to save to PDB file.
    filepath : Path
        Path to outut PDB file.

    Returns
    -------
    None

    """
    if not filepath.parent.is_dir():
        filepath.parent.mkdir(exist_ok=True, parents=True)

    strucio.save_structure(filepath, arr)


def mask_from_res_list(
    atoms: _AtomArrayOrStack, res_list: list[int]
) -> NDArray[np.bool_]:
    mask = np.isin(atoms.res_id, res_list)
    return mask


def apply_mask(atoms: _AtomArrayOrStack, mask: NDArray[np.bool_]) -> _AtomArrayOrStack:
    """Apply a boolean mask to an AtomArray or AtomArrayStack to filter atoms.

    Parameters
    ----------
    atoms : (AtomArray | AtomArrayStack)
        The atoms to be filtered.
    mask : NDArray[np.bool\_]
        The boolean mask that specifies which atoms to keep.

    Returns
    -------
        AtomArray | AtomArrayStack:
            The filtered atoms.
    """
    if isinstance(atoms, AtomArray):
        return atoms[mask]
    elif isinstance(atoms, AtomArrayStack):
        return atoms[..., mask]
    else:
        raise TypeError("atoms must be an AtomArray or AtomArrayStack")


def align_sequences(
    ref_seq: str,
    subject_seq: str,
    ref_numbering: list[int] | None = None,
    subject_numbering: list[int] | None = None,
) -> tuple[str, str, list[int], list[int]]:
    """Aligns two sequences and returns a tuple containing the aligned sequences
    and their respective numbering.

    Parameters
    ----------
    ref_seq : (str)
        The reference sequence to align to.
    subject_seq : (str)
        The subject sequence to be aligned.
    ref_numbering : (list[int], optional)
        List of residue numbers for the reference sequence.
    subject_numbering : (list[int], optional)
        List of residue numbers for the subject sequence.

    Returns:
        tuple: A tuple containing:
            - Aligned reference sequence (str)
            - Aligned subject sequence (str)
            - Numbering of the aligned reference sequence (list[int])
            - Numbering of the aligned subject sequence (list[int])

    Raises
    ------
        ValueError if the sequences cannot be aligned.
    """
    alignments = get_seq_alignments(ref_seq, subject_seq)
    get_seq_identity(alignments=alignments)
    aln = alignments[0]
    s = align.alignment.get_symbols(aln)

    if ref_numbering is None:
        ref_numbering = list(range(1, len(ref_seq) + 1))

    if subject_numbering is None:
        subject_numbering = list(range(1, len(subject_seq) + 1))

    # assigning reference numbering to all residues in a given sequence
    # (e.g. PDB chain)
    # p = subject protein sequence, u = ref sequence
    ref_numbering_mapped = []
    subject_numbering_mapped = []
    subject_sequence_mapped = ""
    ref_sequence_mapped = ""
    ui = -1
    pj = -1
    for p, u in zip(s[0], s[1]):
        if u:
            ui += 1
        if p:
            pj += 1
        if u and p and u != "X" and p != "X":
            # ref_numbering_mapped.append(ui + 1)  # 1-based
            ref_numbering_mapped.append(ref_numbering[ui])
            subject_numbering_mapped.append(subject_numbering[pj])
            ref_sequence_mapped += u
            subject_sequence_mapped += p
    return (
        ref_sequence_mapped,
        subject_sequence_mapped,
        ref_numbering_mapped,
        subject_numbering_mapped,
    )


def _get_structure_and_res_info(
    structure: Path | _AtomArrayOrStack,
) -> tuple[_AtomArrayOrStack, list[int], list[str]]:
    structure = atom_array_from_pdb_file(structure)
    numbering, resn = struc.get_residues(structure)
    return structure, numbering, resn


def _convert_resn_to_sequence_and_numbering(
    structure_info: tuple[_AtomArrayOrStack, list[int], list[str]],
) -> tuple[str, list[int]]:
    """Converts residue names to a sequence and extracts numbering."""
    _, numbering, resn = structure_info
    seq = resn2seq(resn)
    return seq, numbering
