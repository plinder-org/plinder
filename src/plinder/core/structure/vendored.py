from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Union

import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
from biotite import TextFile
from biotite.sequence.align import Alignment
from biotite.structure.atoms import Atom, AtomArray, AtomArrayStack, stack
from numpy.typing import NDArray

import plinder.core.utils.constants as pc
from plinder.core.structure.models import BackboneDefinition
from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)

DOCKQ_BACKBONE_ATOMS = ["C", "CA", "N", "O"]


_AtomArrayOrStack = Union[AtomArray, AtomArrayStack]


def biotite_pdbxfile() -> TextFile:
    from biotite.structure.io.pdbx import CIFFile

    return CIFFile


def biotite_pdbfile() -> TextFile:
    from biotite.structure.io.pdb import PDBFile

    return PDBFile


@lru_cache(maxsize=None)
def rust_pdbfile() -> TextFile:
    # TODO: is this used anywhere?
    try:
        import fastpdb

        return fastpdb.PDBFile
    except ImportError:
        log.warning(
            "Requested fastpdb engine, but its not installed. "
            "Falling back to biotite"
        )
        return biotite_pdbfile()


def atom_array_from_pdb_file(
    structure: Path | AtomArray,
    extra_fields: list[str] | None = None,
) -> _AtomArrayOrStack:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        try:
            arr = strucio.load_structure(
                structure.as_posix(), extra_fields=extra_fields
            )
            assert isinstance(arr, (AtomArray, AtomArrayStack))
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
            raise
    return structure


def atom_vdw_radius(at: Atom) -> float:
    try:
        if at.hetero:
            rad = struc.info.vdw_radius_single(at.element)
        else:
            rad = struc.info.vdw_radius_protor(at.res_name, at.atom_name)
            if rad is None:
                rad = struc.info.vdw_radius_single(at.element)
    except Exception:
        rad = struc.info.vdw_radius_single("N")
    assert isinstance(rad, float)
    return rad


def backbone_mask(
    atoms: _AtomArrayOrStack, backbone_definition: str | BackboneDefinition
) -> NDArray[np.bool_]:
    if backbone_definition == "dockq":
        mask = np.isin(atoms.atom_name, DOCKQ_BACKBONE_ATOMS)
    else:
        mask = struc.filter_peptide_backbone(atoms)
    return mask


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


def filter_atoms(
    atoms: _AtomArrayOrStack,
    calpha_only: bool = False,
    backbone_only: bool = False,
    heavy_only: bool = True,
    backbone_definition: str | BackboneDefinition = "dockq",
) -> _AtomArrayOrStack:
    if calpha_only:
        atoms = apply_mask(atoms, atoms.atom_name == "CA")
    if backbone_only:
        atoms = apply_mask(atoms, backbone_mask(atoms, backbone_definition))
    if heavy_only:
        atoms = apply_mask(atoms, atoms.element != "H")
    return atoms


def get_backbone_atom_masks(
    native: AtomArray,
    decoy_stack: AtomArrayStack,
    backbone_only: bool,
    calpha_only: bool,
    backbone_definition: str | BackboneDefinition = "dockq",
) -> tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    # superimpose
    if backbone_only:
        model_mask = backbone_mask(decoy_stack, backbone_definition)
        native_mask = backbone_mask(native, backbone_definition)
    elif calpha_only:
        model_mask = decoy_stack.atom_name == "CA"
        native_mask = native.atom_name == "CA"
    else:
        model_mask = np.repeat(True, decoy_stack.shape[1])
        native_mask = np.repeat(True, native.shape[0])
    return native_mask, model_mask


def renumber_res_ids(arr: _AtomArrayOrStack) -> _AtomArrayOrStack:
    try:
        # biotite for py39 and below does not have this method
        arr.res_id = struc.create_continuous_res_ids(arr)
    except AttributeError:
        # Deprecated in biotite 0.41.0, and py39 no longer supported
        arr = struc.renumber_res_ids(arr, start=1)
    return arr


def standardize_atom_array(
    arr: _AtomArrayOrStack, standardize: bool = True
) -> _AtomArrayOrStack:
    if not standardize:
        return arr
    std_order = struc.info.standardize_order(arr)
    arr = apply_mask(arr, std_order).copy()
    arr = renumber_res_ids(arr)
    return arr


def stack_filter_intersect(
    arr_list: list[AtomArray],
    remove_elements: bool = False,
    standardize_order: bool = False,
) -> tuple[AtomArrayStack, bool]:
    # Ensure only matching atoms and annotations persisted
    # This way we can construct an AtomArrayStack
    common_decoys: list[AtomArray] = []
    largest_decoy = max(arr.shape[0] for arr in arr_list)

    idx = 0
    while len(common_decoys) < len(arr_list):
        decoy = arr_list[idx].copy()
        if hasattr(decoy, "atom_id"):
            # All annotation categories must be equal. If there are differing
            # atom counts in the decoys, the filter_intersection will return
            # a seemingly invalid, smaller number of intersecting atoms if the
            # atom_id (an atom index) is not equal.
            decoy.set_annotation("atom_id", np.repeat(0, decoy.shape[0]))
        if len(arr_list) == 1:
            # Nothing to filter intersect on, only one decoy
            common_decoys.append(standardize_atom_array(decoy, standardize_order))
            break

        decoy2 = arr_list[idx + 1].copy()
        if hasattr(decoy2, "atom_id"):
            decoy2.set_annotation("atom_id", np.repeat(0, decoy2.shape[0]))

        if remove_elements:
            decoy.set_annotation("element", np.repeat("", decoy.shape[0]))
            decoy2.set_annotation("element", np.repeat("", decoy2.shape[0]))
        decoy_intersect = decoy[struc.filter_intersection(decoy, decoy2)]
        if len(common_decoys):
            decoy_prev = common_decoys[-1].copy()
            decoy_intersect = decoy_intersect[
                struc.filter_intersection(decoy_intersect, decoy_prev)
            ]
        common_decoys.append(standardize_atom_array(decoy_intersect, standardize_order))
        idx += 1
        # We need the final element in the list
        if (idx + 1) == len(arr_list):
            decoy2 = decoy2[struc.filter_intersection(decoy2, decoy_intersect)]
            common_decoys.append(standardize_atom_array(decoy2, standardize_order))

    in_common = min(common.shape[0] for common in common_decoys)
    max_in_common = max(common.shape[0] for common in common_decoys)
    smallest_decoy = sorted(common_decoys, key=lambda x: x.shape[0])[0]  # type: ignore
    mismatch = largest_decoy - in_common
    if in_common < max_in_common:
        # Make sure all of the common atoms that exceed smallest are filtered out
        for i, decoy in enumerate(common_decoys):
            if decoy.shape[0] > in_common:
                decoy_intersect = decoy[
                    struc.filter_intersection(decoy, smallest_decoy)
                ]
                common_decoys[i] = decoy_intersect.copy()

    in_common = min(common.shape[0] for common in common_decoys)
    mismatch = largest_decoy - in_common
    if mismatch > 0:
        log.warning(f"Large atom mismatch detected between decoys: {mismatch} atoms.")
        lost_chains = any(len(set(arr.chain_id)) < 2 for arr in common_decoys)
        if lost_chains and not remove_elements:
            log.warning("Attempting to remove element annotations")
            return stack_filter_intersect(
                arr_list,
                standardize_order=False,
                remove_elements=True,
            )

    try:
        return stack(common_decoys), standardize_order
    except Exception:
        if standardize_order and remove_elements:
            # log.warning(
            #     "Caution: results will only represent the atoms in common! "
            #     f"keeping {in_common} / {largest_decoy} atoms in common"
            # )
            # Escape hatch to prevent infinite recursion
            raise ValueError("Unable to create valid decoy stack! Mismatched!")
        # First see if stack can be made by removing element annotations
        if not remove_elements:
            log.warning("Attempting to remove element annotations")
            return stack_filter_intersect(
                arr_list,
                standardize_order=False,
                remove_elements=True,
            )
        # Finally, see if renumbering + standardizing can fix
        if not standardize_order:
            log.warning("Attempting to renumber and standardize order!")
            return stack_filter_intersect(
                arr_list,
                standardize_order=True,
                remove_elements=True,
            )
        # Escape hatch to prevent infinite recursion
        raise ValueError("Unable to create valid decoy stack! Mismatched!")


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


def calc_num_mismatches(alignments: list[Alignment]) -> tuple[int, int]:
    """Calculate number of mismatches in sequence alignment

    Parameters
    ----------
    alignments : list[Alignment]
        sequence alignment(s)

    Returns
    -------
        tuple[int, int]:
            A tuple containing number of mismatches and matches in sequence alignment, respectively.
    """
    n_matches = 0
    n_mismatches = 0
    for alignment in alignments:
        s = align.alignment.get_symbols(alignment)
        # p = subject protein sequence, u = ref sequence
        ui = -1
        pj = -1
        for p, u in zip(s[0], s[1]):
            if u:
                ui += 1
            if p:
                pj += 1
            if u and p and u != "X" and p != "X":
                n_matches += 1
            else:
                n_mismatches += 1
    return n_mismatches, n_matches


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


def _get_seq_aligned_structures(
    ref: Path | _AtomArrayOrStack,
    subject: Path | _AtomArrayOrStack,
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


def get_buried_sasa(
    a: struc.AtomArray,
    b: struc.AtomArray,
    point_number: int = 200,
    ignore_ions: bool = True,
) -> int:
    """Calculate the rounded dSASA value for a pair of interacting substructures (chains).

    Parameters
    ----------
    a : struc.AtomArray
        The first substructure (chain).
    b : struc.AtomArray
        The second substructure (chain).
    point_number : int, optional
        The number of points per atom, by default 200.
        This parameter affects the speed and the accuracy of the calculation.

    Returns
    -------
    int
        The rounded dSASA value.

    Examples
    --------
    >>> atom1a = struc.atoms.Atom([0,0,0], chain_id="A")
    >>> atom2a = struc.atoms.Atom([1,1,1], chain_id="A")
    >>> atom1b = struc.atoms.Atom([2,2,2], chain_id="B")
    >>> atom2b = struc.atoms.Atom([3,3,3], chain_id="B")
    >>> a = struc.atoms.array([atom1a, atom2a])
    >>> b = struc.atoms.array([atom1b, atom2b])
    >>> get_buried_sasa(a, b, ignore_ions=False)
    0
    """
    dsasa: int = 0
    try:
        sasa_ab = np.sum(
            np.nan_to_num(
                struc.sasa(a + b, point_number=point_number, ignore_ions=ignore_ions)
            )
        )
        sasa_a = np.sum(
            np.nan_to_num(
                struc.sasa(a, point_number=point_number, ignore_ions=ignore_ions)
            )
        )
        sasa_b = np.sum(
            np.nan_to_num(
                struc.sasa(b, point_number=point_number, ignore_ions=ignore_ions)
            )
        )
        dsasa = sasa_a + sasa_b - sasa_ab
        dsasa = int(np.round(dsasa))
    except (KeyError, ValueError) as e:
        log.error(f"Failed to calculate buried sasa using ProtOr vdw radius: {e}")
        sasa_ab = np.sum(
            np.nan_to_num(
                struc.sasa(
                    a + b,
                    point_number=point_number,
                    ignore_ions=ignore_ions,
                    vdw_radii="Single",
                )
            )
        )
        sasa_a = np.sum(
            np.nan_to_num(
                struc.sasa(
                    a,
                    point_number=point_number,
                    ignore_ions=ignore_ions,
                    vdw_radii="Single",
                )
            )
        )
        sasa_b = np.sum(
            np.nan_to_num(
                struc.sasa(
                    b,
                    point_number=point_number,
                    ignore_ions=ignore_ions,
                    vdw_radii="Single",
                )
            )
        )
        dsasa = sasa_a + sasa_b - sasa_ab
        dsasa = int(np.round(dsasa))
    except Exception as e:
        log.error(
            f"Failed to calculate buried sasa even with vdw_radii=Single, returning 0! {e}"
        )
    return dsasa


def get_resolved_resi_from_atom_array(model: _AtomArrayOrStack) -> dict[str, list[int]]:
    """Returns a list of resolved resi (only checking the CA atom) in the PDB file.

    Parameters
    ----------
    model : _AtomArrayOrStack
        AtomArray or Stack

    Returns
    -------
    List[int]
        A list of resolved resi in the PDB file.
    """
    model_ch: set[str] = set(model.chain_id)
    per_chain_resolved: dict[str, list[int]] = {}
    for ch in model_ch:
        per_chain_resolved[ch] = [
            int(res_id)
            for res_id in model.res_id[
                (model.atom_name == "CA") & (model.chain_id == ch)
            ].tolist()
        ]
    return per_chain_resolved


def get_resolved_resi(pdb_path: str) -> dict[str, list[int]]:
    """Returns a list of resolved resi (only checking the CA atom) in the PDB file.

    Parameters
    ----------
    pdb_path : str
        The path to the PDB file.

    Returns
    -------
    List[int]
        A list of resolved resi in the PDB file.
    """
    model: _AtomArrayOrStack = atom_array_from_pdb_file(Path(pdb_path))
    resolved_resi = get_resolved_resi_from_atom_array(model)
    return resolved_resi


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


def cif_to_pdb(
    cif_path: str, pdb_path: str, chains: dict[str, str] | None = None
) -> None:
    """Convert CIF file to PDB file

    Parameters
    ----------
    cif_path : str
        The path to the CIF file to be converted.
    pdb_path : str
        The path where the converted PDB file will be saved.
    chains : Optional[dict], default is None
        If provided, only the chains specified in the dictionary will be kept
        and they will be renamed according to the dictionary.
        Note: The new chain names must be single characters to be PDB compatible.
    """
    # temp workaround for biotite API
    from biotite.structure.io import pdbx

    f = biotite_pdbxfile()
    pdbx_file = f.read(cif_path)
    model = pdbx.get_structure(
        pdbx_file, model=1, use_author_fields=False, extra_fields=["b_factor", "charge"]
    )
    # print(f"DEBUG: {cif_path} has {sorted(np.unique(model.chain_id))} chains")
    if chains is not None:
        mask = np.isin(model.chain_id, list(chains.keys()))
        if np.sum(mask) == 0:
            print(f"WARNING: chains {list(chains.keys())} not found in {cif_path}")
            return
        model = model[mask]
        model_orig = model.copy()
        for old_chain_id, new_chain_id in chains.items():
            mask = model_orig.chain_id == old_chain_id
            model.chain_id[mask] = new_chain_id
    write_pdb(model, Path(pdb_path))
