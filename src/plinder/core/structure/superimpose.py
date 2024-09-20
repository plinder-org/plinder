from __future__ import annotations

from typing import Union

import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.structure as struc
import numpy as np
from numpy.typing import NDArray

from plinder.core.structure.vendored import resn2seq

__all__ = ["superimpose_chain"]


Atoms = Union[struc.AtomArray, struc.AtomArrayStack]


def superimpose_chain(
    fixed: Atoms,
    mobile: Atoms,
    substitution_matrix: str | align.SubstitutionMatrix = "BLOSUM62",
    gap_penalty: tuple[int, int] = (-10, -1),
    min_anchors: int = 3,
    outlier_threshold: float = 1.5,
    max_iterations: int = 10,
) -> tuple[
    Atoms,
    struc.AffineTransformation,
    NDArray[np.int_],
    NDArray[np.int_],
]:
    """Superimpose one protein chain onto another one, considering sequence differences and
    conformational outliers.

    The method finds corresponding residues by sequence alignment and selects their
    :math:`C_{\\alpha}` atoms as superimposition *anchors*.
    Then iteratively the anchor atoms are superimposed and outliers
    (in terms of high squared distance between matching anchor atoms) are removed, until
    either no outliers are left or the maximum number of iterations is reached.

    Parameters
    ----------
    fixed : AtomArray, shape(n,) or AtomArrayStack, shape(m,n) or ndarray, shape(n,), dtype=float or ndarray, shape(m,n), dtype=float
        The fixed structure(s).
        Alternatively coordinates can be given.
    mobile : AtomArray, shape(n,) or AtomArrayStack, shape(m,n) or ndarray, shape(n,), dtype=float or ndarray, shape(m,n), dtype=float
        The structure(s) which is/are superimposed on the `fixed` structure.
        Alternatively coordinates can be given.
    substitution_matrix : str or SubstitutionMatrix,
        The (name of the) substitution matrix used for sequence alignment.
    gap_penalty : int or tuple of int, optional
        The gap penalty for sequence alignment.
        A single value indicates a linear penalty, while a tuple indicates an affine penalty.
    min_anchors: int, optional
        If less than `min_anchors` anchors are found by sequence alignment, the method ditches the
        alignment and matches all :math:`C_{\\alpha}` atoms.
        If the number of :math:`C_{\\alpha}` atoms does not match in this fallback case, an exception
        is raised.
        Furthermore, the outlier removal is stopped, if less than `min_anchors` anchors would be
        left.
    outlier_threshold : float, optional
        The threshold for considering a conformational outlier.
        The threshold is given in units of IQR.
        This means that a residue is considered an outlier, if its squared distance is
        `outlier_threshold` x IQR above the upper quartile.
    max_iterations : int, optional
        The maximum number of iterations for removing conformational outliers.

    Returns
    -------
    fitted : AtomArray or AtomArrayStack or ndarray, shape(n,), dtype=float or ndarray, shape(m,n), dtype=float
        A copy of the `mobile` structure(s), superimposed on the fixed structure.
        Only coordinates are returned, if coordinates were given in `mobile`.
    transform : AffineTransformation
        This object contains the affine transformation(s) that were applied on `mobile`.
        :meth:`AffineTransformation.apply()` can be used to transform another AtomArray
        in the same way.
    fixed_anchor_indices, mobile_anchor_indices : ndarray, shape(k,), dtype=int
        The indices of the anchor :math:`C_{\\alpha}` atoms in the fixed and mobile
        structure, respectively.
        These atoms were used for the superimposition.

    Note:
        As this method relies on amino acid sequence alignment, it works only for proteins
        with decent homology.

    """
    fixed_ca_indices = _get_ca_indices(fixed)
    mobile_ca_indices = _get_ca_indices(mobile)
    if len(fixed_ca_indices) < min_anchors or len(mobile_ca_indices) < min_anchors:
        raise ValueError(
            "Structures have too few CA atoms for required number of anchors"
        )

    anchor_indices = _find_matching_anchors(
        fixed[fixed_ca_indices],
        mobile[mobile_ca_indices],
        substitution_matrix,
        gap_penalty,
    )
    if len(anchor_indices) < min_anchors:
        # Fallback: Match all CA atoms
        if len(fixed_ca_indices) != len(mobile_ca_indices):
            raise ValueError("Tried fallback, but number of CA atoms does not match")
        fixed_anchor_indices = fixed_ca_indices
        mobile_anchor_indices = mobile_ca_indices
    else:
        # The anchor indices point to the CA atoms
        # -> get the corresponding indices for the whole structure
        fixed_anchor_indices = fixed_ca_indices[anchor_indices[:, 0]]
        mobile_anchor_indices = mobile_ca_indices[anchor_indices[:, 1]]

    if max_iterations < 1:
        raise ValueError("Maximum number of iterations must be at least 1")
    # For the first iteration, all anchors are included
    not_outlier_mask = np.ones(len(fixed_anchor_indices), dtype=bool)
    for _ in range(max_iterations):
        # Remove outliers from previous refinement iteration
        fixed_anchor_indices = fixed_anchor_indices[not_outlier_mask]
        mobile_anchor_indices = mobile_anchor_indices[not_outlier_mask]
        # Run superimposition
        fixed_coord = fixed.coord[..., fixed_anchor_indices, :]
        mobile_coord = mobile.coord[..., mobile_anchor_indices, :]
        superimposed_coord, transform = struc.superimpose(fixed_coord, mobile_coord)
        # Find outliers
        sq_dist = struc.distance(fixed_coord, superimposed_coord) ** 2
        lower_quartile, upper_quartile = np.quantile(sq_dist, [0.25, 0.75])
        interquartile_range = upper_quartile - lower_quartile
        not_outlier_mask = (
            sq_dist <= upper_quartile + outlier_threshold * interquartile_range
        )
        if np.all(not_outlier_mask):
            # No outliers anymore -> early termination
            break
        if np.count_nonzero(not_outlier_mask) < min_anchors:
            # Less than min_anchors anchors would be left -> early termination
            break

    return (
        transform.apply(mobile),
        transform,
        fixed_anchor_indices,
        mobile_anchor_indices,
    )


def _get_ca_indices(atoms: Atoms) -> NDArray[np.int_]:
    return np.where((struc.filter_amino_acids(atoms)) & (atoms.atom_name == "CA"))[0]


def _find_matching_anchors(
    fixed_ca: Atoms,
    mobile_ca: Atoms,
    substitution_matrix: str | align.SubstitutionMatrix,
    gap_penalty: tuple[int, int] | int,
) -> NDArray[np.int_]:
    fixed_seq = _get_sequence(fixed_ca)
    mobile_seq = _get_sequence(mobile_ca)

    alphabet = seq.ProteinSequence.alphabet
    if isinstance(substitution_matrix, str):
        substitution_matrix = align.SubstitutionMatrix(
            alphabet, alphabet, substitution_matrix
        )
    score_matrix = substitution_matrix.score_matrix()

    alignment = align.align_optimal(
        fixed_seq,
        mobile_seq,
        substitution_matrix,
        gap_penalty,
        terminal_penalty=False,
        max_number=1,
    )[0]
    alignment_codes = align.get_codes(alignment)
    anchor_mask = (
        # Anchors must be similar amino acids
        # X and X are not similar (at least in built-in matrices)
        (score_matrix[alignment_codes[0], alignment_codes[1]] > 0)
        # Cannot anchor gaps
        & (alignment_codes[0] != -1)
        & (alignment_codes[1] != -1)
    )
    anchors: NDArray[np.int_] = alignment.trace[anchor_mask]
    return anchors


def _get_sequence(ca_atoms: Atoms) -> seq.ProteinSequence:
    # We already know that the atoms are amino acids
    # -> no further check necessary
    return seq.ProteinSequence(resn2seq(ca_atoms.res_name))
