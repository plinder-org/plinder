from __future__ import annotations

import biotite.structure as struc
import numpy as np
from biotite.structure.atoms import AtomArray, AtomArrayStack
from numpy.typing import NDArray

from plinder.core.structure.models import ChainConfig
from plinder.core.structure.vendored import apply_mask, get_seq_aligned_structures
from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)


def remove_annotations(
    structure: AtomArray | AtomArrayStack,
    categories: list[str] = ["element", "ins_code"],
) -> AtomArray | AtomArrayStack:
    if isinstance(structure, AtomArrayStack):
        shape_idx = 1
    else:
        shape_idx = 0

    for annotation in categories:
        val: float | str = 0.0 if annotation == "b_factor" else ""
        annotation_arr: NDArray[np.double | np.str_] = np.repeat(
            val, structure.shape[shape_idx]
        )
        structure.set_annotation(annotation, annotation_arr)
    return structure


def fix_annotation_mismatch(
    ref: AtomArray,
    decoys: AtomArrayStack,
    categories: list[str] = ["element", "ins_code"],
) -> tuple[AtomArray, AtomArrayStack]:
    for annot in ref.get_annotation_categories():
        ref_annot = ref.get_annotation(annot)
        decoy_annot = decoys.get_annotation(annot)
        if not np.array_equal(ref_annot, decoy_annot):
            log.debug(f"Decoy and ref have differing {annot} categories!")
            if annot not in categories:
                continue

            decoys = remove_annotations(decoys, categories=[annot])
            ref = remove_annotations(ref, categories=[annot])
    return ref, decoys


def fix_mismatched_atoms(
    native: AtomArray, decoy_stack: AtomArrayStack, max_atom_delta: int
) -> tuple[AtomArray, AtomArrayStack]:
    identical = np.array_equal(decoy_stack.res_id, native.res_id)
    if identical:
        # Both shape and res_id elements are identical
        return native, decoy_stack

    log.debug("Detected mismatch between native and decoy stack!")
    # Sequence based structural alignment
    # Make the decoy_stack match numbering of native
    log.debug("Attempting sequence-based structural alignment")
    native, decoy_stack = get_seq_aligned_structures(native, decoy_stack)
    native_intersect = native[struc.filter_intersection(native, decoy_stack)]
    in_common = native_intersect.shape[0]
    native_atoms = native.shape[0]
    mismatch = native_atoms - in_common
    if mismatch > max_atom_delta:
        log.debug(
            f"Large atom mismatch detected between native and models: {mismatch} atoms."
        )
        log.debug("Attempting to fix mismatch")

        # In test case there are missing element annotations
        native = remove_annotations(native)
        decoy_stack = remove_annotations(decoy_stack)
        native_intersect = native[struc.filter_intersection(native, decoy_stack)]
        in_common = native_intersect.shape[0]
        log.warning(
            "Caution: results will only represent the atoms in common! "
            f"keeping {in_common} / {native_atoms} atoms in common"
        )

    native = native_intersect.copy()
    decoy_stack = decoy_stack[..., struc.filter_intersection(decoy_stack, native)]
    return native, decoy_stack


def set_canonical_chain_order(
    structure: AtomArray | AtomArrayStack | list[AtomArray],
    chains: "ChainConfig",
    subject: str,
) -> AtomArray | AtomArrayStack | list[AtomArray]:
    # Create set of residues in interface split into receptor and ligand
    # Conflict between residue numbers in different chains is handled by
    # logical mask on array.chain_id and array.res_id
    raise NotImplementedError
    lig_chains = getattr(chains, f"{subject}_ligand")
    rec_chains = getattr(chains, f"{subject}_receptor")

    if isinstance(structure, list):
        for i, arr in enumerate(structure):
            R_mask = np.isin(arr.chain_id, rec_chains)
            L_mask = np.isin(arr.chain_id, lig_chains)
            R = arr[R_mask].copy()
            L = arr[L_mask].copy()
            structure[i] = R + L
        return structure
    else:
        R_mask = np.isin(structure.chain_id, rec_chains)
        L_mask = np.isin(structure.chain_id, lig_chains)
        R = apply_mask(structure, R_mask)
        L = apply_mask(structure, L_mask)
        return R + L


def remove_duplicate_calpha(atoms: AtomArray) -> AtomArray:
    unique_mask = []
    unique = set()
    for at in atoms:
        at_id = f"{at.chain_id}-{at.res_id}"
        mask = not (at_id in unique)
        if mask:
            unique.add(at_id)
        else:
            log.warning(f"{at_id} is duplicated!")
        unique_mask.append(mask)
    atoms = atoms[np.array(unique_mask)].copy()
    return atoms
