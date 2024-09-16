# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from itertools import combinations
from typing import Any, List, NamedTuple, Set, Tuple, Union

import biotite.structure as struc
import numpy as np
import pandas as pd
from biotite.structure.atoms import AtomArray
from numpy.typing import NDArray

from plinder.core.utils.log import setup_logger

_Contacts = Set[Tuple[str, str, int, int]]
_StackContacts = List[_Contacts]
_AtomResContacts = Union[Tuple[_Contacts, _Contacts], _Contacts]
_StackAtomResContacts = Union[Tuple[_StackContacts, _StackContacts], _StackContacts]


class ContactPairs(NamedTuple):
    residue_contacts: _Contacts
    atom_contacts: _Contacts


log = setup_logger(__name__)


def get_atoms_within_coordinate(
    atom_array: AtomArray | NDArray[np.double],
    coords: NDArray[np.double],
    radius: float,
    cell_size: float | None = None,
    as_mask: bool = False,
) -> NDArray[np.int_]:
    """
    Find atom indices (contact matrix) within distance of coords.

    Parameters
    ----------
    atom_array : AtomArray | NDArray[np.double]
        The array of atoms representing the protein structure.
    coords : NDArray[np.double]
        The coordinates to find atoms within.
    radius : float
        The radius within which to consider atoms as being in contact.
    cell_size : float
        The size of the cell for the cell list.
    as_mask : bool, optional
        Whether to return the contacts as a mask, by default False.

    Returns
    -------
    NDArray[np.int_]
        The contact indices.
    """
    cell_list = struc.CellList(atom_array, cell_size=cell_size)
    contacts: NDArray[np.int_] = cell_list.get_atoms(
        coords.reshape(-1, coords.shape[-1]), radius=radius, as_mask=False
    )
    return contacts


def get_atom_neighbors(
    atom_array: AtomArray,
    query_array: AtomArray,
    radius: float,
    cell_size: float | None = None,
    as_mask: bool = False,
) -> AtomArray:
    """
    Find atoms within a specified distance of coordinates.

    Parameters
    ----------
    atom_array : AtomArray
        The array of atoms.
    query_array : AtomArray
        The array of query atoms.
    radius : float
        The radius within which to find neighboring atoms.
    cell_size : float or None, optional
        The cell size for the cell list, by default None. If None, the radius is used as the cell size.
    as_mask : bool, optional
        Whether to return the contacts as a mask, by default False.

    Returns
    -------
    AtomArray
        The array of neighboring atoms.
    """
    if not cell_size:
        cell_size = radius
    contacts = get_atoms_within_coordinate(
        atom_array=query_array,
        coords=atom_array.coord,
        radius=radius,
        cell_size=cell_size,
        as_mask=as_mask,
    )
    neighbor_atoms = set()
    for at_idx, j in np.argwhere(contacts >= 0.0):
        neighbor_atoms.add(at_idx)

    return atom_array[list(neighbor_atoms)]


def get_binder_residues_atom_array(
    structure: AtomArray, binder_chains_resids: list[tuple[str, int]]
) -> AtomArray:
    struc_info = [(atm.chain_id, atm.res_id) for atm in structure]
    return structure[
        np.array([(sinfo in binder_chains_resids) for sinfo in struc_info])
    ]


def get_contact_df(
    prot_arr: AtomArray,
    ch_combos: List[Tuple[Tuple[Any, Any], Any]],
    lig_arr: AtomArray | None = None,
    radius: float = 4.5,
    cell_size: float | None = None,
) -> pd.DataFrame:
    all_contacts = []

    # Chain1 is set to ligand
    for ch1, ch2 in ch_combos:
        if lig_arr is not None:
            arr1 = lig_arr[(lig_arr.chain_id == ch1[0]) & (lig_arr.res_id == ch1[1])]
        else:
            arr1 = prot_arr[prot_arr.chain_id == ch1[0]]
        arr2 = prot_arr[prot_arr.chain_id == ch2[0]]
        contacts = get_atoms_within_coordinate(
            atom_array=arr1,
            coords=arr2.coord,
            radius=radius,
            cell_size=cell_size,
            as_mask=False,
        )
        contact_pairs = set()
        for a2_idx, a1_idx in np.argwhere(contacts != -1):
            # Chain1 residue
            a1_at = arr1.res_id[contacts[a2_idx, a1_idx]]
            # Chain2 residue
            a2_at = arr2.res_id[a2_idx]
            # Residue IDs in contact
            contact_pairs.add((a1_at, a2_at))
        if len(contact_pairs):
            # Note: Chain1 is always ligand
            all_contacts.append(
                {
                    "chain1": ch1[0],
                    "resid1": ch1[1],
                    "chain2": ch2[0],
                    "resid2": ch2[1],
                    "contacts": contact_pairs,
                }
            )
    if all_contacts:
        return pd.DataFrame(all_contacts)
    return pd.DataFrame()


def pairwise_chain_contacts(
    atom_array: AtomArray, radius: float = 4.5, cell_size: float | None = None
) -> pd.DataFrame:
    """
    Calculate pairwise contacts between chains in a protein structure.

    Parameters
    ----------
    atom_array : AtomArray
        The array of atoms representing the protein structure.
    radius : float, optional
        The radius within which to consider atoms as being in contact, by default 4.5.
    cell_size : float | None, optional
        The size of the cell used for the contact calculation. If None, it is set to the value of radius.

    Returns
    -------
    pd.DataFrame | None
        A DataFrame containing the chain pairs and their contacts, or None if no contacts are found.

    Examples
    --------
    >>> from plinder.core.structure.atoms import atom_array_from_pdb_file
    >>> atom_array = atom_array_from_pdb_file('1abc.pdb') # doctest: +SKIP
    >>> pairwise_chain_contacts(atom_array) # doctest: +SKIP
    """
    if not cell_size:
        cell_size = radius

    # Prot array
    prot_arr = atom_array[struc.filter_amino_acids(atom_array)]
    prot_chains_resids = sorted(set(prot_arr.chain_id))
    # Add dummy res_id for protein
    prot_chains_resids_set = {(i, None) for i in prot_chains_resids}

    lig_arr = atom_array[~struc.filter_amino_acids(atom_array) & (atom_array.hetero)]
    lig_chains_resids = sorted(set(zip(lig_arr.chain_id, lig_arr.res_id)))

    pli_combos = [(l, p) for l in lig_chains_resids for p in prot_chains_resids_set]
    ppi_combos = list(combinations(prot_chains_resids_set, 2))

    if len(lig_arr) == 0:
        pli_contacts = pd.DataFrame()
    else:
        pli_contacts = get_contact_df(
            prot_arr, pli_combos, lig_arr, radius=radius, cell_size=cell_size
        )
    ppi_contacts = get_contact_df(
        prot_arr, ppi_combos, radius=radius, cell_size=cell_size
    )

    if isinstance(ppi_contacts, pd.DataFrame) | isinstance(pli_contacts, pd.DataFrame):
        ppi_contacts = pd.DataFrame(ppi_contacts)
        return ppi_contacts, pli_contacts

    return pd.DataFrame(), pd.DataFrame()


def get_ligand_filter(
    complex_array: AtomArray, binder_chain_ids_resids: List[Tuple[str, int]]
) -> np.ndarray[bool, Any]:
    filter = np.array(
        [
            ((atm.chain_id, atm.res_id) in binder_chain_ids_resids)
            for atm in complex_array
        ]
    )
    return filter


def label_structure_gaps(atom_array: AtomArray) -> pd.DataFrame | None:
    """
    Find gaps in residue numbering between C-alpha atoms.

    Parameters
    ----------
    atom_array : AtomArray
        The array of atoms for which to find gaps.

    Returns
    -------
    chains : list
        A sorted list of unique chain IDs from the atom array.
    """
    chains = sorted(set(atom_array.chain_id))
    gaps = []
    for ch in chains:
        prot = atom_array[atom_array.chain_id == ch]
        prot_ca = prot[prot.atom_name == "CA"]
        prot_ca = prot_ca[prot_ca.element != "CA"]
        for i in range(prot_ca.shape[0] - 1):
            if prot_ca[i + 1].res_id != (prot_ca[i].res_id + 1):
                if prot_ca[i + 1].res_id == (prot_ca[i].res_id):
                    log.warning(
                        f"found two Calpha atoms for residue {prot_ca[i].res_id}"
                    )
                    continue
                gaps.append(
                    {
                        "chain": ch,
                        "gap_start": prot_ca[i].res_id,
                        "gap_end": prot_ca[i + 1].res_id,
                    }
                )
    if gaps:
        gaps = pd.DataFrame(gaps)
        return gaps
    return None
