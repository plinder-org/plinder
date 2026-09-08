# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Protein-ligand interaction detection.

Computes interaction fingerprints (hydrogen bonds, salt bridges, pi-stacking,
pi-cation, halogen bonds) via peppr's ``ContactMeasurement``, plus water- and
metal-bridged interactions, covalent connections from ``_struct_conn``, and
crystallographic symmetry-mate contacts.
"""

from __future__ import annotations

import multiprocessing as mp
from collections import OrderedDict, defaultdict
from multiprocessing.connection import Connection
from pathlib import Path
from threading import RLock
from typing import Any, cast

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
import peppr.contacts as peppr_contacts
from biotite.structure import filter_heavy
from numpy.typing import NDArray
from rdkit import Chem

from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)

# PEPPR caches interchangeable tautomers by canonical SMILES inside one
# ContactMeasurement. Plinder creates one measurement per ligand occurrence,
# which otherwise repeats the same expensive protein-residue enumeration many
# times in one ingest worker. Extend the identical cache key across measurements
# while bounding memory across a multi-entry batch.
_MAX_TAUTOMER_CACHE_SIZE = 4096
_TautomerCacheKey = tuple[
    str,
    tuple[tuple[int, int, int, int, bool, int | None], ...],
    tuple[tuple[int, int, float, bool, bool, int], ...],
]
_TAUTOMER_CACHE: OrderedDict[_TautomerCacheKey, tuple[Chem.Mol, ...]] = OrderedDict()
_TAUTOMER_CACHE_LOCK = RLock()
_ORIGINAL_GET_INTERCHANGEABLE_TAUTOMERS = peppr_contacts.get_interchangeable_tautomers
_ORIGINAL_FIND_RESONANCE_CHARGES = peppr_contacts.find_resonance_charges

# Normal PEPPR behavior remains unchanged.  Only chemically large, charged
# residues run resonance enumeration in an interruptible child process.  This
# avoids an unkillable RDKit combinatorial tail while preserving the exact
# result whenever enumeration completes within the generous bound.
_RESONANCE_GUARD_MIN_ATOMS = 64
_RESONANCE_TIMEOUT_SECONDS = 10.0
_MAX_RESONANCE_CACHE_SIZE = 4096
_ResonanceResult = tuple[
    NDArray[np.bool_],
    NDArray[np.bool_],
    NDArray[np.int_],
]
_RESONANCE_CACHE: OrderedDict[_TautomerCacheKey, _ResonanceResult] = OrderedDict()
_RESONANCE_CACHE_LOCK = RLock()


def _tautomer_cache_key(molecule: Chem.Mol) -> _TautomerCacheKey:
    """Include atom order and PEPPR split-residue state in the cache key."""
    heavy_neighbors_property = peppr_contacts._ORIG_NUM_HEAVY_NEIGHS
    atom_signature = tuple(
        (
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            int(atom.GetHybridization()),
            atom.GetTotalNumHs(),
            atom.GetIsAromatic(),
            atom.GetIntProp(heavy_neighbors_property)
            if atom.HasProp(heavy_neighbors_property)
            else None,
        )
        for atom in molecule.GetAtoms()
    )
    bond_signature = tuple(
        sorted(
            (
                min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                bond.GetBondTypeAsDouble(),
                bond.GetIsAromatic(),
                bond.GetIsConjugated(),
                int(bond.GetStereo()),
            )
            for bond in molecule.GetBonds()
        )
    )
    return Chem.MolToSmiles(molecule), atom_signature, bond_signature


def _cached_get_interchangeable_tautomers(molecule: Chem.Mol) -> list[Chem.Mol]:
    """Reuse PEPPR tautomer enumeration across ligand contact measurements."""
    key = _tautomer_cache_key(molecule)
    with _TAUTOMER_CACHE_LOCK:
        cached = _TAUTOMER_CACHE.get(key)
        if cached is not None:
            _TAUTOMER_CACHE.move_to_end(key)
            return [Chem.Mol(tautomer) for tautomer in cached]

    generated = tuple(
        Chem.Mol(tautomer)
        for tautomer in _ORIGINAL_GET_INTERCHANGEABLE_TAUTOMERS(molecule)
    )
    with _TAUTOMER_CACHE_LOCK:
        _TAUTOMER_CACHE[key] = generated
        _TAUTOMER_CACHE.move_to_end(key)
        while len(_TAUTOMER_CACHE) > _MAX_TAUTOMER_CACHE_SIZE:
            _TAUTOMER_CACHE.popitem(last=False)
    return [Chem.Mol(tautomer) for tautomer in generated]


def _copy_resonance_result(result: _ResonanceResult) -> _ResonanceResult:
    return tuple(array.copy() for array in result)  # type: ignore[return-value]


def _run_resonance_worker(molecule: Chem.Mol, connection: Connection) -> None:
    """Run PEPPR's exact RDKit resonance enumeration out of process."""
    try:
        connection.send(("ok", _ORIGINAL_FIND_RESONANCE_CHARGES(molecule)))
    except BaseException as exc:
        connection.send(("error", f"{type(exc).__name__}: {exc}"))
    finally:
        connection.close()


def _deposited_charge_fallback(molecule: Chem.Mol) -> _ResonanceResult:
    """Use the estimated deposited formal-charge assignment without resonance."""
    charges = np.fromiter(
        (atom.GetFormalCharge() for atom in molecule.GetAtoms()),
        dtype=int,
        count=molecule.GetNumAtoms(),
    )
    return (
        charges > 0,
        charges < 0,
        np.arange(molecule.GetNumAtoms(), dtype=int),
    )


def _bounded_find_resonance_charges(molecule: Chem.Mol) -> _ResonanceResult:
    """Preserve exact PEPPR resonance behavior with a rare timeout fallback."""
    is_large_and_charged = molecule.GetNumAtoms() >= _RESONANCE_GUARD_MIN_ATOMS and any(
        atom.GetFormalCharge() != 0 for atom in molecule.GetAtoms()
    )
    if not is_large_and_charged:
        return cast(_ResonanceResult, _ORIGINAL_FIND_RESONANCE_CHARGES(molecule))

    key = _tautomer_cache_key(molecule)
    with _RESONANCE_CACHE_LOCK:
        cached = _RESONANCE_CACHE.get(key)
        if cached is not None:
            _RESONANCE_CACHE.move_to_end(key)
            return _copy_resonance_result(cached)

    # Data ingest is Linux-only and fork lets the child reuse the already
    # constructed RDKit molecule without serializing the full contact model.
    context = mp.get_context("fork")
    receive, send = context.Pipe(duplex=False)
    process = context.Process(
        target=_run_resonance_worker,
        args=(molecule, send),
        daemon=True,
    )
    process.start()
    send.close()
    status = "timeout"
    payload: Any = None
    try:
        if receive.poll(_RESONANCE_TIMEOUT_SECONDS):
            status, payload = receive.recv()
    finally:
        receive.close()
        if process.is_alive():
            process.terminate()
        process.join(timeout=1.0)
        if process.is_alive():
            process.kill()
            process.join()

    if status == "ok":
        result = cast(_ResonanceResult, payload)
    elif status == "error":
        raise struc.BadStructureError(
            f"PEPPR resonance enumeration failed in worker: {payload}"
        )
    else:
        result = _deposited_charge_fallback(molecule)
        log.warning(
            "PEPPR resonance enumeration exceeded %.1fs for a charged "
            "%d-atom residue; using its estimated deposited formal charges",
            _RESONANCE_TIMEOUT_SECONDS,
            molecule.GetNumAtoms(),
        )

    with _RESONANCE_CACHE_LOCK:
        _RESONANCE_CACHE[key] = _copy_resonance_result(result)
        _RESONANCE_CACHE.move_to_end(key)
        while len(_RESONANCE_CACHE) > _MAX_RESONANCE_CACHE_SIZE:
            _RESONANCE_CACHE.popitem(last=False)
    return _copy_resonance_result(result)


peppr_contacts.get_interchangeable_tautomers = cast(
    Any, _cached_get_interchangeable_tautomers
)
peppr_contacts.find_resonance_charges = cast(Any, _bounded_find_resonance_charges)
ContactMeasurement = peppr_contacts.ContactMeasurement


def get_symmetry_mate_contacts(
    mmcif: Path | pdbx.CIFFile, contact_threshold: float = 5.0
) -> dict[tuple[str, int], dict[tuple[str, int], dict[int, set[int]]]]:
    """
    Find inter-residue contacts generated by crystallographic symmetry.

    Only contacts involving symmetry mates (not identity) are returned.

    Parameters
    ----------
    mmcif : Path or CIFFile
        Parsed mmCIF or a path to one (supports .gz). Passing a parsed file
        avoids reading a large source entry again during ingest.
    contact_threshold : float, optional
        Distance cutoff in Angstrom, by default 5.0.

    Returns
    -------
    dict[tuple[str, int], dict[tuple[str, int], dict[int, set[int]]]]
        Mapping of (chain_id, residue_id) to partner residues,
        with atom serials mapped to the symmetry image indices.
    """
    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        get_unit_cell_with_altloc,
        read_mmcif_file,
    )

    cif_file = mmcif if isinstance(mmcif, pdbx.CIFFile) else read_mmcif_file(mmcif)

    # Build the full unit cell (all symmetry copies)
    try:
        unit_cell = get_unit_cell_with_altloc(
            cif_file, model=1, use_author_fields=False
        )
    except Exception:
        # No symmetry information (NMR, computational models)
        return {}
    unit_cell = unit_cell[~struc.filter_solvent(unit_cell)]
    unit_cell = unit_cell[filter_heavy(unit_cell)]

    if unit_cell.box is None:
        return {}
    box_lengths = np.linalg.norm(unit_cell.box, axis=1)
    if not np.all(np.isfinite(box_lengths)) or np.any(
        box_lengths <= 2 * contact_threshold
    ):
        # Non-crystallographic entries can carry a placeholder 1 Å box plus
        # assembly operators (e.g. EM entry 6hbg).  A periodic search radius
        # spanning half a box or more is physically ambiguous and makes the
        # CellList repeat a huge number of meaningless images.
        log.warning(
            "Skipping symmetry-mate contacts for invalid/too-small unit cell %s",
            box_lengths.tolist(),
        )
        return {}

    # Get ASU to determine atoms per symmetry copy
    asu = get_structure_with_altloc(cif_file, model=1, use_author_fields=False)
    asu = asu[~struc.filter_solvent(asu)]
    asu = asu[filter_heavy(asu)]
    n_asu = len(asu)
    n_total = len(unit_cell)
    if n_total == n_asu:
        return {}
    n_copies = n_total // n_asu

    # Label each atom with its symmetry image index
    image_idx = np.repeat(np.arange(n_copies), n_asu)

    # Use periodic CellList to find contacts across unit cell boundaries
    cell_list = struc.CellList(
        unit_cell,
        cell_size=contact_threshold,
        periodic=True,
        box=unit_cell.box,
    )

    results: dict[
        tuple[str, int], dict[tuple[str, int], dict[int, set[int]]]
    ] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

    # For each atom in the ASU (image 0), find contacts with symmetry mates
    for i in range(n_asu):
        neighbors = cell_list.get_atoms(unit_cell.coord[i], radius=contact_threshold)
        for j in neighbors:
            if j == i or image_idx[j] == 0:
                continue
            c1 = (
                unit_cell.label_asym_id[i]
                if hasattr(unit_cell, "label_asym_id")
                else unit_cell.chain_id[i]
            )
            c2 = (
                unit_cell.label_asym_id[j]
                if hasattr(unit_cell, "label_asym_id")
                else unit_cell.chain_id[j]
            )
            r1 = int(unit_cell.res_id[i]) if unit_cell.res_id[i] else 1
            r2 = int(unit_cell.res_id[j]) if unit_cell.res_id[j] else 1
            atom_serial = i + 1
            results[(c1, r1)][(c2, r2)][atom_serial].add(int(image_idx[j]))

    return results


def get_covalent_connections(
    cif_data: pdbx.CIFBlock,
) -> dict[str, list[tuple[str, str]]]:
    """
    Extract covalent connections from CIF block.

    Parameters
    ----------
    cif_data : pdbx.CIFBlock
        biotite CIF block

    Returns
    -------
    dict[str, list[tuple[str, str]]]
        All covalent links as defined by mmcif annotations
    """
    from plinder.data.annotations.cif_utils import parse_struct_conn

    nucleobase_list = {"A", "C", "U", "G", "DA", "DC", "DG", "DT", "PSU"}
    valid_types = {"covale", "metalc", "hydrog"}

    cov_dict: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for c in parse_struct_conn(cif_data):
        if c["conn_type"] not in valid_types:
            continue
        if c["conn_type"] == "hydrog":
            if c["comp1"].strip() not in nucleobase_list:
                continue
        link1 = ":".join(
            [
                c["auth_seq1"],
                c["comp1"],
                c["chain1"],
                c["seq1"],
                c["atom1"],
            ]
        )
        link2 = ":".join(
            [
                c["auth_seq2"],
                c["comp2"],
                c["chain2"],
                c["seq2"],
                c["atom2"],
            ]
        )
        cov_dict[c["conn_type"]].append((link1, link2))
    return cov_dict


def extract_ligand_links_to_neighbouring_chains(
    all_covalent_dict: dict[str, list[tuple[str, str]]],
    ligand_asym_id: str,
    neighboring_asym_ids: set[str],
    link_type: str = "covale",
) -> set[str]:
    """
    Parse covalant dictionary for a given ligand and its neighbours.

    Parameters
    ----------
    all_covalent_dict : dict[str, list[tuple[str, str]]]
        All covalent links as defined by mmcif annotations
    ligand_asym_id : str
        ligand asymmetric identification string
    neighboring_asym_ids : set[str]
        set of neighbour asymmetric identification strings
    link_type : str, optional
        covalent linkage type in dictionary, by default "covale",
        options include:
            "covale": actual covalent linkage
            "metalc": other dative bond, eg. metal-ligand dative bond
            "hydrog": strong hydrogen bonding of nucleic acid

    Returns
    -------
    set[str]
        set of covalent linkages in the entry between the ligand and its neighbours

    Notes
    -----
    For the purpose of covalent annotations, we only consider "covale".
    """
    covalent_linkages = set()
    if link_type in all_covalent_dict:
        for link1, link2 in all_covalent_dict[link_type]:
            chain1, chain2 = link1.split(":")[2], link2.split(":")[2]
            chains = {chain1, chain2}
            if len(chains) == 1:
                # remove linkages that are to the same chain!
                continue
            if ligand_asym_id not in chains:
                # only if the ligand is involved
                continue
            # One chain is the ligand, the other a neighbouring (receptor)
            # chain; emit as receptor__ligand. Match whole asym ids, which can be
            # multi-character (e.g. "AA") once an entry has more than 26 chains.
            if chain1 in neighboring_asym_ids and chain2 == ligand_asym_id:
                covalent_linkages.add(f"{link1}__{link2}")
            elif chain2 in neighboring_asym_ids and chain1 == ligand_asym_id:
                covalent_linkages.add(f"{link2}__{link1}")
    return covalent_linkages


# ---------------------------------------------------------------------------
# Bridged interaction detection
# TODO: remove once a new peppr is released with these methods.
# ---------------------------------------------------------------------------

# Water bridge lower bound: 0.75 * VdW_sum (~2.28 A for O-O)
# avoids clashes but allows short water-mediated H-bonds.
# Upper bound: 1.15 * VdW_sum (~3.50 A for O-O), standard H-bond max.
_WATER_BRIDGE_DISTANCE_SCALING = (0.75, 1.15)

# Metals that form coordination bonds (not spectator ions like Na/Cl/K)
_COORDINATION_METALS = frozenset(
    {
        "MG",
        "CA",
        "ZN",
        "FE",
        "FE2",  # Fe(II)
        "MN",
        "CO",
        "CU",
        "CU1",  # Cu(I)
        "NI",
        "CD",
        "MO",
        "4MO",  # Mo(IV)
        "6MO",  # Mo(VI)
        "W",
        "V",
    }
)
_METAL_ACCEPTOR_PATTERN = (
    "["
    "$([O]),"
    "$([#7;!$([nX3]);!$([NX3]-*=[!#6]);!$([NX3]-[a]);!$([NX4])]),"
    "$([#16]),"
    "$([*;-{1-};!+{1-}])"
    "]"
)


def _find_bridged_interactions(
    receptor: "struc.AtomArray",
    ligand: "struc.AtomArray",
    bridge_atoms: "struc.AtomArray",
    receptor_pattern: str,
    ligand_pattern: str,
    distance_scaling: tuple[float, float],
) -> list[tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]]:
    """Find interactions bridged by intermediary atoms (water or metal).

    TODO: remove once peppr has ContactMeasurement.find_bridged_interactions.
    """
    import biotite.structure.info as info
    from peppr.contacts import ContactMeasurement, find_atoms_by_pattern

    if bridge_atoms.array_length() == 0:
        return []

    try:
        cm = ContactMeasurement(receptor, ligand)
    except Exception as e:
        log.warning(f"ContactMeasurement setup failed: {e}")
        return []

    receptor_matched = find_atoms_by_pattern(cm._binding_site_mol, receptor_pattern)
    ligand_matched = find_atoms_by_pattern(cm._ligand_mol, ligand_pattern)
    if len(receptor_matched) == 0 or len(ligand_matched) == 0:
        return []

    receptor_coords = cm._binding_site.coord[receptor_matched]
    ligand_coords = cm._ligand.coord[ligand_matched]
    lo, hi = sorted(distance_scaling)

    r_vdw = np.array(
        [info.vdw_radius_single(e) for e in cm._binding_site.element[receptor_matched]]
    )
    l_vdw = np.array(
        [info.vdw_radius_single(e) for e in cm._ligand.element[ligand_matched]]
    )

    bridges: list[tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]] = []
    for bi in range(bridge_atoms.array_length()):
        b_coord = bridge_atoms.coord[bi]
        b_vdw = info.vdw_radius_single(bridge_atoms.element[bi])

        r_dists = np.linalg.norm(receptor_coords - b_coord, axis=1)
        r_thresholds = r_vdw + b_vdw
        r_contacts = receptor_matched[
            (r_dists >= lo * r_thresholds) & (r_dists <= hi * r_thresholds)
        ]
        if len(r_contacts) == 0:
            continue

        l_dists = np.linalg.norm(ligand_coords - b_coord, axis=1)
        l_thresholds = l_vdw + b_vdw
        l_contacts = ligand_matched[
            (l_dists >= lo * l_thresholds) & (l_dists <= hi * l_thresholds)
        ]
        if len(l_contacts) == 0:
            continue

        for ri in r_contacts:
            for li in l_contacts:
                bridges.append(
                    (
                        cm._binding_site_indices[ri : ri + 1],
                        np.array([li], dtype=int),
                        np.array([bi], dtype=int),
                    )
                )

    return bridges


def find_water_bridges(
    receptor: "struc.AtomArray",
    ligand: "struc.AtomArray",
    waters: "struc.AtomArray",
    distance_scaling: tuple[float, float] = _WATER_BRIDGE_DISTANCE_SCALING,
) -> list[tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]]:
    """Find water-mediated hydrogen bonds between receptor and ligand."""
    from peppr.common import ACCEPTOR_PATTERN, DONOR_PATTERN

    water_oxygens = waters[waters.element == "O"]
    hbond_pattern = "[" + DONOR_PATTERN[1:-1] + "," + ACCEPTOR_PATTERN[1:-1] + "]"
    return _find_bridged_interactions(
        receptor,
        ligand,
        water_oxygens,
        hbond_pattern,
        hbond_pattern,
        distance_scaling,
    )


def find_metal_bridges(
    receptor: "struc.AtomArray",
    ligand: "struc.AtomArray",
    metals: "struc.AtomArray",
    cutoff: float = 3.0,
) -> list[tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]]:
    """Find metal-mediated coordination between receptor and ligand."""
    from peppr.contacts import ContactMeasurement, find_atoms_by_pattern

    coord_mask = np.isin(metals.res_name, list(_COORDINATION_METALS))
    if not np.any(coord_mask):
        return []
    coord_metals = metals[coord_mask]

    try:
        cm = ContactMeasurement(receptor, ligand)
    except Exception as e:
        log.warning(f"ContactMeasurement setup failed for metal bridges: {e}")
        return []

    receptor_matched = find_atoms_by_pattern(
        cm._binding_site_mol, _METAL_ACCEPTOR_PATTERN
    )
    ligand_matched = find_atoms_by_pattern(cm._ligand_mol, _METAL_ACCEPTOR_PATTERN)
    if len(receptor_matched) == 0 or len(ligand_matched) == 0:
        return []

    bridges: list[tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]] = []
    for bi in range(coord_metals.array_length()):
        b_coord = coord_metals.coord[bi]
        r_dists = np.linalg.norm(
            cm._binding_site.coord[receptor_matched] - b_coord, axis=1
        )
        r_contacts = receptor_matched[r_dists < cutoff]
        if len(r_contacts) == 0:
            continue
        l_dists = np.linalg.norm(cm._ligand.coord[ligand_matched] - b_coord, axis=1)
        l_contacts = ligand_matched[l_dists < cutoff]
        if len(l_contacts) == 0:
            continue
        for ri in r_contacts:
            for li in l_contacts:
                bridges.append(
                    (
                        cm._binding_site_indices[ri : ri + 1],
                        np.array([li], dtype=int),
                        np.array([bi], dtype=int),
                    )
                )
    return bridges


# Canonical interaction-type names recorded when a peppr detector fails, so a
# missing interaction type can be told apart from a genuinely absent one.
PEPPR_INTERACTION_TYPES = (
    "hydrogen_bonds",
    "salt_bridges",
    "pi_stacks",
    "pi_cation",
    "halogen_bonds",
    "water_bridges",
    "metal_complexes",
)


def run_peppr_interactions(
    receptor: struc.AtomArray,
    ligand: struc.AtomArray,
    waters: struc.AtomArray,
    metals: struc.AtomArray,
    ligand_chain: str,
    chain_mapping: dict[str, str],
) -> tuple[dict[str, dict[int, list[str]]], set[tuple[str, int]]]:
    """Compute interaction hash using peppr ContactMeasurement.

    Parameters
    ----------
    receptor : AtomArray
        Receptor heavy atoms.
    ligand : AtomArray
        Ligand heavy atoms.
    waters : AtomArray
        Water heavy atoms.
    metals : AtomArray
        Metal ion heavy atoms (used for metal bridge detection).
    ligand_chain : str
        Ligand chain identifier ({instance}.{chain}).
    chain_mapping : dict[str, str]
        Identity mapping over the chain IDs already in
        ``{instance}.{chain}`` form (kept as a parameter for legacy
        reasons; callers pass ``{c: c for c in np.unique(...)}``).

    Returns
    -------
    interaction_hashes : dict
        {instance.chain: {residue_number: [interaction_strings]}}
    water_set : set
        {(instance.chain, residue_number)} of bridging waters.
    failed_interaction_types : list[str]
        Interaction types (from :data:`PEPPR_INTERACTION_TYPES`) whose peppr
        detector raised, so an empty result for them means "not computed",
        not "none found". A failure to build the complex records all types.
    """
    interaction_hashes: dict[str, dict[int, list[str]]] = {}
    water_set: set[tuple[str, int]] = set()
    failed_interaction_types: list[str] = []

    try:
        cm = ContactMeasurement(receptor, ligand)
    except Exception as e:
        components = sorted({str(name) for name in ligand.res_name})
        log.warning(
            "run_peppr_interactions: peppr could not build the ligand %s + "
            "binding-site complex for contact analysis; skipping its "
            "interactions. Components: %s. Underlying error: %s",
            ligand_chain,
            components,
            e,
        )
        # No detector could run — every type is unknown, not absent.
        return interaction_hashes, water_set, list(PEPPR_INTERACTION_TYPES)

    def _add(chain: str, resnr: int, attr: str) -> None:
        if chain == ligand_chain:
            return
        if chain not in interaction_hashes:
            interaction_hashes[chain] = {}
        if resnr not in interaction_hashes[chain]:
            interaction_hashes[chain][resnr] = []
        interaction_hashes[chain][resnr].append(attr)

    _PROTEIN_MAINCHAIN = {"N", "CA", "C", "O"}
    _NA_MAINCHAIN = {"P", "O5'", "C5'", "C4'", "C3'", "O3'"}
    _mainchain_mask = (
        np.isin(receptor.atom_name, list(_PROTEIN_MAINCHAIN))
        & struc.filter_amino_acids(receptor)
    ) | (
        np.isin(receptor.atom_name, list(_NA_MAINCHAIN))
        & struc.filter_nucleotides(receptor)
    )

    def _is_sidechain(atom_idx: int) -> bool:
        return not _mainchain_mask[atom_idx]

    # H-bonds
    try:
        rec_donates, lig_donates = cm.find_hbonds()
        for ri, _li in rec_donates:
            c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
            sc = _is_sidechain(ri)
            _add(
                c,
                int(receptor.res_id[ri]),
                f"type:hydrogen_bonds__protisdon:True__sidechain:{sc}",
            )
        for ri, _li in lig_donates:
            c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
            sc = _is_sidechain(ri)
            _add(
                c,
                int(receptor.res_id[ri]),
                f"type:hydrogen_bonds__protisdon:False__sidechain:{sc}",
            )
    except Exception as e:
        log.warning(f"run_peppr_interactions: find_hbonds failed: {e}")
        failed_interaction_types.append("hydrogen_bonds")

    # Salt bridges
    try:
        salt_bridges = cm.find_salt_bridges()
        for ri, _li in salt_bridges:
            c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
            _add(c, int(receptor.res_id[ri]), "type:salt_bridges__protispos:True")
    except Exception as e:
        log.warning(f"run_peppr_interactions: find_salt_bridges failed: {e}")
        failed_interaction_types.append("salt_bridges")

    # Pi-stacking (deduplicate per residue)
    try:
        from biotite.structure import PiStacking

        stacking = cm.find_stacking_interactions()
        seen_stacking: set[tuple[str, int, str]] = set()
        for rec_idx, _lig_idx, stack_type in stacking:
            ri = rec_idx[0]
            c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
            stype = "T" if stack_type == PiStacking.PERPENDICULAR else "P"
            key = (c, int(receptor.res_id[ri]), stype)
            if key not in seen_stacking:
                seen_stacking.add(key)
                _add(c, int(receptor.res_id[ri]), f"type:pi_stacks__stack_type:{stype}")
    except Exception as e:
        log.warning(f"run_peppr_interactions: find_stacking_interactions failed: {e}")
        failed_interaction_types.append("pi_stacks")

    # Pi-cation
    try:
        pi_cation = cm.find_pi_cation_interactions()
        for rec_idx, _lig_idx, cation_in_receptor in pi_cation:
            ri = rec_idx[0]
            c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
            if cation_in_receptor:
                _add(
                    c,
                    int(receptor.res_id[ri]),
                    "type:pi_cation__lig_group:Aromatic__protcharged:True",
                )
            else:
                _add(
                    c,
                    int(receptor.res_id[ri]),
                    "type:pi_cation__lig_group:Cation__protcharged:False",
                )
    except Exception as e:
        log.warning(f"run_peppr_interactions: find_pi_cation_interactions failed: {e}")
        failed_interaction_types.append("pi_cation")

    # Halogen bonds
    try:
        from peppr.common import (
            ACCEPTOR_PATTERN,
            HALOGEN_DISTANCE_SCALING,
            HALOGEN_PATTERN,
        )

        halogen_bonds = cm.find_contacts_by_pattern(
            ACCEPTOR_PATTERN,
            HALOGEN_PATTERN,
            HALOGEN_DISTANCE_SCALING,
        )
        for ri, _li in halogen_bonds:
            c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
            sc = _is_sidechain(ri)
            _add(c, int(receptor.res_id[ri]), f"type:halogen_bonds__sidechain:{sc}")
    except Exception as e:
        log.warning(f"run_peppr_interactions: halogen_bonds failed: {e}")
        failed_interaction_types.append("halogen_bonds")

    # Water bridges (via plinder patch — peppr public doesn't have this yet)
    try:
        if waters.array_length() > 0:
            from peppr.common import DONOR_PATTERN
            from peppr.contacts import find_atoms_by_pattern

            receptor_donors = set(
                find_atoms_by_pattern(cm._binding_site_mol, DONOR_PATTERN)
            )
            w_bridges = find_water_bridges(receptor, ligand, waters)
            for rec_idx, _lig_idx, water_idx in w_bridges:
                ri = rec_idx[0]
                wi = water_idx[0]
                # Check if receptor atom is a donor
                bs_idx = None
                for j, orig_idx in enumerate(cm._binding_site_indices):
                    if orig_idx == ri:
                        bs_idx = j
                        break
                protisdon = bs_idx in receptor_donors if bs_idx is not None else True
                c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
                _add(
                    c,
                    int(receptor.res_id[ri]),
                    f"type:water_bridges__protisdon:{protisdon}",
                )
                w_chain = chain_mapping.get(waters.chain_id[wi], waters.chain_id[wi])
                water_set.add((w_chain, int(waters.res_id[wi])))
    except Exception as e:
        log.warning(f"run_peppr_interactions: find_water_bridges failed: {e}")
        failed_interaction_types.append("water_bridges")

    # Metal bridges (via plinder patch — peppr public doesn't have this yet)
    try:
        if metals.array_length() > 0:
            m_bridges = find_metal_bridges(receptor, ligand, metals)
            for rec_idx, _lig_idx, metal_idx in m_bridges:
                ri = rec_idx[0]
                mi = metal_idx[0]
                c = chain_mapping.get(receptor.chain_id[ri], receptor.chain_id[ri])
                metal_elem = metals.element[mi]
                _add(
                    c,
                    int(receptor.res_id[ri]),
                    f"type:metal_complexes__metal_type:{metal_elem}",
                )
    except Exception as e:
        log.warning(f"run_peppr_interactions: find_metal_bridges failed: {e}")
        failed_interaction_types.append("metal_complexes")

    return interaction_hashes, water_set, failed_interaction_types
