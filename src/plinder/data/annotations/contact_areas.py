# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Chain-pair contact areas from a Voronota-LT radical tessellation.

One tessellation of a biological assembly yields the exact contact area between
every pair of chains in a single near-linear pass.  Areas are Voronoi contact
areas between van der Waals spheres expanded by a solvent probe, not buried
SASA, so a pair of chains that never touch has no entry rather than zero.

Cost, measured on the complete icosahedral assembly of 7kmx (1.44 million heavy
non-water atoms, 840 chains): about 13 s and 3.2 GB of peak memory, i.e. roughly
9 µs and 2.2 KB per atom, almost all of it inside the C++ tessellation.  Callers
that process arbitrary assemblies should cap the atom count (see
``tessellation_atom_limit`` on :class:`plinder.data.annotations.aggregate_annotations.Entry`).
"""

from __future__ import annotations

from collections.abc import Mapping

import biotite.structure as struc
import numpy as np
from numpy.typing import NDArray

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

DEFAULT_PROBE_RADIUS = 1.4
"""Solvent probe radius (Å) added to every atom before tessellating."""

ChainPairAreas = Mapping[tuple[str, str], float]


def tessellation_atom_mask(atoms: struc.AtomArray) -> NDArray[np.bool_]:
    """Select the atoms that take part in a tessellation: heavy, non-water."""
    return np.asarray(
        struc.filter_heavy(atoms) & ~struc.filter_solvent(atoms), dtype=np.bool_
    )


def chain_pair_contact_areas(
    atoms: struc.AtomArray,
    *,
    probe: float = DEFAULT_PROBE_RADIUS,
) -> dict[tuple[str, str], float]:
    """Return the Voronota-LT contact area (Å²) between every pair of chains.

    Hydrogens and waters are excluded before tessellating, so waters count as
    solvent rather than as a chain.  Radii follow Voronota-LT's molecular
    assignment rules keyed on residue and atom names.

    Parameters
    ----------
    atoms : struc.AtomArray
        Atoms of one biological assembly (or any set of chains) carrying the
        ``ins_code`` annotation biotite reads from mmCIF.
    probe : float
        Solvent probe radius in Å added to every atom.

    Returns
    -------
    dict[tuple[str, str], float]
        Contact area per unordered chain pair, keyed by the lexicographically
        sorted ``(chain_a, chain_b)`` tuple.  Pairs without contact are absent.
    """
    from voronotalt import voronotalt_python as vlt

    selected = atoms[tessellation_atom_mask(atoms)]
    if selected.array_length() == 0:
        return {}
    # biotite recreates a deleted ``ins_code`` category as an empty array when
    # slicing, so check the length rather than the presence of the annotation.
    insertion_codes = (
        selected.ins_code
        if "ins_code" in selected.get_annotation_categories()
        else None
    )
    if insertion_codes is None or len(insertion_codes) != selected.array_length():
        insertion_codes = np.full(selected.array_length(), "", dtype=str)
    # Fill the C++ vector in place: a Python list of SWIG proxies would cost
    # about 1 KB per atom on top of the tessellation itself.
    balls = vlt.VectorMolecularAtomBall()
    for chain_id, res_id, ins_code, res_name, atom_name, (x, y, z) in zip(
        selected.chain_id,
        selected.res_id,
        insertion_codes,
        selected.res_name,
        selected.atom_name,
        selected.coord,
    ):
        balls.push_back(
            vlt.MolecularAtomBall(
                str(chain_id),
                int(res_id),
                str(ins_code),
                str(res_name),
                str(atom_name),
                float(x),
                float(y),
                float(z),
            )
        )
    if balls.size() != selected.array_length():
        raise ValueError(
            "atom annotations have inconsistent lengths: built "
            f"{balls.size()} balls for {selected.array_length()} atoms"
        )
    tessellation = vlt.MolecularRadicalTessellation.from_atoms(
        balls,
        probe=probe,
        compute_only_inter_chain_contacts=True,
        record_inter_chain_contact_summaries=True,
        record_everything_possible=False,
    )
    areas: dict[tuple[str, str], float] = {}
    for summary in tessellation.inter_chain_contact_summaries:
        chain_a, chain_b = str(summary.ID1_chain), str(summary.ID2_chain)
        if chain_a == chain_b:
            continue
        pair = (chain_a, chain_b) if chain_a < chain_b else (chain_b, chain_a)
        areas[pair] = areas.get(pair, 0.0) + float(summary.area)
    return areas


def partner_contact_areas(
    chain_pair_areas: ChainPairAreas,
    member_chains: set[str],
) -> dict[str, float]:
    """Sum the contact area of one molecule's chains with every other chain.

    ``member_chains`` are the instance chains forming one molecule (several for
    a covalently linked multi-chain ligand); contacts between member chains are
    internal and ignored.  The result maps each partner chain to its total
    contact area with the molecule, sorted by chain ID.
    """
    partners: dict[str, float] = {}
    for (chain_a, chain_b), area in chain_pair_areas.items():
        if chain_a in member_chains and chain_b not in member_chains:
            partners[chain_b] = partners.get(chain_b, 0.0) + area
        elif chain_b in member_chains and chain_a not in member_chains:
            partners[chain_a] = partners.get(chain_a, 0.0) + area
    return dict(sorted(partners.items()))
