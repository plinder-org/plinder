# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Protein-interface detection in deposited biological assemblies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Mapping

import biotite.structure as struc
import numpy as np
import pyarrow as pa

if TYPE_CHECKING:
    from plinder.data.annotations.ligand_utils import BiounitSpatialIndex
    from plinder.data.annotations.protein_utils import Chain

PROTEIN_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})

INTERFACE_ANNOTATION_SCHEMA = pa.schema(
    [
        ("entry_pdb_id", pa.string()),
        ("system_id", pa.string()),
        ("system_biounit_id", pa.string()),
        ("interface_chain_1", pa.string()),
        ("interface_chain_2", pa.string()),
        ("interface_chain_1_residue_numbers", pa.list_(pa.int32())),
        ("interface_chain_1_residue_indices", pa.list_(pa.int32())),
        ("interface_chain_2_residue_numbers", pa.list_(pa.int32())),
        ("interface_chain_2_residue_indices", pa.list_(pa.int32())),
        ("interface_num_contact_residue_pairs", pa.int64()),
    ]
)


def interface_system_id(
    pdb_id: str,
    biounit_id: str,
    chain_1: str,
    chain_2: str,
) -> str:
    """Return the system ID for an unordered assembly-chain pair."""
    if chain_1 == chain_2:
        raise ValueError("an interface requires two distinct chain instances")
    first, second = sorted((str(chain_1), str(chain_2)))
    return f"{pdb_id.lower()}__{biounit_id}__{first}--{second}"


@dataclass(frozen=True)
class ProteinInterface:
    """One unordered protein-chain interface in one biological assembly."""

    pdb_id: str
    biounit_id: str
    chain_1: str
    chain_2: str
    chain_1_residue_numbers: tuple[int, ...]
    chain_1_residue_indices: tuple[int, ...]
    chain_2_residue_numbers: tuple[int, ...]
    chain_2_residue_indices: tuple[int, ...]
    num_contact_residue_pairs: int

    def __post_init__(self) -> None:
        if self.chain_1 >= self.chain_2:
            raise ValueError("interface chains must use canonical sorted ordering")
        if len(self.chain_1_residue_numbers) != len(self.chain_1_residue_indices):
            raise ValueError("chain 1 interface residue mappings have unequal lengths")
        if len(self.chain_2_residue_numbers) != len(self.chain_2_residue_indices):
            raise ValueError("chain 2 interface residue mappings have unequal lengths")
        if self.num_contact_residue_pairs < 1:
            raise ValueError("an interface must contain at least one residue contact")

    @property
    def system_id(self) -> str:
        return interface_system_id(
            self.pdb_id,
            self.biounit_id,
            self.chain_1,
            self.chain_2,
        )

    def to_row(self) -> dict[str, object]:
        """Return one Arrow-compatible annotation row."""
        return {
            "entry_pdb_id": self.pdb_id,
            "system_id": self.system_id,
            "system_biounit_id": self.biounit_id,
            "interface_chain_1": self.chain_1,
            "interface_chain_2": self.chain_2,
            "interface_chain_1_residue_numbers": list(self.chain_1_residue_numbers),
            "interface_chain_1_residue_indices": list(self.chain_1_residue_indices),
            "interface_chain_2_residue_numbers": list(self.chain_2_residue_numbers),
            "interface_chain_2_residue_indices": list(self.chain_2_residue_indices),
            "interface_num_contact_residue_pairs": self.num_contact_residue_pairs,
        }


def _asym_id(instance_chain: str) -> str:
    """Extract a label asym ID from a ``copy.asym`` assembly-chain ID."""
    try:
        _, asym_id = instance_chain.split(".", maxsplit=1)
    except ValueError as exc:
        raise ValueError(
            f"biological-assembly chain has no copy prefix: {instance_chain!r}"
        ) from exc
    return asym_id


def _residue_mappings(
    *,
    instance_chain: str,
    residue_numbers: Iterable[int],
    chains: Mapping[str, Chain],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Map resolved label sequence numbers to zero-based structure indices."""
    asym_id = _asym_id(instance_chain)
    chain = chains[asym_id]
    numbers = tuple(sorted(set(int(value) for value in residue_numbers)))
    try:
        indices = tuple(int(chain.residues[number].index) for number in numbers)
    except KeyError as exc:
        raise ValueError(
            f"assembly interface residue is absent from ASU chain {asym_id}: "
            f"{exc.args[0]}"
        ) from exc
    return numbers, indices


def detect_protein_interfaces(
    atoms: struc.AtomArray,
    *,
    pdb_id: str,
    biounit_id: str,
    chains: Mapping[str, Chain],
    contact_radius: float = 10.0,
    min_chain_length: int = 12,
    min_interface_residues: int = 3,
    spatial_index: BiounitSpatialIndex | None = None,
) -> list[ProteinInterface]:
    """Detect protein-chain interfaces from backbone contacts.

    Every eligible polypeptide chain has a SEQRES length of at least
    ``min_chain_length``.  A residue belongs to an interface when one of its
    ``N``, ``CA``, ``C`` or ``O`` atoms lies within ``contact_radius`` of a
    backbone atom in another eligible assembly-chain instance.  Interfaces
    are retained only when both sides contain at least
    ``min_interface_residues`` resolved residues.
    """
    if contact_radius <= 0:
        raise ValueError("interface contact radius must be positive")
    if min_chain_length < 1 or min_interface_residues < 1:
        raise ValueError("interface chain and residue limits must be positive")
    if spatial_index is None:
        from plinder.data.annotations.ligand_utils import BiounitSpatialIndex

        spatial_index = BiounitSpatialIndex.from_atoms(atoms, contact_radius)
    elif spatial_index.max_radius < contact_radius:
        raise ValueError(
            "interface contact radius exceeds the biological-assembly index radius"
        )

    eligible_chains = {
        instance_chain
        for instance_chain in spatial_index.chain_ids
        if (chain := chains.get(_asym_id(instance_chain))) is not None
        and "polypeptide" in chain.chain_type_str.lower()
        and chain.length >= min_chain_length
    }
    if len(eligible_chains) < 2:
        return []

    eligible_backbone = np.fromiter(
        (
            str(instance_chain) in eligible_chains
            and str(atom_name) in PROTEIN_BACKBONE_ATOMS
            for instance_chain, atom_name in zip(atoms.chain_id, atoms.atom_name)
        ),
        dtype=bool,
        count=len(atoms),
    )
    backbone_indices = np.flatnonzero(eligible_backbone)
    if not len(backbone_indices):
        return []

    # Store residue-pair contacts in canonical chain order. Iterating each
    # eligible backbone atom once avoids an O(number_of_chains^2) chain-pair
    # scan for very large assemblies.
    contacts: dict[tuple[str, str], set[tuple[int, int]]] = {}
    for atom_index in backbone_indices:
        neighbors = np.asarray(
            spatial_index.cell_list.get_atoms(
                atoms.coord[atom_index], radius=contact_radius
            ),
            dtype=int,
        ).reshape(-1)
        neighbors = neighbors[(neighbors > atom_index) & (neighbors < len(atoms))]
        neighbors = neighbors[eligible_backbone[neighbors]]
        query_chain = str(atoms.chain_id[atom_index])
        query_residue = int(atoms.res_id[atom_index])
        for target_index in neighbors:
            target_chain = str(atoms.chain_id[target_index])
            if target_chain == query_chain:
                continue
            target_residue = int(atoms.res_id[target_index])
            if query_chain < target_chain:
                key = (query_chain, target_chain)
                residue_pair = (query_residue, target_residue)
            else:
                key = (target_chain, query_chain)
                residue_pair = (target_residue, query_residue)
            contacts.setdefault(key, set()).add(residue_pair)

    interfaces: list[ProteinInterface] = []
    for (chain_1, chain_2), residue_pairs in sorted(contacts.items()):
        chain_1_numbers, chain_1_indices = _residue_mappings(
            instance_chain=chain_1,
            residue_numbers=(pair[0] for pair in residue_pairs),
            chains=chains,
        )
        chain_2_numbers, chain_2_indices = _residue_mappings(
            instance_chain=chain_2,
            residue_numbers=(pair[1] for pair in residue_pairs),
            chains=chains,
        )
        if (
            len(chain_1_numbers) < min_interface_residues
            or len(chain_2_numbers) < min_interface_residues
        ):
            continue
        interfaces.append(
            ProteinInterface(
                pdb_id=pdb_id.lower(),
                biounit_id=str(biounit_id),
                chain_1=chain_1,
                chain_2=chain_2,
                chain_1_residue_numbers=chain_1_numbers,
                chain_1_residue_indices=chain_1_indices,
                chain_2_residue_numbers=chain_2_numbers,
                chain_2_residue_indices=chain_2_indices,
                num_contact_residue_pairs=len(residue_pairs),
            )
        )
    return interfaces


def protein_interfaces_to_table(
    interfaces: Iterable[ProteinInterface],
) -> pa.Table:
    """Return a typed table, including for entries with no interfaces."""
    return pa.Table.from_pylist(
        [interface.to_row() for interface in interfaces],
        schema=INTERFACE_ANNOTATION_SCHEMA,
    )
