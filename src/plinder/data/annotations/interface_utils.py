# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Protein-interface detection in deposited biological assemblies."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Mapping

import biotite.structure as struc
import numpy as np
import pyarrow as pa

if TYPE_CHECKING:
    from plinder.data.annotations.ligand_utils import BiounitSpatialIndex
    from plinder.data.annotations.protein_utils import Chain

PROTEIN_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})
DEFAULT_MIN_INTERFACE_RESIDUES = 7
MIN_INTERFACE_RESIDUES_METADATA_KEY = b"plinder.interface.min_interface_residues"
PRODIGY_CONTACT_RADIUS = 5.0
PRODIGY_FEATURE_NAMES = (
    "CP",
    "AC",
    "AP",
    "AA",
    "ALA",
    "CYS",
    "GLU",
    "ASP",
    "GLY",
    "PHE",
    "ILE",
    "HIS",
    "MET",
    "LEU",
    "GLN",
    "PRO",
    "SER",
    "ARG",
    "THR",
    "VAL",
    "TYR",
    "link_density",
)
PRODIGY_AMINO_ACID_CLASSES = {
    "ALA": "A",
    "CYS": "A",
    "GLU": "C",
    "ASP": "C",
    "GLY": "A",
    "PHE": "A",
    "ILE": "A",
    "HIS": "C",
    "LYS": "C",
    "MET": "A",
    "LEU": "A",
    "ASN": "P",
    "GLN": "P",
    "PRO": "A",
    "SER": "P",
    "ARG": "C",
    "THR": "P",
    "TRP": "A",
    "VAL": "A",
    "TYR": "A",
}
# The model and feature ordering are pinned to official PRODIGY-cryst:
# https://github.com/haddocking/prodigy-cryst/tree/29d99d535cae3014608c9394b8809619a60fcfa2
# Its legacy scikit-learn pickle is represented as numeric arrays so the same
# trees can be evaluated on every Python version supported by Plinder.
_PRODIGY_MODEL_PATH = Path(__file__).parent / "static_files" / "prodigy_classifier.npz"

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
        ("prodigy_is_annotated", pa.bool_()),
        ("prodigy_label", pa.string()),
        ("prodigy_probability_bio", pa.float32()),
        ("prodigy_link_density", pa.float32()),
        ("prodigy_intermolecular_contacts", pa.int32()),
        ("prodigy_charged_charged_contacts", pa.int32()),
        ("prodigy_charged_polar_contacts", pa.int32()),
        ("prodigy_charged_apolar_contacts", pa.int32()),
        ("prodigy_polar_polar_contacts", pa.int32()),
        ("prodigy_apolar_polar_contacts", pa.int32()),
        ("prodigy_apolar_apolar_contacts", pa.int32()),
    ]
)


@dataclass(frozen=True)
class ProdigyCrystalAnnotation:
    """PRODIGY-cryst contact features and BIO/XTAL prediction."""

    label: str
    probability_bio: float
    link_density: float
    intermolecular_contacts: int
    charged_charged_contacts: int
    charged_polar_contacts: int
    charged_apolar_contacts: int
    polar_polar_contacts: int
    apolar_polar_contacts: int
    apolar_apolar_contacts: int


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
    prodigy: ProdigyCrystalAnnotation | None = None

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
        row: dict[str, object] = {
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
            "prodigy_is_annotated": self.prodigy is not None,
        }
        prodigy_fields = {
            "prodigy_label": "label",
            "prodigy_probability_bio": "probability_bio",
            "prodigy_link_density": "link_density",
            "prodigy_intermolecular_contacts": "intermolecular_contacts",
            "prodigy_charged_charged_contacts": "charged_charged_contacts",
            "prodigy_charged_polar_contacts": "charged_polar_contacts",
            "prodigy_charged_apolar_contacts": "charged_apolar_contacts",
            "prodigy_polar_polar_contacts": "polar_polar_contacts",
            "prodigy_apolar_polar_contacts": "apolar_polar_contacts",
            "prodigy_apolar_apolar_contacts": "apolar_apolar_contacts",
        }
        row.update(
            {
                column: (
                    getattr(self.prodigy, attribute)
                    if self.prodigy is not None
                    else None
                )
                for column, attribute in prodigy_fields.items()
            }
        )
        return row


@cache
def _load_prodigy_classifier() -> dict[str, np.ndarray]:
    """Load the official PRODIGY-cryst random forest once per process."""
    if not _PRODIGY_MODEL_PATH.is_file():
        raise FileNotFoundError(
            f"missing packaged PRODIGY-cryst classifier: {_PRODIGY_MODEL_PATH}"
        )
    with np.load(_PRODIGY_MODEL_PATH, allow_pickle=False) as model:
        arrays = {name: model[name] for name in model.files}
    if tuple(str(value) for value in arrays["feature_names"]) != PRODIGY_FEATURE_NAMES:
        raise ValueError("packaged PRODIGY-cryst classifier has incompatible features")
    if tuple(str(value) for value in arrays["classes"]) != ("BIO", "XTAL"):
        raise ValueError("packaged PRODIGY-cryst classifier has incompatible classes")
    return arrays


def _predict_prodigy_probability(features: np.ndarray) -> tuple[str, float]:
    """Evaluate the exact trees from official PRODIGY-cryst."""
    model = _load_prodigy_classifier()
    offsets = model["tree_offsets"]
    left = model["children_left"]
    right = model["children_right"]
    node_features = model["feature"]
    thresholds = model["threshold"]
    leaf_probabilities = model["probability"]
    probabilities = np.zeros(2, dtype=np.float64)
    for tree_index in range(len(offsets) - 1):
        node = int(offsets[tree_index])
        while left[node] >= 0:
            if features[node_features[node]] <= thresholds[node]:
                node = int(left[node])
            else:
                node = int(right[node])
        probabilities += leaf_probabilities[node]
    probabilities /= len(offsets) - 1
    return ("BIO" if probabilities[0] >= probabilities[1] else "XTAL"), float(
        probabilities[0]
    )


def annotate_prodigy_crystal_interface(
    atoms: struc.AtomArray,
    *,
    chain_1: str,
    chain_2: str,
    spatial_index: BiounitSpatialIndex,
) -> ProdigyCrystalAnnotation | None:
    """Classify one chain interface using official PRODIGY-cryst.

    Its contact definitions, feature ordering and fixed random forest are
    preserved.  Biotite supplies the already reconstructed biological
    assembly, and NumPy evaluates the legacy scikit-learn trees.
    Interfaces containing unsupported modified residues remain in the dataset
    with null PRODIGY fields.
    """
    chain_1_mask = atoms.chain_id == chain_1
    chain_2_mask = atoms.chain_id == chain_2
    if not np.any(chain_1_mask) or not np.any(chain_2_mask):
        raise ValueError(
            f"interface chains are absent from assembly: {chain_1}, {chain_2}"
        )
    chain_2_atom_mask = np.asarray(chain_2_mask, dtype=bool)
    contact_pairs: set[tuple[tuple[int, str], tuple[int, str]]] = set()
    chain_1_indices = np.flatnonzero(chain_1_mask)
    for atom_index in chain_1_indices:
        neighbors = np.asarray(
            spatial_index.cell_list.get_atoms(
                atoms.coord[atom_index], radius=PRODIGY_CONTACT_RADIUS
            ),
            dtype=int,
        ).reshape(-1)
        neighbors = neighbors[(neighbors >= 0) & (neighbors < len(atoms))]
        neighbors = neighbors[chain_2_atom_mask[neighbors]]
        residue_1 = (int(atoms.res_id[atom_index]), str(atoms.ins_code[atom_index]))
        for target_index in neighbors:
            residue_2 = (
                int(atoms.res_id[target_index]),
                str(atoms.ins_code[target_index]),
            )
            contact_pairs.add((residue_1, residue_2))
    if not contact_pairs:
        return None

    residue_name_1 = {
        (int(atoms.res_id[index]), str(atoms.ins_code[index])): str(
            atoms.res_name[index]
        )
        for index in np.flatnonzero(chain_1_mask)
    }
    residue_name_2 = {
        (int(atoms.res_id[index]), str(atoms.ins_code[index])): str(
            atoms.res_name[index]
        )
        for index in np.flatnonzero(chain_2_mask)
    }
    contacted_residue_names = {
        residue_name_1[residue_1] for residue_1, _ in contact_pairs
    } | {residue_name_2[residue_2] for _, residue_2 in contact_pairs}
    if not contacted_residue_names.issubset(PRODIGY_AMINO_ACID_CLASSES):
        return None
    bins = {name: 0 for name in ("AA", "PP", "CC", "AP", "CP", "AC")}
    bins.update({name: 0 for name in PRODIGY_AMINO_ACID_CLASSES})
    for residue_1, residue_2 in contact_pairs:
        name_1 = residue_name_1[residue_1]
        name_2 = residue_name_2[residue_2]
        contact_type = "".join(
            sorted(
                (
                    PRODIGY_AMINO_ACID_CLASSES[name_1],
                    PRODIGY_AMINO_ACID_CLASSES[name_2],
                )
            )
        )
        bins[contact_type] += 1
        bins[name_1] += 1
        bins[name_2] += 1
    contacted_1 = {pair[0] for pair in contact_pairs}
    contacted_2 = {pair[1] for pair in contact_pairs}
    link_density = len(contact_pairs) / (len(contacted_1) * len(contacted_2))
    features = np.asarray(
        [float(bins[name]) for name in PRODIGY_FEATURE_NAMES[:-1]] + [link_density],
        dtype=np.float64,
    )
    label, probability_bio = _predict_prodigy_probability(features)
    return ProdigyCrystalAnnotation(
        label=label,
        probability_bio=probability_bio,
        link_density=link_density,
        intermolecular_contacts=len(contact_pairs),
        charged_charged_contacts=bins["CC"],
        charged_polar_contacts=bins["CP"],
        charged_apolar_contacts=bins["AC"],
        polar_polar_contacts=bins["PP"],
        apolar_polar_contacts=bins["AP"],
        apolar_apolar_contacts=bins["AA"],
    )


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
    min_interface_residues: int = DEFAULT_MIN_INTERFACE_RESIDUES,
    annotate_prodigy: bool = True,
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
                prodigy=(
                    annotate_prodigy_crystal_interface(
                        atoms,
                        chain_1=chain_1,
                        chain_2=chain_2,
                        spatial_index=spatial_index,
                    )
                    if annotate_prodigy
                    else None
                ),
            )
        )
    return interfaces


def protein_interfaces_to_table(
    interfaces: Iterable[ProteinInterface],
    *,
    min_interface_residues: int = DEFAULT_MIN_INTERFACE_RESIDUES,
) -> pa.Table:
    """Return a typed table with the ingest threshold in schema metadata."""
    if min_interface_residues < 1:
        raise ValueError("minimum interface residues must be positive")
    schema = INTERFACE_ANNOTATION_SCHEMA.with_metadata(
        {MIN_INTERFACE_RESIDUES_METADATA_KEY: str(min_interface_residues).encode()}
    )
    return pa.Table.from_pylist(
        [interface.to_row() for interface in interfaces],
        schema=schema,
    )


def min_interface_residues_from_schema(schema: pa.Schema) -> int:
    """Read the frozen interface threshold from one V3 interface table."""
    metadata = schema.metadata or {}
    value = metadata.get(MIN_INTERFACE_RESIDUES_METADATA_KEY)
    if value is None:
        raise ValueError("interface table does not record its minimum residue count")
    try:
        threshold = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid minimum interface residue metadata") from exc
    if threshold < 1:
        raise ValueError("minimum interface residues must be positive")
    return threshold
