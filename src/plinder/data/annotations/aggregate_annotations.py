# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import re
import shutil
import typing as ty
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
from biotite.file import DeserializationError, InvalidFileError
from biotite.structure import filter_heavy
from pydantic import BeforeValidator, Field, PrivateAttr
from rdkit import RDLogger

from plinder.core.utils.log import setup_logger
from plinder.data.annotations.cif_utils import (
    build_biounit,
    get_chain_external_mappings,
    get_entry_info,
    get_entry_taxonomy,
    get_label_asym_sequences,
    get_model_count,
    get_structure_with_altloc,
    remove_nonphysical_bonds,
)
from plinder.data.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
    ResidueValidationThresholds,
)
from plinder.data.annotations.interaction_utils import (
    get_covalent_connections,
    get_symmetry_mate_contacts,
)
from plinder.data.annotations.interface_utils import (
    DEFAULT_MIN_INTERFACE_RESIDUES,
    ProteinInterface,
    detect_protein_interfaces,
)
from plinder.data.annotations.ligand_utils import (
    BiounitSpatialIndex,
    Ligand,
    get_artifact_codes,
    get_water_chain_ids,
    is_known_artifact_ligand,
    validate_chain_residue,
)
from plinder.data.annotations.protein_utils import (
    Chain,
    _is_polynucleotide,
    _is_polypeptide,
    detect_ligand_chains,
    detect_ligand_chains_from_cif,
    get_atom_site_author_ids,
    get_receptor_type,
)
from plinder.data.annotations.save_utils import save_ligands
from plinder.data.annotations.utils import (
    DocBaseModel,
    description_excluded_from_flat_export,
)

LOG = setup_logger(__name__)
RDLogger.DisableLog("rdApp.*")
SymmetryMateContacts = ty.Annotated[
    dict[tuple[str, int], dict[tuple[str, int], dict[int, set[int]]]],
    BeforeValidator(validate_chain_residue),
    Field(default_factory=dict),
]
CUSTOM_STRUCTURE_MODES = ("as_is", "pdb")


def _require_mmcif_path(path: Path) -> Path:
    """Reject legacy PDB files and unrelated structure formats."""
    path = Path(path)
    name = path.name.lower()
    if not name.endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        raise ValueError(
            f"custom structures must be mmCIF files (.cif/.mmcif, optionally "
            f"gzip-compressed), not {path.name!r}"
        )
    return path


def _chain_type_from_coordinates(atoms: struc.AtomArray) -> str:
    """Identify standard polymer chains when entity metadata is absent."""
    if np.all(struc.filter_solvent(atoms)):
        return "water"
    residue_starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    if len(residue_starts) < 2:
        return "non-polymer"
    protein_residues = 0
    nucleotide_residues = 0
    dna_residues = 0
    rna_residues = 0
    ambiguous_nucleotide_residues = 0
    for start, stop in zip(residue_starts[:-1], residue_starts[1:]):
        residue = atoms[start:stop]
        atom_names = set(residue.atom_name.astype(str))
        if np.any(struc.filter_amino_acids(residue)) or {
            "N",
            "CA",
            "C",
        }.issubset(atom_names):
            protein_residues += 1
        if np.any(struc.filter_nucleotides(residue)):
            nucleotide_residues += 1
            residue_name = str(residue.res_name[0]).upper()
            if residue_name in {"DA", "DC", "DG", "DT", "DU", "DI"}:
                dna_residues += 1
            elif residue_name in {"A", "C", "G", "U", "I"} or atom_names.intersection(
                {"O2'", "O2*"}
            ):
                rna_residues += 1
            else:
                ambiguous_nucleotide_residues += 1
    residue_count = len(residue_starts) - 1
    if protein_residues == residue_count:
        return "polypeptide(L)"
    if nucleotide_residues == residue_count:
        if dna_residues and not rna_residues and not ambiguous_nucleotide_residues:
            return "polydeoxyribonucleotide"
        if rna_residues and not dna_residues and not ambiguous_nucleotide_residues:
            return "polyribonucleotide"
        # Mixed or modified nucleotides without a decisive sugar marker must
        # not be mislabeled as pure RNA. The hybrid type keeps both receptor
        # classes visible downstream.
        return "polydeoxyribonucleotide/polyribonucleotide hybrid"
    return "non-polymer"


def _sequence_from_coordinates(atoms: struc.AtomArray) -> str:
    """Return the resolved polymer sequence, or an empty string on ambiguity."""
    try:
        sequences, _ = struc.to_sequence(atoms, allow_hetero=True)
    except (IndexError, TypeError, ValueError, struc.BadStructureError):
        return ""
    return str(sequences[0]) if len(sequences) == 1 else ""


class _LigandChainClasses(ty.NamedTuple):
    """Ligand-like asym IDs split by how ingest treats them.

    Monoatomic ions and known artifacts are deferred: they only become
    ligands when they connect to a proper primary ligand's pocket.
    """

    monoatomic_ion_asym_ids: set[str]
    known_artifact_asym_ids: set[str]
    primary_asym_ids: set[str]


def _selected_assembly_ids(
    available: ty.Iterable[str],
    selected: ty.Iterable[str] | None,
) -> list[str]:
    """Return a validated assembly subset in caller-provided order."""
    available_ids = list(dict.fromkeys(str(value) for value in available))
    if selected is None:
        return available_ids
    selected_values = [selected] if isinstance(selected, str) else selected
    selected_ids = list(dict.fromkeys(str(value) for value in selected_values))
    if not selected_ids:
        raise ValueError("assembly_ids must not be empty when provided")
    missing = sorted(set(selected_ids).difference(available_ids))
    if missing:
        raise ValueError(
            f"requested assembly IDs are absent from the mmCIF: {missing}; "
            f"available={available_ids}"
        )
    return selected_ids


def remove_alphabets(x: str) -> int:
    """Turn alphanumeric to integer

    Parameters
    ----------
    x : str
        Alphanumeric string

    Returns
    -------
    int
    """

    return int(re.sub("[^0-9]", "", str(x)))


@dataclass
class QualityCriteria:
    max_entry_resolution: float = 3.5
    max_entry_r: float = 0.4
    max_entry_rfree: float = 0.45
    max_entry_r_minus_rfree: float = 0.05
    ligand_max_num_unresolved_heavy_atoms: int = 0
    ligand_max_alt_count: int = 1
    ligand_min_average_occupancy: float = 0.8
    ligand_min_average_rscc: float = 0.8
    ligand_max_average_rsr: float = 0.3
    ligand_max_percent_outliers_clashes: float = 0
    ligand_max_fraction_atoms_with_crystal_contacts: float = 0
    pocket_max_num_unresolved_heavy_atoms: int = 0
    pocket_max_alt_count: int = 1
    pocket_min_average_occupancy: float = 0.8
    pocket_min_average_rscc: float = 0.8
    pocket_max_average_rsr: float = 0.3
    pocket_max_percent_outliers_clashes: int = 100


class System(DocBaseModel):
    pdb_id: str = Field(description="[EXCLUDE] PDB ID")
    biounit_id: str = Field(description="Biounit ID")
    id_legacy: str = Field(
        default="",
        description="Historical system ID using global assembly-operation chain instances",
    )
    ligands: list[Ligand] = Field(description="[EXCLUDE] List of Ligands in a systems")
    receptor_type: str = Field(
        description=(
            "Receptor polymer composition: protein, DNA, RNA, other, or a "
            "+-joined combination"
        )
    )
    ligand_validation: ResidueListValidation | None = Field(
        default=None,
        description="[EXCLUDE] Validation object for the ligand residues in the system",
    )
    pocket_validation: ResidueListValidation | None = Field(
        default=None,
        description="[EXCLUDE] Validation object for the system's pocket residues",
    )
    pass_criteria: bool | None = Field(
        default=None, description="Whether the system passes validation criteria"
    )  # TODO: remove as attribute and have as function

    """
    This class defines a system which includes a protein-ligand complex
    and its neighboring ligands and receptor residues

    """

    @classmethod
    def document_properties(
        cls, prefix: str
    ) -> ty.Generator[tuple[str, str | None, str], ty.Any, ty.Any]:
        """Describe model fields plus the flat columns emitted by ``format()``."""
        yield from super().document_properties(prefix)
        yield (
            f"{prefix}_water_residues",
            "list[str]",
            "Interacting water residues encoded as "
            "<instance>.<asym>_<residue_number>",
        )
        for mapping_name in ("CATH", "Pfam", "SCOP2", "SCOP2B", "UniProt"):
            yield (
                f"{prefix}_pocket_{mapping_name}",
                "str",
                f"Most frequent {mapping_name} mapping among pocket residues",
            )

    def proper_ligands(self) -> list[Ligand]:
        return [ligand for ligand in self.ligands if ligand.is_proper]

    @cached_property
    def protein_chains_asym_id(self) -> list[str]:
        """
        Interacting receptor chains (protein or nucleic acid) of the system.

        The historical property name is retained in the annotation schema.
        """
        return sorted(
            set(
                chain
                for ligand in self.ligands
                for chain in ligand.protein_chains_asym_id
            )
        )

    @cached_property
    def id_no_biounit(self) -> str:
        """
        ID of the system without the biounit
        """
        return "__".join(
            [
                self.pdb_id,
                "_".join(
                    x.split(".", maxsplit=1)[1] for x in self.protein_chains_asym_id
                ),
                "_".join(x.split(".", maxsplit=1)[1] for x in self.ligand_chains),
            ]
        )

    @cached_property
    def ligand_chains(self) -> list[str]:
        """
        Ligand chains of the system
        """
        return [f"{ligand.instance}.{ligand.asym_id}" for ligand in self.ligands]

    @cached_property
    def num_pocket_residues(self) -> int:
        """
        Number of pocket residues of the system
        """
        pocket_residues = set()
        for ligand in self.ligands:
            ligand_pocket_residues = ligand.pocket_residues
            for chain in ligand_pocket_residues:
                pocket_residues |= set(
                    (chain, residue) for residue in ligand_pocket_residues[chain]
                )
        return len(pocket_residues)

    @cached_property
    def proper_num_pocket_residues(self) -> int:
        """
        Number of pocket residues of the system excluding ions and artifacts
        """
        pocket_residues = set()
        for chain in self.pocket_residues:
            for residue in self.pocket_residues[chain]:
                pocket_residues.add((chain, residue))
        return len(pocket_residues)

    @cached_property
    def num_interactions(self) -> int:
        """
        Number of interactions of the system
        """
        return sum(l.num_interactions for l in self.ligands)

    @cached_property
    def proper_num_interactions(self) -> int:
        """
        Number of interactions of the system excluding ions and artifacts
        """
        return sum(l.num_interactions for l in self.proper_ligands())

    @cached_property
    def num_unique_interactions(self) -> int:
        """
        Number of unique interactions of the system
        """
        return sum(l.num_unique_interactions for l in self.ligands)

    @cached_property
    def proper_num_unique_interactions(self) -> int:
        """
        Number of unique interactions of the system excluding ions and artifacts
        """
        return sum(l.num_unique_interactions for l in self.proper_ligands())

    @cached_property
    def num_covalent_ligands(self) -> int:
        """
        Number of covalent ligands of the system
        """
        return sum(ligand.is_covalent for ligand in self.ligands)

    @cached_property
    def proper_num_covalent_ligands(self) -> int:
        """
        Number of covalent ligands of the system excluding ions and artifacts
        """
        return sum(ligand.is_covalent for ligand in self.proper_ligands())

    @cached_property
    def id(self) -> str:
        """
        ID of the system
        """
        return "__".join(
            [
                self.pdb_id,
                self.biounit_id,
                "_".join(self.protein_chains_asym_id),
                "_".join(self.ligand_chains),
            ]
        )

    @cached_property
    def system_type(self) -> str:
        """
        Type of the system (one of: holo, ion, artifact)
        """
        if any(not l.is_ion and not l.is_artifact for l in self.ligands):
            return "holo"
        elif any(l.is_ion for l in self.ligands):
            return "ion"
        else:
            return "artifact"

    @cached_property
    def has_binding_affinity(self) -> bool:
        """
        Whether any ligand in the system has a binding affinity from BindingDB
        """
        return any(l.binding_affinity is not None for l in self.ligands)

    @cached_property
    def pocket_residues(self) -> dict[str, dict[int, str]]:
        """[EXCLUDE] Pockets residues of the system"""
        all_residues: dict[str, dict[int, str]] = defaultdict(dict)
        for ligand in self.ligands:
            if not ligand.is_proper:
                continue
            ligand_pocket_residues = ligand.pocket_residues
            for chain in ligand_pocket_residues:
                all_residues[chain].update(ligand_pocket_residues[chain])
        return all_residues

    @cached_property
    def interactions(self) -> dict[str, dict[int, list[str]]]:
        """[EXCLUDE] Interactions of the system"""
        all_interactions: dict[str, dict[int, list[str]]] = defaultdict(
            lambda: defaultdict(list)
        )
        for ligand in self.ligands:
            if not ligand.is_proper:
                continue
            for chain in ligand.interactions:
                for residue in ligand.interactions[chain]:
                    all_interactions[chain][residue].extend(
                        ligand.interactions[chain][residue]
                    )
        return all_interactions

    @cached_property
    def interactions_counter(self) -> dict[str, dict[int, ty.Counter[str]]]:
        """[EXCLUDE] Counter of interactions of the system"""
        interactions_counter: dict[str, dict[int, ty.Counter[str]]] = {}
        for chain in self.interactions:
            interactions_counter[chain] = {}
            for residue in self.interactions[chain]:
                interactions_counter[chain][residue] = Counter(
                    self.interactions[chain][residue]
                )
        return interactions_counter

    def format_chains(
        self,
        chain_type: str,
        chains: dict[str, Chain],
        sub_chains: list[str] | None = None,
    ) -> dict[str, ty.Any]:
        if sub_chains is None:
            if chain_type == "protein":
                sub_chains = self.protein_chains_asym_id
            elif chain_type == "ligand":
                sub_chains = self.ligand_chains
            else:
                raise ValueError(f"chain_type={chain_type} requires sub_chains")
        sub_chain_list = [ch.split(".", maxsplit=1) for ch in sub_chains]

        sub_chains_data = [
            chains[c].format(int(instance)) for instance, c in sub_chain_list
        ]
        if chain_type == "other":
            for sub_chain, (_, asym_id) in zip(sub_chains_data, sub_chain_list):
                sub_chain["chain_type"] = chains[asym_id].chain_type_str

        if len(sub_chains_data) == 0:
            return {}
        data: dict[str, list[ty.Any]] = defaultdict(list)
        for sub_chain in sub_chains_data:
            for key in sub_chain:
                data[f"system_{chain_type}_chains_{key}"].append(sub_chain[key])
        return data

    @cached_property
    def num_protein_chains(self) -> int:
        """
        Number of interacting receptor chains of the system.
        """
        return len(self.protein_chains_asym_id)

    @cached_property
    def proper_num_protein_chains(self) -> int:
        """
        Number of interacting receptor chains excluding ions and artifacts.
        """
        return len(
            set(
                chain
                for ligand in self.proper_ligands()
                for chain in ligand.protein_chains_asym_id
            )
        )

    @cached_property
    def num_ligand_chains(self) -> int:
        """
        Number of ligand chains of the system
        """
        return len(self.ligands)

    @cached_property
    def proper_num_ligand_chains(self) -> int:
        """
        Number of ligand chains of the system excluding ions and artifacts
        """
        return len(self.proper_ligands())

    def format_validation(
        self, entry_pass_criteria: bool | None, criteria: QualityCriteria
    ) -> dict[str, ty.Any]:
        data = {}
        if self.ligand_validation:
            ligand_validation = self.ligand_validation.format()
            data.update({f"system_ligand_{k}": v for k, v in ligand_validation.items()})
        if self.pocket_validation:
            pocket_validation = self.pocket_validation.format()
            data.update({f"system_pocket_{k}": v for k, v in pocket_validation.items()})
        if (
            self.pocket_validation is None
            or self.ligand_validation is None
            or entry_pass_criteria is None
            or not entry_pass_criteria
            or data["system_ligand_validation_max_alt_count"] is None
            or data["system_ligand_validation_average_occupancy"] is None
            or data["system_ligand_validation_average_rscc"] is None
            or data["system_ligand_validation_average_rsr"] is None
            or data["system_ligand_validation_percent_outliers_clashes"] is None
            or data["system_pocket_validation_num_unresolved_heavy_atoms"] is None
            or data["system_pocket_validation_max_alt_count"] is None
            or data["system_pocket_validation_average_occupancy"] is None
            or data["system_pocket_validation_average_rscc"] is None
            or data["system_pocket_validation_average_rsr"] is None
            or data["system_pocket_validation_percent_outliers_clashes"] is None
        ):
            self.pass_criteria = False
        else:
            """
            Quality criteria for the system
            """
            quality = [
                # LIGAND
                self.num_unresolved_heavy_atoms is not None
                and self.num_unresolved_heavy_atoms
                <= self.num_covalent_ligands
                + criteria.ligand_max_num_unresolved_heavy_atoms,
                data["system_ligand_validation_max_alt_count"]
                <= criteria.ligand_max_alt_count,
                data["system_ligand_validation_average_occupancy"]
                >= criteria.ligand_min_average_occupancy,
                data["system_ligand_validation_average_rscc"]
                >= criteria.ligand_min_average_rscc,
                data["system_ligand_validation_average_rsr"]
                <= criteria.ligand_max_average_rsr,
                data["system_ligand_validation_percent_outliers_clashes"]
                <= criteria.ligand_max_percent_outliers_clashes,
                self.fraction_atoms_with_crystal_contacts is not None
                and self.fraction_atoms_with_crystal_contacts
                <= criteria.ligand_max_fraction_atoms_with_crystal_contacts,
                # POCKET
                data["system_pocket_validation_num_unresolved_heavy_atoms"]
                <= criteria.pocket_max_num_unresolved_heavy_atoms,
                data["system_pocket_validation_max_alt_count"]
                <= criteria.pocket_max_alt_count,
                data["system_pocket_validation_average_occupancy"]
                >= criteria.pocket_min_average_occupancy,
                data["system_pocket_validation_average_rscc"]
                >= criteria.pocket_min_average_rscc,
                data["system_pocket_validation_average_rsr"]
                <= criteria.pocket_max_average_rsr,
                data["system_pocket_validation_percent_outliers_clashes"]
                <= criteria.pocket_max_percent_outliers_clashes,
            ]
            self.pass_criteria = all(quality)
        data["system_pass_validation_criteria"] = self.pass_criteria
        return data

    def format(
        self,
        chains: dict[str, Chain],
        entry_pass_criteria: bool | None,
        criteria: QualityCriteria = QualityCriteria(),
    ) -> dict[str, ty.Any]:
        data: dict[str, ty.Any] = defaultdict(str)
        for field, (description, _) in self.get_descriptions_and_types().items():
            if description_excluded_from_flat_export(description):
                continue
            if not field.startswith("system_"):
                name = f"system_{field}"
            else:
                name = field
            data[name] = getattr(self, field, None)

        pocket_mapping = self.get_pocket_domains(chains)
        for mapping in pocket_mapping:
            data[f"system_pocket_{mapping}"] = pocket_mapping[mapping]
        data.update(self.format_validation(entry_pass_criteria, criteria))
        for chain_type in ["protein", "ligand"]:
            data.update(self.format_chains(chain_type, chains))
        data["system_water_residues"] = sorted(
            f"{chain_id}_{residue_number}"
            for chain_id, residue_numbers in self.waters.items()
            for residue_number in set(residue_numbers)
        )
        return data

    @cached_property
    def waters(self) -> dict[str, list[int]]:
        """[EXCLUDE] Waters interacting with any of the ligands in the system"""
        waters: dict[str, list[int]] = defaultdict(list)
        for ligand in self.ligands:
            for chain in ligand.waters:
                waters[chain] += ligand.waters[chain]
        return waters

    def select_waters(self) -> str:
        query = []
        for chain in self.waters:
            chain_query = " or ".join(f"rnum={resnum}" for resnum in self.waters[chain])
            query.append(f"(chain='{chain}' and ({chain_query}))")
        return " or ".join(query)

    @cached_property
    def num_crystal_contacted_residues(self) -> int:
        """
        Number of residues from other symmetry mates which are in contact with any ligand in the system.
        """
        residues = set()
        for ligand in self.ligands:
            residues |= set(ligand.crystal_contacts.keys())
        return len(residues)

    @cached_property
    def num_atoms_with_crystal_contacts(self) -> int:
        """
        Number of atoms in the system ligands which are in contact with residues from other symmetry mates.
        """
        return sum(ligand.num_atoms_with_crystal_contacts for ligand in self.ligands)

    @cached_property
    def num_heavy_atoms(self) -> int | None:
        """
        Number of heavy atoms in the system ligands
        """
        if any(ligand.num_heavy_atoms is None for ligand in self.ligands):
            return None
        return sum(
            ligand.num_heavy_atoms
            for ligand in self.ligands
            if ligand.num_heavy_atoms is not None
        )

    @cached_property
    def num_resolved_heavy_atoms(self) -> int | None:
        """
        Number of resolved heavy atoms in the system ligands
        """
        if any(ligand.num_resolved_heavy_atoms is None for ligand in self.ligands):
            return None
        return sum(
            ligand.num_resolved_heavy_atoms
            for ligand in self.ligands
            if ligand.num_resolved_heavy_atoms is not None
        )

    @cached_property
    def ligand_max_qed(self) -> float:
        """
        Maximum QED of the system ligands
        """
        return max(
            ligand.qed if ligand.qed is not None else -1.0 for ligand in self.ligands
        )

    @cached_property
    def ligand_max_molecular_weight(self) -> float:
        """
        Maximum molecular weight of the system ligands
        """
        max_vals = [
            ligand.molecular_weight if ligand.molecular_weight is not None else -1.0
            for ligand in self.ligands
            if ligand.molecular_weight is not None
        ]
        if len(max_vals):
            return max(max_vals)
        else:
            return 1

    @cached_property
    def proper_ligand_max_molecular_weight(self) -> float:
        """
        Maximum molecular weight of the system ligands excluding ions and artifacts
        """
        ligands = self.proper_ligands()
        if len(ligands) == 0:
            return -1.0
        weights = [
            ligand.molecular_weight if ligand.molecular_weight is not None else -1.0
            for ligand in ligands
        ]
        return max(weights)

    @cached_property
    def num_unresolved_heavy_atoms(self) -> int | None:
        """
        Number of unresolved heavy atoms in the system ligands
        """
        if any(ligand.num_unresolved_heavy_atoms is None for ligand in self.ligands):
            return None
        return sum(
            ligand.num_unresolved_heavy_atoms
            for ligand in self.ligands
            if ligand.num_unresolved_heavy_atoms is not None
        )

    @cached_property
    def fraction_atoms_with_crystal_contacts(self) -> float | None:
        """
        Fraction of atoms in the system ligands which are in contact with residues from other symmetry mates.
        """
        if not self.num_heavy_atoms:
            return None
        return self.num_atoms_with_crystal_contacts / self.num_heavy_atoms

    def selection(self, include_waters: bool = True) -> str:
        ligand_selection = " or ".join(
            f"({ligand.selection})" for ligand in self.ligands
        )
        protein_selection = " or ".join(
            f"(cname='{chain}')" for chain in self.protein_chains_asym_id
        )
        selection = f"({ligand_selection}) or ({protein_selection})"
        if include_waters and len(self.waters):
            selection += f" or {self.select_waters()}"
        return selection

    def set_validation(
        self,
        chains: dict[str, Chain],
        thresholds: ResidueValidationThresholds = ResidueValidationThresholds(),
    ) -> None:
        self.ligand_validation = ResidueListValidation.from_residues(
            [
                chains[c.split(".", maxsplit=1)[1]].residues[r].validation  # type: ignore
                for c in self.ligand_chains
                for r in chains[c.split(".", maxsplit=1)[1]].residues
            ],
            thresholds,
        )
        self.pocket_validation = ResidueListValidation.from_residues(
            [
                chains[c.split(".", maxsplit=1)[1]].residues[r].validation  # type: ignore
                for c in self.pocket_residues
                for r in self.pocket_residues[c]
            ],
            thresholds,
        )

    def get_pocket_domains(self, chains_dict: dict[str, Chain]) -> dict[str, str]:
        pocket_mapping: dict[str, dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )
        for (
            neighboring_chain,
            neighboring_residues_list,
        ) in self.pocket_residues.items():
            neighboring_chain = neighboring_chain.split(".", maxsplit=1)[-1]
            neighboring_residues_set = {int(i) for i in neighboring_residues_list}
            for mapping_name in chains_dict[neighboring_chain].mappings:
                if mapping_name == "BIRD":
                    continue
                for domain in chains_dict[neighboring_chain].mappings[mapping_name]:
                    unrolled_v = {
                        tuple(range(*[remove_alphabets(j) for j in i]))
                        for i in chains_dict[neighboring_chain].mappings[mapping_name][
                            domain
                        ]
                        if i is not None
                    }
                    unrolled_v_2 = {j for i in unrolled_v for j in i}
                    intersection = neighboring_residues_set.intersection(unrolled_v_2)
                    if len(intersection) > 0:
                        pocket_mapping[mapping_name][domain] += len(intersection)
        result = {}
        for k, values in pocket_mapping.items():
            result[k] = sorted(values.items(), key=lambda x: x[1], reverse=True)[0][0]
        return result


class Entry(DocBaseModel):
    _ligand_contacts_requested: bool = PrivateAttr(default=False)

    pdb_id: str = Field(
        default_factory=str,
        description="RCSB PDB ID. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entry.id.html",
    )
    release_date: str = Field(
        default_factory=str,
        description="RCSB structure release date. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_database_PDB_rev.date_original.html",
    )
    oligomeric_state: str | None = Field(
        default_factory=str,
        description="Author's provided description of quaternary structure in RCSB. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_pdbx_struct_assembly.oligomeric_details.html",
    )
    determination_method: str | None = Field(
        default_factory=str,
        description="RCSB method of structure determination. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl.method.html",
    )
    keywords: str | None = Field(
        default_factory=str,
        description="RCSB keywords describing the structure. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_keywords.pdbx_keywords.html",
    )
    pH: str | None = Field(
        default_factory=str,
        description="pH at which structure is solved. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl_crystal_grow.pH.html",
    )
    resolution: float | None = Field(
        default_factory=float,
        description="RCSB structure resolution. See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_refine.ls_d_res_high.html",
    )
    source_taxonomy_ids: list[int] = Field(
        default_factory=list,
        description="Distinct NCBI taxonomy IDs of the deposited entities' source organisms",
    )
    source_organism_names: list[str] = Field(
        default_factory=list,
        description="Distinct scientific names of the deposited entities' source organisms",
    )
    host_taxonomy_ids: list[int] = Field(
        default_factory=list,
        description="Distinct NCBI taxonomy IDs of recombinant expression hosts",
    )
    host_organism_names: list[str] = Field(
        default_factory=list,
        description="Distinct scientific names of recombinant expression hosts",
    )
    chains: dict[str, Chain] = Field(
        default_factory=dict,
        description="[EXCLUDE] Chains dictionary with chain name mapped to chain object",
    )
    ligand_like_chains: dict[str, str] = Field(
        default_factory=dict,
        description="[EXCLUDE] Chain: chain type for other ligand-like chains in the entry",
    )
    systems: dict[str, System] = Field(
        default_factory=dict,
        description="[EXCLUDE] System dictionary with system id mapped to system object",
    )
    interfaces: list[ProteinInterface] = Field(
        default_factory=list,
        description="[EXCLUDE] Protein-chain interfaces across deposited assemblies",
    )
    covalent_bonds: dict[str, list[tuple[str, str]]] = Field(
        default_factory=dict,
        description="[EXCLUDE] All covalent interactions in the entry as defined by mmcif annotations. They types are separated by dictionary key and they include: "
        + "covale: actual covalent linkage, metalc: other dative bond interactions like metal-ligand dative bond, "
        + "hydrog: strong hydrogen bonding of nucleic acid. For the purpose of covalent annotations, we use only covale for downstream processing.",
    )
    chain_to_seqres: dict[str, str] = Field(
        default_factory=dict,
        description="[EXCLUDE] Chain to sequence mapping",
    )
    validation: EntryValidation | None = Field(
        default=None, description="[EXCLUDE] Entry validation"
    )
    pass_criteria: bool | None = Field(
        default=None, description="Whether the entry passes validation criteria"
    )
    water_chains: list[str] = Field(
        default_factory=list, description="[EXCLUDE] Water chains in the entry"
    )
    biounit_chain_ids: dict[str, list[str]] = Field(
        default_factory=dict,
        description="[EXCLUDE] Resolved biological-assembly chain instances by assembly ID",
    )
    biounit_legacy_chain_ids: dict[str, dict[str, str]] = Field(
        default_factory=dict,
        description="[EXCLUDE] Canonical-to-historical chain instance IDs by assembly ID",
    )
    biounit_ligand_contact_counts: dict[str, dict[str, dict[str, int]]] = Field(
        default_factory=dict,
        description=(
            "[EXCLUDE] Counts of ion, artifact, and other ligand chains near "
            "each biological-assembly receptor chain instance"
        ),
    )
    # TODO: consider surfacing this as an exported plindex column (drop the
    # ``__`` prefix + wire column_descriptions/schema) so a dropped assembly is
    # visible downstream, not only in this in-memory metadata + the error log.
    failed_assembly_ids: list[str] = Field(
        default_factory=list,
        description="[EXCLUDE] Biological-assembly IDs biotite could not build; their "
        "systems are absent from this entry.",
    )
    symmetry_mate_contacts: SymmetryMateContacts = Field(
        default_factory=dict,
        description="[EXCLUDE] Symmetry mate contacts in the entry",
    )

    def prune(
        self,
        *,
        clear_non_pocket_residues: bool = True,
        load_for_scoring: bool = True,
        max_protein_chains: int = 5,
        max_ligand_chains: int = 5,
    ) -> Entry:
        """
        Update an entry in place to reduce memory footprint and discard
        unuseful systems for downstream processing.

        Parameters
        ----------
        clear_non_pocket_residues : bool, default=True
            remove residues not present in pockets
        load_for_scoring : bool, default=True
            remove systems that are not holo or have too many protein/ligand chains
        max_protein_chains : int, default=5
            only keep systems with leq max_protein_chains
        max_ligand_chains : int, default=5
            only keep systems with leq max_ligand_chains

        Returns
        -------
        self : Entry
            the updated entry
        """
        if clear_non_pocket_residues:
            self.clear_non_pocket_residues()
        if load_for_scoring:
            n_before = len(self.systems)
            self.systems = {
                s.id: s
                for s in self.systems.values()
                if s.system_type == "holo"
                and len(s.protein_chains_asym_id) <= max_protein_chains
                and len(s.ligand_chains) <= max_ligand_chains
            }
            n_dropped = n_before - len(self.systems)
            if n_dropped:
                LOG.info(
                    f"{self.pdb_id}: prune(load_for_scoring) dropped {n_dropped} "
                    f"of {n_before} systems (non-holo or exceeding "
                    f"{max_protein_chains} protein / {max_ligand_chains} ligand "
                    "chains)"
                )
        return self

    def _populate_chains(
        self,
        atoms: struc.AtomArray,
        block: pdbx.CIFBlock,
    ) -> None:
        """Set entry.chains and entry.water_chains from biotite data."""
        chain_starts = struc.get_chain_starts(atoms, add_exclusive_stop=True)
        chain_segments: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for start, stop in zip(chain_starts[:-1], chain_starts[1:]):
            chain_segments[str(atoms.chain_id[start])].append((int(start), int(stop)))
        solvent_mask = struc.filter_solvent(atoms)
        water_chains = {
            chain_id
            for chain_id, segments in chain_segments.items()
            if all(bool(solvent_mask[start:stop].all()) for start, stop in segments)
        }
        entity_by_asym: dict[str, str] = {}
        if "struct_asym" in block:
            struct_asym = block["struct_asym"]
            entity_by_asym = {
                str(asym_id): str(entity_id)
                for asym_id, entity_id in zip(
                    struct_asym["id"].as_array(),
                    struct_asym["entity_id"].as_array(),
                )
            }
        type_by_entity: dict[str, str] = {}
        if "entity" in block:
            entity = block["entity"]
            type_by_entity.update(
                {
                    str(entity_id): str(entity_type)
                    for entity_id, entity_type in zip(
                        entity["id"].as_array(), entity["type"].as_array()
                    )
                }
            )
        if "entity_poly" in block:
            entity_poly = block["entity_poly"]
            type_by_entity.update(
                {
                    str(entity_id): str(entity_type)
                    for entity_id, entity_type in zip(
                        entity_poly["entity_id"].as_array(),
                        entity_poly["type"].as_array(),
                    )
                }
            )
        auth_id_by_asym, residue_author_ids_by_asym = get_atom_site_author_ids(block)
        self.chains = {}
        # Chain metadata does not use bonds.  Temporarily detaching the global
        # BondList prevents every small chain slice from scanning and
        # reindexing the full entry bond graph.
        bonds = atoms.bonds
        atoms.bonds = None
        try:
            for chain_id in sorted(set(chain_segments) - water_chains):
                segments = chain_segments[chain_id]
                if len(segments) == 1:
                    start, stop = segments[0]
                    chain_atoms = atoms[start:stop]
                else:
                    chain_atoms = struc.concatenate(
                        [atoms[start:stop] for start, stop in segments]
                    )
                entity_id = entity_by_asym.get(chain_id, "")
                chain_type = type_by_entity.get(entity_id, "unknown")
                if chain_type == "unknown":
                    chain_type = _chain_type_from_coordinates(chain_atoms)
                if chain_id not in self.chain_to_seqres and (
                    _is_polypeptide(chain_type) or _is_polynucleotide(chain_type)
                ):
                    sequence = _sequence_from_coordinates(chain_atoms)
                    if sequence:
                        self.chain_to_seqres[chain_id] = sequence
                self.chains[chain_id] = Chain.from_cif_data(
                    chain_id,
                    block,
                    chain_atoms,
                    len(self.chain_to_seqres.get(chain_id, "")),
                    entity_id=entity_id,
                    auth_id=auth_id_by_asym.get(chain_id),
                    chain_type_str=chain_type,
                    residue_author_ids=residue_author_ids_by_asym.get(chain_id, {}),
                )
        finally:
            atoms.bonds = bonds

        self.water_chains = sorted(water_chains)

    @staticmethod
    def _clear_ligand_files(save_folder: Path | None, pdb_id: str) -> None:
        """Remove ligand SDFs left by a prior annotation of this entry.

        Re-ingest must never merge newly retained ligands with SDFs left in
        ``<save_folder>/<pdb_id>/ligand_files`` by an earlier run; clear the
        directory up front so even an early no-system return leaves it empty.
        """
        if save_folder is None:
            return
        ligand_dir = Path(save_folder) / pdb_id / "ligand_files"
        if ligand_dir.exists():
            shutil.rmtree(ligand_dir)

    @staticmethod
    def _load_clean_atoms(
        cif_file_obj: pdbx.CIFFile,
        *,
        source: str,
        no_bonds_error: str,
        multimodel_note: str = "",
        require_bonds: bool = True,
    ) -> struc.AtomArray:
        """Parse model 1 as a bond-clean, heavy-atom ``AtomArray``.

        Shared by both ingest paths so their structure loading and bond
        handling cannot drift apart. Loads model 1 with altlocs and
        ``include_bonds=True``, drops hydrogens, and prunes non-physical
        bonds.

        Parameters
        ----------
        cif_file_obj : pdbx.CIFFile
            Parsed CIF (already bond-enriched for custom structures).
        source : str
            Structure label for the multi-model warning, e.g.
            ``"PDB '1abc'"`` or ``"Custom CIF"``.
        no_bonds_error : str
            Message for the ``ValueError`` raised when biotite derives no
            bonds despite ``include_bonds=True``.
        multimodel_note : str, optional
            Extra guidance appended to the multi-model warning.
        require_bonds : bool, optional
            Raise when biotite derives no bonds. Ligand ingest needs the bond
            graph; interface-only and protein-only ingest tolerate its absence.

        Returns
        -------
        struc.AtomArray
            Heavy atoms of model 1 with a physical bond graph.

        Raises
        ------
        ValueError
            If ``require_bonds`` and biotite returns no bonds (``atoms.bonds is None``).

        Notes
        -----
        ``use_author_fields=False`` keeps ligand-chain detection keyed on
        label_seq_id, while ``auth_seq_id`` is loaded separately so author
        residue numbers survive into :meth:`Chain.from_cif_data`. Only
        model 1 is used; NMR ensembles and multi-sample predictions must be
        processed one model at a time.
        """
        n_models = get_model_count(cif_file_obj)
        if n_models > 1:
            LOG.warning(
                f"{source} has {n_models} models — using model 1 only."
                f"{multimodel_note}"
            )
        atoms = get_structure_with_altloc(
            cif_file_obj,
            model=1,
            use_author_fields=False,
            include_bonds=True,
            extra_fields=["auth_seq_id"],
        )
        atoms = atoms[filter_heavy(atoms)]
        if atoms.bonds is None and require_bonds:
            raise ValueError(no_bonds_error)
        remove_nonphysical_bonds(atoms)
        return atoms

    def _finalize(
        self,
        ligands: dict[str, Ligand],
        atoms: struc.AtomArray,
        *,
        save_folder: Path | None = None,
        min_shared_pocket_members: int = 3,
        cif_file_obj: pdbx.CIFFile | None = None,
        symmetry_mate_contact_threshold: float = 5.0,
        include_ligands: bool = True,
    ) -> None:
        """Set systems, label crystal contacts and chains, and save ligand SDFs.

        The shared tail of both ingest paths.

        Parameters
        ----------
        ligands : dict[str, Ligand]
            Candidate ligands collected across all biounits.
        atoms : struc.AtomArray
            Heavy-atom model used to write the ligand SDFs.
        save_folder : Path | None, optional
            Root for canonical ASU ligand SDFs; ``None`` skips saving.
        cif_file_obj : pdbx.CIFFile | None, optional
            Parsed CIF used for crystal-contact detection; ``None`` skips it.
        symmetry_mate_contact_threshold : float, optional
            Distance (Å) for symmetry-mate contacts.
        include_ligands : bool, optional
            Build ligand systems and write their SDFs. Interface-only and
            protein-only ingest skip both but still label the chains.

        Notes
        -----
        Crystal contacts are computed only when ``cif_file_obj`` is given
        (deposited structures with crystallographic symmetry) and only for
        ligands that ended up in a system; custom/predicted CIFs pass
        ``None`` to skip that step. One SDF is written per retained ligand,
        keyed by its primary asym and spanning every member chain, so
        covalently-linked ligand chains are saved as a single molecule.
        """
        if include_ligands:
            self.set_systems(
                ligands, min_shared_pocket_members=min_shared_pocket_members
            )
        if self.systems and cif_file_obj is not None:
            self.symmetry_mate_contacts = get_symmetry_mate_contacts(
                cif_file_obj,
                symmetry_mate_contact_threshold,
            )
            if self.symmetry_mate_contacts:
                for system in self.systems.values():
                    for ligand in system.ligands:
                        ligand.label_crystal_contacts(self.symmetry_mate_contacts)
        self.label_chains()
        retained_ligand_chain_groups = {
            ligand.asym_id: ligand.member_asym_ids
            for system in self.systems.values()
            for ligand in system.ligands
        }
        if include_ligands and save_folder is not None and retained_ligand_chain_groups:
            save_ligands(
                atoms,
                retained_ligand_chain_groups,
                save_folder / self.pdb_id / "ligand_files",
            )

    def _collect_ligands_from_biounit(
        self,
        biounit: struc.AtomArray,
        biounit_id: str,
        interaction_search_threshold: float,
        neighboring_residue_threshold: float,
        neighboring_ligand_threshold: float,
        data_dir: Path | None,
        ligand_smiles_dict: dict[str, str] | None = None,
        ligand_ccd_code_dict: dict[str, str] | None = None,
        ligand_asym_ids: set[str] | None = None,
        ligand_instance_chains: set[str] | None = None,
        water_chains: set[str] | None = None,
        spatial_index: BiounitSpatialIndex | None = None,
    ) -> dict[str, "Ligand"]:
        """Create Ligand objects for every ligand chain in a single biounit.

        ``ligand_smiles_dict`` is passed through to :meth:`Ligand.from_pli`
        and is only set by :meth:`Entry.from_custom_cif_file` — it lets
        user-supplied SMILES act as the CCD fallback for stereo
        validation and SMILES assignment on custom residues.

        ``ligand_ccd_code_dict`` records explicit custom-component to CCD
        references so the output ligand annotation reports the supplied code.

        ``ligand_asym_ids`` optionally limits work to selected ligand chains.
        PDB ingest uses this to probe non-ion ligands before calculating
        contacts for potentially thousands of monoatomic ions.

        ``ligand_instance_chains`` optionally selects exact biological-
        assembly copies, rather than every copy of an asymmetric-unit chain.

        ``water_chains`` may be precomputed once per biological assembly and
        shared across calls that process different ligand subsets.

        ``spatial_index`` shares the assembly CellList and atom hierarchy
        across every ligand and across the non-ion/ion passes.
        """
        ligands: dict[str, Ligand] = {}
        if spatial_index is None:
            spatial_index = BiounitSpatialIndex.from_atoms(
                biounit,
                max(
                    interaction_search_threshold,
                    neighboring_residue_threshold,
                    neighboring_ligand_threshold,
                ),
            )
        # Find ligand chains: chain_id format is "{instance}.{asym_id}"
        biounit_ligand_chains = [
            c
            for c in spatial_index.chain_ids
            if "." in c
            and c.split(".", maxsplit=1)[1] in self.ligand_like_chains
            and (
                ligand_asym_ids is None
                or c.split(".", maxsplit=1)[1] in ligand_asym_ids
            )
            and (ligand_instance_chains is None or c in ligand_instance_chains)
        ]
        if not biounit_ligand_chains:
            return ligands
        if water_chains is None:
            water_chains = get_water_chain_ids(biounit)
        # Covalently-linked ligand chains are one molecule (e.g. a macrocycle
        # whose non-polymer parts are deposited as separate chains) and are
        # built as a single ligand spanning all member chains.
        ligand_groups = self._covalent_ligand_groups(biounit_ligand_chains)
        group_count = len(ligand_groups)
        if group_count >= 100:
            LOG.info(
                "PDB %s assembly %s: processing %d ligand groups",
                self.pdb_id,
                biounit_id,
                group_count,
            )
        for ligand_index, group in enumerate(ligand_groups, start=1):
            if group_count >= 100 and ligand_index % 100 == 0:
                LOG.info(
                    "PDB %s assembly %s: processed %d/%d ligand groups",
                    self.pdb_id,
                    biounit_id,
                    ligand_index,
                    group_count,
                )
            # Residue numbers per member instance-chain in the group.
            member_residue_numbers: dict[str, list[int]] = {}
            for member_chain in sorted(group):
                member_atoms = spatial_index.take_atoms(
                    biounit,
                    spatial_index.atom_indices_for_chain(member_chain),
                    include_bonds=False,
                )
                member_residue_numbers[member_chain] = list(
                    dict.fromkeys(int(r) for r in member_atoms.res_id)
                )
            # Primary chain (deterministic): the first sorted member. The
            # ligand is keyed on it, but its atoms span every member chain.
            primary_chain = sorted(group)[0]
            primary_instance, primary_asym_id = primary_chain.split(".", maxsplit=1)
            ligand = Ligand.from_pli(
                pdb_id=self.pdb_id,
                biounit_id=biounit_id,
                biounit=biounit,
                ligand_instance=int(primary_instance),
                ligand_chain=self.chains[primary_asym_id],
                residue_numbers=member_residue_numbers[primary_chain],
                ligand_like_chains=self.ligand_like_chains,
                all_covalent_dict=self.covalent_bonds,
                interaction_search_threshold=interaction_search_threshold,
                neighboring_residue_threshold=neighboring_residue_threshold,
                neighboring_ligand_threshold=neighboring_ligand_threshold,
                data_dir=data_dir,
                chain_to_seqres=self.chain_to_seqres,
                ligand_smiles_dict=ligand_smiles_dict,
                ligand_ccd_code_dict=ligand_ccd_code_dict,
                water_chains=water_chains,
                spatial_index=spatial_index,
                member_residue_numbers=member_residue_numbers,
            )
            if ligand is not None:
                ligands[ligand.id] = ligand
        return ligands

    def _covalent_ligand_groups(
        self, biounit_ligand_chains: list[str]
    ) -> list[set[str]]:
        """Group ligand instance-chains linked by covalent bonds.

        Two ligand chains join the same group when a ``covale`` bond in
        ``self.covalent_bonds`` connects their asym ids within the same
        biounit instance. Covalently-linked ligand chains are the same
        physical molecule and must be reported as one ligand rather than
        several. Chains with no ligand-ligand covalent partner form
        singleton groups (the previous one-ligand-per-chain behaviour).
        """
        covale_edges: set[frozenset[str]] = set()
        for link1, link2 in self.covalent_bonds.get("covale", []):
            asym1, asym2 = link1.split(":")[2], link2.split(":")[2]
            if asym1 != asym2:
                covale_edges.add(frozenset((asym1, asym2)))

        parent = {chain: chain for chain in biounit_ligand_chains}

        def find(node: str) -> str:
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node

        def union(a: str, b: str) -> None:
            parent[find(a)] = find(b)

        # Only union chains that share an instance: biological-assembly
        # copies reuse asym labels but are physically independent.
        by_instance: dict[str, dict[str, str]] = defaultdict(dict)
        for chain in biounit_ligand_chains:
            instance, asym = chain.split(".", maxsplit=1)
            by_instance[instance][asym] = chain
        for asym_to_chain in by_instance.values():
            for edge in covale_edges:
                asym1, asym2 = tuple(edge)
                if asym1 in asym_to_chain and asym2 in asym_to_chain:
                    union(asym_to_chain[asym1], asym_to_chain[asym2])

        groups: dict[str, set[str]] = defaultdict(set)
        for chain in biounit_ligand_chains:
            groups[find(chain)].add(chain)
        return sorted(groups.values(), key=lambda g: sorted(g))

    @staticmethod
    def _ligand_pocket_members(ligand: Ligand) -> set[str]:
        """Return the exact members used to connect a non-artifact ligand."""
        members = {
            f"res:{chain}:{residue}"
            for chain, residues in ligand.neighboring_residues.items()
            for residue in residues
        }
        members.update(
            f"lig:{chain}"
            for chain in ligand.neighboring_ligands + ligand.interacting_ligands
        )
        return members

    def _upper_pocket_members(
        self,
        biounit: struc.AtomArray,
        spatial_index: BiounitSpatialIndex,
        instance_chain: str,
        *,
        neighboring_residue_threshold: float,
        interaction_search_threshold: float,
    ) -> set[str]:
        """Cheap superset of the pocket members for one deferred ligand.

        Every actual receptor neighbor is included at the normal pocket
        radius.  Every ligand-like chain within the interaction search radius
        is included as a possible ligand interaction.  The result may retain
        extra deferred ligands, but cannot discard an edge later created by
        :meth:`set_systems`.
        """
        chain_indices = spatial_index.atom_indices_for_chain(instance_chain)
        if chain_indices.size == 0:
            return set()
        coordinates = biounit.coord[chain_indices]
        receptor_asym_ids = {
            asym_id
            for asym_id, chain in self.chains.items()
            if asym_id not in self.ligand_like_chains
            and (
                _is_polypeptide(chain.chain_type_str)
                or _is_polynucleotide(chain.chain_type_str)
            )
        }
        members: set[str] = set()
        neighbor_indices = spatial_index.atom_indices_near(
            coordinates, neighboring_residue_threshold
        )
        for index in neighbor_indices:
            chain = str(biounit.chain_id[index])
            asym_id = chain.split(".", maxsplit=1)[-1]
            if asym_id in receptor_asym_ids:
                members.add(f"res:{chain}:{int(biounit.res_id[index])}")

        possible_interactions = spatial_index.atom_indices_near(
            coordinates, interaction_search_threshold
        )
        for chain in np.unique(biounit.chain_id[possible_interactions]):
            chain = str(chain)
            if chain == instance_chain or "." not in chain:
                continue
            if chain.split(".", maxsplit=1)[1] in self.ligand_like_chains:
                members.add(f"lig:{chain}")
        return members

    def _connected_deferred_ligand_chains(
        self,
        biounit: struc.AtomArray,
        spatial_index: BiounitSpatialIndex,
        primary_ligands: ty.Iterable[Ligand],
        deferred_instance_chains: set[str],
        *,
        min_shared_pocket_members: int,
        neighboring_residue_threshold: float,
        interaction_search_threshold: float,
    ) -> set[str]:
        """Keep deferred non-artifacts that may connect to a primary ligand."""
        import networkit as nk

        if not deferred_instance_chains:
            return set()
        if min_shared_pocket_members <= 0:
            return deferred_instance_chains

        primary_members = {
            ligand.instance_chain: self._ligand_pocket_members(ligand)
            for ligand in primary_ligands
            if not ligand.is_artifact
        }
        if not primary_members:
            return set()
        pocket_members = dict(primary_members)
        for chain in deferred_instance_chains:
            pocket_members[chain] = self._upper_pocket_members(
                biounit,
                spatial_index,
                chain,
                neighboring_residue_threshold=neighboring_residue_threshold,
                interaction_search_threshold=interaction_search_threshold,
            )

        chain_ids = list(pocket_members)
        graph = nk.Graph(len(chain_ids))
        member_chains: dict[str, list[int]] = defaultdict(list)
        for chain_index, members in enumerate(pocket_members.values()):
            for member in members:
                member_chains[member].append(chain_index)
        shared_counts: dict[tuple[int, int], int] = defaultdict(int)
        for member_indices in member_chains.values():
            for left, right in combinations(member_indices, 2):
                pair = (min(left, right), max(left, right))
                shared_counts[pair] += 1
                if shared_counts[pair] == min_shared_pocket_members:
                    graph.addEdge(*pair)

        primary_chains = set(primary_members)
        retained: set[str] = set()
        components = nk.components.ConnectedComponents(graph).run().getComponents()
        for component in components:
            component_chains = {chain_ids[index] for index in component}
            if component_chains & primary_chains:
                retained.update(component_chains & deferred_instance_chains)
        return retained

    def _record_biounit_ligand_contact_counts(
        self,
        biounit: struc.AtomArray,
        biounit_id: str,
        spatial_index: BiounitSpatialIndex,
        *,
        monoatomic_ion_asym_ids: set[str],
        known_artifact_asym_ids: set[str],
        contact_threshold: float,
    ) -> None:
        """Count every ligand-like chain contacting each receptor instance."""
        self._ligand_contacts_requested = True
        receptor_asym_ids = {
            asym_id
            for asym_id, chain in self.chains.items()
            if asym_id not in self.ligand_like_chains
            and (
                _is_polypeptide(chain.chain_type_str)
                or _is_polynucleotide(chain.chain_type_str)
            )
        }
        counts: dict[str, Counter[str]] = defaultdict(Counter)
        for ligand_instance in spatial_index.chain_ids:
            if "." not in ligand_instance:
                continue
            ligand_asym_id = ligand_instance.split(".", maxsplit=1)[1]
            if ligand_asym_id not in self.ligand_like_chains:
                continue
            ligand_indices = spatial_index.atom_indices_for_chain(ligand_instance)
            if ligand_indices.size == 0:
                continue
            nearby = spatial_index.atom_indices_near(
                biounit.coord[ligand_indices], contact_threshold
            )
            receptor_instances = {
                str(chain_instance)
                for chain_instance in np.unique(biounit.chain_id[nearby])
                if "." in str(chain_instance)
                and str(chain_instance).split(".", maxsplit=1)[1] in receptor_asym_ids
            }
            if ligand_asym_id in monoatomic_ion_asym_ids:
                contact_type = "ions"
            elif ligand_asym_id in known_artifact_asym_ids:
                contact_type = "artifacts"
            else:
                contact_type = "other_ligands"
            for receptor_instance in receptor_instances:
                counts[receptor_instance][contact_type] += 1
        self.biounit_ligand_contact_counts[str(biounit_id)] = {
            chain_instance: {
                "ions": int(chain_counts["ions"]),
                "artifacts": int(chain_counts["artifacts"]),
                "other_ligands": int(chain_counts["other_ligands"]),
            }
            for chain_instance, chain_counts in counts.items()
        }

    @classmethod
    def _from_cif_block(cls, cif_data: pdbx.CIFBlock, *, pdb_id: str) -> Entry:
        """Create an entry carrying the block's deposition metadata.

        Every extractor tolerates absent categories, so predicted or otherwise
        minimal custom mmCIFs simply leave the corresponding fields unset.
        """
        from plinder.data.annotations.cif_utils import _cif_scalar

        entry_info = get_entry_info(cif_data)
        entry_taxonomy = get_entry_taxonomy(cif_data)
        release_date = _cif_scalar(
            cif_data, "pdbx_audit_revision_history", "revision_date"
        )
        resolution = entry_info.get("entry_resolution")
        r = None
        if resolution is not None:
            try:
                r = float(resolution)
            except ValueError:
                r = None

        entry = cls(
            pdb_id=pdb_id,
            release_date=release_date or "",
            oligomeric_state=str(entry_info.get("entry_oligomeric_state"))
            if entry_info.get("entry_oligomeric_state") is not None
            else None,
            determination_method=str(entry_info.get("entry_determination_method"))
            if entry_info.get("entry_determination_method") is not None
            else None,
            keywords=str(entry_info.get("entry_keywords"))
            if entry_info.get("entry_keywords", None) is not None
            else None,
            pH=str(entry_info.get("entry_pH"))
            if entry_info.get("entry_pH") is not None
            else None,
            resolution=r,
            **entry_taxonomy,
        )
        entry.covalent_bonds = get_covalent_connections(cif_data)
        entry.chain_to_seqres = get_label_asym_sequences(cif_data)
        return entry

    def _attach_chains(
        self,
        atoms: struc.AtomArray,
        cif_data: pdbx.CIFBlock,
        *,
        min_polymer_size: int,
        data_dir: Path | None,
        include_ligands: bool,
        include_interfaces: bool,
    ) -> _LigandChainClasses:
        """Populate chains and mappings, then classify the ligand-like chains.

        Returns the ion / known-artifact / primary split that every biounit
        pass uses to defer cheap-to-miss ligands until a proper primary ligand
        connects to them.  Artifact detection needs ``data_dir``; without it
        no chain is treated as a known artifact.
        """
        self._populate_chains(atoms, cif_data)
        per_chain = get_chain_external_mappings(cif_data)
        for chain in per_chain:
            # External databases may annotate an asym ID that is absent from
            # the parsed model (for example, a chain omitted from model 1).
            # Such mappings cannot be attached to an Entry chain.
            if chain in self.chains:
                self.chains[chain].mappings = per_chain[chain]
        self.ligand_like_chains = detect_ligand_chains(self, min_polymer_size)
        if not include_ligands:
            return _LigandChainClasses(set(), set(), set())
        monoatomic_ion_mask = struc.filter_monoatomic_ions(atoms)
        ion_only_chains = set(
            str(chain) for chain in atoms.chain_id[monoatomic_ion_mask]
        )
        ion_only_chains.difference_update(
            str(chain) for chain in atoms.chain_id[~monoatomic_ion_mask]
        )
        monoatomic_ion_asym_ids = set(self.ligand_like_chains) & ion_only_chains
        known_artifact_asym_ids: set[str] = set()
        if data_dir is not None:
            artifact_codes = get_artifact_codes(data_dir)
            known_artifact_asym_ids = {
                asym_id
                for asym_id in self.ligand_like_chains
                if asym_id not in monoatomic_ion_asym_ids
                and is_known_artifact_ligand(
                    (
                        residue.name
                        for residue in self.chains[asym_id].residues.values()
                    ),
                    artifact_codes,
                )
            }
        primary_asym_ids = (
            set(self.ligand_like_chains)
            - monoatomic_ion_asym_ids
            - known_artifact_asym_ids
        )
        if not primary_asym_ids and include_interfaces:
            LOG.info(
                "PDB %r has no primary ligand chains; processing only "
                "protein interfaces",
                self.pdb_id,
            )
        elif not primary_asym_ids:
            LOG.info("PDB %r has no primary ligand chains", self.pdb_id)
        if monoatomic_ion_asym_ids or known_artifact_asym_ids:
            LOG.info(
                "PDB %s: deferring %d ion and %d known-artifact ASU ligand chains",
                self.pdb_id,
                len(monoatomic_ion_asym_ids),
                len(known_artifact_asym_ids),
            )
        return _LigandChainClasses(
            monoatomic_ion_asym_ids, known_artifact_asym_ids, primary_asym_ids
        )

    def _annotate_biounit(
        self,
        biounit: struc.AtomArray,
        biounit_id: str,
        *,
        ligand_classes: _LigandChainClasses,
        include_ligands: bool,
        include_interfaces: bool,
        data_dir: Path | None,
        interaction_search_threshold: float,
        neighboring_residue_threshold: float,
        neighboring_ligand_threshold: float,
        min_shared_pocket_members: int,
        interface_contact_radius: float,
        interface_min_chain_length: int,
        interface_min_residues: int,
        interface_annotate_prodigy: bool,
        ligand_smiles_dict: dict[str, str] | None = None,
        ligand_ccd_code_dict: dict[str, str] | None = None,
    ) -> dict[str, Ligand]:
        """Annotate one biological assembly and return its candidate ligands.

        The shared body of both ingest paths: PDB ingest calls it once per
        deposited assembly, custom ``as_is`` ingest once for the supplied
        coordinates.  It records the instance chain IDs, the per-receptor
        ligand contact counts, and the protein interfaces, then collects the
        primary ligands first and only those deferred ions and known artifacts
        that connect to a proper primary ligand.  The biounit must carry a
        ``legacy_chain_id`` annotation.
        """
        self.biounit_chain_ids[biounit_id] = sorted(
            str(chain_id) for chain_id in np.unique(biounit.chain_id)
        )
        self.biounit_legacy_chain_ids[biounit_id] = {
            str(chain_id): str(legacy_chain_id)
            for chain_id, legacy_chain_id in zip(
                biounit.chain_id, biounit.legacy_chain_id
            )
        }
        spatial_radii: list[float] = []
        if include_ligands:
            spatial_radii.extend(
                [
                    interaction_search_threshold,
                    neighboring_residue_threshold,
                    neighboring_ligand_threshold,
                ]
            )
        if include_interfaces:
            spatial_radii.append(interface_contact_radius)
        spatial_index = (
            BiounitSpatialIndex.from_atoms(biounit, max(spatial_radii))
            if spatial_radii
            else None
        )
        if include_ligands:
            assert spatial_index is not None
            self._record_biounit_ligand_contact_counts(
                biounit,
                biounit_id,
                spatial_index,
                monoatomic_ion_asym_ids=ligand_classes.monoatomic_ion_asym_ids,
                known_artifact_asym_ids=ligand_classes.known_artifact_asym_ids,
                contact_threshold=neighboring_residue_threshold,
            )
        if include_interfaces:
            self.interfaces.extend(
                detect_protein_interfaces(
                    biounit,
                    pdb_id=self.pdb_id,
                    biounit_id=biounit_id,
                    chains=self.chains,
                    contact_radius=interface_contact_radius,
                    min_chain_length=interface_min_chain_length,
                    min_interface_residues=interface_min_residues,
                    annotate_prodigy=interface_annotate_prodigy,
                    spatial_index=spatial_index,
                )
            )
        if not include_ligands or not ligand_classes.primary_asym_ids:
            return {}
        assert spatial_index is not None
        water_chains = get_water_chain_ids(biounit)

        def collect(**selection: set[str] | None) -> dict[str, Ligand]:
            return self._collect_ligands_from_biounit(
                biounit,
                biounit_id,
                interaction_search_threshold,
                neighboring_residue_threshold,
                neighboring_ligand_threshold,
                data_dir,
                ligand_smiles_dict=ligand_smiles_dict,
                ligand_ccd_code_dict=ligand_ccd_code_dict,
                water_chains=water_chains,
                spatial_index=spatial_index,
                **selection,
            )

        primary_ligands = collect(ligand_asym_ids=ligand_classes.primary_asym_ids)
        proper_primary_ligands = [
            ligand for ligand in primary_ligands.values() if ligand.is_proper
        ]
        if not proper_primary_ligands:
            return primary_ligands

        monoatomic_ion_asym_ids = ligand_classes.monoatomic_ion_asym_ids
        known_artifact_asym_ids = ligand_classes.known_artifact_asym_ids
        ion_instance_chains = {
            chain
            for chain in spatial_index.chain_ids
            if "." in chain
            and chain.split(".", maxsplit=1)[1] in monoatomic_ion_asym_ids
        }
        retained_ion_chains = self._connected_deferred_ligand_chains(
            biounit,
            spatial_index,
            proper_primary_ligands,
            ion_instance_chains,
            min_shared_pocket_members=min_shared_pocket_members,
            neighboring_residue_threshold=neighboring_residue_threshold,
            interaction_search_threshold=interaction_search_threshold,
        )
        ion_ligands = collect(ligand_instance_chains=retained_ion_chains)
        non_artifact_ligands = [
            ligand
            for ligand in (*primary_ligands.values(), *ion_ligands.values())
            if not ligand.is_artifact
        ]
        referenced_artifact_chains = {
            chain
            for ligand in non_artifact_ligands
            for chain in (ligand.neighboring_ligands + ligand.interacting_ligands)
            if "." in chain
            and chain.split(".", maxsplit=1)[1] in known_artifact_asym_ids
        }
        # Ligand interaction detection is not guaranteed to be symmetric
        # for short polymer ligands.  Include every known artifact inside
        # the interaction search radius as a conservative superset; the
        # exact set_systems() graph drops any false-positive candidates.
        for ligand in non_artifact_ligands:
            ligand_indices = spatial_index.atom_indices_for_chain(ligand.instance_chain)
            nearby_indices = spatial_index.atom_indices_near(
                biounit.coord[ligand_indices], interaction_search_threshold
            )
            referenced_artifact_chains.update(
                str(chain)
                for chain in np.unique(biounit.chain_id[nearby_indices])
                if "." in str(chain)
                and str(chain).split(".", maxsplit=1)[1] in known_artifact_asym_ids
            )
        artifact_ligands = collect(ligand_instance_chains=referenced_artifact_chains)
        return {**primary_ligands, **ion_ligands, **artifact_ligands}

    @classmethod
    def from_cif_file(
        cls,
        cif_file: Path,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        min_polymer_size: int = 12,
        data_dir: Path | None = None,
        save_folder: Path | None = None,
        interaction_search_threshold: float = 10.0,
        symmetry_mate_contact_threshold: float = 5.0,
        min_shared_pocket_members: int = 3,
        interface_contact_radius: float = 10.0,
        interface_min_chain_length: int = 12,
        interface_min_residues: int = DEFAULT_MIN_INTERFACE_RESIDUES,
        interface_annotate_prodigy: bool = True,
        include_ligands: bool = True,
        include_interfaces: bool = True,
        assembly_ids: ty.Iterable[str] | None = None,
        protein_only: bool = False,
    ) -> Entry:
        """
        Load an entry object from mmCIF files in the pipeline

        Parameters
        ----------
        cif_file : Path
            mmCIF file of interest
        neighboring_residue_threshold : float
            Distance from ligand for protein residues to be considered a ligand
        neighboring_ligand_threshold : float
            Distance from ligand for other ligands to be considered a ligand
        min_polymer_size : int = 12
            Minimum residue count for a polymer chain to be receptor.
            Shorter polymers are classified as ligands.  Set to 12 as
            the minimum length for meaningful MMseqs2/Foldseek searches.
        data_dir : Path | None
            Optional PLINDER data root used only for reference datasets such
            as CCD synonyms, cofactors, artifacts, and binding affinities.
            This is independent of ``save_folder``.
        save_folder : Path
            Root directory for one canonical ASU ligand SDF per chain.
            Files are written to ``<save_folder>/<pdb_id>/ligand_files``.
        interaction_search_threshold : float
            Receptor residues within this distance (Å) of a ligand form the
            complex handed to peppr for interaction detection. Keep it at or
            above peppr's 8 Å contact cutoff, which plinder leaves at its
            default; a smaller value silently truncates the contact search.
        min_shared_pocket_members : int
            Minimum shared pocket residues to group non-artifact ligands.
        interface_contact_radius : float
            Maximum backbone-atom distance defining a protein interface.
        interface_min_chain_length : int
            Minimum SEQRES length for an interface chain.
        interface_min_residues : int
            Minimum number of contacting residues required on each side.
        interface_annotate_prodigy : bool
            Whether to add PRODIGY-cryst contact and BIO/XTAL annotations.
        include_ligands : bool
            Whether to derive ligand systems and canonical ligand SDFs.
        include_interfaces : bool
            Whether to derive protein-protein interfaces and PRODIGY annotations.
        assembly_ids : Iterable[str] | None
            Optional subset of deposited biological assemblies. By default all
            assemblies listed by the mmCIF are processed.
        protein_only : bool
            Permit chain and assembly extraction without ligand or interface
            annotation. Used by receptor-only custom scoring.

        Returns
        -------
        Entry
            Entry object for the given pdbid
        """
        if not include_ligands and not include_interfaces and not protein_only:
            raise ValueError("entry ingest must include ligands, interfaces, or both")

        cif_file = Path(cif_file)
        data_dir = Path(data_dir) if data_dir is not None else None
        save_folder = Path(save_folder) if save_folder is not None else None

        from plinder.data.annotations.cif_utils import (
            _cif_scalar,
            read_mmcif_file,
        )

        cif_file_obj = read_mmcif_file(cif_file)
        cif_data = list(cif_file_obj.values())[0]
        # Extract metadata from CIF block
        pdb_id = (_cif_scalar(cif_data, "entry", "id") or "").lower()
        if include_ligands:
            cls._clear_ligand_files(save_folder, pdb_id)
        entry = cls._from_cif_block(cif_data, pdb_id=pdb_id)
        entry._ligand_contacts_requested = include_ligands
        if include_ligands and not include_interfaces and not protein_only:
            # Ligand-only ingest can skip loading the bonded structure when the
            # entry has no ligand-like chains at all; interface and protein-only
            # ingest still need the chains and assemblies.
            ligand_preflight = detect_ligand_chains_from_cif(
                cif_data,
                min_polymer_size,
            )
            if ligand_preflight == {}:
                LOG.info(
                    f"PDB {pdb_id!r} has no ligand-like chains; skipping bonded "
                    "structure loading"
                )
                return entry

        atoms = cls._load_clean_atoms(
            cif_file_obj,
            source=f"PDB {pdb_id!r}",
            no_bonds_error=(
                f"{pdb_id}: biotite returned no bonds despite include_bonds=True"
            ),
            require_bonds=include_ligands,
        )
        ligand_classes = entry._attach_chains(
            atoms,
            cif_data,
            min_polymer_size=min_polymer_size,
            data_dir=data_dir,
            include_ligands=include_ligands,
            include_interfaces=include_interfaces,
        )
        ligands: dict[str, Ligand] = {}
        selected_assemblies = _selected_assembly_ids(
            pdbx.list_assemblies(cif_file_obj), assembly_ids
        )
        for assembly_id in selected_assemblies:
            try:
                biounit = build_biounit(cif_file_obj, assembly_id)
            except Exception as e:
                # Skip this assembly but record it: its systems are silently
                # missing from the entry otherwise. Other assemblies proceed.
                LOG.error(
                    f"Could not build assembly {assembly_id} for "
                    f"{entry.pdb_id!r}: {e}"
                )
                entry.failed_assembly_ids.append(assembly_id)
                continue
            ligands.update(
                entry._annotate_biounit(
                    biounit,
                    str(assembly_id),
                    ligand_classes=ligand_classes,
                    include_ligands=include_ligands,
                    include_interfaces=include_interfaces,
                    data_dir=data_dir,
                    interaction_search_threshold=interaction_search_threshold,
                    neighboring_residue_threshold=neighboring_residue_threshold,
                    neighboring_ligand_threshold=neighboring_ligand_threshold,
                    min_shared_pocket_members=min_shared_pocket_members,
                    interface_contact_radius=interface_contact_radius,
                    interface_min_chain_length=interface_min_chain_length,
                    interface_min_residues=interface_min_residues,
                    interface_annotate_prodigy=interface_annotate_prodigy,
                )
            )
        entry._finalize(
            ligands,
            atoms,
            save_folder=save_folder,
            min_shared_pocket_members=min_shared_pocket_members,
            cif_file_obj=cif_file_obj,
            symmetry_mate_contact_threshold=symmetry_mate_contact_threshold,
            include_ligands=include_ligands,
        )
        return entry

    @classmethod
    def from_custom_cif_file(
        cls,
        pdb_id: str | None,
        cif_file: Path,
        ligand_smiles_dict: dict[str, str] | None = None,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        min_polymer_size: int = 12,
        interaction_search_threshold: float = 10.0,
        save_folder: Path | None = None,
        min_shared_pocket_members: int = 3,
        save_fixed_cif: Path | None = None,
        structure_mode: ty.Literal["as_is", "pdb"] = "as_is",
        assembly_ids: ty.Iterable[str] | None = None,
        symmetry_mate_contact_threshold: float = 5.0,
        interface_contact_radius: float = 10.0,
        interface_min_chain_length: int = 12,
        interface_min_residues: int = DEFAULT_MIN_INTERFACE_RESIDUES,
        interface_annotate_prodigy: bool = True,
        include_ligands: bool = True,
        include_interfaces: bool = True,
        data_dir: Path | None = None,
        ligand_ccd_code_dict: dict[str, str] | None = None,
    ) -> Entry:
        """
        Create an entry from an already assembled or deposited PDB mmCIF.

        Parameters
        ----------
        pdb_id : str | None
            Identifier used for an ``as_is`` structure. Must be ``None`` in
            ``pdb`` mode, where the deposited ``_entry.id`` is authoritative.
        cif_file : Path
            mmcif files of interest
        ligand_smiles_dict : dict[str, str] | None, optional
            Mapping of component ID (e.g. ``LIG``) to SMILES.
            Required for unknown ligands without ``_chem_comp_bond``
            (typical of cofolding outputs). Known CCD compounds
            are handled automatically.
        ligand_ccd_code_dict : dict[str, str] | None, optional
            Mapping of a custom component ID to the CCD component whose atom
            names, bonds, and canonical SMILES should be used, for example
            ``{"LIG": "ATP"}``. A component may be present in either this
            mapping or ``ligand_smiles_dict``, but not both.
        neighboring_residue_threshold : float, optional
            Max distance (Å) for neighboring receptor residues, by default 6.0
        neighboring_ligand_threshold : float, optional
            Max distance (Å) for neighboring ligands, by default 4.0
        min_polymer_size : int, optional
            Minimum residue count for a chain to be polymer (not ligand), by default 10
        save_folder : Path | None, optional
            Root directory for canonical ASU ligand SDFs, by default None.
            This path is never interpreted as a PLINDER data root.
        data_dir : Path | None, optional
            Optional PLINDER data root used for reference annotations. This is
            independent of ``save_folder`` and is not inferred from it.
        save_fixed_cif : Path | None, optional
            If provided and the CIF needed bond-order enrichment, write
            the enriched copy to this path. The input CIF at ``cif_file``
            is never mutated. Raises ``FileExistsError`` if the target
            already exists and ``ValueError`` if it resolves to the same
            path as ``cif_file``. By default (``None``) no file is written.
        structure_mode : {"as_is", "pdb"}
            ``as_is`` treats model 1 as one already assembled structure without
            symmetry expansion. ``pdb`` uses deposited biological-assembly
            operators through the production PDB ingest path.
        assembly_ids : Iterable[str] | None
            Biological assemblies selected in ``pdb`` mode. By default all are
            processed. Assembly selection is invalid in ``as_is`` mode.
        include_ligands : bool
            Whether to derive ligand systems and canonical ligand SDFs.
        include_interfaces : bool
            Whether to derive protein-protein interfaces.

        Returns
        -------
        Entry
            Entry object for the given pdbid

        Notes
        -----
        ``as_is`` mode runs the same entry, chain, and per-assembly annotation
        as :meth:`from_cif_file`, with the supplied coordinates standing in for
        the single biological assembly ``1``: deposition metadata, covalent
        links, external chain mappings, ion and known-artifact deferral (the
        latter needs ``data_dir``), ligand contact counts, interfaces, and
        ligand collection are identical.  Only the ligand-only preflight and
        the crystal-contact labelling of deposited entries are skipped.

        Raises
        ------
        MissingBondOrderError
            If the CIF contains unknown ligands and neither a SMILES nor CCD
            code is provided for them.
        FileExistsError
            If ``save_fixed_cif`` already exists.
        ValueError
            If ``save_fixed_cif`` points at the input ``cif_file``;
            if biotite returns no bonds at all (corrupted / missing
            bond information CIF); or any error propagated from
            :func:`~plinder.data.annotations.cif_utils.enrich_cif_with_smiles_bonds`
            (invalid SMILES, atom-count / element-order mismatch in the
            positional path, sanitize / template-match failure in the
            opt-in substructure path, or multi-instance comp_id
            divergence).
        """
        from plinder.data.annotations.cif_utils import (
            MissingBondOrderError,
            check_cif_bond_orders,
            check_custom_mmcif_fields,
            enrich_cif_with_ccd_bonds,
            enrich_cif_with_smiles_bonds,
            get_unknown_ligand_ids,
            read_mmcif_file,
        )

        cif_file = _require_mmcif_path(cif_file)
        data_dir = Path(data_dir) if data_dir is not None else None
        save_folder = Path(save_folder) if save_folder is not None else None
        ligand_ccd_code_dict = (
            {
                str(component_id).strip(): str(ccd_code).strip().upper()
                for component_id, ccd_code in ligand_ccd_code_dict.items()
            }
            if ligand_ccd_code_dict is not None
            else None
        )
        if structure_mode not in CUSTOM_STRUCTURE_MODES:
            raise ValueError(
                f"invalid structure_mode {structure_mode!r}; "
                f"expected one of {CUSTOM_STRUCTURE_MODES}"
            )
        if not include_ligands and (
            ligand_smiles_dict is not None or ligand_ccd_code_dict is not None
        ):
            raise ValueError("ligand chemistry overrides require include_ligands=True")
        overlapping_chemistry = set(ligand_smiles_dict or {}).intersection(
            ligand_ccd_code_dict or {}
        )
        if overlapping_chemistry:
            raise ValueError(
                "provide either SMILES or a CCD code for each ligand component, "
                f"not both: {sorted(overlapping_chemistry)}"
            )
        try:
            cif_file_obj = read_mmcif_file(cif_file)
            cif_data = list(cif_file_obj.values())[0]
        except (DeserializationError, IndexError, InvalidFileError, OSError) as exc:
            raise ValueError(f"cannot parse custom mmCIF {cif_file}: {exc}") from exc
        if structure_mode == "pdb" and pdb_id is not None:
            raise ValueError(
                "pdb_id must be None in pdb mode; the deposited _entry.id "
                "is authoritative"
            )
        check_custom_mmcif_fields(
            cif_data,
            source=cif_file,
            structure_mode=structure_mode,
            require_label_ids=True,
        )
        if structure_mode == "pdb":
            if (
                "entry" not in cif_data
                or "id" not in cif_data["entry"]
                or cif_data["entry"].row_count == 0
            ):
                raise ValueError(f"custom mmCIF {cif_file} needs _entry.id in pdb mode")
            if (
                ligand_smiles_dict is not None
                or ligand_ccd_code_dict is not None
                or save_fixed_cif is not None
            ):
                raise ValueError(
                    "pdb mode uses deposited PDB chemistry and does not accept "
                    "ligand chemistry overrides or save_fixed_cif"
                )
            requested_assemblies = (
                tuple(str(value) for value in assembly_ids)
                if assembly_ids is not None and not isinstance(assembly_ids, str)
                else assembly_ids
            )
            entry = cls.from_cif_file(
                cif_file,
                neighboring_residue_threshold=neighboring_residue_threshold,
                neighboring_ligand_threshold=neighboring_ligand_threshold,
                min_polymer_size=min_polymer_size,
                data_dir=data_dir,
                save_folder=save_folder,
                interaction_search_threshold=interaction_search_threshold,
                symmetry_mate_contact_threshold=symmetry_mate_contact_threshold,
                min_shared_pocket_members=min_shared_pocket_members,
                interface_contact_radius=interface_contact_radius,
                interface_min_chain_length=interface_min_chain_length,
                interface_min_residues=interface_min_residues,
                interface_annotate_prodigy=interface_annotate_prodigy,
                include_ligands=include_ligands,
                include_interfaces=include_interfaces,
                assembly_ids=requested_assemblies,
                protein_only=not include_ligands and not include_interfaces,
            )
            if not entry.biounit_chain_ids:
                raise ValueError(
                    f"pdb mode requires deposited biological assemblies: {cif_file}"
                )
            if requested_assemblies is not None:
                requested = (
                    {requested_assemblies}
                    if isinstance(requested_assemblies, str)
                    else set(requested_assemblies)
                )
                missing = sorted(requested.difference(entry.biounit_chain_ids))
                if missing:
                    raise ValueError(
                        f"failed to construct requested biological assemblies: {missing}"
                    )
            return entry
        if pdb_id is None or not str(pdb_id).strip():
            raise ValueError("pdb_id is required in as_is mode")
        if assembly_ids is not None:
            raise ValueError("assembly_ids are only valid in pdb mode")
        pdb_id = str(pdb_id).strip().lower()

        if include_ligands:
            cls._clear_ligand_files(save_folder, pdb_id)

        # Resolve explicit CCD references before checking which components still
        # need a user-provided SMILES. Both paths write _chem_comp_bond into the
        # in-memory copy without touching the caller's file.
        effective_smiles = dict(ligand_smiles_dict or {})
        enrichment_applied = False
        if include_ligands and ligand_ccd_code_dict:
            enrich_cif_with_ccd_bonds(
                cif_file_obj,
                ligand_ccd_codes=ligand_ccd_code_dict,
            )
            enrichment_applied = True

        unknown_ids = get_unknown_ligand_ids(cif_file_obj) if include_ligands else []
        if unknown_ids:
            missing_chemistry = set(unknown_ids).difference(effective_smiles)
            if missing_chemistry:
                raise MissingBondOrderError(
                    f"CIF contains unknown ligands {sorted(missing_chemistry)} with no "
                    "_chem_comp_bond and no CCD match. "
                    "Provide a SMILES in ligand_smiles_dict or a CCD code in "
                    "ligand_ccd_code_dict to assign bond orders."
                )
            enrich_cif_with_smiles_bonds(
                cif_file_obj,
                ligand_smiles={
                    comp_id: effective_smiles[comp_id] for comp_id in unknown_ids
                },
            )
            # Postcondition: enrichment must have resolved every flagged ligand.
            # Fail early/clearly if any unknown-without-bonds remain, rather than
            # feeding a partially-bonded CIF into structure parsing.
            check_cif_bond_orders(cif_file_obj)
            enrichment_applied = True
        ligand_smiles_dict = effective_smiles or None

        # Optionally persist the enriched CIF. Guard against overwriting
        # the caller's input or an existing file.
        if save_fixed_cif is not None and enrichment_applied:
            save_fixed_cif = Path(save_fixed_cif)
            if save_fixed_cif.resolve() == Path(cif_file).resolve():
                raise ValueError(
                    "save_fixed_cif must not point at the input cif_file — "
                    "the input file is never overwritten."
                )
            if save_fixed_cif.exists():
                raise FileExistsError(
                    f"save_fixed_cif target already exists: {save_fixed_cif}"
                )
            cif_file_obj.write(str(save_fixed_cif))

        atoms = cls._load_clean_atoms(
            cif_file_obj,
            source=f"Custom CIF {cif_file}",
            no_bonds_error=(
                f"Custom CIF {cif_file}: biotite returned no bonds at all "
                "after enrichment — the CIF is corrupted or missing all "
                "bond information (_chem_comp_bond, _struct_conn, and CCD "
                "coverage are all absent)."
            ),
            multimodel_note=(
                " Call from_custom_cif_file once per model for ensemble analysis."
            ),
            require_bonds=include_ligands,
        )
        entry = cls._from_cif_block(cif_data, pdb_id=pdb_id)
        entry._ligand_contacts_requested = include_ligands
        ligand_classes = entry._attach_chains(
            atoms,
            cif_data,
            min_polymer_size=min_polymer_size,
            data_dir=data_dir,
            include_ligands=include_ligands,
            include_interfaces=include_interfaces,
        )
        # Treat model 1 as one already assembled structure: a single biological
        # assembly "1" without symmetry expansion that keeps the supplied chain
        # IDs, then annotate it exactly like a deposited assembly.
        biounit = atoms.copy()
        biounit.chain_id = np.array([f"1.{c}" for c in biounit.chain_id])
        biounit.set_annotation("legacy_chain_id", biounit.chain_id.copy())
        ligands = entry._annotate_biounit(
            biounit,
            "1",
            ligand_classes=ligand_classes,
            include_ligands=include_ligands,
            include_interfaces=include_interfaces,
            data_dir=data_dir,
            interaction_search_threshold=interaction_search_threshold,
            neighboring_residue_threshold=neighboring_residue_threshold,
            neighboring_ligand_threshold=neighboring_ligand_threshold,
            min_shared_pocket_members=min_shared_pocket_members,
            interface_contact_radius=interface_contact_radius,
            interface_min_chain_length=interface_min_chain_length,
            interface_min_residues=interface_min_residues,
            interface_annotate_prodigy=interface_annotate_prodigy,
            ligand_smiles_dict=ligand_smiles_dict,
            ligand_ccd_code_dict=ligand_ccd_code_dict,
        )
        # Crystal contacts are skipped: as_is coordinates carry no
        # crystallographic symmetry to expand.
        entry._finalize(
            ligands,
            atoms,
            save_folder=save_folder,
            min_shared_pocket_members=min_shared_pocket_members,
            include_ligands=include_ligands,
        )
        return entry

    def set_systems(
        self,
        ligands: dict[str, Ligand],
        min_shared_pocket_members: int = 3,
    ) -> None:
        """Group ligands into systems by shared pocket and proximity.

        Non-artifact ligands (drug-like, cofactors, ions) are grouped
        if they share at least *min_shared_pocket_members* pocket
        members (receptor residues + neighboring ligand chains).
        Pocket members use chain instance IDs (e.g. ``1.A``) so
        ligands in different subunits only merge when they genuinely
        share residues on the same chain copy. Biological assemblies
        are always independent: their instance-chain labels share a
        namespace but do not describe the same physical chains.

        Artifacts (GOL, PEG, etc.) are only attached to a system if
        they are within 4 Å of a non-artifact ligand.

        Parameters
        ----------
        ligands : dict[str, Ligand]
            All ligands in the entry keyed by ligand ID.
        min_shared_pocket_members : int
            Minimum shared pocket members to group non-artifact ligands.
        """
        import networkit as nk

        ligand_ids = list(ligands.keys())
        G = nk.Graph(len(ligand_ids))

        # Step 1: group non-artifact ligands by shared pocket residues
        pocket_members: dict[int, set[str]] = {}
        for i, lid in enumerate(ligand_ids):
            lig = ligands[lid]
            if lig.is_artifact:
                continue
            members: set[str] = set()
            for chain, resnums in lig.neighboring_residues.items():
                for rn in resnums:
                    members.add(f"res:{chain}:{rn}")
            for lc in lig.neighboring_ligands + lig.interacting_ligands:
                members.add(f"lig:{lc}")
            pocket_members[i] = members

        groupable = list(pocket_members)
        if min_shared_pocket_members <= 0:
            for i, j in combinations(groupable, 2):
                if (
                    ligands[ligand_ids[i]].biounit_id
                    == ligands[ligand_ids[j]].biounit_id
                ):
                    G.addEdge(i, j)
        else:
            member_ligands: dict[str, list[int]] = defaultdict(list)
            for ligand_index, members in pocket_members.items():
                biounit_id = ligands[ligand_ids[ligand_index]].biounit_id
                for member in members:
                    member_ligands[f"{biounit_id}:{member}"].append(ligand_index)
            shared_member_counts: dict[tuple[int, int], int] = defaultdict(int)
            for member_indices in member_ligands.values():
                for i, j in combinations(member_indices, 2):
                    pair = (min(i, j), max(i, j))
                    shared_member_counts[pair] += 1
                    if shared_member_counts[pair] == min_shared_pocket_members:
                        G.addEdge(*pair)

        # Step 2: attach artifacts within 4A of a non-artifact ligand
        ligand_id_to_index = {ligand_id: i for i, ligand_id in enumerate(ligand_ids)}
        for i, lid in enumerate(ligand_ids):
            lig = ligands[lid]
            if not lig.is_artifact:
                continue
            for neighbor_chain in lig.neighboring_ligands + lig.interacting_ligands:
                neighbor_id = "__".join([self.pdb_id, lig.biounit_id, neighbor_chain])
                j_idx = ligand_id_to_index.get(neighbor_id)
                if j_idx is not None and not ligands[ligand_ids[j_idx]].is_artifact:
                    G.addEdge(i, j_idx)
        cc = nk.components.ConnectedComponents(G)
        cc.run()
        components = cc.getComponents()
        system_ligands: dict[int, list[Ligand]] = {}
        for idx, component in enumerate(sorted(components, key=len, reverse=True)):
            system_ligands[idx + 1] = []
            for node_idx in component:
                system_ligands[idx + 1].append(ligands[ligand_ids[node_idx]])
        self.systems: dict[str, System] = {}
        for ligs in system_ligands.values():
            if not ligs:
                continue
            biounit_ids = {ligand.biounit_id for ligand in ligs}
            if len(biounit_ids) != 1:
                raise RuntimeError(
                    f"system ligands span biological assemblies: {sorted(biounit_ids)}"
                )
            instance_chains = [ligand.instance_chain for ligand in ligs]
            if len(instance_chains) != len(set(instance_chains)):
                raise RuntimeError(
                    "system contains repeated ligand instance chains: "
                    f"{instance_chains}"
                )
            receptor_asym_ids = sorted(
                {
                    instance_chain.split(".", maxsplit=1)[1]
                    for ligand in ligs
                    for instance_chain in ligand.protein_chains_asym_id
                }
            )
            system = System(
                pdb_id=self.pdb_id,
                biounit_id=next(iter(biounit_ids)),
                ligands=sorted(ligs, key=lambda x: x.id),
                receptor_type=get_receptor_type(
                    self.chains[asym_id].chain_type_str for asym_id in receptor_asym_ids
                ),
            )
            if system.proper_ligands() and len(system.protein_chains_asym_id):
                self.systems[system.id] = system

    @cached_property
    def author_to_asym(self) -> dict[str, str]:
        """[EXCLUDE] Map author chain id to asym id"""
        return {
            c.auth_id: c.asym_id
            for c in self.chains.values()
            if "polypeptide" in c.chain_type_str.lower()
        }

    def chains_for_alignment(self, chain_type: str, aln_type: str) -> list[str]:
        """
        Get chains for foldseek/mmseqs alignment
        Parameters
        ----------
        self : Entry
            Entry object
        chain_type : str
            Chain type (apo/holo)
        aln_type : str
            Alignment type (folseek, mmses, etc)
        Returns
        -------
        list[str]
        """
        assert chain_type in (
            "apo",
            "holo",
            "pred",
        ), "chain_type must be 'apo', 'holo', or 'pred'"
        if chain_type == "holo":
            receptor_asym_ids = {
                i_c.split(".", maxsplit=1)[1]
                for system in self.systems.values()
                if system.system_type == "holo"
                for i_c in system.protein_chains_asym_id
            }
            na_chains = sorted(
                asym
                for asym in receptor_asym_ids
                if _is_polynucleotide(self.chains[asym].chain_type_str)
            )
            if na_chains:
                LOG.warning(
                    f"PDB {self.pdb_id!r}: nucleic acid receptor chains "
                    f"{na_chains} are excluded from {aln_type} alignment "
                    "because similarity databases and scores are protein-only."
                )
            chains = {
                self.chains[asym].auth_id
                for asym in receptor_asym_ids
                if _is_polypeptide(self.chains[asym].chain_type_str)
            }
        elif chain_type == "apo":
            holo_entities = set(
                self.chains[c].entity_id for c in self.chains if self.chains[c].holo
            )
            chains = set(
                self.chains[c].auth_id
                for c in self.chains
                if not self.chains[c].holo
                and self.chains[c].entity_id not in holo_entities
                and _is_polypeptide(self.chains[c].chain_type_str)
            )
        elif chain_type == "pred":
            chains = set()
            for c in self.chains:
                if self.chains[c].holo and _is_polypeptide(
                    self.chains[c].chain_type_str
                ):
                    chains |= set(self.chains[c].mappings.get("UniProt", {}))
        if aln_type == "foldseek":
            if chain_type == "pred":
                return [f"AF-{c}-F1-model_v4_A" for c in chains]
            else:
                return [f"pdb_0000{self.pdb_id}_xyz-enrich_{c}" for c in chains]
        elif aln_type == "mmseqs":
            if chain_type == "pred":
                return list(chains)
            else:
                return [f"{self.pdb_id}_{c}" for c in chains]
        return []

    def label_chains(self) -> None:
        """
        Label chains as apo/holo/ligand
        Parameters
        ----------
        self : Entry
            Entry object

        Returns
        -------
        None
        """
        holo_chains = set()
        for system in self.systems.values():
            if system.system_type == "holo":
                holo_chains.update(
                    [c.split(".", maxsplit=1)[1] for c in system.protein_chains_asym_id]
                )
        for interface in self.interfaces:
            holo_chains.update(
                {
                    interface.chain_1.split(".", maxsplit=1)[-1],
                    interface.chain_2.split(".", maxsplit=1)[-1],
                }
            )
        for chain in self.chains:
            self.chains[chain].holo = chain in holo_chains

    def format_validation(
        self, criteria: QualityCriteria = QualityCriteria()
    ) -> dict[str, ty.Any]:
        assert self.validation is not None
        data = self.validation.format()
        data = {f"entry_{k}": v for k, v in data.items()}
        if data["entry_validation_r"] is None:
            self.pass_criteria = False
        else:
            quality = [
                # ENTRY
                data["entry_validation_resolution"] is not None
                and data["entry_validation_resolution"]
                <= criteria.max_entry_resolution,
                data["entry_validation_r"] is not None
                and data["entry_validation_r"] <= criteria.max_entry_r,
                data["entry_validation_rfree"] is not None
                and data["entry_validation_rfree"] <= criteria.max_entry_rfree,
                data["entry_validation_r_minus_rfree"] is not None
                and data["entry_validation_r_minus_rfree"]
                <= criteria.max_entry_r_minus_rfree,
            ]
            self.pass_criteria = all(quality)
        data["entry_pass_validation_criteria"] = self.pass_criteria
        return data

    def format(
        self, criteria: QualityCriteria = QualityCriteria()
    ) -> dict[str, ty.Any]:
        """
        Format label for entry-level annotations by prepending label with "entry_"

        Parameters
        ----------
        self : Entry
            Entry object

        Returns
        -------
        dict[str, ty.Any]
        """
        data: dict[str, ty.Any] = defaultdict(str)
        columns = [
            "pdb_id",
            "release_date",
            "oligomeric_state",
            "determination_method",
            "keywords",
            "pH",
            "resolution",
            "source_taxonomy_ids",
            "source_organism_names",
            "host_taxonomy_ids",
            "host_organism_names",
        ]
        for field in columns:
            name = f"entry_{field}"
            data[name] = getattr(self, field, None)

        if self.validation:
            data.update(self.format_validation(criteria))
            data["entry_pass_validation_criteria"] = self.pass_criteria
        return data

    def chains_to_df(self) -> pd.DataFrame:
        """Return one metadata row for each receptor polymer chain."""
        columns = [
            "entry_pdb_id",
            "chain_asym_id",
            "chain_auth_id",
            "chain_entity_id",
            "chain_type",
            "chain_receptor_type",
            "chain_sequence",
            "chain_length",
            "chain_num_unresolved_residues",
            "chain_is_holo",
            "chain_is_ligand_like",
            "chain_uniprot_ids",
        ]
        rows = []
        for chain_id in sorted(self.chains):
            chain = self.chains[chain_id]
            if not (
                _is_polypeptide(chain.chain_type_str)
                or _is_polynucleotide(chain.chain_type_str)
            ):
                continue
            rows.append(
                {
                    "entry_pdb_id": self.pdb_id,
                    "chain_asym_id": chain.asym_id,
                    "chain_auth_id": chain.auth_id,
                    "chain_entity_id": chain.entity_id,
                    "chain_type": chain.chain_type_str,
                    "chain_receptor_type": get_receptor_type([chain.chain_type_str]),
                    "chain_sequence": self.chain_to_seqres.get(chain.asym_id, ""),
                    "chain_length": chain.length,
                    "chain_num_unresolved_residues": chain.num_unresolved_residues,
                    "chain_is_holo": chain.holo,
                    "chain_is_ligand_like": chain_id in self.ligand_like_chains,
                    "chain_uniprot_ids": sorted(chain.mappings.get("UniProt", {})),
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def metadata_to_df(self) -> pd.DataFrame:
        """Return one row of entry-level annotations."""
        return pd.DataFrame([self.format()])

    def biounit_chains_to_df(self) -> pd.DataFrame:
        """Return biological-assembly membership once per chain instance."""
        columns = [
            "entry_pdb_id",
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
        ]
        contact_columns = [
            "chain_num_contacting_ions",
            "chain_num_contacting_artifacts",
            "chain_num_contacting_other_ligands",
        ]
        contacts_computed = self._ligand_contacts_requested or bool(
            self.biounit_ligand_contact_counts
        )
        if contacts_computed:
            columns.extend(contact_columns)
        rows = []
        water_chains = set(self.water_chains)
        ligand_chains = set(self.ligand_like_chains)
        for biounit_id, chain_instances in sorted(self.biounit_chain_ids.items()):
            contact_counts = self.biounit_ligand_contact_counts.get(str(biounit_id), {})
            for chain_instance in sorted(set(chain_instances)):
                asym_id = chain_instance.split(".", maxsplit=1)[-1]
                if asym_id in water_chains:
                    role = "water"
                elif asym_id in ligand_chains:
                    role = "ligand"
                else:
                    role = "receptor"
                chain_counts = contact_counts.get(chain_instance, {})
                row: dict[str, object] = {
                    "entry_pdb_id": self.pdb_id,
                    "biounit_id": str(biounit_id),
                    "chain_instance": chain_instance,
                    "chain_asym_id": asym_id,
                    "chain_role": role,
                }
                if contacts_computed:
                    row.update(
                        {
                            "chain_num_contacting_ions": int(
                                chain_counts.get("ions", 0)
                            ),
                            "chain_num_contacting_artifacts": int(
                                chain_counts.get("artifacts", 0)
                            ),
                            "chain_num_contacting_other_ligands": int(
                                chain_counts.get("other_ligands", 0)
                            ),
                        }
                    )
                rows.append(row)
        return pd.DataFrame(rows, columns=columns)

    def to_df(self) -> pd.DataFrame:
        """
        Convert entry data object to pd.DataFrame
        Parameters
        ----------
        self : Entry
            Entry object

        Returns
        -------
        pd.DataFrame
        """
        if self.validation is not None:
            self.format_validation()
        rows = []
        entry_data = {"entry_pdb_id": self.pdb_id}
        for system in self.systems:
            annotation = self.systems[system]
            legacy_mapping = self.biounit_legacy_chain_ids.get(
                annotation.biounit_id, {}
            )
            legacy_protein_chains = sorted(
                legacy_mapping.get(chain_id, chain_id)
                for chain_id in annotation.protein_chains_asym_id
            )
            legacy_ligand_chains = sorted(
                legacy_mapping.get(chain_id, chain_id)
                for chain_id in annotation.ligand_chains
            )
            annotation.id_legacy = "__".join(
                [
                    annotation.pdb_id,
                    annotation.biounit_id,
                    "_".join(legacy_protein_chains),
                    "_".join(legacy_ligand_chains),
                ]
            )
            system_data = annotation.format(
                self.chains,
                self.pass_criteria,
            )
            for ligand in self.systems[system].ligands:
                legacy_instance_chain = legacy_mapping.get(
                    ligand.instance_chain, ligand.instance_chain
                )
                ligand.id_legacy = "__".join(
                    [ligand.pdb_id, ligand.biounit_id, legacy_instance_chain]
                )
                ligand_data = ligand.format(self.chains)
                rows.append({**entry_data, **system_data, **ligand_data})
        return pd.DataFrame(rows)

    def clear_non_pocket_residues(self) -> None:
        """
        Remove non-pocket residues from chains
        """
        all_pocket_residues: dict[str, set[int]] = defaultdict(set)
        for system in self.systems.values():
            for ligand in system.ligands:
                for chain in ligand.pocket_residues:
                    all_pocket_residues[chain.split(".", maxsplit=1)[1]].update(
                        ligand.pocket_residues[chain].keys()
                    )
        n_before = sum(len(c.residues) for c in self.chains.values())
        for chain in self.chains:
            if chain in all_pocket_residues:
                self.chains[chain].residues = {
                    r: self.chains[chain].residues[r]
                    for r in all_pocket_residues[chain]
                }
            else:
                self.chains[chain].residues = {}
        n_dropped = n_before - sum(len(c.residues) for c in self.chains.values())
        if n_dropped:
            LOG.info(
                f"{self.pdb_id}: clear_non_pocket_residues dropped {n_dropped} of "
                f"{n_before} chain residues (kept pocket residues only)"
            )

    def set_validation(
        self,
        validation_file: Path,
        cif_file: Path,
        thresholds: ResidueValidationThresholds = ResidueValidationThresholds(),
    ) -> None:
        from PDBValidation.ValidationFactory import ValidationFactory

        if self.determination_method != "X-RAY DIFFRACTION":
            LOG.warning(
                f"set_validation: Skipping validation for {self.pdb_id} as method is not X-RAY DIFFRACTION"
            )
            return
        if not validation_file.exists():
            LOG.error(f"set_validation: Validation file not found {validation_file}")
            return
        try:
            doc = ValidationFactory(
                str(validation_file), mmcif_path=str(cif_file)
            ).getValidation()
            self.validation = EntryValidation.from_entry(doc)
            if self.validation and self.validation.r is not None:
                system_chain_ids = {
                    instance_chain.split(".", maxsplit=1)[-1]
                    for system in self.systems.values()
                    for instance_chain in system.protein_chains_asym_id
                }
                system_chain_ids.update(
                    ligand.asym_id
                    for system in self.systems.values()
                    for ligand in system.ligands
                )
                # TODO: consider per-chain / per-system try/except here. This
                # entry-level catch is coarse: one bad chain or system abandons
                # validation for the whole entry, leaving partial state. Per-
                # residue validation already degrades gracefully.
                for chain in system_chain_ids:
                    self.chains[chain].set_validation(doc, thresholds)
                for system in self.systems:
                    self.systems[system].set_validation(self.chains, thresholds)
        except Exception as e:
            LOG.error(
                f"set_validation: Error setting validation for {self.pdb_id}: {e}"
            )


def document(output_dir: Path) -> None:
    """
    Document the columns in the annotation files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    Entry.document_properties_to_tsv(prefix="entry", filename=output_dir / "entry.tsv")
    EntryValidation.document_properties_to_tsv(
        prefix="entry_validation", filename=output_dir / "entry_validation.tsv"
    )
    System.document_properties_to_tsv(
        prefix="system", filename=output_dir / "system.tsv"
    )
    ResidueListValidation.document_properties_to_tsv(
        prefix="system_pocket", filename=output_dir / "system_pocket_validation.tsv"
    )
    ResidueListValidation.document_properties_to_tsv(
        prefix="system_ligand", filename=output_dir / "system_ligand_validation.tsv"
    )
    Chain.document_properties_to_tsv(
        prefix="system_protein_chains",
        filename=output_dir / "system_protein_chains.tsv",
    )
    Chain.document_properties_to_tsv(
        prefix="system_ligand_chains", filename=output_dir / "system_ligand_chains.tsv"
    )
    Ligand.document_properties_to_tsv(
        prefix="ligand", filename=output_dir / "ligands.tsv"
    )
