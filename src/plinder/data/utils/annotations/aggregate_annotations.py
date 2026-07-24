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
import networkit as nk
import numpy as np
import pandas as pd
from PDBValidation.ValidationFactory import ValidationFactory
from pydantic import BeforeValidator, Field
from rdkit import RDLogger

from plinder.core.structure.atoms import is_hydrogen_isotope
from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.cif_utils import (
    apply_struct_conn_bonds,
    build_biounit,
    get_chain_external_mappings,
    get_entry_info,
    get_label_asym_sequences,
    get_model_count,
    get_structure_with_altloc,
)
from plinder.data.utils.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
    ResidueValidationThresholds,
)
from plinder.data.utils.annotations.interaction_utils import (
    get_covalent_connections,
    get_symmetry_mate_contacts,
)
from plinder.data.utils.annotations.ligand_utils import (
    BiounitSpatialIndex,
    Ligand,
    get_artifact_codes,
    get_water_chain_ids,
    is_known_artifact_ligand,
    validate_chain_residue,
)
from plinder.data.utils.annotations.protein_utils import (
    Chain,
    _is_polynucleotide,
    _is_polypeptide,
    detect_ligand_chains,
    detect_ligand_chains_from_cif,
    get_receptor_type,
)
from plinder.data.utils.annotations.save_utils import save_ligands
from plinder.data.utils.annotations.utils import DocBaseModel

LOG = setup_logger(__name__)
RDLogger.DisableLog("rdApp.*")
SymmetryMateContacts = ty.Annotated[
    dict[tuple[str, int], dict[tuple[str, int], dict[int, set[int]]]],
    BeforeValidator(validate_chain_residue),
    Field(default_factory=dict),
]


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
    pdb_id: str = Field(description="__PDB ID")
    biounit_id: str = Field(description="Biounit ID")
    ligands: list[Ligand] = Field(description="__List of Ligands in a systems")
    receptor_type: str = Field(
        description=(
            "Receptor polymer composition: protein, DNA, RNA, other, or a "
            "+-joined combination"
        )
    )
    ligand_validation: ResidueListValidation | None = Field(
        default=None,
        description="__Validation object for the ligand residues in the system",
    )
    pocket_validation: ResidueListValidation | None = Field(
        default=None, description="__Validation object for the system's pocket residues"
    )
    pass_criteria: bool | None = Field(
        default=None, description="__Passes quality criteria"
    )  # TODO: remove as attribute and have as function

    """
    This class defines a system which includes a protein-ligand complex
    and its neighboring ligands and receptor residues

    """

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
                "_".join(x.split(".")[1] for x in self.protein_chains_asym_id),
                "_".join(x.split(".")[1] for x in self.ligand_chains),
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
        Number of interactions of the system
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
        Number of unique interactions of the system
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
        Number of covalent ligands of the system
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
        """
        __Pockets residues of the system
        """
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
        """
        __Interactions of the system
        """
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
        """
        __Counter of interactions of the system
        """
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
        for field, desc_type in self.get_descriptions_and_types().items():
            # blacklist fields that will be added with custom formatters below or that we don't want to add to the plindex
            descr = str(desc_type[0]).lstrip().replace("\n", " ")
            if descr.startswith("__"):
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
        """
        __Waters interacting with any of the ligands in the system
        """
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
                chains[c.split(".")[1]].residues[r].validation  # type: ignore
                for c in self.ligand_chains
                for r in chains[c.split(".")[1]].residues
            ],
            thresholds,
        )
        self.pocket_validation = ResidueListValidation.from_residues(
            [
                chains[c.split(".")[1]].residues[r].validation  # type: ignore
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
            neighboring_chain = neighboring_chain.split(".")[-1]
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
    chains: dict[str, Chain] = Field(
        default_factory=dict,
        description="__Chains dictionary with chain name mapped to chain object",
    )
    ligand_like_chains: dict[str, str] = Field(
        default_factory=dict,
        description="__Chain: chain type for other ligand-like chains in the entry",
    )
    systems: dict[str, System] = Field(
        default_factory=dict,
        description="__System dictionary with system id mapped to system object",
    )
    covalent_bonds: dict[str, list[tuple[str, str]]] = Field(
        default_factory=dict,
        description="__All covalent interactions in the entry as defined by mmcif annotations. They types are separated by dictionary key and they include: "
        + "covale: actual covalent linkage, metalc: other dative bond interactions like metal-ligand dative bond, "
        + "hydrogc: strong hydorogen bonding of nucleic acid. For the purpose of covalent annotations, we use only covale for downstream processing.",
    )
    chain_to_seqres: dict[str, str] = Field(
        default_factory=dict, description="__Chain to sequence mapping"
    )
    validation: EntryValidation | None = Field(
        default=None, description="__Entry validation"
    )
    pass_criteria: bool | None = Field(
        default=None, description="__Entry pass validation criteria"
    )
    water_chains: list[str] = Field(
        default_factory=list, description="__Water chains in the entry"
    )
    biounit_chain_ids: dict[str, list[str]] = Field(
        default_factory=dict,
        description="__Resolved biological-assembly chain instances by assembly ID",
    )
    symmetry_mate_contacts: SymmetryMateContacts = Field(
        default_factory=dict, description="__Symmetry mate contacts in the entry"
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
            self.systems = {
                s.id: s
                for s in self.systems.values()
                if s.system_type == "holo"
                and len(s.protein_chains_asym_id) <= max_protein_chains
                and len(s.ligand_chains) <= max_ligand_chains
            }
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
                self.chains[chain_id] = Chain.from_cif_data(
                    chain_id,
                    block,
                    chain_atoms,
                    len(self.chain_to_seqres.get(chain_id, "")),
                    entity_id=entity_id,
                    chain_type_str=type_by_entity.get(entity_id, "unknown"),
                )
        finally:
            atoms.bonds = bonds

        self.water_chains = sorted(water_chains)

    def _finalize(
        self,
        ligands: dict[str, Ligand],
        min_shared_pocket_members: int = 3,
    ) -> None:
        """Label crystal contacts and set systems."""
        if self.symmetry_mate_contacts:
            for ligand in ligands.values():
                ligand.label_crystal_contacts(self.symmetry_mate_contacts)
        self.set_systems(ligands, min_shared_pocket_members=min_shared_pocket_members)
        self.label_chains()

    def _collect_ligands_from_biounit(
        self,
        biounit: struc.AtomArray,
        biounit_id: str,
        plip_complex_threshold: float,
        neighboring_residue_threshold: float,
        neighboring_ligand_threshold: float,
        data_dir: Path | None,
        ligand_smiles_dict: dict[str, str] | None = None,
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
                    plip_complex_threshold,
                    neighboring_residue_threshold,
                    neighboring_ligand_threshold,
                ),
            )
        # Find ligand chains: chain_id format is "{instance}.{asym_id}"
        biounit_ligand_chains = [
            c
            for c in spatial_index.chain_ids
            if "." in c
            and c.split(".")[1] in self.ligand_like_chains
            and (ligand_asym_ids is None or c.split(".")[1] in ligand_asym_ids)
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
            primary_instance, primary_asym_id = primary_chain.split(".")
            ligand = Ligand.from_pli(
                pdb_id=self.pdb_id,
                biounit_id=biounit_id,
                biounit=biounit,
                ligand_instance=int(primary_instance),
                ligand_chain=self.chains[primary_asym_id],
                residue_numbers=member_residue_numbers[primary_chain],
                ligand_like_chains=self.ligand_like_chains,
                all_covalent_dict=self.covalent_bonds,
                plip_complex_threshold=plip_complex_threshold,
                neighboring_residue_threshold=neighboring_residue_threshold,
                neighboring_ligand_threshold=neighboring_ligand_threshold,
                data_dir=data_dir,
                chain_to_seqres=self.chain_to_seqres,
                ligand_smiles_dict=ligand_smiles_dict,
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
            instance, asym = chain.split(".")
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

    @classmethod
    def from_cif_file(
        cls,
        cif_file: Path,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        min_polymer_size: int = 12,
        data_dir: Path | None = None,
        save_folder: Path | None = None,
        plip_complex_threshold: float = 10.0,
        symmetry_mate_contact_threshold: float = 5.0,
        min_shared_pocket_members: int = 3,
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
        save_folder : Path
            Root directory for one canonical ASU ligand SDF per chain.
            Files are written to ``<save_folder>/<pdb_id>/ligand_files``.
        plip_complex_threshold : float
            Maximum distance (Å) from ligand for interaction analysis
        min_shared_pocket_members : int
            Minimum shared pocket residues to group non-artifact ligands.

        Returns
        -------
        Entry
            Entry object for the given pdbid
        """
        from plinder.data.utils.annotations.cif_utils import (
            _cif_scalar,
            read_mmcif_file,
        )

        cif_file_obj = read_mmcif_file(cif_file)
        cif_data = list(cif_file_obj.values())[0]
        entry_info = get_entry_info(cif_data)

        # Extract metadata from CIF block
        pdb_id = (_cif_scalar(cif_data, "entry", "id") or "").lower()
        if save_folder is not None:
            # Re-ingest must never merge newly retained ligands with SDFs from
            # a previous annotation of the same entry.  Clear this derived
            # directory before any early no-system return as well.
            ligand_dir = Path(save_folder) / pdb_id / "ligand_files"
            if ligand_dir.exists():
                shutil.rmtree(ligand_dir)
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
        )
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

        # Load structure with biotite
        # Multi-model PDBs (e.g. NMR ensembles) silently use model 1 here;
        # warn so callers know other models are dropped.
        n_models = get_model_count(cif_file_obj)
        if n_models > 1:
            LOG.warning(f"PDB {pdb_id!r} has {n_models} models — using model 1 only.")
        atoms = get_structure_with_altloc(
            cif_file_obj, model=1, use_author_fields=False, include_bonds=True
        )
        atoms = atoms[~is_hydrogen_isotope(atoms.element)]
        if atoms.bonds is None:
            raise ValueError(
                f"{pdb_id}: biotite returned no bonds despite include_bonds=True"
            )
        apply_struct_conn_bonds(atoms, cif_data)
        chain_to_seqres = get_label_asym_sequences(cif_data)

        entry.covalent_bonds = get_covalent_connections(cif_data)
        entry.chain_to_seqres = chain_to_seqres
        entry._populate_chains(atoms, cif_data)

        if save_folder is not None and data_dir is None:
            data_dir = save_folder.parent.parent
        per_chain = get_chain_external_mappings(cif_data)
        for chain in per_chain:
            # External databases may annotate an asym ID that is absent from
            # the parsed model (for example, a chain omitted from model 1).
            # Such mappings cannot be attached to an Entry chain.
            if chain in entry.chains:
                entry.chains[chain].mappings = per_chain[chain]
        entry.ligand_like_chains = detect_ligand_chains(entry, min_polymer_size)
        if not entry.ligand_like_chains:
            # There can be no systems without ligand-like chains.  In
            # particular, avoid building biological assemblies for large
            # receptor-only entries that cannot contribute annotation rows.
            return entry
        monoatomic_ion_mask = struc.filter_monoatomic_ions(atoms)
        ion_only_chains = set(
            str(chain) for chain in atoms.chain_id[monoatomic_ion_mask]
        )
        ion_only_chains.difference_update(
            str(chain) for chain in atoms.chain_id[~monoatomic_ion_mask]
        )
        monoatomic_ion_asym_ids = set(entry.ligand_like_chains) & ion_only_chains
        known_artifact_asym_ids: set[str] = set()
        if data_dir is not None:
            artifact_codes = get_artifact_codes(data_dir)
            known_artifact_asym_ids = {
                asym_id
                for asym_id in entry.ligand_like_chains
                if asym_id not in monoatomic_ion_asym_ids
                and is_known_artifact_ligand(
                    (
                        residue.name
                        for residue in entry.chains[asym_id].residues.values()
                    ),
                    artifact_codes,
                )
            }
        primary_asym_ids = (
            set(entry.ligand_like_chains)
            - monoatomic_ion_asym_ids
            - known_artifact_asym_ids
        )
        if not primary_asym_ids:
            LOG.info(
                f"PDB {entry.pdb_id!r} has only known artifact or "
                "monoatomic-ion ligand chains; "
                "skipping biological assemblies"
            )
            return entry
        if monoatomic_ion_asym_ids or known_artifact_asym_ids:
            LOG.info(
                "PDB %s: deferring %d ion and %d known-artifact ASU ligand chains",
                entry.pdb_id,
                len(monoatomic_ion_asym_ids),
                len(known_artifact_asym_ids),
            )
        ligands: dict[str, Ligand] = {}

        assembly_ids = pdbx.list_assemblies(cif_file_obj)
        for assembly_id in assembly_ids:
            try:
                biounit = build_biounit(cif_file_obj, assembly_id)
            except Exception as e:
                LOG.warning(f"Could not build assembly {assembly_id}: {e}")
                continue
            entry.biounit_chain_ids[assembly_id] = sorted(
                str(chain_id) for chain_id in np.unique(biounit.chain_id)
            )
            spatial_index = BiounitSpatialIndex.from_atoms(
                biounit,
                max(
                    plip_complex_threshold,
                    neighboring_residue_threshold,
                    neighboring_ligand_threshold,
                ),
            )
            water_chains = get_water_chain_ids(biounit)
            primary_ligands = entry._collect_ligands_from_biounit(
                biounit,
                assembly_id,
                plip_complex_threshold,
                neighboring_residue_threshold,
                neighboring_ligand_threshold,
                data_dir,
                ligand_asym_ids=primary_asym_ids,
                water_chains=water_chains,
                spatial_index=spatial_index,
            )
            proper_primary_ligands = [
                ligand for ligand in primary_ligands.values() if ligand.is_proper
            ]
            if not proper_primary_ligands:
                ligands.update(primary_ligands)
                continue

            ion_instance_chains = {
                chain
                for chain in spatial_index.chain_ids
                if "." in chain
                and chain.split(".", maxsplit=1)[1] in monoatomic_ion_asym_ids
            }
            retained_ion_chains = entry._connected_deferred_ligand_chains(
                biounit,
                spatial_index,
                proper_primary_ligands,
                ion_instance_chains,
                min_shared_pocket_members=min_shared_pocket_members,
                neighboring_residue_threshold=neighboring_residue_threshold,
                interaction_search_threshold=plip_complex_threshold,
            )
            ion_ligands = entry._collect_ligands_from_biounit(
                biounit,
                assembly_id,
                plip_complex_threshold,
                neighboring_residue_threshold,
                neighboring_ligand_threshold,
                data_dir,
                ligand_instance_chains=retained_ion_chains,
                water_chains=water_chains,
                spatial_index=spatial_index,
            )
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
                ligand_indices = spatial_index.atom_indices_for_chain(
                    ligand.instance_chain
                )
                nearby_indices = spatial_index.atom_indices_near(
                    biounit.coord[ligand_indices], plip_complex_threshold
                )
                referenced_artifact_chains.update(
                    str(chain)
                    for chain in np.unique(biounit.chain_id[nearby_indices])
                    if "." in str(chain)
                    and str(chain).split(".", maxsplit=1)[1] in known_artifact_asym_ids
                )
            artifact_ligands = entry._collect_ligands_from_biounit(
                biounit,
                assembly_id,
                plip_complex_threshold,
                neighboring_residue_threshold,
                neighboring_ligand_threshold,
                data_dir,
                ligand_instance_chains=referenced_artifact_chains,
                water_chains=water_chains,
                spatial_index=spatial_index,
            )
            ligands.update(primary_ligands)
            ligands.update(ion_ligands)
            ligands.update(artifact_ligands)
        entry.set_systems(
            ligands,
            min_shared_pocket_members=min_shared_pocket_members,
        )
        if entry.systems:
            entry.symmetry_mate_contacts = get_symmetry_mate_contacts(
                cif_file_obj,
                symmetry_mate_contact_threshold,
            )
            if entry.symmetry_mate_contacts:
                for system in entry.systems.values():
                    for ligand in system.ligands:
                        ligand.label_crystal_contacts(entry.symmetry_mate_contacts)
        entry.label_chains()
        # One SDF per ligand, keyed by primary asym, spanning every member
        # chain so covalently-linked ligand chains are saved as one molecule.
        retained_ligand_chain_groups = {
            ligand.asym_id: ligand.member_asym_ids
            for system in entry.systems.values()
            for ligand in system.ligands
        }
        if save_folder is not None and retained_ligand_chain_groups:
            save_ligands(
                atoms,
                retained_ligand_chain_groups,
                save_folder / entry.pdb_id / "ligand_files",
            )
        return entry

    @classmethod
    def from_custom_cif_file(
        cls,
        pdb_id: str,
        cif_file: Path,
        ligand_smiles_dict: dict[str, str] | None = None,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        min_polymer_size: int = 12,
        plip_complex_threshold: float = 10.0,
        save_folder: Path | None = None,
        min_shared_pocket_members: int = 3,
        save_fixed_cif: Path | None = None,
    ) -> Entry:
        """
        Creates entry from an extrernal (non-PDB) mmCIF file

        Parameters
        ----------
        pdb_id : str
            annotation be used in PDB ID column
        cif_file : Path
            mmcif files of interest
        ligand_smiles_dict : dict[str, str] | None, optional
            Mapping of component ID (e.g. ``LIG``) to SMILES.
            Required for unknown ligands without ``_chem_comp_bond``
            (typical of cofolding outputs). Known CCD compounds
            are handled automatically.
        neighboring_residue_threshold : float, optional
            Max distance (Å) for neighboring receptor residues, by default 6.0
        neighboring_ligand_threshold : float, optional
            Max distance (Å) for neighboring ligands, by default 4.0
        min_polymer_size : int, optional
            Minimum residue count for a chain to be polymer (not ligand), by default 10
        save_folder : Path | None, optional
            Root directory for canonical ASU ligand SDFs, by default None.
        save_fixed_cif : Path | None, optional
            If provided and the CIF needed bond-order enrichment, write
            the enriched copy to this path. The input CIF at ``cif_file``
            is never mutated. Raises ``FileExistsError`` if the target
            already exists and ``ValueError`` if it resolves to the same
            path as ``cif_file``. By default (``None``) no file is written.

        Returns
        -------
        Entry
            Entry object for the given pdbid

        Raises
        ------
        MissingBondOrderError
            If the CIF contains unknown ligands and no ``ligand_smiles_dict``
            is provided.
        FileExistsError
            If ``save_fixed_cif`` already exists.
        ValueError
            If ``save_fixed_cif`` points at the input ``cif_file``;
            if biotite returns no bonds at all (corrupted / missing
            bond information CIF); or any error propagated from
            :func:`~plinder.data.utils.annotations.cif_utils.enrich_cif_with_smiles_bonds`
            (invalid SMILES, atom-count / element-order mismatch in the
            positional path, sanitize / template-match failure in the
            opt-in substructure path, or multi-instance comp_id
            divergence).
        """
        from plinder.data.utils.annotations.cif_utils import (
            MissingBondOrderError,
            enrich_cif_with_smiles_bonds,
            get_unknown_ligand_ids,
            read_mmcif_file,
        )

        if save_folder is not None:
            ligand_dir = Path(save_folder) / pdb_id / "ligand_files"
            if ligand_dir.exists():
                shutil.rmtree(ligand_dir)

        # Read CIF once into memory — we mutate this copy only, never the file on disk.
        cif_file_obj = read_mmcif_file(cif_file)

        # Multi-model CIFs (NMR ensembles, Boltz multi-sample, PyMOL
        # states) are processed using model 1 only — surface a warning
        # so users know other models were dropped and can call this
        # function per-model if they need ensemble analysis.
        n_models = get_model_count(cif_file_obj)
        if n_models > 1:
            LOG.warning(
                f"Custom CIF has {n_models} models — using model 1 only. "
                "Call from_custom_cif_file once per model for ensemble analysis."
            )

        # Check for missing bond orders and enrich CIF in-memory if needed
        unknown_ids = get_unknown_ligand_ids(cif_file_obj)
        enrichment_applied = False
        if unknown_ids:
            if ligand_smiles_dict is None:
                raise MissingBondOrderError(
                    f"CIF contains unknown ligands {unknown_ids} with no "
                    "_chem_comp_bond and no CCD match. "
                    "Provide ligand_smiles_dict to assign bond orders."
                )
            enrich_cif_with_smiles_bonds(
                cif_file_obj,
                ligand_smiles=ligand_smiles_dict,
            )
            enrichment_applied = True

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

        cif_data = list(cif_file_obj.values())[0]
        atoms = get_structure_with_altloc(
            cif_file_obj, model=1, use_author_fields=False, include_bonds=True
        )
        atoms = atoms[~is_hydrogen_isotope(atoms.element)]
        if atoms.bonds is None:
            # ``include_bonds=True`` returning ``None`` means biotite
            # derived **no bonds at all** for the structure — every
            # residue lookup failed. This is a fundamentally broken or
            # corrupted CIF (post-enrichment, no _chem_comp_bond, no
            # _struct_conn, and no CCD coverage for any residue).
            raise ValueError(
                f"Custom CIF {cif_file}: biotite returned no bonds at all "
                "after enrichment — the CIF is corrupted or missing all "
                "bond information (_chem_comp_bond, _struct_conn, and CCD "
                "coverage are all absent)."
            )
        apply_struct_conn_bonds(atoms, cif_data)
        chain_to_seqres = get_label_asym_sequences(cif_data)

        entry = cls(
            pdb_id=pdb_id,
            chain_to_seqres=chain_to_seqres,
        )
        entry._populate_chains(atoms, cif_data)
        entry.ligand_like_chains = detect_ligand_chains(entry, min_polymer_size)
        # Create single biounit with "1." prefix on chain IDs
        biounit = atoms.copy()
        biounit.chain_id = np.array([f"1.{c}" for c in biounit.chain_id])
        entry.biounit_chain_ids["1"] = sorted(
            str(chain_id) for chain_id in np.unique(biounit.chain_id)
        )
        spatial_index = BiounitSpatialIndex.from_atoms(
            biounit,
            max(
                plip_complex_threshold,
                neighboring_residue_threshold,
                neighboring_ligand_threshold,
            ),
        )
        water_chains = get_water_chain_ids(biounit)
        ligands = entry._collect_ligands_from_biounit(
            biounit,
            "1",  # single assembly; custom CIFs lack _pdbx_struct_assembly
            plip_complex_threshold,
            neighboring_residue_threshold,
            neighboring_ligand_threshold,
            data_dir=None,
            ligand_smiles_dict=ligand_smiles_dict,
            water_chains=water_chains,
            spatial_index=spatial_index,
        )
        entry._finalize(
            ligands,
            min_shared_pocket_members=min_shared_pocket_members,
        )
        # One SDF per ligand, keyed by primary asym, spanning every member
        # chain so covalently-linked ligand chains are saved as one molecule.
        retained_ligand_chain_groups = {
            ligand.asym_id: ligand.member_asym_ids
            for system in entry.systems.values()
            for ligand in system.ligands
        }
        if save_folder is not None and retained_ligand_chain_groups:
            save_ligands(
                atoms,
                retained_ligand_chain_groups,
                save_folder / entry.pdb_id / "ligand_files",
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
        """
        __Map author chain id to asym id
        """
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
                i_c.split(".")[1]
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
        ligand_chains = set()
        for system in self.systems.values():
            if system.system_type == "holo":
                holo_chains.update(
                    [c.split(".")[1] for c in system.protein_chains_asym_id]
                )
            ligand_chains.update([l.asym_id for l in system.ligands])
        for chain in self.chains:
            if chain not in ligand_chains and chain not in holo_chains:
                self.chains[chain].holo = False
            elif chain in holo_chains:
                self.chains[chain].holo = True

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
        ]
        for field in columns:
            name = f"entry_{field}"
            data[name] = getattr(self, field, None)

        if self.validation:
            data.update(self.format_validation(criteria))
            data["entry_pass_validation_criteria"] = self.pass_criteria
        return data

    def chains_to_df(self) -> pd.DataFrame:
        """Return one normalized metadata row for each receptor polymer chain."""
        columns = [
            "entry_pdb_id",
            "chain_asym_id",
            "chain_auth_id",
            "chain_entity_id",
            "chain_type",
            "chain_receptor_type",
            "chain_length",
            "chain_num_unresolved_residues",
            "chain_is_holo",
            "chain_uniprot_ids",
        ]
        rows = []
        for chain_id in sorted(self.chains):
            chain = self.chains[chain_id]
            if chain_id in self.ligand_like_chains:
                continue
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
                    "chain_length": chain.length,
                    "chain_num_unresolved_residues": chain.num_unresolved_residues,
                    "chain_is_holo": chain.holo,
                    "chain_uniprot_ids": sorted(chain.mappings.get("UniProt", {})),
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def biounit_chains_to_df(self) -> pd.DataFrame:
        """Return biological-assembly membership once per chain instance."""
        columns = [
            "entry_pdb_id",
            "biounit_id",
            "chain_instance",
            "chain_asym_id",
            "chain_role",
        ]
        rows = []
        water_chains = set(self.water_chains)
        ligand_chains = set(self.ligand_like_chains)
        for biounit_id, chain_instances in sorted(self.biounit_chain_ids.items()):
            for chain_instance in sorted(set(chain_instances)):
                asym_id = chain_instance.split(".", maxsplit=1)[-1]
                if asym_id in water_chains:
                    role = "water"
                elif asym_id in ligand_chains:
                    role = "ligand"
                else:
                    role = "receptor"
                rows.append(
                    {
                        "entry_pdb_id": self.pdb_id,
                        "biounit_id": str(biounit_id),
                        "chain_instance": chain_instance,
                        "chain_asym_id": asym_id,
                        "chain_role": role,
                    }
                )
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
        rows = []
        entry_data = self.format()
        for system in self.systems:
            annotation = self.systems[system]
            system_data = annotation.format(
                self.chains,
                self.pass_criteria,
            )
            for ligand in self.systems[system].ligands:
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
                    all_pocket_residues[chain.split(".")[1]].update(
                        ligand.pocket_residues[chain].keys()
                    )
        for chain in self.chains:
            if chain in all_pocket_residues:
                self.chains[chain].residues = {
                    r: self.chains[chain].residues[r]
                    for r in all_pocket_residues[chain]
                }
            else:
                self.chains[chain].residues = {}

    def set_validation(
        self,
        validation_file: Path,
        cif_file: Path,
        thresholds: ResidueValidationThresholds = ResidueValidationThresholds(),
    ) -> None:
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
