# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import re
import typing as ty
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import networkx as nx
import pandas as pd
from ost import io, mol
from PDBValidation.ValidationFactory import ValidationFactory
from plip.basic import config
from posebusters import PoseBusters
from pydantic import BeforeValidator, Field
from rdkit import RDLogger

from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
    ResidueValidationThresholds,
)
from plinder.data.utils.annotations.interaction_utils import (
    get_covalent_connections,
    get_symmetry_mate_contacts,
)
from plinder.data.utils.annotations.interface_gap import annotate_interface_gaps
from plinder.data.utils.annotations.ligand_utils import Ligand, validate_chain_residue
from plinder.data.utils.annotations.protein_utils import (
    Chain,
    detect_ligand_chains,
    get_chain_external_mappings,
    get_entry_info,
    read_mmcif_container,
)
from plinder.data.utils.annotations.save_utils import (
    save_cif_file,
    save_ligands,
    save_pdb_file,
)
from plinder.data.utils.annotations.utils import DocBaseModel

LOG = setup_logger(__name__)
RDLogger.DisableLog("rdApp.*")
ECOD_DATA = None

# Ignore Biolip artifacts
config.biolip_list = []


SymmetryMateContacts = ty.Annotated[
    dict[tuple[str, int], dict[tuple[str, int], set[int]]],
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

    return int(re.sub("[^0-9]", "", x))


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
    and it's neighboring ligands and protein residues

    """

    def proper_ligands(self) -> list[Ligand]:
        return [ligand for ligand in self.ligands if ligand.is_proper]

    @cached_property
    def protein_chains_asym_id(self) -> list[str]:
        """
        Interacting protein chains of the system
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
        return sum(l.num_pocket_residues for l in self.ligands)

    @cached_property
    def proper_num_pocket_residues(self) -> int:
        """
        Number of pocket residues of the system excluding ions and artifacts
        """
        return sum(l.num_pocket_residues for l in self.proper_ligands())

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
    def has_kinase_inhibitor(self) -> bool:
        """
        Whether the system has a kinase inhibitor
        """
        return any(l.is_kinase_inhibitor for l in self.ligands)

    @cached_property
    def has_binding_affinity(self) -> bool:
        """
        Whether any ligand in the system has a binding affinity from BindingDB
        """
        return any(l.binding_affinity is not None for l in self.ligands)

    @cached_property  # TODO: change this to exclude residues only interacting with artifacts or ions
    def pocket_residues(self) -> dict[str, dict[int, str]]:
        """
        __Pockets residues of the system
        """
        all_residues: dict[str, dict[int, str]] = defaultdict(dict)
        for ligand in self.ligands:
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
            if ligand.is_artifact:
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
    ) -> dict[str, ty.Any]:
        if chain_type == "protein":
            sub_chains = self.protein_chains_asym_id
        elif chain_type == "ligand":
            sub_chains = self.ligand_chains
        else:
            raise ValueError(f"chain_type={chain_type} not understood")
        sub_chain_list = [ch.split(".") for ch in sub_chains]

        sub_chains_data = [
            chains[c].format(int(instance)) for instance, c in sub_chain_list
        ]

        if len(sub_chains_data) == 0:
            return {}
        data: dict[str, list[str]] = defaultdict(list)
        for sub_chain in sub_chains_data:
            for key in sub_chain:
                data[f"system_{chain_type}_chains_{key}"].append(sub_chain[key])
        return data

    @cached_property
    def num_protein_chains(self) -> int:
        """
        Number of interacting protein chains of the system
        """
        return len(self.protein_chains_asym_id)

    @cached_property
    def proper_num_protein_chains(self) -> int:
        """
        Number of interacting protein chains of the system excluding ions and artifacts
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
        return data

    @cached_property
    def waters(self) -> dict[str, list[int]]:
        """
        __Waters interacting (as detected by PLIP) with any of the ligands in the system
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
        return max(
            ligand.molecular_weight if ligand.molecular_weight is not None else -1.0
            for ligand in self.ligands
            if ligand.molecular_weight is not None
        )

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
        if self.num_heavy_atoms is None:
            return None
        return self.num_atoms_with_crystal_contacts / self.num_heavy_atoms

    def selection(self, include_waters: bool = True) -> str:
        ligand_selection = " or ".join(
            f"({ligand.selection})" for ligand in self.ligands
        )
        protein_selection = " or ".join(
            f"(cname={mol.QueryQuoteName(chain)})"
            for chain in self.protein_chains_asym_id
        )
        selection = f"({ligand_selection}) or ({protein_selection})"
        if include_waters and len(self.waters):
            selection += f" or {self.select_waters()}"
        return selection

    def save_system(
        self,
        chain_to_seqres: dict[str, str],
        biounit: mol.EntityHandle,
        info: io.MMCifInfoBioUnit,
        system_folder: Path,
        include_waters: bool = True,
    ) -> None:
        system_folder.mkdir(exist_ok=True)
        with open(system_folder / "sequences.fasta", "w") as f:
            for i_c in self.protein_chains_asym_id:
                c = i_c.split(".")[1]
                if c in chain_to_seqres:
                    f.write(f">{i_c}\n")
                    f.write(chain_to_seqres[c] + "\n")
        selection = self.selection(include_waters=include_waters)
        ent_system = mol.CreateEntityFromView(
            biounit.Select(selection),
            True,
        )
        (system_folder / "ligand_files").mkdir(exist_ok=True)
        save_ligands(
            ent_system,
            [ligand.selection for ligand in self.ligands],
            self.ligand_chains,
            [l.smiles for l in self.ligands],
            [l.num_unresolved_heavy_atoms for l in self.ligands],
            system_folder / "ligand_files",
        )
        save_cif_file(ent_system, info, self.id, system_folder / "system.cif")
        selection = " or ".join(f"chain='{c}'" for c in self.protein_chains_asym_id)
        if include_waters and len(self.waters):
            selection += f" or {self.select_waters()}"
        save_cif_file(
            ent_system.Select(selection),
            info,
            self.id,
            system_folder / "receptor.cif",
        )
        try:
            # TODO: move out and add a flag instead
            save_pdb_file(
                biounit,
                mol.CreateEntityFromView(ent_system.Select(selection), True),
                self.protein_chains_asym_id,
                [],
                system_folder / "receptor.pdb",
                system_folder / "chain_mapping.json",
                self.waters if include_waters else {},
                system_folder / "water_mapping.json",
            )
        except Exception as e:
            LOG.error(f"save_system: Error saving system in PDB format {self.id}: {e}")

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

    def run_posebusters_on_system(self, system_folder: Path) -> None:
        """
        Run posebusters on the system
        """
        pb = PoseBusters(config="redock")
        receptor_file = system_folder / "receptor.pdb"
        if not receptor_file.exists():
            return
        for ligand in self.ligands:
            ligand_file = (
                system_folder / "ligand_files" / f"{ligand.instance_chain}.sdf"
            )
            if not ligand_file.exists():
                continue
            try:
                result_dict = pb.bust(
                    mol_pred=str(ligand_file),
                    mol_true=str(ligand_file),
                    mol_cond=str(receptor_file),
                    full_report=True,
                ).to_dict()
            except Exception as e:
                LOG.error(
                    f"run_posebusters: Error running posebusters on {ligand.id}: {e}"
                )
                continue
            key = (str(ligand_file), ligand.instance_chain)
            ligand.posebusters_result = {
                k: v.get(key) for k, v in result_dict.items() if v.get(key)
            }

    def get_pocket_domains(self, chains_dict: dict[str, Chain]) -> dict[str, str]:
        global ECOD_DATA
        if ECOD_DATA is None:
            data_dir = Path(get_config().data.plinder_dir)
            ECOD_DATA = pd.read_parquet(data_dir / "dbs" / "ecod" / "ecod.parquet")
        ecod_df = ECOD_DATA[ECOD_DATA["pdb"] == self.pdb_id]
        ecod_mapping = dict(zip(ecod_df["domainid"], ecod_df["domain"]))
        pocket_mapping: dict[str, dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )
        for (
            neighboring_chain,
            neighboring_residues_list,
        ) in self.pocket_residues.items():
            neighboring_chain = neighboring_chain.split(".")[-1]
            neighboring_residues_set = {int(i) for i in neighboring_residues_list}
            neighboring_residue_numbers_dict = {
                int(chains_dict[neighboring_chain].residues[i].auth_number): i
                for i in neighboring_residues_list
            }
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
                    if mapping_name == "ECOD":
                        intersection = set(
                            neighboring_residue_numbers_dict[i]
                            for i in neighboring_residue_numbers_dict
                            if i in unrolled_v_2
                        )
                    else:
                        intersection = neighboring_residues_set.intersection(
                            unrolled_v_2
                        )
                    if len(intersection) > 0:
                        pocket_mapping[mapping_name][domain] += len(intersection)
                        if mapping_name == "ECOD":
                            pocket_mapping[f"{mapping_name}_t_name"][
                                ecod_mapping[domain]
                            ] += len(intersection)
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

    @classmethod
    def from_json(
        cls,
        entry_json: Path,
        clear_non_pocket_residues: bool = False,
        load_for_scoring: bool = False,
        max_protein_chains: int = 5,
        max_ligand_chains: int = 5,
    ) -> Entry:
        """
        Load an entry from a JSON file and prune it. See prune
        for more details.
        """
        with entry_json.open("rb") as f:
            entry = cls.model_validate_json(f.read())
        return entry.prune(
            clear_non_pocket_residues=clear_non_pocket_residues,
            load_for_scoring=load_for_scoring,
            max_protein_chains=max_protein_chains,
            max_ligand_chains=max_ligand_chains,
        )

    @classmethod
    def from_cif_file(
        cls,
        cif_file: Path,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        min_polymer_size: int = 10,  # TODO: this used to be max_non_small_mol_ligand_length
        max_non_small_mol_ligand_length: int = 20,  # TODO: review and make consistent
        save_folder: Path | None = None,
        max_protein_chains_to_save: int = 5,
        max_ligand_chains_to_save: int = 5,
        plip_complex_threshold: float = 10.0,
        skip_save_systems: bool = False,
        skip_posebusters: bool = False,
        symmetry_mate_contact_threshold: float = 5.0,
    ) -> Entry:
        """
        Load an entry object from mmcif files

        Parameters
        ----------
        cif_file : Path
            mmcif files of interest
        neighboring_residue_threshold : float
            Distance from ligand for protein \
                residues to be considered a ligand
        neighboring_ligand_threshold : float
            Distance from ligand for other ligans \
                to be considered a ligand
        min_polymer_size : int = 10
            Minimum number of residues for chain to be seen as a \
                polymer, or Maximum number of residues for chain to be seen as a ligand \
        max_non_small_mol_ligand_length: int = 20
            Maximum length of polymer that should be assessed for potentially being ligand
        save_folder : Path
            Path to save files
        max_protein_chains_to_save : int
            Maximum number of protein chains to save
        max_ligand_chains_to_save : int
            Maximum number of protein chains to save
        plip_complex_threshold=10
            Maximum distance from ligand to residues to be
            included for plip calculations.
        skip_save_systems: bool = False
            skips saving system files
        skip_posebusters: bool = False
            skips running posebusters analysis

        Returns
        -------
        Entry
            Entry object for the given pdbid
        """
        ent, seqres, info = io.LoadMMCIF(
            str(cif_file), seqres=True, info=True, remote=False
        )
        cif_data = read_mmcif_container(cif_file)
        symmetry_mate_contacts = get_symmetry_mate_contacts(
            cif_file, symmetry_mate_contact_threshold
        )
        entry_info = get_entry_info(cif_data)
        per_chain = get_chain_external_mappings(cif_data)
        # TODO: annotate_interface_gaps does not use the same ligand chain definitions as the rest
        # move this to later after protein/ligand chain assignment?
        interface_proximal_gaps = annotate_interface_gaps(cif_file)
        resolution = entry_info.get("entry_resolution")
        r = None
        if resolution is not None:
            try:
                r = float(resolution)
            except ValueError:
                r = None
        entry = cls(
            pdb_id=info.struct_details.entry_id.lower(),
            release_date=info.revisions.GetDateOriginal(),
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
            covalent_bonds=get_covalent_connections(cif_data),
            chain_to_seqres={c.name: c.string for c in seqres},
            symmetry_mate_contacts=symmetry_mate_contacts,
        )
        entry.chains = {
            chain.name: Chain.from_ost_chain(
                chain, info, len(entry.chain_to_seqres.get(chain.name, ""))
            )
            for chain in ent.chains
            if chain.type != mol.CHAINTYPE_WATER
        }
        entry.water_chains = [
            chain.name for chain in ent.chains if chain.type == mol.CHAINTYPE_WATER
        ]

        data_dir = None
        if save_folder is not None:
            data_dir = save_folder.parent.parent
        for chain in per_chain:
            entry.chains[chain].mappings = per_chain[chain]
        if data_dir is not None:
            entry.add_ecod()
            entry.add_panther(data_dir / "dbs" / "panther")
            entry.add_kinase(data_dir / "dbs" / "kinase" / "kinase_uniprotac.parquet")
        entry.ligand_like_chains = detect_ligand_chains(
            ent, entry, min_polymer_size, max_non_small_mol_ligand_length
        )
        ligands = {}
        biounits = {}
        for biounit_info in info.biounits:
            biounit = mol.alg.CreateBU(ent, biounit_info)
            biounit_ligand_chains = [
                chain.name
                for chain in biounit.chains
                if chain.name.split(".")[1] in entry.ligand_like_chains
            ]
            for ligand_chain in biounit_ligand_chains:
                ligand_instance, ligand_asym_id = ligand_chain.split(".")
                data_dir = None
                if save_folder is not None:
                    data_dir = save_folder.parent.parent
                residue_numbers = [
                    residue.number.num
                    for residue in biounit.FindChain(ligand_chain).residues
                ]
                ligand = Ligand.from_pli(
                    pdb_id=entry.pdb_id,
                    biounit_id=biounit_info.id,
                    biounit=biounit,
                    ligand_instance=int(ligand_instance),
                    ligand_chain=entry.chains[ligand_asym_id],
                    residue_numbers=residue_numbers,
                    ligand_like_chains=entry.ligand_like_chains,
                    interface_proximal_gaps=interface_proximal_gaps,
                    all_covalent_dict=entry.covalent_bonds,
                    plip_complex_threshold=plip_complex_threshold,
                    neighboring_residue_threshold=neighboring_residue_threshold,
                    neighboring_ligand_threshold=neighboring_ligand_threshold,
                    data_dir=data_dir,
                )
                if ligand is not None:
                    ligands[ligand.id] = ligand
            biounits[biounit_info.id] = biounit
        entry.set_systems(ligands)
        entry.label_chains()
        if save_folder is not None and not skip_save_systems:
            entry.save_systems(
                info,
                biounits,
                save_folder,
                max_protein_chains_to_save,
                max_ligand_chains_to_save,
            )
        # TODO: this is backwards because it assumes save_systems
        #       has already run but will fail if it hadn't run previously
        #       so we just check if save_folder is None (which it's not in the pipeline)
        # VO: added option to skip to speed up testing!
        if not skip_posebusters:
            entry.run_posebusters(
                save_folder,
                max_protein_chains_to_save,
                max_ligand_chains_to_save,
            )
        return entry

    def set_systems(self, ligands: dict[str, Ligand]) -> None:
        """
        Setter method for system ids for ligands

        Parameters
        ----------
        ligands : dict[str, Ligand]

        Returns
        -------
        None
        """

        G = nx.Graph()
        for ligand_id in ligands:
            G.add_node(ligand_id)
            for neighboring_ligand_instance_chain in (
                ligands[ligand_id].neighboring_ligands
                + ligands[ligand_id].interacting_ligands
            ):
                neighboring_ligand_id = "__".join(
                    [
                        self.pdb_id,
                        ligands[ligand_id].biounit_id,
                        f"{neighboring_ligand_instance_chain}",
                    ]
                )
                if neighboring_ligand_id in ligands:
                    G.add_edge(ligand_id, neighboring_ligand_id)
        system_ligands: dict[int, list[Ligand]] = {}
        for idx, component in enumerate(
            sorted(nx.connected_components(G), key=len, reverse=True)
        ):
            system_ligands[idx + 1] = []
            for ligand_id in component:
                system_ligands[idx + 1].append(ligands[ligand_id])
        self.systems: dict[str, System] = {}
        for ligs in system_ligands.values():
            system = System(
                pdb_id=self.pdb_id,
                biounit_id=ligs[0].biounit_id,
                ligands=sorted(ligs, key=lambda x: x.id),
            )
            if len(system.protein_chains_asym_id):
                self.systems[system.id] = system

    @cached_property
    def author_to_asym(self) -> dict[str, str]:
        """
        __Map author chain id to asym id
        """
        return {
            c.auth_id: c.asym_id
            for c in self.chains.values()
            if c.chain_type == mol.CHAINTYPE_POLY_PEPTIDE_L
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
            chains = set(
                self.chains[i_c.split(".")[1]].auth_id
                for system in self.systems.values()
                for i_c in system.protein_chains_asym_id
                if system.system_type == "holo"
            )
        elif chain_type == "apo":
            holo_entities = set(
                self.chains[c].entity_id for c in self.chains if self.chains[c].holo
            )
            chains = set(
                self.chains[c].auth_id
                for c in self.chains
                if not self.chains[c].holo
                and not self.chains[c].entity_id in holo_entities
            )
        elif chain_type == "pred":
            chains = set()
            for c in self.chains:
                if self.chains[c].holo:
                    chains |= set(self.chains[c].mappings.get("UniProt", {}))
        if aln_type == "foldseek":
            if chain_type == "pred":
                return [f"AF-{c}-F1-model_v4_A" for c in chains]
            else:
                return [f"pdb_0000{self.pdb_id}_xyz-enrich.cif_{c}" for c in chains]
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
        Format label for entry-level annotations by prepending \
            label with "entry_"
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
            system_data = self.systems[system].format(self.chains, self.pass_criteria)
            for ligand in self.systems[system].ligands:
                ligand_data = ligand.format(self.chains)
                rows.append({**entry_data, **system_data, **ligand_data})
        return pd.DataFrame(rows)

    def iter_systems(
        self, max_protein_chains: int, max_ligand_chains: int
    ) -> ty.Iterator[tuple[str, System]]:
        for system_id, system in self.systems.items():
            if (
                system.system_type == "holo"
                and len(system.protein_chains_asym_id) <= max_protein_chains
                and len(system.ligand_chains) <= max_ligand_chains
            ):
                yield system_id, system

    def run_posebusters(
        self,
        save_folder: Path | None,
        max_protein_chains: int,
        max_ligand_chains: int,
    ) -> None:
        if save_folder is None:
            LOG.warn("run_posebusters got save_folder=None so skipping")
            return
        for system_id, system in self.iter_systems(
            max_protein_chains, max_ligand_chains
        ):
            save_folder_system = save_folder / system.id
            self.systems[system_id].run_posebusters_on_system(save_folder_system)

    def save_systems(
        self,
        info: io.MMCifInfoBioUnit,
        biounits: mol.EntityHandle,
        save_folder: Path,
        max_protein_chains: int = 5,
        max_ligand_chains: int = 5,
    ) -> None:
        """
        Save system files
        Parameters
        ----------
        self : Entry
            Entry object

        Returns
        -------
        pd.DataFrame
        """
        for _, system in self.iter_systems(max_protein_chains, max_ligand_chains):
            save_folder_system = save_folder / system.id
            system.save_system(
                self.chain_to_seqres,
                biounits[system.biounit_id],
                info,
                save_folder_system,
            )

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
        self.label_crystal_contacts()
        if not validation_file.exists():
            LOG.error(f"set_validation: Validation file not found {validation_file}")
            return
        try:
            doc = ValidationFactory(
                str(validation_file), mmcif_path=str(cif_file)
            ).getValidation()
            self.validation = EntryValidation.from_entry(doc)
            if self.validation and self.validation.r is not None:
                for chain in self.chains:
                    self.chains[chain].set_validation(doc, thresholds)
                for system in self.systems:
                    self.systems[system].set_validation(self.chains, thresholds)
        except Exception as e:
            LOG.error(
                f"set_validation: Error setting validation for {self.pdb_id}: {e}"
            )

    def label_crystal_contacts(self) -> None:
        """
        Label contacts of ligand residues to other symmetry mates
        Excludes neighboring residues (i.e same biounit)
        """
        for system in self.systems:
            for ligand in self.systems[system].ligands:
                crystal_contacts: dict[tuple[str, int], set[int]] = defaultdict(set)
                for residue_number in ligand.residue_numbers:
                    # get all contacts with chains in other asymmetric units
                    contacts = self.symmetry_mate_contacts.get(
                        (ligand.asym_id, residue_number), dict()
                    )
                    for x, y in contacts.items():
                        # keep only contacts with receptor
                        if x[0] not in self.ligand_like_chains:
                            crystal_contacts[x] |= y
                ligand.set_crystal_contacts(crystal_contacts)

    def add_ecod(self) -> None:
        """
        Add ECOD annotations to chains
        """
        global ECOD_DATA
        if ECOD_DATA is None:
            data_dir = Path(get_config().data.plinder_dir)
            ECOD_DATA = pd.read_parquet(data_dir / "dbs" / "ecod" / "ecod.parquet")
        ecod_df = ECOD_DATA[ECOD_DATA["pdb"] == self.pdb_id]
        if not len(ecod_df):
            return
        for chain_id, chain in self.chains.items():
            ecod_df_chain = ecod_df[ecod_df["chain"] == chain.auth_id]
            if not len(ecod_df_chain):
                continue
            mappings: dict[str, set[tuple[str, str]]] = defaultdict(set)
            for row in ecod_df_chain.T.to_dict().values():
                mappings[row["domainid"]].add((row["pdb_from"], row["pdb_to"]))
            self.chains[chain_id].mappings["ECOD"] = {
                k: list(v) for k, v in mappings.items()
            }

    def add_panther(self, panther_data_dir: Path) -> None:
        """
        Add panther annotations to chains
        """
        dfs = {}
        for chain_id, chain in self.chains.items():
            uniprots = chain.mappings.get("UniProt", {})
            mappings: dict[str, set[tuple[str, str]]] = defaultdict(set)
            for uniprot in uniprots:
                shard = uniprot[-1]
                if shard not in dfs:
                    dfs[shard] = pd.read_parquet(
                        panther_data_dir / f"panther_{shard}.parquet"
                    ).set_index("uniprotac")
                if uniprot in dfs[shard].index:
                    mappings[dfs[shard].loc[uniprot]["panther"]] |= set(
                        uniprots[uniprot]  # type: ignore
                    )
            if len(mappings):
                self.chains[chain_id].mappings["PANTHER"] = {
                    k: list(v) for k, v in mappings.items()
                }

    def add_kinase(self, kinase_data_path: Path) -> None:
        """
        Add kinase annotations to chains
        """
        kinase_df = pd.read_parquet(kinase_data_path)
        kinase_df = kinase_df[kinase_df["pdb"] == self.pdb_id]
        if not len(kinase_df):
            return
        for chain_id, chain in self.chains.items():
            kinase_df_chain = kinase_df[kinase_df["chain"] == chain.auth_id]
            if not len(kinase_df_chain):
                continue
            mappings: dict[str, set[tuple[str, str]]] = defaultdict(set)
            uniprots = chain.mappings.get("UniProt", {})
            for row in kinase_df_chain.T.to_dict().values():
                uniprot = uniprots.get(row["uniprot"])
                if uniprot is not None:
                    mappings[row["kinase"]] |= set(uniprot)  # type: ignore
            if len(mappings):
                self.chains[chain_id].mappings["kinase_name"] = {
                    k: list(v) for k, v in mappings.items()
                }


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
