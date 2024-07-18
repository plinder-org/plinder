# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import re
import typing as ty
from collections import Counter, defaultdict
from functools import cached_property
from pathlib import Path

import networkx as nx
import pandas as pd
from ost import conop, io, mol
from PDBValidation.ValidationFactory import ValidationFactory
from plip.basic import config
from posebusters import PoseBusters
from pydantic import BaseModel, Field
from rdkit import RDLogger

from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.get_ligand_validation import (
    EntryValidation,
    ResidueListValidation,
    SystemValidationThresholds,
)
from plinder.data.utils.annotations.interaction_utils import (
    get_covalent_connections,
    run_plip_on_split_structure,
)
from plinder.data.utils.annotations.interface_gap import annotate_interface_gaps
from plinder.data.utils.annotations.ligand_utils import Ligand
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

LOG = setup_logger(__name__)
RDLogger.DisableLog("rdApp.*")
COMPOUND_LIB = conop.GetDefaultLib()
ECOD_DATA = None

# Ignore Biolip artifacts
config.biolip_list = []


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


class System(BaseModel):
    pdb_id: str
    biounit_id: str
    ligands: list[Ligand]
    ligand_validation: ResidueListValidation | None = None
    pocket_validation: ResidueListValidation | None = None
    pass_criteria: bool | None = None

    """
    This dataclass defines as system which included a protein-ligand complex
    and it's neighboring ligands and protein residues

    Parameters
    ----------
    pdb_id : str
        4-letter pdbid of interest
    biounit_id: str
        biounit id
    ligands: list[Ligand]
        List of Ligands in a systems

    Attributes
    ----------
    interacting_protein_chains : list[tuple[int, str]]

    neighboring_protein_chains: list[tuple[int, str]]
        biounit id
    ligand_chains: list[tuple[int, str]]
        List of Ligands in a systems
    num_pocket_residues : int
        Number of protein residues in pocket
    num_interactions : int
        Number of interactions
    id : str
        system id
    system_type : str
        system type, whether real ligand, ion, artifacts, etc

    Methods
    -------
    format_chains : dict[str, str]
        Add descriptive-formating to chain, where neighboring, \
            interacting residues or ligand
    format_system : dict[str, ty.Any]
        Add descriptive-formating to systems
    save_system : None
        Save systems into cif, odb and sdf files

    Examples
    --------
    TODO: Add example

    """

    @cached_property
    def interacting_protein_chains(self) -> list[str]:
        return sorted(
            set(
                chain
                for ligand in self.ligands
                for chain in ligand.interacting_protein_chains
            )
        )

    @cached_property
    def neighboring_protein_chains(self) -> list[str]:
        return sorted(
            set(
                chain
                for ligand in self.ligands
                for chain in ligand.neighboring_protein_chains
            )
        )

    @cached_property
    def ligand_chains(self) -> list[str]:
        return [f"{ligand.instance}.{ligand.asym_id}" for ligand in self.ligands]

    @cached_property
    def num_pocket_residues(self) -> int:
        return sum(l.num_pocket_residues for l in self.ligands)

    @cached_property
    def num_interactions(self) -> int:
        return sum(l.num_interactions for l in self.ligands)

    @cached_property
    def id(self) -> str:
        return "__".join(
            [
                self.pdb_id,
                self.biounit_id,
                "_".join(self.interacting_protein_chains),
                "_".join(self.ligand_chains),
            ]
        )

    @cached_property
    def system_type(self) -> str:
        if any(not l.is_ion and not l.is_artifact for l in self.ligands):
            return "holo"
        elif any(l.is_ion for l in self.ligands):
            return "ion"
        else:
            return "artifact"

    @cached_property
    def has_kinase_inhibitor(self) -> bool:
        return any(l.is_kinase_inhibitor for l in self.ligands)

    @cached_property
    def pocket_residues(self) -> dict[str, dict[int, str]]:
        all_residues: dict[str, dict[int, str]] = defaultdict(dict)
        for ligand in self.ligands:
            ligand_pocket_residues = ligand.pocket_residues
            for chain in ligand_pocket_residues:
                all_residues[chain].update(ligand_pocket_residues[chain])
        return all_residues

    @cached_property
    def interactions(self) -> dict[str, dict[int, list[str]]]:
        all_interactions: dict[str, dict[int, list[str]]] = defaultdict(
            lambda: defaultdict(list)
        )
        for ligand in self.ligands:
            for chain in ligand.interactions:
                for residue in ligand.interactions[chain]:
                    all_interactions[chain][residue].extend(
                        ligand.interactions[chain][residue]
                    )
        return all_interactions

    @cached_property
    def interactions_counter(self) -> dict[str, dict[int, ty.Counter[str]]]:
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
    ) -> dict[str, str]:
        if chain_type == "interacting_protein":
            sub_chains = self.interacting_protein_chains
        elif chain_type == "neighboring_protein":
            sub_chains = self.neighboring_protein_chains
        elif chain_type == "ligand":
            sub_chains = self.ligand_chains
        else:
            raise ValueError(f"chain_type={chain_type} not understood")
        sub_chain_list = [ch.split(".") for ch in sub_chains]

        sub_chain_list = [
            chains[c].to_dict(int(instance))  # type: ignore
            for instance, c in sub_chain_list
        ]

        if len(sub_chain_list) == 0:
            return {}
        data: dict[str, list[str]] = defaultdict(list)
        for sub_chain in sub_chain_list:
            for key in sub_chain:
                data[f"system_{chain_type}_chains{key}"].append(str(sub_chain[key]))  # type: ignore
        return {k: ";".join(v) for k, v in data.items()}

    @property
    def num_interacting_protein_chains(self) -> int:
        return len(self.interacting_protein_chains)

    @property
    def num_ligand_chains(self) -> int:
        return len(self.ligands)

    def format_system(self, chains: dict[str, Chain]) -> dict[str, ty.Any]:
        data = {
            "system_id": self.id,
            "system_pdb_id": self.pdb_id,
            "system_biounit_id": self.biounit_id,
            "system_type": self.system_type,
            "system_num_pocket_residues": self.num_pocket_residues,
            "system_num_interactions": self.num_interactions,
            "system_num_interacting_protein_chains": len(
                self.interacting_protein_chains
            ),
            "system_num_neighboring_protein_chains": len(
                self.neighboring_protein_chains
            ),
            "system_num_ligand_chains": len(self.ligand_chains),
            "system_has_kinase_inhibitor": self.has_kinase_inhibitor,
        }
        pocket_mapping = self.get_pocket_domains(chains)
        for mapping in pocket_mapping:
            data[f"system_pocket_{mapping}"] = pocket_mapping[mapping]
        data["suitable_for_ml_training"] = self.suitable_for_ml_training()
        if self.ligand_validation:
            ligand_validation = self.ligand_validation.to_dict()
            data.update({f"system_ligand_{k}": v for k, v in ligand_validation.items()})
        if self.pocket_validation:
            pocket_validation = self.pocket_validation.to_dict()
            data.update({f"system_pocket_{k}": v for k, v in pocket_validation.items()})
        data["system_pass_validation_criteria"] = self.pass_criteria
        for chain_type in ["interacting_protein", "neighboring_protein", "ligand"]:
            data.update(self.format_chains(chain_type, chains))
        return data

    @property
    def waters(self) -> dict[str, list[int]]:
        waters: dict[str, list[int]] = defaultdict(list)
        for ligand in self.ligands:
            for chain in ligand.waters:
                waters[chain] += ligand.waters[chain]
        return waters

    @property
    def select_waters(self) -> str:
        query = []
        for chain in self.waters:
            chain_query = " or ".join(f"rnum={resnum}" for resnum in self.waters[chain])
            query.append(f"(chain='{chain}' and ({chain_query}))")
        return " or ".join(query)

    def save_system(
        self,
        chain_to_seqres: dict[str, str],
        biounit: mol.EntityHandle,
        info: io.MMCifInfoBioUnit,
        system_folder: Path,
    ) -> None:
        system_folder.mkdir(exist_ok=True)
        with open(system_folder / "sequences.fasta", "w") as f:
            for i_c in self.interacting_protein_chains:
                c = i_c.split(".")[1]
                if c in chain_to_seqres:
                    f.write(f">{i_c}\n")
                    f.write(chain_to_seqres[c] + "\n")
        selection = " or ".join(
            f"chain='{c}'" for c in self.interacting_protein_chains + self.ligand_chains
        )
        if len(self.waters):
            selection += f" or {self.select_waters}"
        ent_system = mol.CreateEntityFromView(
            biounit.Select(selection),
            True,
        )
        (system_folder / "ligand_files").mkdir(exist_ok=True)
        save_ligands(
            ent_system,
            self.ligand_chains,
            [l.smiles for l in self.ligands],
            [l.num_unresolved_heavy_atoms for l in self.ligands],
            system_folder / "ligand_files",
        )
        save_cif_file(ent_system, info, self.id, system_folder / "system.cif")
        selection = " or ".join(f"chain='{c}'" for c in self.interacting_protein_chains)
        if len(self.waters):
            selection += f" or {self.select_waters}"
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
                ent_system.Copy(),
                self.interacting_protein_chains,
                self.ligand_chains,
                system_folder / "system.pdb",
                system_folder / "chain_mapping.json",
                self.waters,
                system_folder / "water_mapping.json",
            )
            with open(system_folder / "chain_mapping.json") as f:
                chain_mapping = json.load(f)
            system_pdb = io.LoadPDB(str(system_folder / "system.pdb"))
            io.SavePDB(
                system_pdb.Select(
                    " or ".join(
                        f"chain='{chain_mapping[c]}'"
                        for c in self.interacting_protein_chains
                    )
                    + " or chain='_'"
                ),
                str(system_folder / "receptor.pdb"),
            )
        except Exception as e:
            LOG.error(f"save_system: Error saving system in PDB format {self.id}: {e}")

    def set_validation(
        self,
        chains: dict[str, Chain],
        thresholds: SystemValidationThresholds,
        entry_pass_criteria: bool,
    ) -> None:
        self.ligand_validation = ResidueListValidation.from_residues(
            [
                chains[c.split(".")[1]].residues[r].validation  # type: ignore
                for c in self.ligand_chains
                for r in chains[c.split(".")[1]].residues
            ],
            thresholds.residue_thresholds,
        )
        self.pocket_validation = ResidueListValidation.from_residues(
            [
                chains[c.split(".")[1]].residues[r].validation  # type: ignore
                for c in self.pocket_residues
                for r in self.pocket_residues[c]
            ],
            thresholds.residue_thresholds,
        )
        if self.ligand_validation is not None and self.pocket_validation is not None:
            self.pass_criteria = (
                entry_pass_criteria
                and self.ligand_validation.pass_criteria(
                    thresholds.residue_list_thresholds["ligand"]
                )
                and self.pocket_validation.pass_criteria(
                    thresholds.residue_list_thresholds["pocket"]
                )
            )
        else:
            self.pass_criteria = None

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
            from plinder.data.pipeline.config import IngestConfig

            data_dir = Path(IngestConfig().plinder_dir)
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

    def suitable_for_ml_training(self) -> bool:
        """Defines what suitable for ml training

        Returns
        -------
        bool
            True if it should be included for ml training
        """
        return bool(
            (self.num_interactions > 2)
            and (not any(ligand.is_artifact for ligand in self.ligands))
            and (not any(ligand.is_invalid for ligand in self.ligands))
            and (not any(ligand.is_ion for ligand in self.ligands))
            and (
                all(
                    ligand.posebusters_result.get("mol_pred_loaded", False)
                    for ligand in self.ligands
                )
            )
            and (
                all(
                    ligand.posebusters_result.get("passes_valence_checks", False)
                    for ligand in self.ligands
                )
            )
            and (
                all(
                    ligand.posebusters_result.get("sanitization", False)
                    for ligand in self.ligands
                )
            )
            and (
                all(
                    ligand.posebusters_result.get("passes_kekulization", False)
                    for ligand in self.ligands
                )
            )
        )


class Entry(BaseModel):
    pdb_id: str
    release_date: str
    oligomeric_state: str | None
    determination_method: str | None
    keywords: str | None
    pH: str | None
    resolution: float | None
    chains: dict[str, Chain] = Field(default_factory=dict)
    ligand_like_chains: dict[str, str] = Field(default_factory=dict)
    systems: dict[str, System] = Field(default_factory=dict)
    covalent_bonds: dict[str, list[dict[str, str]]] = Field(default_factory=dict)
    chain_to_seqres: dict[str, str] = Field(default_factory=dict)
    validation: EntryValidation | None = None
    pass_criteria: bool | None = None
    water_chains: list[str] = Field(default_factory=list)

    """
    This dataclass defines as system which included a protein-ligand complex
    and it's neighboring ligands and protein residues. For access to all
    cached properties when loading from a stored JSON file, please use Entry.from_json.

    Parameters
    ----------
    pdb_id : str
        4-letter pdb code
    release_date: str
        Entry release date
    oligomeric_state: str
        Oligomeric state decription
    determination_method: str
        Structure determination method, CRYSTAL, NMR, etc
    keywords: str
        Other keywords like "SUGAR-BINDING PROTEIN"
    pH: str
        pH used in structure determination
    resolution: float
        Strucure resolution
    chains: dict[str, Chain] = field(default_factory=dict)
        Chains dictionary with chain name mapped to chain object
    systems: dict[str, System] = field(default_factory=dict)
        System dictionary with system id mapped to system object
    covalent_bonds: dict[str, set[str]] = field(default_factory=dict)
        Covalent interation type mapped to all linkages of that type
    chain_to_seqres: dict[str, str] = field(default_factory=dict)
        Chain to sequence mapping

    Attributes
    ----------
    author_to_asym : dict[str, str]
        Mapping of author to asym chain ids

    Methods
    -------
    from_cif_file
    set_systems
    author_to_asym
    chains_for_alignment
    label_artifacts
    label_chains
    format_entry
    to_df
    save_systems

    Examples
    --------
    TODO: Add example
    """

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
                and len(s.interacting_protein_chains) <= max_protein_chains
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
        artifact_within_entry_threshold: int = 15,
        artifact_interacting_residue_count_threshold: int = 2,
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
        entry_info = get_entry_info(cif_data)
        per_chain = get_chain_external_mappings(cif_data)
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
            oligomeric_state=entry_info.get("entry_oligomeric_state", None),
            determination_method=entry_info.get("entry_determination_method", None),
            keywords=entry_info.get("entry_keywords", None),
            pH=entry_info.get("entry_pH", None),
            resolution=r,
            covalent_bonds=get_covalent_connections(cif_data),
            chain_to_seqres={c.name: c.string for c in seqres},
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
            data_dir = Path(save_folder).parent.parent
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
            entry_pdb_id = entry.pdb_id
            for ligand_chain in biounit_ligand_chains:
                ligand_instance, ligand_asym_id = ligand_chain.split(".")
                plip_output = run_plip_on_split_structure(
                    entry_pdb_id,
                    biounit,
                    ligand_chain,  # see func definition - consider adding entry.ligand_like_chains
                    plip_complex_threshold,
                )
                if plip_output is None:
                    continue
                (interactions, plip_chain_mapping) = plip_output
                data_dir = None
                if save_folder is not None:
                    data_dir = Path(save_folder).parent.parent
                ligand = Ligand.from_pli(
                    entry.pdb_id,
                    biounit_info.id,
                    biounit,
                    int(ligand_instance),
                    entry.chains[ligand_asym_id],
                    entry.ligand_like_chains,
                    interactions,
                    plip_chain_mapping,
                    interface_proximal_gaps,
                    entry.covalent_bonds,  # type: ignore
                    neighboring_residue_threshold,
                    neighboring_ligand_threshold,
                    data_dir=data_dir,
                )
                ligands[ligand.id] = ligand
            biounits[biounit_info.id] = biounit
        entry.set_systems(ligands)
        # NOTE: this is done at init now
        # entry.label_artifacts(
        #     # within_entry_threshold=artifact_within_entry_threshold,
        #     # interacting_residue_count_threshold=artifact_interacting_residue_count_threshold,
        # )
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
            if len(system.interacting_protein_chains):
                self.systems[system.id] = system

    @cached_property
    def author_to_asym(self) -> dict[str, str]:
        """
        Map autho chain id to asym id
        Parameters
        ----------
        self : Entry

        Returns
        -------
        dict[str, str]
        """
        return {
            c.auth_id: c.asym_id
            for c in self.chains.values()
            if c.chain_type == mol.CHAINTYPE_POLY_PEPTIDE_L
        }

    def chains_for_alignment(self, chain_type: str, aln_type: str) -> list[str]:
        """
        Map autho chain id to asym id
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
                for i_c in system.interacting_protein_chains
                if system.system_type == "holo"
            )
        elif chain_type == "apo":
            # ignore chains from systems that have more than one interaction
            # ignore_chains_from_systems = set(
            #     self.chains[i_c.split(".")[1]].auth_id
            #     for system in self.systems.values()
            #     for i_c in system.interacting_protein_chains
            #     if system.system_type == "artifact" and system.num_interactions > 1
            # )
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

    # def label_artifacts(
    #     self,
    #     within_entry_threshold: int = 15,
    #     interacting_residue_count_threshold: int = 2,
    # ) -> None:
    #     """
    #     Label artifacts

    #     Parameters
    #     ----------
    #     self : Entry
    #         Entry object
    #     within_entry_threshold : int
    #         Exclusion criteria for maximum \
    #         count of a particular ligand in entry
    #     interacting_residue_count_threshold : int
    #         Interacting residue count

    #     Returns
    #     -------
    #     TODO:Check to be sure what this function should not be returning anything
    #     None
    #     """
    #     artifact_ligand_chains = {}
    #     for system in self.systems:
    #         for ligand in self.systems[system].ligands:
    #             if ligand.in_artifact_list:
    #                 artifact_ligand_chains[ligand.asym_id] = ligand.ccd_code
    #     entry_ccd_to_count: dict[str, int] = defaultdict(int)
    #     for c in self.chains:
    #         if c in artifact_ligand_chains:
    #             entry_ccd_to_count[artifact_ligand_chains[c]] += 1

    #     for s in self.systems:
    #         for ligand in self.systems[s].ligands:
    #             ligand.identify_artifact()

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
                    [c.split(".")[1] for c in system.interacting_protein_chains]
                )
            ligand_chains.update([l.asym_id for l in system.ligands])
        for chain in self.chains:
            if chain not in ligand_chains and chain not in holo_chains:
                self.chains[chain].holo = False
            elif chain in holo_chains:
                self.chains[chain].holo = True

    def format_entry(self) -> dict[str, ty.Any]:
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

        data = {
            "entry_pdb_id": self.pdb_id,
            "entry_release_date": self.release_date,
            "entry_oligomeric_state": self.oligomeric_state,
            "entry_determination_method": self.determination_method,
            "entry_keywords": self.keywords,
            "entry_pH": self.pH,
            "entry_resolution": self.resolution,
        }
        if self.validation:
            validation = self.validation.to_dict()
            data.update(
                {
                    f"entry_{k}": v
                    for k, v in validation.items()
                    if f"entry_{k}" not in data
                }
            )
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
        entry_data = self.format_entry()
        for system in self.systems:
            system_data = self.systems[system].format_system(self.chains)
            for ligand in self.systems[system].ligands:
                ligand_data = ligand.format_ligand(self.chains)
                rows.append({**entry_data, **system_data, **ligand_data})
        return pd.DataFrame(rows)

    def iter_systems(
        self, max_protein_chains: int, max_ligand_chains: int
    ) -> ty.Iterator[tuple[str, System]]:
        for system_id, system in self.systems.items():
            if (
                system.system_type == "holo"
                and len(system.interacting_protein_chains) <= max_protein_chains
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
            save_folder_system = Path(save_folder) / system.id
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
            save_folder_system = Path(save_folder) / system.id
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
        thresholds: SystemValidationThresholds = SystemValidationThresholds(),
    ) -> None:
        if not validation_file.exists():
            LOG.error(f"set_validation: Validation file not found {validation_file}")
            return

        if self.determination_method != "X-RAY DIFFRACTION":
            LOG.warning(
                f"set_validation: Skipping validation for {self.pdb_id} as method is not X-RAY DIFFRACTION"
            )
            return
        try:
            doc = ValidationFactory(
                str(validation_file), mmcif_path=str(cif_file)
            ).getValidation()
            self.validation = EntryValidation.from_entry(doc)
            if self.validation:
                self.pass_criteria = self.validation.pass_criteria(
                    thresholds.entry_thresholds
                )
                if self.pass_criteria is None:
                    return
                for chain in self.chains:
                    self.chains[chain].set_validation(
                        doc, thresholds.residue_thresholds
                    )
                for system in self.systems:
                    self.systems[system].set_validation(
                        self.chains, thresholds, self.pass_criteria
                    )
        except Exception as e:
            LOG.error(
                f"set_validation: Error setting validation for {self.pdb_id}: {e}"
            )

    def add_ecod(self) -> None:
        """
        Add ECOD annotations to chains
        """
        global ECOD_DATA
        if ECOD_DATA is None:
            from plinder.data.pipeline.config import IngestConfig

            data_dir = Path(IngestConfig().plinder_dir)
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
                self.chains[chain_id].mappings["Kinase name"] = {
                    k: list(v) for k, v in mappings.items()
                }
