# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
from collections import defaultdict
from functools import cached_property
from pathlib import Path
from typing import Any

from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxReader import PdbxReader
from ost import conop, io, mol
from PDBValidation.Validation import PDBValidation
from pydantic import ConfigDict, Field

from plinder.data.utils.annotations.get_ligand_validation import (
    ResidueListValidation,
    ResidueValidation,
    ResidueValidationThresholds,
)
from plinder.data.utils.annotations.utils import DocBaseModel

NON_SMALL_MOL_LIG_TYPES = [
    mol.CHAINTYPE_POLY,
    mol.CHAINTYPE_POLY_PEPTIDE_D,
    mol.CHAINTYPE_POLY_PEPTIDE_L,
    mol.CHAINTYPE_POLY_PEPTIDE_D,
    mol.CHAINTYPE_POLY_PEPTIDE_L,
    mol.CHAINTYPE_POLY_DN,
    mol.CHAINTYPE_POLY_RN,
    mol.CHAINTYPE_POLY_SAC_D,
    mol.CHAINTYPE_POLY_SAC_L,
    mol.CHAINTYPE_POLY_DN_RN,
    mol.CHAINTYPE_MACROLIDE,
    mol.CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE,
    mol.CHAINTYPE_POLY_PEPTIDE_DN_RN,
    mol.CHAINTYPE_BRANCHED,
    mol.CHAINTYPE_OLIGOSACCHARIDE,
    mol.CHAINTYPE_N_CHAINTYPES,
]


def read_mmcif_container(mmcif_filename: Path) -> DataContainer:
    """Parse mmcif file with PDBxReader

    Parameters
    ----------
    mmcif_filename : Path
    Returns
    -------
    DataContainer

    """
    data: list[DataContainer] = []
    if mmcif_filename.suffix == ".gz":
        with gzip.open(str(mmcif_filename), "rt", encoding="utf-8") as f:
            prd = PdbxReader(f)
            prd.read(data)
    else:
        prd = PdbxReader(mmcif_filename)
        prd.read(data)
    return data[0]


def get_entry_info(data: DataContainer) -> dict[str, str | float | None]:
    """Get entry-level information from DataContainer

    Parameters
    ----------
    data : DataContainer
        Data container fot mmcif attributes
    Returns
    -------
    dict[str, str | float | None]
        Dictionary of entry-level information

    """
    entry_info = {}
    mappings = [
        ("entry_oligomeric_state", "pdbx_struct_assembly", "oligomeric_details"),
        ("entry_determination_method", "exptl", "method"),
        ("entry_keywords", "struct_keywords", "pdbx_keywords"),
        ("entry_pH", "exptl_crystal_grow", "pH"),
    ]
    for key, obj_name, attr_name in mappings:
        x = data.getObj(obj_name)
        if x is not None:
            entry_info[key] = x.getValueOrDefault(attr_name)
    resolution_options = [
        ("refine", "ls_d_res_high"),
        # ("em_3d_reconstruction", "resolution"), # TODO: add this back for next annotation rerun
    ]
    resolution = None
    for obj_name, attr_name in resolution_options:
        x = data.getObj(obj_name)
        if x is not None:
            r = x.getValueOrDefault(attr_name)
            if r is not None:
                resolution = r
                break
    entry_info["entry_resolution"] = resolution
    return entry_info


def get_chain_external_mappings(
    data: DataContainer
) -> dict[str, dict[str, dict[str, list[tuple[str, str] | None]]]]:
    """Get additional metadata directory from nextgen mmcif

    Parameters
    ----------
    cif_file : Path
        Next-gen mmcif file

    Returns
    -------
    Tuple[Dict[Any, Any], Dict[Any, Any]]
    """
    per_chain: dict[str, dict[str, dict[str, set[tuple[str, str] | None]]]] = {}

    # SIFTS mapping
    xref = data.getObj("pdbx_sifts_xref_db_segments")
    if xref is not None:
        columns = xref.getAttributeList()
        for a in xref:
            a = dict(zip(columns, a))
            if a["asym_id"] not in per_chain:
                per_chain[a["asym_id"]] = defaultdict(lambda: defaultdict(set))
            per_chain[a["asym_id"]][a["xref_db"]][a["xref_db_acc"]].add(
                (a["seq_id_start"], a["seq_id_end"])
            )

    # UniProt mapping
    uniprot = data.getObj("pdbx_sifts_unp_segments")
    if uniprot is not None:
        columns = uniprot.getAttributeList()
        for a in uniprot:
            a = dict(zip(columns, a))
            if a["asym_id"] not in per_chain:
                per_chain[a["asym_id"]] = defaultdict(lambda: defaultdict(set))
            per_chain[a["asym_id"]]["UniProt"][a["unp_acc"]].add(
                (a["seq_id_start"], a["seq_id_end"])
            )

    # BIRD entries with PRD codes: https://www.wwpdb.org/data/bird
    pdbx_molecule = data.getObj("pdbx_molecule")
    # see: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/pdbx_molecule.html
    if pdbx_molecule is not None:
        columns = pdbx_molecule.getAttributeList()
        for a in pdbx_molecule:
            a = dict(zip(columns, a))
            if a["asym_id"] not in per_chain:
                per_chain[a["asym_id"]] = defaultdict(lambda: defaultdict(set))
            per_chain[a["asym_id"]]["BIRD"][f"{a['asym_id']}"].add(None)
    per_chain_list: dict[str, dict[str, dict[str, list[tuple[str, str] | None]]]] = {}
    for chain in per_chain:
        per_chain_list[chain] = {}
        for mapping in per_chain[chain]:
            per_chain_list[chain][mapping] = {
                k: list(v) for k, v in per_chain[chain][mapping].items()
            }
    return per_chain_list


def detect_ligand_chains(
    entity: Any,
    entry: Any,
    min_polymer_size: int = 10,
    max_non_small_mol_ligand_length: int = 20,
) -> dict[str, str]:
    """
    Note
    ----
    entity is the first element of the tuple returned by ost.io.LoadMMCIF
    entry is an Entry object that contains appropriate mappings
    """
    ligand_chains = dict()
    for chain in entity.chains:
        if chain.type == mol.CHAINTYPE_WATER:
            continue
        # classifying by polymer length and annotations
        chain_length = len(chain.residues)
        bird_id = list(entry.chains[chain.name].mappings.get("BIRD", {"": None}))[0]
        uniprot_id = list(entry.chains[chain.name].mappings.get("UniProt", {"": None}))[
            0
        ]
        # TODO: Let's revisit this at some point, but I think this logic is too
        # complicated, could be simplified.
        if (
            # chain has PRD id based on BIRD annotation:
            # https://www.wwpdb.org/data/bird
            # thus can be considered as ligand!
            bird_id
        ) or (
            # short/medium synthetic peptides that do not map to UniProt are considered ligand
            chain.is_polypeptide
            and chain_length <= max_non_small_mol_ligand_length
            and not uniprot_id
        ):
            ligand_chains[chain.name] = str(chain.type)

        elif (
            # Polymers that do not fall for the exception above and
            # are longer that certain length to be considered as ligands
            # TODO: these are all polymer chains, we might want to separate into protein chains and other
            (chain.is_polypeptide and chain_length >= min_polymer_size)
            or (chain.is_polynucleotide and chain_length >= min_polymer_size)
            or (
                (chain.is_oligosaccharide or chain.is_polysaccharide)
                and chain_length >= min_polymer_size
            )
            or (chain.type == mol.CHAINTYPE_POLY and chain_length >= min_polymer_size)
        ):
            # these are excluded!
            continue
        else:
            # the rest is guessed as being ligand
            ligand_chains[chain.name] = str(chain.type)
    return ligand_chains


class Residue(DocBaseModel):
    chain: str
    index: int
    number: int
    auth_number: str
    one_letter_code: str
    name: str
    chem_type: str
    validation: ResidueValidation | None = None
    """
    This dataclass defines as system which included a protein-ligand complex
    and it's neighboring ligands and protein residues

    Parameters
    ----------
    chain : str
        chain name
    index : int
        residue index
    number : int
        residue number
    one_letter_code : str
        residue one-letter code
    name : str
        residue name
    chem_type: mol.ChemType
        residue chemical type

    Attributes
    ----------
    is_ptm: bool
        Does residue have post translational modification

    Examples
    --------
    TODO: Add example
    """

    @cached_property
    def is_ptm(self) -> bool:
        """
        Does the residue have a post translational modification.
        """
        return (
            mol.ChemType(self.chem_type).IsAminoAcid()
            and self.name not in conop.STANDARD_AMINOACIDS
        )


class Chain(DocBaseModel):
    asym_id: str = Field(description="Chain asymmetric id")
    auth_id: str = Field(description="Chain author id")
    entity_id: str = Field(description="Chain entity id")
    chain_type_str: str = Field(
        description="__Chain type string representation as defined https://openstructure.org/docs/2.8/mol/base/entity/#ost.mol.ChainType"
    )
    residues: dict[int, Residue] = Field(
        description="__Dictionary of residues in chain with keys as residue number"
    )
    length: int = Field(description="SEQRES length")
    num_unresolved_residues: int = Field(
        description="Number of unresolved residues (SEQRES length - len(residues))"
    )
    mappings: dict[str, dict[str, list[tuple[str, str] | None]]] = Field(
        default_factory=dict,
        description="__Mapping of metadata associated with chain with keys as chain asym id",
    )
    holo: bool = Field(
        default=True, description="__Is the chain part of a holo system or not"
    )
    validation: ResidueListValidation | None = Field(
        default=None,
        description="__Crystal validation information for the residues in the chain",
    )

    # Added this because pydantic doesn't know how to validate mol.ChainType
    model_config = ConfigDict(
        arbitrary_types_allowed=True,
    )

    @classmethod
    def from_ost_chain(
        cls, chain: mol.ChainHandle, info: io.MMCifInfo, length: int
    ) -> Chain:
        """Load Chain from ost Chain.

        Parameters
        ----------
        cls : Chain
            Chain class
        chain: mol.ChainHandle :
            Openstructure mol.ChainHandle
        info: io.MMCifInfo
            Openstructure io.MMCifInfo
        length: int
            SEQRES length

        Returns
        -------
        Chain
        """
        residues = {
            residue.number.num: Residue(
                chain=chain.name,
                index=residue_index,
                number=residue.number.num,
                auth_number=residue.GetStringProp("pdb_auth_resnum"),
                one_letter_code=residue.one_letter_code,
                name=residue.name,
                chem_type=str(residue.chem_type),
            )
            for residue_index, residue in enumerate(chain.residues)
        }
        auth_id = ""
        if "." not in chain.name:
            auth_id = chain.GetStringProp("pdb_auth_chain_name")
        return cls(
            asym_id=chain.name,
            auth_id=auth_id,
            entity_id=info.GetMMCifEntityIdTr(chain.name),
            chain_type_str=str(mol.StringFromChainType(chain.type)),
            residues=residues,
            length=length,
            num_unresolved_residues=length - len(residues),
        )

    @cached_property
    def chain_type(self) -> mol.ChainType:
        """
        __Chain type as defined https://openstructure.org/docs/2.8/mol/base/entity/#ost.mol.ChainType
        """
        return mol.ChainTypeFromString(self.chain_type_str)

    @cached_property
    def residue_index_to_number(self) -> dict[int, int]:
        """
        __Dictionary of residue index to residue number
        """
        return {self.residues[r].index: r for r in self.residues}

    def format(self, instance: int) -> dict[str, Any]:
        """
        Format chain as a dictionary
        """
        data: dict[str, Any] = {
            "asym_id": f"{instance}.{self.asym_id}",
            "auth_id": self.auth_id,
            "entity_id": self.entity_id,
            "length": self.length,
            "num_unresolved_residues": self.num_unresolved_residues,
        }
        # data.update(self.mappings)
        if self.validation is not None:
            data.update(self.validation.format())
        return data

    def set_validation(
        self, doc: PDBValidation, residue_thresholds: ResidueValidationThresholds
    ) -> None:
        for residue in self.residues:
            self.residues[residue].validation = ResidueValidation.from_residue(
                self.asym_id, residue, self.entity_id, doc
            )
        validations = []
        for r in self.residues.values():
            if r.validation is not None:
                validations.append(r.validation)
        self.validation = ResidueListValidation.from_residues(
            validations, residue_thresholds
        )
