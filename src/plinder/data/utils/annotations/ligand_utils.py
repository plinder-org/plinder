# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import itertools
import logging
import typing as ty
from collections import Counter, defaultdict
from functools import cache, cached_property
from pathlib import Path

import numpy as np
import pandas as pd
from mmcif.api.PdbxContainers import DataContainer
from openbabel import pybel
from ost import io, mol
from ost.conop import GetDefaultLib
from pydantic import BaseModel, BeforeValidator, Field
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, AllChem, Crippen, rdMolDescriptors
from rdkit.Chem.rdchem import Mol, RWMol

from plinder.core.utils.config import get_config
from plinder.data.common.constants import BASE_DIR
from plinder.data.utils.annotations.extras import (
    get_ccd_smiles_dict,
    sort_ccd_codes,
)
from plinder.data.utils.annotations.interaction_utils import (
    extract_ligand_links_to_neighbouring_chains,
    get_plip_hash,
    pdbize,
    run_plip_on_split_structure,
)
from plinder.data.utils.annotations.protein_utils import Chain
from plinder.data.utils.annotations.rdkit_utils import set_smiles_from_ligand_ost

COMPOUND_LIB = GetDefaultLib()
PEPTIDE_TYPES = [
    mol.CHAINTYPE_POLY,
    mol.CHAINTYPE_POLY_PEPTIDE_D,
    mol.CHAINTYPE_POLY_PEPTIDE_L,
]

DNA_TYPES = [mol.CHAINTYPE_POLY_DN]

RNA_TYPES = [mol.CHAINTYPE_POLY_RN]

MIXED_NUCLEIC_ACID_TYPES = [mol.CHAINTYPE_POLY_DN_RN, mol.CHAINTYPE_POLY_PEPTIDE_DN_RN]

OLIGOSACCHARIDE_TYPES = [
    mol.CHAINTYPE_POLY_SAC_D,
    mol.CHAINTYPE_POLY_SAC_L,
    mol.CHAINTYPE_OLIGOSACCHARIDE,
    mol.CHAINTYPE_BRANCHED,
]

MACROCYCLE_TYPES = [mol.CHAINTYPE_MACROLIDE, mol.CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE]

NON_SMALL_MOL_LIG_TYPES = (
    PEPTIDE_TYPES
    + DNA_TYPES
    + RNA_TYPES
    + MIXED_NUCLEIC_ACID_TYPES
    + OLIGOSACCHARIDE_TYPES
    + MACROCYCLE_TYPES
)


def validate_chain_residue(obj: dict[str, ty.Any]) -> dict[str, ty.Any]:
    clean = {}
    for k, v in obj.items():
        key = tuple(k.split(",")) if isinstance(k, str) else k
        if isinstance(v, dict):
            clean[key] = validate_chain_residue(v)
        else:
            clean[key] = v
    return clean  # type: ignore


def lig_has_dummies(
    ligand_code: str,
    dummy_lig_list: list[str] = [
        "DUM",
        "UNX",
        "ASX",
        "GLX",
        "UNL",
        "UNK",
        "UPL",
        "DN",
        "N",
    ],
) -> bool:
    """Check for ccd codes containing dummy/unknown entries

    Args:
        ligand_code str: ligand CCD code
        dummy_lig_list (list, optional): list of ccd codes for unknown or dummy entries.
        Defaults to ['DUM', 'UNX', 'UNL', 'UNK', 'UPL', 'DN', 'N'].

    Returns:
        bool: if ligand considered as dummy and treated as artifact
    """
    # check for dummy list including composites, too!
    return len(set(ligand_code.split("-")).intersection(dummy_lig_list)) > 0


@cache
def get_ccd_synonyms(data_dir: Path) -> tuple[list[set[str]], dict[str, str]]:
    """Get Synonym dictonary for CCD SMILES
    CCD smiles dict from download_components_cif
    and get_ccd_smiles_dict

    Returns:
        Dict[str, str]: dictonary mapping synonymous CCD code to preferred one
    """
    from plinder.data.pipeline.io import download_components_cif

    ccd_lib_cifpath = download_components_cif(data_dir=data_dir)
    smidict = get_ccd_smiles_dict(ccd_lib_cifpath)
    ccd_df = pd.DataFrame.from_dict(smidict, orient="index").reset_index()
    # note: SMILES assumed to be CANONICALIZED by OE read by get_ccd_smiles_dict()
    ccd_df.columns = ["ccd_code", "SMILES"]
    # remove dummies
    ccd_df = ccd_df[~ccd_df["ccd_code"].apply(lig_has_dummies)]
    ccd_sets = ccd_df.groupby("SMILES").aggregate(list).reset_index()
    # ccd_sets
    ccd_sets["ccd_synonym_count"] = ccd_sets["ccd_code"].apply(lambda x: len(x))
    ccd_dups = ccd_sets[ccd_sets.ccd_synonym_count > 1].copy()
    list_of_synonym_sets = [set(x) for x in ccd_dups["ccd_code"].to_list()]
    # keep unique_ccd as sorted first entry
    ccd_dups["unique_ccd"] = ccd_dups["ccd_code"].apply(lambda x: sort_ccd_codes(x)[0])
    ccd_dups_exp = ccd_dups.explode("ccd_code")
    ccd_dups_exp.index = ccd_dups_exp["ccd_code"]
    ccd_synonym_dict = ccd_dups_exp["unique_ccd"].to_dict()
    return (list_of_synonym_sets, ccd_synonym_dict)


# lazy evaluate data fetches referenced as module globals
# TODO : clean this up and deduplicate extras with pipeline.io
COFACTORS = None
LIST_OF_CCD_SYNONYMS = None
CCD_SYNONYMS_DICT = None
# instantiate artifact list once and reuse variable
ARTIFACTS = None
KINASE_INHIBITORS = None
BINDING_AFFINITY = None


def add_missed_synonyms(current_set: set[str]) -> set[str]:
    assert LIST_OF_CCD_SYNONYMS is not None
    missed_synonyms = [
        x.difference(current_set)
        for x in LIST_OF_CCD_SYNONYMS
        if len(x.intersection(current_set)) > 0
    ]
    return set(itertools.chain(*missed_synonyms)).union(current_set)


def get_unique_ccd_longname(longname: str) -> str:
    assert CCD_SYNONYMS_DICT is not None

    if longname.startswith("PRD_"):
        return longname
    else:
        return "-".join([CCD_SYNONYMS_DICT.get(s, s) for s in longname.split("-")])


def get_ligand_chainid_comp_id_map(data: DataContainer) -> dict[str, set[str]]:
    atom_sites = data.getObj("atom_site")
    atom_site_columns = ["group_PDB", "label_comp_id", "label_asym_id"]
    if atom_sites is None:
        return {}

    chain_comp_id_map = defaultdict(set)
    for atom in atom_sites.getCombinationCountsWithConditions(
        atom_site_columns, [("group_PDB", "eq", "HETATM")]
    ):
        chain_comp_id_map[atom[2]].add(atom[1])
    return chain_comp_id_map


def get_bond_info(
    data: DataContainer, comp_ids: set[str]
) -> dict[str, list[tuple[str, str, str]]]:
    comp_bond_info = data.getObj("chem_comp_bond")
    if comp_bond_info is None:
        return {}
    comp_bond_cols = ["comp_id", "atom_id_1", "atom_id_2", "value_order"]
    bonds_dict = defaultdict(list)
    for bond in comp_bond_info.getCombinationCountsWithConditions(
        comp_bond_cols, [("comp_id", "in", comp_ids)]
    ):
        if bond[0] == "HOH":
            continue
        bonds_dict[bond[0]].append((bond[1], bond[2], bond[3]))
    return bonds_dict


def get_all_indices(lst: list[str], item: ty.Any) -> list[int]:
    """
    Get all indices of an item in a list

    Parameters
    ----------
    lst : lst
        The first parameter.

    Returns
    -------
    list
       list of indices of the item in the list
    """
    my_list = np.array(lst)
    return [int(i) for i in np.where(my_list == item)[0]]


def bond_pdb_order(value_order: str) -> Chem.rdchem.BondType:
    """
    Get rdkit bond type from pdb order

    Parameters
    ----------
    value_order : str

    Returns
    -------
    Chem.rdchem.BondType
    """
    if value_order.casefold() == "sing":
        return Chem.rdchem.BondType(1)
    if value_order.casefold() == "doub":
        return Chem.rdchem.BondType(2)
    if value_order.casefold() == "trip":
        return Chem.rdchem.BondType(3)
    return None


def get_rdkit_mol_from_pdb_block(
    pdb_block: str, bonds_dict: dict[str, list[tuple[str, str, str]]]
) -> str:
    mol = AllChem.MolFromPDBBlock(pdb_block)
    atoms_ids = [
        f"{atm.GetPDBResidueInfo().GetResidueName().strip()}"
        + f":{atm.GetPDBResidueInfo().GetName().strip()}"
        for atm in mol.GetAtoms()
    ]

    rw_mol = RWMol(mol)
    for comp_id, bonds in bonds_dict.items():
        for row in bonds:
            atom_1 = row[0]
            atom_2 = row[1]

            if (f"{comp_id}:{atom_1}" not in atoms_ids) | (
                f"{comp_id}:{atom_2}" not in atoms_ids
            ):
                # extra atom
                pass

            if atom_1.startswith("H") | atom_2.startswith("H"):
                # skip hydrogens
                pass
            else:
                try:
                    atom_1_ids = get_all_indices(atoms_ids, f"{comp_id}:{atom_1}")
                    atom_2_ids = get_all_indices(atoms_ids, f"{comp_id}:{atom_2}")
                    for atom_1_id, atom_2_id, order in zip(
                        atom_1_ids, atom_2_ids, np.repeat(row[2], len(atom_1_ids))
                    ):
                        bond_order = bond_pdb_order(order)
                        rw_mol.RemoveBond(int(atom_1_id), int(atom_2_id))
                        rw_mol.AddBond(int(atom_1_id), int(atom_2_id), bond_order)
                except ValueError:
                    print(
                        f"Error perceiving {atom_1} - {atom_2} bond in _chem_comp_bond"
                    )
                except RuntimeError:
                    print(f"Duplicit bond {atom_1} - {atom_2}")

    bonded_mol = rw_mol.GetMol()
    return str(Chem.MolToSmiles(bonded_mol))


def get_smiles_from_cif(
    data: DataContainer, ent: io.EntityHandle, polymer_cutoff: int = 20
) -> dict[str, str]:
    rdk_mols = {}
    chain_id_comp_id_map = get_ligand_chainid_comp_id_map(data)
    for chain_id, list_of_comp_ids in chain_id_comp_id_map.items():
        bonds_dict = get_bond_info(data, list_of_comp_ids)
        mol_ent = mol.CreateEntityFromView(
            ent.Select(f"chain='{chain_id}'"),
            True,
        )
        # If number of residue is withing cutoff range
        if len(mol_ent.residues) < polymer_cutoff:
            pdb_block = io.EntityToPDBStr(pdbize(ent, mol_ent)[0])
            rdk_mols[chain_id] = get_rdkit_mol_from_pdb_block(pdb_block, bonds_dict)
        # Skip water
        elif sum([res.name == "HOH" for res in mol_ent.residues]) > 0:
            continue
    return rdk_mols


def get_rdkit_mol_with_bond_order_from_cif(
    rdk_smiles_dict: dict[str, str], chain_id: str
) -> str:
    return rdk_smiles_dict.get(chain_id, "")


def get_chain_type(chain_type: str) -> str:
    """
    Get chain type

    Parameter
    ---------
    chain_type : str,
        ost chain type

    Return
    ------
    str
        ligand type
    """
    if chain_type == mol.CHAINTYPE_NON_POLY:
        return "SMALLMOLECULE"
    if chain_type in PEPTIDE_TYPES:
        return "PEPTIDE"
    elif chain_type in DNA_TYPES:
        return "DNA"
    elif chain_type in RNA_TYPES:
        return "RNA"
    elif chain_type in MIXED_NUCLEIC_ACID_TYPES:
        return "MIXED"
    elif chain_type in OLIGOSACCHARIDE_TYPES:
        return "SACCHARIDE"
    elif chain_type in MACROCYCLE_TYPES:
        return "MACROCYCLES"
    else:
        return "UNKNOWN"


@cache
def parse_cofactors(data_dir: Path) -> set[str]:
    """Download and parse cofactors.

    Returns
    -------
    Set[str]
        Set of cofactors

    """
    from plinder.data.pipeline.io import download_cofactors

    cofactors_json = download_cofactors(data_dir=data_dir)
    extra = {
        "Ascorbic acid": ["UU3"],
        "Coenzyme F420": ["6J4", "F42"],
        "Factor F430": ["F43", "M43"],
        "Pantetheine": ["PNY"],
        "Pantothenic acids": ["66S", "8Q1", "PAU"],
        "Nicotinamide": ["NCA"],
        "Adenosine nucleotides": [
            "A",
            "AMP",
            "ATP",
            "ADP",
        ],  # + ["ANP"],  # ANP is a mimic-inhibitor
        "Guanosine nucleotides": [
            "G",
            "GTP",
            "GDP",
            "GMP",
            "CPG",
            "G25",
            "5GP",
        ],  # + ["GNP", "GTN"],  # GNP/GTN is inhibitor
        "Cytidine nucleotides": ["C", "C5P", "C25", "CDP", "CTP"],
        "Thymidine nucleotides": ["T", "TMP", "DT", "TTP", "THM", "TYD"],
        "Uridine nucleotides": ["U", "DU", "U5P", "U25", "UMP", "UDP", "UTP"],
        "MIO": ["CRW"],
        "NAD": ["NAH"],
        "Glutathione": ["CYP"],
        "Biopterin": ["HBL", "BH4", "THB"],
        "Tetrahydrofolic acid": ["MEF"],
        "Lumazine": ["DLZ"],
        "Menaquinone": ["MQ8", "MQ9", "MQE", "7MQ"],
        "Heme": ["1CP", "CP3", "MMP", "UP2", "UP3"],
        "Methanopterin": ["H4M", "H4Z"],
        "Lipoamide": ["LPM"],
        "Ubiquinone": ["DCQ", "HQE", "PLQ"],
        "Pyridoxal": ["PXL", "UEG"],
        "Siderophores": ["488", "EB4", "SE8"],
        "Methanofuran": ["MFN"],
        "Vitamin A": ["BCR", "ECH", "EQ3", "RAW"],
        "Vitamin K1": ["PQN"],
        "CHLOROPHYLL and similar": [
            "CLA",
            "CHL",
            "CL0",
            "CL1",
            "CL2",
            "CL7",
            "BCB",
            "BCL",
            "07D",
            "G9R",
            "PEB",
            "PUB",
            "CYC",
            "BPH",
        ],
        # "Lipids": ["SPH"], # TODO: ?
        # "Sugars": ["NAG", "BCG", "GLC"], # TODO: more?
        #
    }
    cofactors = set()
    for c in cofactors_json:
        for c_list in cofactors_json[c]:
            cofactors |= set(c_list.get("cofactors", []))
    for c in extra:
        cofactors |= set(extra[c])

    # add missed synonyms
    cofactors = add_missed_synonyms(cofactors)

    return cofactors


@cache
def parse_artifacts() -> set[str]:
    """Get and parse artifacts
    Returns:
        set[str]: set[str]
    """
    artifact_log = BASE_DIR / "utils/annotations/static_files/artifacts_badlist.csv"
    with open(artifact_log, "r") as f:
        lines = f.readlines()
    artifacts = {l.strip() for l in lines if not l.startswith("#")}
    # add missed synonyms
    artifacts = add_missed_synonyms(artifacts)
    return artifacts


@cache
def parse_kinase_inhibitors(data_dir: Path) -> set[str]:
    from plinder.data.pipeline.io import download_kinase_data

    kinase_ligand_path = download_kinase_data(data_dir=data_dir)
    kinase_ligand_path = kinase_ligand_path.with_name("kinase_ligand_ccd_codes.parquet")
    kinase_ligand_df = pd.read_parquet(kinase_ligand_path)
    return set(kinase_ligand_df["PDB-code"])


@cache
def get_binding_affinity(data_dir: Path) -> ty.Any:
    from plinder.data.pipeline.io import download_affinity_data

    return download_affinity_data(data_dir=data_dir)


def get_num_resolved_heavy_atoms(resolved_smiles: str) -> int:
    obmol = pybel.readstring("smi", resolved_smiles)
    obmol.removeh()
    return len(obmol.atoms)


def get_len_of_longest_linear_hydrocarbon_linker(
    mol: Mol,
    max_count: int = 50,
    link_unit_smarts: str = "[#6D2R0]",
) -> int:
    """Estimate maximum linker length defined by link_unit_smarts, eg.
    unbranched hydrocarbons (default)

    Args:
        mol (Mol): RDKit molecule
        max_count (int, optional):
            Max count for linker. Defaults to 50.
        link_unit_smarts (str, optional):
            Linker unit defined by SMARTS. Defaults to "[#6D2R0]".

    Returns:
        int: maximum linker length defined by link_unit_smarts (default: unbranched hydrocarbon)
    """
    try:
        # needs ring info!
        Chem.SanitizeMol(
            mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
        )
        # length of longest hydrocarbon chain (excludes the ends and rings)
        for i in range(max_count):
            # chain_smarts = "[#6D2R0,#6D1R0]" * (i+1) # includes the ends
            chain_smarts = "~".join([link_unit_smarts] * (i + 1))
            if len(mol.GetSubstructMatches(Chem.MolFromSmarts(chain_smarts))) == 0:
                return i
        # TODO: what to do if fails or not found? now returns -1
        return -1
    except:
        return -1


def is_excluded_mol(
    smiles: str,
    min_C_threshold: int = 2,
    min_HA_threshold: int = 5,
    max_charge: int = 2,
    max_linear_hydrocarbon_linker: int = 12,
) -> bool:
    """Exclude some molecules by default as useless for druglikeness
    Uses OR logic for violating rules:
        - less than 2 carbon atoms
        - less than 5 non-hydrogen atoms
        - charge larger than +/- 2
        - unbranched hydrocarbon linker no longer than 12

    Args:
        smiles (str): molecule SMILES
        min_C_threshold (int, optional):
            Minimum carbon atom count. Defaults to 2.
        min_HA_threshold (int, optional):
            Minimum non-hydrogen atom count. Defaults to 5.
        max_charge (int, optional):
            Maximum allowed absolute charge. Defaults to 2.
        max_linear_hydrocarbon_linker (int, optional):
            Maximum allowed unbranched hydrocarbon linker. Defaults to 12.

    Returns:
        bool: should molecule be considered as artifact
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)

    # get heavy atom and carbon counts
    carbon = Chem.MolFromSmarts("[#6]")
    numC = len(mol.GetSubstructMatches(carbon))
    numHA = mol.GetNumHeavyAtoms()

    if numHA < min_HA_threshold or numC < min_C_threshold:
        return True

    # get formal charge
    charge = Chem.rdmolops.GetFormalCharge(mol)
    if abs(charge) > max_charge:
        return True
    elif (
        get_len_of_longest_linear_hydrocarbon_linker(mol)
        > max_linear_hydrocarbon_linker
    ):
        return True
    else:
        return False


def is_single_atom_or_ion(mol: Mol) -> bool:
    numHA = mol.GetNumHeavyAtoms()
    skip_single_elems = Chem.MolFromSmarts("[#6,#1,#0,#7,#8,#15,#16,#34,#52]")
    numCHNOPSetc = len(mol.GetSubstructMatches(skip_single_elems))
    return numHA == 1 and numCHNOPSetc == 0


def annotate_interface_gaps_per_chain(
    interface_proximal_gaps: dict[str, dict[tuple[str, str], dict[str, int]]],
    asym_id: str,
) -> tuple[int | None, ...]:
    try:
        ppi_atoms_within_4A_of_gap = sum(
            [
                v["interface_atom_gaps_4A"]
                for k, v in interface_proximal_gaps[
                    "ppi_interface_gap_annotation"
                ].items()
                if asym_id in k
            ]
        )
    except TypeError:
        ppi_atoms_within_4A_of_gap = None

    try:
        ppi_atoms_within_8A_of_gap = sum(
            [
                v["interface_atom_gaps_8A"]
                for k, v in interface_proximal_gaps[
                    "ppi_interface_gap_annotation"
                ].items()
                if asym_id in k
            ]
        )
    except TypeError:
        ppi_atoms_within_8A_of_gap = None
    try:
        num_missing_ppi_interface_residues = sum(
            [
                v["missing_interface_residues_4A"]
                for k, v in interface_proximal_gaps[
                    "ppi_interface_gap_annotation"
                ].items()
                if asym_id in k
            ]
        )
    except TypeError:
        num_missing_ppi_interface_residues = None
    try:
        pli_atoms_within_4A_of_gap = sum(
            [
                v["interface_atom_gaps_4A"]
                for k, v in interface_proximal_gaps[
                    "ligand_interface_gap_annotation"
                ].items()
                if asym_id in k
            ]
        )
    except TypeError:
        pli_atoms_within_4A_of_gap = None

    try:
        pli_atoms_within_8A_of_gap = sum(
            [
                v["interface_atom_gaps_8A"]
                for k, v in interface_proximal_gaps[
                    "ligand_interface_gap_annotation"
                ].items()
                if asym_id in k
            ]
        )
    except TypeError:
        pli_atoms_within_8A_of_gap = None
    try:
        num_missing_pli_interface_residues = sum(
            [
                v["missing_interface_residues_4A"]
                for k, v in interface_proximal_gaps[
                    "ligand_interface_gap_annotation"
                ].items()
                if asym_id in k
            ]
        )
    except TypeError:
        num_missing_pli_interface_residues = None

    return (
        ppi_atoms_within_4A_of_gap,
        ppi_atoms_within_8A_of_gap,
        num_missing_ppi_interface_residues,
        pli_atoms_within_4A_of_gap,
        pli_atoms_within_8A_of_gap,
        num_missing_pli_interface_residues,
    )


CrystalContacts = ty.Annotated[
    dict[tuple[str, int], set[int]],
    BeforeValidator(validate_chain_residue),
    Field(default_factory=dict),
]


class Ligand(BaseModel):
    pdb_id: str
    biounit_id: str
    asym_id: str
    instance: int
    ccd_code: str
    plip_type: str
    bird_id: str
    centroid: list[float]
    smiles: str
    resolved_smiles: str
    residue_numbers: list[int] = Field(default_factory=list)
    rdkit_canonical_smiles: str | None = None
    molecular_weight: float | None = None
    crippen_clogp: float | None = None
    num_rot_bonds: int | None = None
    num_hbd: int | None = None
    num_hba: int | None = None
    num_rings: int | None = None
    num_heavy_atoms: int | None = None
    is_covalent: bool = False
    covalent_linkages: set[str] = Field(default_factory=list)
    neighboring_residues: dict[str, list[int]] = Field(default_factory=dict)
    # neighboring_protein_chain_objects: list[Chain] = Field(default_factory=list)
    neighboring_ligands: list[str] = Field(default_factory=list)
    interacting_residues: dict[str, list[int]] = Field(default_factory=dict)
    interacting_ligands: list[str] = Field(default_factory=list)
    interactions: dict[str, dict[int, list[str]]] = Field(default_factory=dict)
    neighboring_residue_threshold: float = 6.0
    neighboring_ligand_threshold: float = 4.0
    num_neighboring_ppi_atoms_within_4A_of_gap: int | None = None
    num_neighboring_ppi_atoms_within_8A_of_gap: int | None = None
    num_missing_ppi_interface_residues: int | None = None
    num_pli_atoms_within_4A_of_gap: int | None = None
    num_pli_atoms_within_8A_of_gap: int | None = None
    num_missing_pli_interface_residues: int | None = None
    num_resolved_heavy_atoms: int | None = None
    num_unresolved_heavy_atoms: int | None = None
    tpsa: float | None = None
    qed: float | None = None
    is_ion: bool = False
    is_lipinski: bool = False
    is_fragment: bool = False
    is_oligo: bool = False
    is_cofactor: bool = False
    in_artifact_list: bool = False
    is_artifact: bool = False
    is_other: bool = False
    is_invalid: bool = False
    posebusters_result: dict[str, ty.Any] = Field(default_factory=dict)
    unique_ccd_code: str | None = None
    waters: dict[str, list[int]] = Field(default_factory=dict)
    crystal_contacts: CrystalContacts

    """
    This dataclass defines as system which included a protein-ligand complex
    and it's neighboring ligands and protein residues

    Attributes
    ----------
    pdb_id: str
        4-letter pdb code
    biounit_id: str
        Biounit id
    asym_id: str
        Ligand chain asym chain id
    instance: int
        Biouni instance ID
    ccd_code: str
        Ligand ccd code
    plip_type: str
        plip ligand type
    bird_id: str
        ligand bird id
    num_rot_bonds: int
        Number of rotatable bonds
    centroid: np.ndarray[float, Any]
        Ligand centroid
    smiles: str
        Ligand openbabel smiles
    resolved_smiles: str
        SMILES of only resolved ligand atoms
    rdkit_canonical_smiles: str
        rdkit smiles
    molecular_weight: float
        Molecular weight
    crippen_clogp: float
        Crippen clogp
    num_hbd: int
        Number of hydrogen bonds donor
    num_hba: int
        Number of hydrogen bind acceptor
    num_rings: int
        Number of rings
    num_heavy_atoms: int
        Number of heavy atoms
    tpsa: float
        Topological polar surface area
    qed: float
        ligand QED score
    is_covalent: bool = False
        Is ligand covalently bound to protein
    covalent_linkages: set[str] = field(default_factory=list)
        Residue tags of covalently linked residues in ligand
    neighboring_residues: dict[str, list[int]] = field(default_factory=dict)
        Dictionary of neighboring residues, with \
            {instance}.{chain} key and residue number value
    neighboring_ligands: list[str] = field(default_factory=list)
        List of neighboring ligands {instance}.{chain}
    interacting_residues: dict[str, list[int]] = field(default_factory=dict)
        Dictionary of interacting residues, with \
            {instance}.{chain} key and residue number value
    interacting_ligands: list[str] = field(default_factory=list)
        List of interacting ligands {instance}.{chain}
    interactions: dict[str, dict[int, list[str]]] = field(
                    default_factory=dict
                )
        Dictionary of {instance}.{chain} to residue number to list of PLIP hashes
    neighboring_residue_threshold: float = 6.0
        Maximum distance to consider protein residues neighboring
    neighboring_ligand_threshold: float = 4.0
        Maximum distance to consider ligands neighboring
    num_resolved_heavy_atoms: int | None = None
        Number of heavy atoms resolved in ligand
    num_unresolved_heavy_atoms: int | None = None
        Number of heavy atoms not resolved in ligand
    is_ion: bool = False
        Is ligand an ion
    is_lipinski: bool = False
        Does ligand satisfy Lipinski Ro5
    is_fragment: bool = False
        Does ligand satisfy fragment Ro3
    is_oligo: bool = False
        Does ligand match smarts for oligonucleotide, oligosacharide or oligopeptide
    is_cofactor: bool = False
        Is ligand in cofactor list
    in_artifact_list: bool = False
        Is ligand in artifact list
    is_artifact: bool = False
        Is ligand an artifact
    is_invalid: bool = False
        Is ligand classified as invalid
    is_other: bool = False
        Is ligand classified as other
    is_invalid: bool = False
        RDKit invalid
    posebusters_result: dict[str, ty.Any]
        Results from running posebusters with "re-dock"
    unique_ccd_code: str | None = None
        de-duplicated CCD code
    waters: dict[str, list[int]]
        residue numbers and chains of interacting waters
    residue_numbers: list[int] | None = None
        residue number(s) of ligand
    """

    def set_rdkit(self) -> None:
        try:
            rdkit_compatible_mol = Chem.MolFromSmiles(self.smiles)
            # TODO: Watch out for round trip issues reported in
            # https://github.com/rdkit/rdkit/issues/1740
            self.rdkit_canonical_smiles = Chem.CanonSmiles(self.smiles)
            self.molecular_weight = rdMolDescriptors.CalcExactMolWt(
                rdkit_compatible_mol
            )
            self.num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(
                rdkit_compatible_mol
            )
            self.num_hba = rdMolDescriptors.CalcNumHBA(rdkit_compatible_mol)
            self.num_hbd = rdMolDescriptors.CalcNumHBD(rdkit_compatible_mol)
            self.crippen_clogp = Crippen.MolLogP(rdkit_compatible_mol)
            self.num_rings = rdMolDescriptors.CalcNumRings(rdkit_compatible_mol)
            self.num_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(
                rdkit_compatible_mol
            )
            self.tpsa = rdMolDescriptors.CalcTPSA(rdkit_compatible_mol)
            self.qed = QED.qed(rdkit_compatible_mol)
            self.num_resolved_heavy_atoms = get_num_resolved_heavy_atoms(
                self.resolved_smiles
            )

            if self.num_heavy_atoms and self.num_resolved_heavy_atoms:
                self.num_unresolved_heavy_atoms = (
                    self.num_heavy_atoms - self.num_resolved_heavy_atoms
                )
            # classify ligand based on above molecule
            self.classify_ligand_type(rdkit_compatible_mol)

        except:
            logging.warning(f"Error in setting rdkit for {self.id}")
            self.is_invalid = True
            pass

    def classify_ligand_type(self, mol: Mol) -> None:
        """Get more granular classification of ligand.
        Use oligo smarts and lipinski rules to assign ligand classifications in
        addition to ligand classification obtained from PLIP

        Note
        ----
        Oligo smarts obtained from https://doi.org/10.1021/acs.jcim.3c01573

        Args:
            mol (Mol): RDKit compatible molecule
        """
        RDLogger.DisableLog("rdApp.*")
        oligo_smarts = {
            "oligopeptide": Chem.MolFromSmarts("C(=O)C[N;D2,D3]C(=[O;D1])CN"),
            "oligosaccharide": Chem.MolFromSmarts(
                "O-[C;R0,R1]-[C;R0,R1]-[O,S;R0;D2]-[C;R1]-[O;R1]"
            ),
            "oligonucleotide": Chem.MolFromSmarts("P(=O)([O-,OH])(OC[C;r5])O[C;r5]"),
        }

        if mol is not None:
            if is_single_atom_or_ion(mol):
                self.is_ion = True
            else:
                for sm in oligo_smarts:
                    try:
                        if mol.HasSubstructMatch(oligo_smarts[sm]):
                            # TODO: review - reduced to one class
                            self.is_oligo = True
                    except RuntimeError:
                        self.is_invalid = True
        else:
            self.is_invalid = True

        # Lipinski like Ro3 and Ro5
        if (
            self.molecular_weight
            and self.crippen_clogp
            and self.num_hbd
            and self.num_hba
        ):
            if (
                self.molecular_weight < 300
                and self.crippen_clogp < 3
                and self.num_hbd <= 3
                and self.num_hba <= 3
            ):
                self.is_fragment = True
                self.is_lipinski = True
            elif (
                self.molecular_weight < 500
                and self.crippen_clogp < 5
                and self.num_hbd <= 5
                and self.num_hba <= 10
            ):
                self.is_lipinski = True
        else:
            self.is_other = True

    def classify_ligand_by_attributes(self) -> None:
        """Get more granular classification of ligand.
        Use Lipinski rules to assign ligand classifications
        """

    @classmethod
    def from_pli(
        cls,
        pdb_id: str,
        biounit_id: str,
        biounit: mol.EntityHandle,
        ligand_instance: int,
        ligand_chain: Chain,
        residue_numbers: list[int],
        ligand_like_chains: dict[str, str],
        interface_proximal_gaps: dict[str, dict[tuple[str, str], dict[str, int]]],
        all_covalent_dict: dict[str, list[tuple[str, str]]],
        plip_complex_threshold: float = 10.0,
        neighboring_residue_threshold: float = 6.0,
        neighboring_ligand_threshold: float = 4.0,
        data_dir: ty.Optional[Path] = None,
    ) -> Ligand | None:
        """
        Load ligand object from protein-ligand interaction complex
        along with other information binding site information
        Parameters
        ----------
        cls : Ligand
            Ligand class
        pdb_id : str
            pdb code
        biounit_id : str
            Biounit id
        biounit : mol.EntityHandle
            Biounit openstructure mol.EntityHandle
        ligand_instance : int
            Ligand biounit instance
        ligand_chain : Chain
            Ligand chain object
        residue_numbers: list[int]
            List of residue numbers to use for ligand
        ligand_like_chains: dict[str, str]
            Chain: chain type for other ligand-like chains in the entry
        interface_proximal_gaps: dict[str, dict[tuple[str, str], dict[str, int]]]
            TODO: document
        all_covalent_dict : dict[str, list[tuple[str, str]]]
            All "covalent" residue in entry as defined by mmcif annotations.
            They types are separated by dictionary key and they include:
                "covale": actual covalent linkage
                "metalc": other dative bond interactions\
                     like metal-ligand dative bond
                "hydrogc": strong hydorogen bonding of nucleic acid
            For the purpose of covalent annotations, we selected "covale" for
            downstream processing.
        plip_complex_threshold: float = 10.0
            Maximum distance from ligand to residues to be
            included for pli calculations.
        neighboring_residue_threshold : float = 6.0
            Maximum distance for protein residue to be considered neighboring
        neighboring_ligand_threshold : float = 4.0
            Maximum distance for protein residue to be considered neighboring
        data_dir : Path, optional
            location of plinder root
        """
        global \
            COFACTORS, \
            ARTIFACTS, \
            LIST_OF_CCD_SYNONYMS, \
            CCD_SYNONYMS_DICT, \
            KINASE_INHIBITORS, \
            BINDING_AFFINITY
        if LIST_OF_CCD_SYNONYMS is None or CCD_SYNONYMS_DICT is None:
            if data_dir is None:
                raise ValueError(
                    "data_dir must be provided if CCD_SYNONYMS_DICT or LIST_OF_CCD_SYNONYMS is None"
                )
            LIST_OF_CCD_SYNONYMS, CCD_SYNONYMS_DICT = get_ccd_synonyms(data_dir)
        if COFACTORS is None:
            COFACTORS = parse_cofactors(data_dir)
        if ARTIFACTS is None:
            ARTIFACTS = parse_artifacts()
        if KINASE_INHIBITORS is None:
            KINASE_INHIBITORS = parse_kinase_inhibitors(data_dir)
        if BINDING_AFFINITY is None:
            BINDING_AFFINITY = get_binding_affinity(data_dir)
        ligand_instance_chain = f"{ligand_instance}.{ligand_chain.asym_id}"
        residue_selection = " or ".join(f"rnum={rnum}" for rnum in residue_numbers)
        ligand_selection = f"cname={mol.QueryQuoteName(ligand_instance_chain)} and ({residue_selection})"
        biounit_selection = mol.CreateEntityFromView(
            biounit.Select(
                f"{plip_complex_threshold} <> [{ligand_selection}]",
                mol.QueryFlag.MATCH_RESIDUES,
            ),
            True,
        )
        plip_output = run_plip_on_split_structure(
            biounit,
            biounit_selection,
            ligand_instance_chain,
        )
        if plip_output is None:
            return None
        (interactions, plip_chain_mapping) = plip_output
        ccd_code = "-".join(
            biounit.FindResidue(ligand_instance_chain, residue_number).name
            for residue_number in residue_numbers
        )
        ligand_ost_ent = mol.CreateEntityFromView(
            biounit.Select(ligand_selection), True
        )
        smiles = set_smiles_from_ligand_ost(ligand_ost_ent)
        ligand = cls(
            pdb_id=pdb_id,
            biounit_id=biounit_id,
            asym_id=ligand_chain.asym_id,
            instance=ligand_instance,
            ccd_code=ccd_code,
            plip_type=get_chain_type(
                ligand_chain.chain_type
            ),  # TODO: rename variable, no longer uses plip
            bird_id=list(ligand_chain.mappings.get("BIRD", {"": None}))[0],  # type: ignore
            centroid=list(ligand_ost_ent.GetCenterOfMass()),
            smiles=smiles,
            neighboring_residue_threshold=neighboring_residue_threshold,
            neighboring_ligand_threshold=neighboring_ligand_threshold,
            resolved_smiles=interactions.ligand.smiles,  # TODO: only thing left that depends on PLIP
            residue_numbers=residue_numbers,
        )

        neighboring_residue_selection = biounit.Select(
            f"{ligand.neighboring_residue_threshold} <> [{ligand_selection}]"
            + " and protein=True"
        )

        (
            ligand.num_neighboring_ppi_atoms_within_4A_of_gap,
            ligand.num_neighboring_ppi_atoms_within_8A_of_gap,
            ligand.num_missing_ppi_interface_residues,
            ligand.num_pli_atoms_within_4A_of_gap,
            ligand.num_pli_atoms_within_8A_of_gap,
            ligand.num_missing_pli_interface_residues,
        ) = annotate_interface_gaps_per_chain(
            interface_proximal_gaps, ligand_chain.asym_id
        )

        for residue in neighboring_residue_selection.residues:
            instance_chain = residue.chain.name
            if instance_chain == ligand.instance_chain:
                # this chain is considered a ligand, thus, skip!
                continue
            if instance_chain not in ligand.neighboring_residues:
                ligand.neighboring_residues[instance_chain] = []
            ligand.neighboring_residues[instance_chain].append(residue.number.num)

        neighboring_asym_ids = {
            ch.name.split(".")[-1]
            for ch in neighboring_residue_selection.chains
            if ch.name != ligand.instance_chain
        }

        # DONE: output should be sufficient for RFAA, eg. [(("A", "74", "ND2"), ("B", "1"), ("CW", "null"))]
        # see: https://github.com/baker-laboratory/RoseTTAFold-All-Atom?tab=readme-ov-file#predicting-covalently-modified-proteins

        ligand.covalent_linkages = extract_ligand_links_to_neighbouring_chains(
            all_covalent_dict, ligand.asym_id, neighboring_asym_ids, link_type="covale"
        )

        ligand.is_covalent = (
            len(ligand.covalent_linkages) > 0
        )  # TODO: Check to make sure we are catching all edge cases

        neighboring_ligand_selection = biounit.Select(
            f"{ligand.neighboring_ligand_threshold} <> [{ligand_selection}]"
        )

        ligand.neighboring_ligands = list(
            set(
                residue.chain.name
                for residue in neighboring_ligand_selection.residues
                if residue.chain.name != ligand.instance_chain
                and residue.chain.name.split(".")[1] in ligand_like_chains
            )
        )
        water_chains = set(
            c.name for c in biounit.chains if c.type == mol.CHAINTYPE_WATER
        )
        ligand.waters = defaultdict(list)
        for residue in interactions.interacting_res:
            residue_number, plip_chain = int(residue[:-1]), residue[-1]
            instance_chain = plip_chain_mapping[plip_chain]
            if instance_chain == ligand.instance_chain:
                continue
            if instance_chain in water_chains:
                ligand.waters[instance_chain].append(int(residue_number))
                continue
            if instance_chain.split(".")[1] in ligand_like_chains:
                ligand.interacting_ligands.append(instance_chain)
            else:
                if instance_chain not in ligand.interacting_residues:
                    ligand.interacting_residues[instance_chain] = []
                ligand.interacting_residues[instance_chain].append(int(residue_number))
        ligand.interactions, waters = get_plip_hash(
            interactions, ligand.instance_chain, plip_chain_mapping
        )
        for plip_chain, resnum in waters:
            ligand.waters[plip_chain_mapping[plip_chain]].append(resnum)
        # add rdkit properties and type assignments
        ligand.set_rdkit()
        # set is_artifact
        ligand.identify_artifact()
        # unique code parsing!
        ligand.unique_ccd_code = get_unique_ccd_longname(ligand.ccd_code)

        return ligand

    @property
    def selection(self) -> str:
        residue_selection = " or ".join(f"rnum={rnum}" for rnum in self.residue_numbers)
        ligand_selection = f"cname={mol.QueryQuoteName(self.instance_chain)}"
        if len(self.residue_numbers):
            ligand_selection += f"and ({residue_selection})"
        return ligand_selection

    @cached_property
    def interacting_protein_chains(self) -> list[str]:
        if self.is_artifact:
            return []
        else:
            return self.neighboring_protein_chains

    @cached_property
    def neighboring_protein_chains(self) -> list[str]:
        return list(sorted(self.neighboring_residues.keys()))

    @cached_property
    def num_interacting_residues(self) -> int:
        return sum(
            len(self.interacting_residues[chain]) for chain in self.interacting_residues
        )

    @cached_property
    def num_neighboring_residues(self) -> int:
        return sum(
            len(self.neighboring_residues[chain]) for chain in self.neighboring_residues
        )

    @cached_property
    def num_interactions(self) -> int:
        return sum(
            sum(len(i) for i in self.interactions[chain].values())
            for chain in self.interactions
        )

    @cached_property
    def num_unique_interactions(self) -> int:
        return sum(
            sum(len(set(i)) for i in self.interactions[chain].values())
            for chain in self.interactions
        )

    @cached_property
    def pocket_residues(self) -> dict[str, dict[int, str]]:
        residues: dict[str, dict[int, str]] = {}
        for chain in self.neighboring_residues:
            if chain not in residues:
                residues[chain] = {}
            for residue in self.neighboring_residues[chain]:
                residues[chain][residue] = "neighboring"
        for chain in self.interacting_residues:
            if chain not in residues:
                residues[chain] = {}
            for residue in self.interacting_residues[chain]:
                residues[chain][residue] = "interacting"
        return residues

    def get_pocket_residues_set(self) -> set[tuple[str, int]]:
        pocket_residues_set = set()
        for chain in self.pocket_residues:
            for residue_number in self.pocket_residues[chain]:
                pocket_residues_set.add((chain.split(".")[1], residue_number))
        return pocket_residues_set

    def set_crystal_contacts(
        self, crystal_contacts: dict[tuple[str, int], set[int]]
    ) -> None:
        # exclude contacts from neighboring residues in same biounit
        pocket_residues = self.get_pocket_residues_set()
        self.crystal_contacts = {
            x: y for x, y in crystal_contacts.items() if x not in pocket_residues
        }

    @cached_property
    def num_crystal_contacted_residues(self) -> int:
        return len(self.crystal_contacts)

    @cached_property
    def num_atoms_with_crystal_contacts(self) -> int:
        all_atoms = set()
        for x in self.crystal_contacts.values():
            all_atoms |= x
        return len(all_atoms)

    @cached_property
    def fraction_atoms_with_crystal_contacts(self) -> float | None:
        if self.num_heavy_atoms is None:
            return None
        return self.num_atoms_with_crystal_contacts / self.num_heavy_atoms

    @cached_property
    def num_pocket_residues(self) -> int:
        return sum([len(self.pocket_residues[chain]) for chain in self.pocket_residues])

    @cached_property
    def id(self) -> str:
        return "__".join([self.pdb_id, self.biounit_id, self.instance_chain])

    @cached_property
    def instance_chain(self) -> str:
        return f"{self.instance}.{self.asym_id}"

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

    @cached_property
    def is_kinase_inhibitor(self) -> bool:
        global KINASE_INHIBITORS
        if KINASE_INHIBITORS is None:
            data_dir = Path(get_config().data.plinder_dir)
            KINASE_INHIBITORS = parse_kinase_inhibitors(data_dir)
        return any(c in KINASE_INHIBITORS for c in self.ccd_code.split("-"))

    @cached_property
    def binding_affinity(self) -> float | None:
        """
        Get binding affinity
        """
        global BINDING_AFFINITY
        pdbid_ligid = f"{self.pdb_id}_{self.ccd_code}".upper()
        if BINDING_AFFINITY is None:
            data_dir = Path(get_config().data.plinder_dir)
            BINDING_AFFINITY = get_binding_affinity(data_dir)
        affinity = BINDING_AFFINITY.get(pdbid_ligid)
        if affinity is not None:
            return float(affinity)
        return None

    def identify_artifact(
        self,
    ) -> None:
        """
        Label artifacts, cofactors

        Parameters
        ----------
        self : Ligand
            Ligand object
        Returns
        -------
        dict[str, str]
        """
        assert COFACTORS is not None
        assert ARTIFACTS is not None
        if self.ccd_code in COFACTORS:
            self.is_cofactor = True
        if self.ccd_code in ARTIFACTS:
            self.in_artifact_list = True

        if self.is_ion:
            self.is_artifact = False
        elif self.in_artifact_list:
            self.is_artifact = True
        elif lig_has_dummies(self.ccd_code):
            # check for dummy list including composites, too!
            self.is_artifact = True
        elif is_excluded_mol(self.smiles):
            self.is_artifact = True
        else:
            self.is_artifact = False

    def format_chains(
        self,
        chain_type: str,
        chains: dict[str, Chain],
    ) -> dict[str, str]:
        """
        Format chains for pd.DataFrame

        Parameters
        ----------
        self : Ligand
            Ligand object
        chain_type: str
            Chain tyoe
        chains: dict[str, Chain]
            Chain id : chain mapping
        Returns
        -------
        dict[str, str]
        """
        if chain_type == "interacting_protein":
            sub_chains = self.interacting_protein_chains
        elif chain_type == "neighboring_protein":
            sub_chains = self.neighboring_protein_chains
        elif chain_type == "interacting_ligand":
            sub_chains = self.interacting_ligands
        elif chain_type == "neighboring_ligand":
            sub_chains = self.neighboring_ligands
        else:
            raise ValueError(f"chain_type={chain_type} not understood")
        sub_chains = [
            chains[instance_chain.split(".")[-1]].to_dict(instance_chain.split(".")[0])  # type: ignore
            for instance_chain in sub_chains
        ]
        data: dict[str, list[ty.Any]] = defaultdict(list)
        if len(sub_chains) == 0:
            return {}
        for sub_chain in sub_chains:
            for key in sub_chain:
                data[f"ligand_{chain_type}_chains{key}"].append(str(sub_chain[key]))  # type: ignore
        return {k: ";".join(v) for k, v in data.items()}

    def format_residues(
        self, residue_type: str, chains: dict[str, Chain]
    ) -> dict[str, str]:
        """
        Format residues for pd.DataFrame

        Parameters
        ----------
        self : Ligand
            Ligand object
        residue_type : str
            Chain tyoe
        chains : dict[str, Chain]
            Chain id : chain mapping

        Returns
        -------
        dict[str, str]
        """
        if residue_type == "interacting":
            residues = self.interacting_residues
        elif residue_type == "neighboring":
            residues = self.neighboring_residues
        res = []
        for instance_chain in residues:
            _, chain = instance_chain.split(".")
            for residue_number in residues[instance_chain]:
                res.append(
                    f"{instance_chain}_{residue_number}_{chains[chain].residues[residue_number].index}"
                )
        return {f"ligand_{residue_type}_residues": ";".join(res)}

    def format_interactions(self) -> dict[str, str]:
        """
        Format interactions for pd.DataFrame

        Parameters
        ----------
        self : Ligand

        Returns
        -------
        dict[str, str]

        """
        interactions: list[str] = []
        for chain in self.interactions:
            for residue in self.interactions[chain]:
                for interaction in self.interactions[chain][int(residue)]:
                    interactions.append(f"{chain}_{residue}_{interaction}")
        return {"ligand_interactions": ";".join(interactions)}

    def format_ligand(self, chains: dict[str, Chain]) -> dict[str, ty.Any]:
        """
        Format interactions for pd.DataFrame

        Parameters
        ----------
        self : Ligand
            Ligand object
        chains : dict[str, Chain]
            Chain dictionary

        Returns
        -------
        dict[str, Any]
        """
        data = {
            "ligand_id": self.id,
            "ligand_instance": self.instance,
            "ligand_asym_id": self.asym_id,
            "ligand_auth_id": chains[self.asym_id].auth_id,
            "ligand_ccd_code": self.ccd_code,
            "ligand_unique_ccd_code": self.unique_ccd_code,
            "ligand_plip_type": self.plip_type,
            "ligand_num_rot_bonds": self.num_rot_bonds,
            "ligand_centroid": ";".join(map(str, self.centroid)),
            "ligand_smiles": self.smiles,
            "ligand_resolved_smiles": self.resolved_smiles,
            "ligand_rdkit_canonical_smiles": self.rdkit_canonical_smiles,
            "ligand_molecular_weight": self.molecular_weight,
            "ligand_crippen_clogp": self.crippen_clogp,
            "ligand_num_hbd": self.num_hbd,
            "ligand_num_hba": self.num_hba,
            "ligand_num_rings": self.num_rings,
            "ligand_num_heavy_atoms": self.num_heavy_atoms,
            "ligand_tpsa": self.tpsa,
            "ligand_qed": self.qed,
            "ligand_is_covalent": self.is_covalent,
            "ligand_covalent_linkages": ";".join(self.covalent_linkages),
            "ligand_is_ion": self.is_ion,
            "ligand_is_lipinski": self.is_lipinski,
            "ligand_is_fragment": self.is_fragment,
            "ligand_is_oligo": self.is_oligo,
            "ligand_is_cofactor": self.is_cofactor,
            "ligand_is_other": self.is_other,
            "ligand_is_invalid": self.is_invalid,
            "ligand_in_artifact_list": self.in_artifact_list,
            "ligand_is_artifact": self.is_artifact,
            "ligand_num_interacting_residues": self.num_interacting_residues,
            "ligand_num_neighboring_residues": self.num_neighboring_residues,
            "ligand_num_interactions": self.num_interactions,
            "ligand_num_missing_pli_interface_residues": self.num_missing_pli_interface_residues,
            "ligand_num_pli_atoms_within_4A_of_gap": self.num_pli_atoms_within_4A_of_gap,
            "ligand_num_pli_atoms_within_8A_of_gap": self.num_pli_atoms_within_8A_of_gap,
            "ligand_num_missing_neighboring_ppi_residues": self.num_missing_ppi_interface_residues,
            "ligand_num_neighboring_ppi_atoms_within_8A_of_gap": self.num_neighboring_ppi_atoms_within_8A_of_gap,
            "ligand_num_neighboring_ppi_atoms_within_4A_of_gap": self.num_neighboring_ppi_atoms_within_4A_of_gap,
            "ligand_num_resolved_heavy_atoms": self.num_resolved_heavy_atoms,
            "ligand_num_unresolved_heavy_atoms": self.num_unresolved_heavy_atoms,
            "ligand_is_kinase_inhibitor": self.is_kinase_inhibitor,
            "ligand_num_atoms_with_crystal_contacts": self.num_atoms_with_crystal_contacts,
            "ligand_fraction_atoms_with_crystal_contacts": self.fraction_atoms_with_crystal_contacts,
            "ligand_num_crystal_contacted_residues": self.num_crystal_contacted_residues,
            "ligand_binding_affinity": self.binding_affinity,
        }
        if self.posebusters_result is not None:
            for k in self.posebusters_result:
                data[f"ligand_{k}"] = self.posebusters_result[k]

        data.update(self.format_interactions())
        for chain_type in [
            "interacting_protein",
            "neighboring_protein",
            "interacting_ligand",
            "neighboring_ligand",
        ]:
            data.update(self.format_chains(chain_type, chains))
        for residue_type in ["interacting", "neighboring"]:
            data.update(self.format_residues(residue_type, chains))
        return data
