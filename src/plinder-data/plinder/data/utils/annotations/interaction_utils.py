from __future__ import annotations

import gzip
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
from biotite.structure.io.pdbx import PDBxFile
from mmcif.api.PdbxContainers import DataContainer
from ost import io, mol
from plip.structure.preparation import PDBComplex, PLInteraction
from plip.basic.supplemental import whichchain, whichresnumber

from plinder.data.common.log import setup_logger

log = setup_logger(__name__)

INTERACTION_TYPES = [
    "hbonds_ldon",
    "hbonds_pdon",
    "hydrophobic_contacts",
    "pication_laro",
    "pication_paro",
    "halogen_bonds",
    "pistacking",
    "water_bridges",
    "saltbridge_lneg",
    "saltbridge_pneg",
    "metal_complexes",
]

# Define available names for chains in PDB format
PDB_AVAILABLE_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
PDB_AVAILABLE_CHAINS += PDB_AVAILABLE_CHAINS.lower() + "0123456789"


def convert_category(category: dict[str, np.ndarray[int, Any]]) -> dict[Any, Any]:
    """
    Convert a PDBx/mmCIF category to a dictionary indexed by sequential ids.
    with keys and values taken from the original value arrays.

    Parameters
    ----------
    category: Dict[str, np.ndarray[int, Any]])
        PDBx/mmCIF category to a dictionary


    Returns
    -------
    Dict[Any, Any]
        Dictionary indexed by sequential ids.
    """
    category_dict: dict[Any, Any] = {}
    if category is not None:
        for i in range(len(category[list(category.keys())[0]])):
            category_dict[i] = {}
            for key, value in category.items():
                category_dict[i][key] = value[i]
    return category_dict


def read_mmcif_file(mmcif_filename: Path) -> PDBxFile:
    """
    Read a PDBx/mmCIF file.
    """
    if mmcif_filename.suffix == ".gz":
        with gzip.open(mmcif_filename, "rt") as mmcif_file:
            pdbx_file = PDBxFile.read(mmcif_file)
    else:
        pdbx_file = PDBxFile.read(mmcif_filename)
    return pdbx_file


def get_covalent_connections(data: DataContainer) -> dict[str, list[dict[str, str]]]:
    """
    Get covalent connections from any mmcif file with
    _struct_conn. attribute

    Parameters
    ----------
    mmcif_file : Path
        mmcif file with _struct_conn. attribute

    Returns
    -------
    Dict[str, List[Set[str]]]
        Mapping of covalent residues
    """

    graph_dict = defaultdict(list)
    nucleobase_list = ["A", "C", "U", "G", "DA", "DC", "DG", "DT", "PSU"]
    node_lookup = defaultdict(list)
    to_extract = [
        "ptnr1_label_asym_id",
        "ptnr2_label_asym_id",
        "ptnr1_label_comp_id",
        "ptnr2_label_comp_id",
        "ptnr1_label_atom_id",
        "ptnr2_label_atom_id",
        "conn_type_id",
    ]
    cons = data.getObj("struct_conn")
    if cons is None:
        return {}
    for con in cons.getCombinationCountsWithConditions(
        to_extract, [("conn_type_id", "in", ["covale", "metalc", "hydrog"])]
    ):
        con = dict(zip(to_extract, con))
        if con["conn_type_id"] == "covale":
            node1 = f'{con["ptnr1_label_comp_id"]}:' + f'{con["ptnr1_label_asym_id"]}'
            node2 = f'{con["ptnr2_label_comp_id"]}:' + f'{con["ptnr2_label_asym_id"]}'
            node_lookup[node1].append(node1 + f':{con["ptnr1_label_atom_id"]}')
            node_lookup[node2].append(node2 + f':{con["ptnr2_label_atom_id"]}')
            graph_dict["covale"].append((node1, node2))
        elif con["conn_type_id"] == "metalc":
            node1 = f'{con["ptnr1_label_comp_id"]}:' + f'{con["ptnr1_label_asym_id"]}'
            node2 = f'{con["ptnr2_label_comp_id"]}:' + f'{con["ptnr2_label_asym_id"]}'
            node_lookup[node1].append(node1 + f':{con["ptnr1_label_atom_id"]}')
            node_lookup[node2].append(node2 + f':{con["ptnr2_label_atom_id"]}')
            graph_dict["metalc"].append((node1, node2))
        elif con["conn_type_id"] == "hydrog":
            node1 = f'{con["ptnr1_label_comp_id"]}:' + f'{con["ptnr1_label_asym_id"]}'
            node2 = f'{con["ptnr2_label_comp_id"]}:' + f'{con["ptnr2_label_asym_id"]}'
            node_lookup[node1].append(node1)
            node_lookup[node2].append(node2)
            if con["ptnr1_label_comp_id"].strip() in nucleobase_list:
                graph_dict["hydrog"].append((node1, node2))

    output_graph = {k: nx.from_edgelist(set(v)) for k, v in graph_dict.items()}
    components = defaultdict(list)
    for k, v in output_graph.items():
        for com in nx.connected_components(v):
            components[k].append(com)
    return {
        cov_type: [{j: i for i in d for j in node_lookup[i]} for d in data]
        for cov_type, data in components.items()
    }


def merge_covalent_ligands_with_gaps(
    list_of_linked_ligands: list[set[str]],
) -> list[set[str]]:
    """
    Merge ligands artifical considered as not linked
    because of gap. The logic is, if it's the same chain,
    then it's linked.

    Parameters
    ----------
    mmcif_file :  list[set[str]]
        List of set of linked residues tags.
        Residue tags formated as <resname>:<chainid>:<atomid>

    Returns
    -------
    list[set[str]]
    """
    merged = defaultdict(set)
    for idx, links1 in enumerate(list_of_linked_ligands):
        links1_chains = [res.split(":")[1] for res in links1]
        links1_chain_dict = dict(Counter(links1_chains))
        links1_main_chain = max(links1_chain_dict, key=links1_chain_dict.get)  # type: ignore
        for links2 in list_of_linked_ligands[: idx + 1]:
            links2_chains = {res.split(":")[1] for res in links2}
            links2_chain_dict = dict(Counter(links1_chains))
            links2_main_chain = max(links2_chain_dict, key=links2_chain_dict.get)  # type: ignore
            if len(set(links1_chains).intersection(set(links2_chains))) > 0:
                merged[links1_main_chain].update(set(links1).union(set(links2)))
            else:
                merged[links1_main_chain].update(set(links1))
                merged[links2_main_chain].update(set(links2))
    return list(merged.values())


def extract_other_covalent_subunit(
    all_linkages_dict: dict[str, list[set[str]]],
    cov_units: set[str],
    link_type: str = "covale",
) -> list[str]:
    """
    Lookup other linked subunit of a multipart ligand given the
    residue name and chain id of one of it's subunits.

    Parameters
    ----------
    all_linkages_dict :  Dict[str, List[Set[str]]]
        Mapping of all covalent linkages
    cov_units: set[str],
        Set of covalently linked residue tags
    link_type="covale : str
        linkage type, could be "covale", "metalc", "hydrog"
        representing covalent, metal adduct or \
        hydrogen bonding (mostly nucleic acid related)

    Returns
    -------
    Dict[str, List[Set[str]]]
        Mapping of covalent residues
    """
    if link_type not in all_linkages_dict:
        return []
    try:
        specific_linkages = all_linkages_dict[link_type]
        specific_linkages = merge_covalent_ligands_with_gaps(specific_linkages)
        linkage_no_atoms = [
            {":".join(link.split(":")[:2]) for link in links}
            for links in specific_linkages
        ]
        fragments = [
            specific_linkages[idx]
            for idx, link in enumerate(linkage_no_atoms)
            if len(cov_units.intersection(link))
        ]
        if len(fragments):
            return list(fragments[0])
        else:
            return []
    except (IndexError, KeyError) as e:
        log.error(
            f"extract_other_covalent_subunit: Could not find {link_type} linkages: {e}"
        )
        return []


def run_plip(biounit_pdbized: mol.EntityHandle) -> PDBComplex:
    """Load pdbized biounit and run plip analysis.

    Parameters
    ----------
    biounit_pdbized : mol.EntityHandle
        pdbized biounit

    Returns
    -------
    PDBComplex
        Complex interaction object with all plip related annotation computed
    """

    complex_obj = PDBComplex()
    complex_obj.load_pdb(io.EntityToPDBStr(biounit_pdbized).strip(), as_string=True)
    complex_obj.analyze()
    return complex_obj


def pdbize(
    full_biounit: mol.EntityHandle, entity: mol.EntityHandle
) -> mol.EntityHandle:
    """PDBize entity chains

    Parameters
    ----------
    entity : mol.EntityHandle
        Entity handle
    Returns
    -------
    mol.EntityHandle
        PDBized entity handle
    """
    # Intermediate renaming step
    intermediate_names = {}
    edi = entity.EditXCS(mol.BUFFERED_EDIT)
    for i, chain in enumerate(entity.GetChainList()):
        intermediate_names[f"T{i}"] = chain.name
        edi.RenameChain(chain, f"T{i}")
    edi.UpdateICS()
    # Final renaming step
    chain_index = 0
    name_mapping = {}
    for chain in entity.GetChainList():
        original_name = intermediate_names[chain.name]
        original_chain = full_biounit.FindChain(original_name)
        if chain_index >= len(PDB_AVAILABLE_CHAINS):
            raise ValueError(f"Too many chains ({chain_index}) in entity")
        final_name = PDB_AVAILABLE_CHAINS[chain_index]
        chain_index += 1
        edi.RenameChain(chain, final_name)
        edi.SetChainDescription(chain, original_chain.description)
        edi.SetChainType(chain, original_chain.type)
        name_mapping[original_name] = final_name
    edi.UpdateICS()
    return entity, name_mapping


def run_plip_on_split_structure(
    entry_pdb_id: str,
    biounit: mol.EntityHandle,
    ligand_chain: str,
    plip_threshold: float = 10,
) -> tuple[PLInteraction, dict[str, str]] | None:
    """Split structure into small PLI complex by ligand

    For every ligand, create a smaller complex for faster plip\
    process and deal with cases where plip ignores small molecule
    ligands in the presence of peptides

    Parameters
    ----------
    entry_pdb_id : str
        pdb entry
    biounit : mol.EntityHandle
        biounit
    ligand_chain : int
        {instance}.{chain} of ligand
    plip_threshold: int = 10
        distance from ligand chain for atoms to be included in plip calculations

    Returns
    -------
    Tuple[PlipLigand, PDBComplex, dict[str, str]] | None
        PLIP ligand, PLIP complex object, mapping of original chain to plip chain
    """
    from plip.basic import config

    config.biolip_list = []
    split_structure, chain_mapping = pdbize(
        biounit,
        mol.CreateEntityFromView(
            biounit.Select(
                f"{plip_threshold} <> [cname='{ligand_chain}']",
                mol.QueryFlag.MATCH_RESIDUES,
            ),
            True,
        ),
    )
    ligand_plip_chain = chain_mapping[ligand_chain]
    config.PEPTIDES = (
        [ligand_plip_chain] if biounit.FindChain(ligand_chain).is_polymer else []
    )
    # TODO: review this - we might be treating some peptidic ligands as protein here
    # Consider passing ligand_like_chains to here, too

    complex_obj = run_plip(split_structure)
    ligand_list = [l for l in complex_obj.ligands if l.chain == ligand_plip_chain]
    if not len(ligand_list):
        log.warning(
            f"Could not find ligand at chain {ligand_plip_chain}, originally {ligand_chain} in {entry_pdb_id}"
        )
        return None
    ligand = ligand_list[0]
    lig_tag = f"{ligand.hetid}:{ligand.chain}:{ligand.position}"
    interactions = complex_obj.interaction_sets[lig_tag]
    chain_mapping = {v: k for k, v in chain_mapping.items()}
    return interactions, chain_mapping


def get_hash(
    interactions: PLInteraction,
    plip_chain_mapping: dict[str, str],
) -> tuple[dict[str, dict[int, list[str]]], set[(tuple[str, int])]]:
    """Get fingerprint hash from plip interaction object

    Parameters
    ----------
    intercation : PLInteraction
        plip interaction object for a given ligand
    plip_chain_mapping : Dict[str, str]
        chain mapping from plip chain ID to instance.asym ID

    Returns
    -------
    str
        plip fingerprint hash
    """

    interaction_hashes: dict[str, dict[int, list[str]]] = dict()
    waters = set()
    for int_type in INTERACTION_TYPES:
        int_objs = getattr(interactions, int_type)
        interaction_attributes = []
        if int_type in ["hbonds_ldon", "hbonds_pdon"]:
            for int_obj in int_objs:
                interaction_attributes.append(
                    "type:hydrogen_bonds"
                    # + f"__donortype:{int_obj.dtype}__acceptortype:{int_obj.atype}"
                    + f"__protisdon:{int_obj.protisdon}__sidechain:{int_obj.sidechain}"
                )
        elif int_type == "water_bridges":
            for int_obj in int_objs:
                interaction_attributes.append(
                    "type:water_bridges"
                    # + f"__donortype:{int_obj.dtype}__acceptortype:{int_obj.atype}"
                    + f"__protisdon:{int_obj.protisdon}"
                )
                waters.add((whichchain(int_obj.water), whichresnumber(int_obj.water)))
        elif int_type == "hydrophobic_contacts":
            for int_obj in int_objs:
                interaction_attributes.append("type:hydrophobic_contacts")
        elif int_type in ["pication_laro", "pication_paro"]:
            for int_obj in int_objs:
                if int_obj.protcharged:
                    group = "Aromatic"
                else:
                    # group = int_obj.charge.fgroup
                    group = "Cation"
                interaction_attributes.append(
                    "type:pi_cation"
                    + f"__lig_group:{group}"
                    + f"__protcharged:{int_obj.protcharged}"
                )
        elif int_type == "halogen_bonds":
            for int_obj in int_objs:
                interaction_attributes.append(
                    "type:halogen_bonds"
                    # + f"__donortype:{int_obj.donortype}"
                    # + f"__acceptortype:{int_obj.acctype}"
                    + f"__sidechain:{int_obj.sidechain}"
                )
        elif int_type == "pistacking":
            for int_obj in int_objs:
                interaction_attributes.append(
                    f"type:pi_stacks__stack_type:{int_obj.type}"
                )
        elif int_type in ["saltbridge_lneg", "saltbridge_pneg"]:
            for int_obj in int_objs:
                interaction_attributes.append(
                    "type:salt_bridges"
                    # + f"__pos_group:{int_obj.positive.fgroup}__neg_group:{int_obj.negative.fgroup}"
                    + f"__protispos:{int_obj.protispos}"
                )
        elif int_type == "metal_complexes":
            for int_obj in int_objs:
                interaction_attributes.append(
                    "type:metal_complexes"
                    + f"__metal_type:{int_obj.metal_type}__target_type:"
                    + f"{int_obj.target_type}__coordination:{int_obj.coordination_num}__geometry:"
                    + f"{int_obj.geometry}__location:{int_obj.location}"
                )
        for int_obj, int_attr in zip(int_objs, interaction_attributes):
            instance_chain, resnr = (
                plip_chain_mapping[int_obj.reschain],
                int(int_obj.resnr),
            )
            if instance_chain not in interaction_hashes:
                interaction_hashes[instance_chain] = dict()
            if resnr not in interaction_hashes[instance_chain]:
                interaction_hashes[instance_chain][resnr] = []
            interaction_hashes[instance_chain][resnr].append(int_attr)
    return interaction_hashes, waters
