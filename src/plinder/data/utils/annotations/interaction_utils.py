# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import gemmi
from mmcif.api.PdbxContainers import DataContainer
from ost import io, mol
from plip.basic.supplemental import whichchain, whichresnumber
from plip.structure.preparation import PDBComplex, PLInteraction

from plinder.core.utils.log import setup_logger

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


def get_symmetry_mate_contacts(
    mmcif_filename: Path, contact_threshold: float = 5.0
) -> dict[tuple[str, int], dict[tuple[str, int], set[int]]]:
    """
    Get all contacts within a given threshold between residues which are not in the same asymmetric unit (symmetry mates)
    """
    cif = gemmi.read_structure(mmcif_filename.__str__(), merge_chain_parts=False)
    cif.setup_entities()
    ns = gemmi.NeighborSearch(cif[0], cif.cell, contact_threshold).populate(
        include_h=False
    )
    cs = gemmi.ContactSearch(contact_threshold)
    cs.ignore = gemmi.ContactSearch.Ignore.SameAsu
    cs.twice = True
    pairs = cs.find_contacts(ns)
    results: dict[tuple[str, int], dict[tuple[str, int], set[int]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for p in pairs:
        if p.partner1.residue.is_water() or p.partner2.residue.is_water():
            continue
        i1, i2 = p.partner1.residue.label_seq, p.partner2.residue.label_seq
        if i1 is None:
            i1 = 1
        if i2 is None:
            i2 = 1
        c1, c2 = p.partner1.residue.subchain, p.partner2.residue.subchain
        results[(c1, i1)][(c2, i2)].add(p.partner1.atom.serial)
    return results


def get_covalent_connections(data: DataContainer) -> dict[str, list[tuple[str, str]]]:
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

    cov_dict = defaultdict(list)
    nucleobase_list = ["A", "C", "U", "G", "DA", "DC", "DG", "DT", "PSU"]

    to_extract = [
        "ptnr1_label_asym_id",
        "ptnr2_label_asym_id",
        "ptnr1_label_seq_id",
        "ptnr2_label_seq_id",
        "ptnr1_auth_seq_id",
        "ptnr2_auth_seq_id",
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

        if con["conn_type_id"] == "hydrog":
            if con["ptnr1_label_comp_id"].strip() not in nucleobase_list:
                continue
        cov_dict[con["conn_type_id"]].append(
            (
                ":".join(
                    [
                        con["ptnr1_auth_seq_id"],
                        con["ptnr1_label_comp_id"],
                        con["ptnr1_label_asym_id"],
                        con["ptnr1_label_seq_id"],
                        con["ptnr1_label_atom_id"],
                    ]
                ),
                ":".join(
                    [
                        con["ptnr2_auth_seq_id"],
                        con["ptnr2_label_comp_id"],
                        con["ptnr2_label_asym_id"],
                        con["ptnr2_label_seq_id"],
                        con["ptnr2_label_atom_id"],
                    ]
                ),
            )
        )
    return cov_dict


def extract_ligand_links_to_neighbouring_chains(
    all_covalent_dict: dict[str, list[tuple[str, str]]],
    ligand_asym_id: str,
    neighboring_asym_ids: set[str],
    link_type: str = "covale",
) -> set[str]:
    covalent_linkages = set()
    if link_type in all_covalent_dict:
        for link1, link2 in all_covalent_dict[link_type]:
            chain1, chain2 = link1.split(":")[2], link2.split(":")[2]
            chains = {chain1, chain2}
            if len(chains) == 1:
                # remove linkages that are to the same chain!
                continue
            if set(chains).intersection(ligand_asym_id):
                # only if ligand is involved
                # now check that one chain is neighbour and the other is ligand
                # enforce receptor_ligand ordering
                if neighboring_asym_ids.intersection([chain1]) and set(
                    [chain2]
                ).intersection(ligand_asym_id):
                    covalent_linkages.add(f"{link1}__{link2}")
                elif neighboring_asym_ids.intersection(chain2) and set(
                    [chain1]
                ).intersection(ligand_asym_id):
                    covalent_linkages.add(f"{link2}__{link1}")
    return covalent_linkages


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
    biounit: mol.EntityHandle,
    biounit_selection: mol.EntityHandle,
    ligand_chain: str,
) -> tuple[PLInteraction, dict[str, str]] | None:
    """Split structure into small PLI complex by ligand

    For every ligand, create a smaller complex for faster plip\
    process and deal with cases where plip ignores small molecule
    ligands in the presence of peptides

    Parameters
    ----------
    biounit : mol.EntityHandle
        biounit
    biounit_selection : mol.EntityHandle
        selection in biounit of threshold around ligand
    ligand_chain: str
        {instance}.{chain} of ligand

    Returns
    -------
    Tuple[PlipLigand, PDBComplex, dict[str, str]] | None
        PLIP ligand, PLIP complex object, mapping of original chain to plip chain
    """
    from plip.basic import config

    config.biolip_list = []
    split_structure, chain_mapping = pdbize(biounit, biounit_selection)
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
            f"Could not find ligand at chain {ligand_plip_chain}, originally {ligand_chain}"  # in {entry_pdb_id}"
        )
        return None
    ligand = ligand_list[0]
    lig_tag = f"{ligand.hetid}:{ligand.chain}:{ligand.position}"
    interactions = complex_obj.interaction_sets[lig_tag]
    chain_mapping = {v: k for k, v in chain_mapping.items()}
    return interactions, chain_mapping


def get_plip_hash(
    interactions: PLInteraction,
    chain: str,
    plip_chain_mapping: dict[str, str],
) -> tuple[dict[str, dict[int, list[str]]], set[(tuple[str, int])]]:
    """Get fingerprint hash from plip interaction object

    Parameters
    ----------
    interactions : PLInteraction
        plip interaction object for a given ligand
    chain: str
        ligand chain
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
            if instance_chain == chain:
                continue
            if instance_chain not in interaction_hashes:
                interaction_hashes[instance_chain] = dict()
            if resnr not in interaction_hashes[instance_chain]:
                interaction_hashes[instance_chain][resnr] = []
            interaction_hashes[instance_chain][resnr].append(int_attr)
    return interaction_hashes, waters
