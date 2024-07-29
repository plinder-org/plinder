# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import typing as ty
from pathlib import Path

from ost import conop, io, mol
from rdkit import Chem

from plinder.data.utils.annotations.rdkit_utils import ligand_ost_ent_to_rdkit_mol

# Define available names for protein and ligand chains in PDB format
PDB_PROTEIN_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
PDB_LIGAND_CHAINS = PDB_PROTEIN_CHAINS.lower() + "0123456789"
WATER_CHAIN_NAME = "_"


def save_ligands(
    ent: mol.EntityHandle,
    ligand_selections: list[str],
    ligand_chains: list[str],
    ligand_smiles: list[str],
    ligand_num_unresolved_heavy_atoms: list[int | None],
    output_folder: str | Path,
) -> None:
    for selection, chain, smiles, num_unresolved_heavy_atoms in zip(
        ligand_selections,
        ligand_chains,
        ligand_smiles,
        ligand_num_unresolved_heavy_atoms,
    ):
        ligand_ost = mol.CreateEntityFromView(ent.Select(selection), True)
        rdkit_mol = ligand_ost_ent_to_rdkit_mol(
            ligand_ost, smiles, num_unresolved_heavy_atoms or 0
        )
        rdkit_mol.SetProp("_Name", chain)
        with Chem.SDWriter(str(Path(output_folder) / f"{chain}.sdf")) as w:
            w.write(rdkit_mol)


def save_pdb_file(
    full_biounit: mol.EntityHandle,
    ent: mol.EntityHandle,
    protein_chains: ty.List[str],
    ligand_chains: ty.List[str],
    output_pdb_file: ty.Union[str, Path],
    output_mapping_file: ty.Union[str, Path],
    waters: ty.Dict[str, ty.List[int]],
    water_mapping_file: ty.Union[str, Path],
) -> None:
    """Renames protein and ligand chains to fit the PDB single-letter chain name convention, saves the entity to a PDB file
    and saves the mapping between original and new chain names to a JSON file

    Args:
        full_biounit (mol.EntityHandle): original biounit entity with all chains
        ent (mol.EntityHandle): selected entity with only the protein and ligand chains of the system
        protein_chains (ty.List[str]): list of protein chains
        ligand_chains (ty.List[str]): list of ligand chains
        waters (ty.Dict[str, ty.List[int]]): dictionary with water chain names as keys and list of water residue indices as values
        output_pdb_file (ty.Union[str, Path]): path to the output PDB file
        output_mapping_file (ty.Union[str, Path]): path to the output JSON file
    """
    if len(waters):
        ent = mol.CreateEntityFromView(ent.Select("water=False"), True)

    # Intermediate renaming step
    intermediate_names = {}
    edi = ent.EditXCS(mol.BUFFERED_EDIT)
    for i, chain in enumerate(ent.GetChainList()):
        intermediate_names[f"T{i}"] = chain.name
        edi.RenameChain(chain, f"T{i}")
    edi.UpdateICS()

    # Final renaming step
    protein_chain_index = 0
    ligand_chain_index = 0
    name_mapping = {}
    water_mapping: dict[str, dict[int, int]] = {}

    for chain in ent.GetChainList():
        original_name = intermediate_names[chain.name]
        original_chain = full_biounit.FindChain(original_name)
        if original_name in protein_chains:
            final_name = PDB_PROTEIN_CHAINS[protein_chain_index]
            protein_chain_index += 1
        elif original_name in ligand_chains:
            final_name = PDB_LIGAND_CHAINS[ligand_chain_index]
            ligand_chain_index += 1
        edi.RenameChain(chain, final_name)
        edi.SetChainDescription(chain, original_chain.description)
        edi.SetChainType(chain, original_chain.type)
        name_mapping[original_name] = final_name

    if len(waters):
        water_chain = edi.InsertChain(WATER_CHAIN_NAME)
        edi.SetChainDescription(water_chain, "Interacting waters")
        edi.SetChainType(water_chain, mol.CHAINTYPE_WATER)
        index = 1
        for chain in waters:
            water_mapping[chain] = {}
            for resnum in waters[chain]:
                new_residue = edi.AppendResidue(
                    water_chain,
                    full_biounit.FindChain(chain).FindResidue(resnum),
                    deep=True,
                )
                edi.SetResidueNumber(new_residue, index)
                water_mapping[chain][int(resnum)] = new_residue.number.num
                index += 1

    edi.UpdateICS()
    io.SavePDB(ent, str(output_pdb_file))
    with open(output_mapping_file, "w") as f:
        json.dump(name_mapping, f)
    if len(waters):
        with open(water_mapping_file, "w") as f:
            json.dump(water_mapping, f)


def save_cif_file(
    ent: mol.EntityHandle,
    info: io.MMCifInfo,
    name: str,
    output_cif_file: ty.Union[str, Path],
) -> None:
    lib = conop.GetDefaultLib()
    entity_info = io.MMCifWriterEntityList()
    entity_ids = set(
        info.GetMMCifEntityIdTr(ch.name.split(".")[-1]) for ch in ent.chains
    )
    for entity_id in info.GetEntityIdsOfType("polymer"):
        if entity_id not in entity_ids:
            continue
        # Get entity description from info object
        entity_desc = info.GetEntityDesc(entity_id)
        e = io.MMCifWriterEntity.FromPolymer(
            entity_desc.entity_poly_type, entity_desc.mon_ids, lib
        )
        entity_info.append(e)
        # search all chains assigned to the entity we just added
        for ch in ent.chains:
            if info.GetMMCifEntityIdTr(ch.name.split(".")[-1]) == entity_id:
                entity_info[-1].asym_ids.append(ch.name)
        # deal with heterogeneities
        for a, b in zip(entity_desc.hetero_num, entity_desc.hetero_ids):
            entity_info[-1].AddHet(a, b)
    writer = io.MMCifWriter()
    writer.SetStructure(ent, lib, entity_info=entity_info)
    writer.Write(name, str(output_cif_file))
