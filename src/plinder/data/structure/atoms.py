# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Set, Tuple, Union

import biotite.structure as struc
import numpy as np
from biotite import TextFile
from biotite.structure import connect_via_residue_names
from biotite.structure.atoms import Atom, AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import PDBxFile, get_structure, set_structure
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Mol
from rdkit.DataStructs import cDataStructs
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

from plinder.core.utils.log import setup_logger
from plinder.data.structure.models import BackboneDefinition

log = setup_logger(__name__)

DOCKQ_BACKBONE_ATOMS = ["C", "CA", "N", "O"]

THREE_LETTER_AA = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "SER",
    "THR",
    "VAL",
    "TYR",
    "TRP",
    "PRO",
    "HIS",
    "PHE",
]

BINDING_SITE_METALS = [
    "MG",
    "K",
    "MN",
    "NA",
    "ZN",
    "MG",
    "CA",
    "CD",
    "FE",
    "CU",
    "4MO",
    "CO",
]

_AtomArrayOrStack = Union[AtomArray, AtomArrayStack]


def biotite_pdbfile() -> TextFile:
    from biotite.structure.io.pdb import PDBFile

    return PDBFile


def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import PDBxFile

    return PDBxFile


def atom_array_from_pdb_file(structure: Path | AtomArray) -> _AtomArrayOrStack | None:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_pdbfile()
        try:
            model = reader.read(str(structure))
            arr = model.get_structure(model=1)  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


def atom_array_from_cif_file(
    structure: Path | AtomArray, use_author_fields: bool = True
) -> AtomArray | None:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_ciffile()
        try:
            if structure.suffix == ".gz":
                with gzip.open(str(structure), "rt", encoding="utf-8") as f:
                    mod = reader.read(f)
            else:
                mod = reader.read(structure)
            arr = get_structure(mod, model=1, use_author_fields=use_author_fields)  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


def atom_vdw_radius(at: Atom) -> float:
    try:
        if at.hetero:
            rad = struc.info.vdw_radius_single(at.element)
        else:
            rad = struc.info.vdw_radius_protor(at.res_name, at.atom_name)
            if rad is None:
                rad = struc.info.vdw_radius_single(at.element)
    except Exception:
        rad = struc.info.vdw_radius_single("N")
    assert isinstance(rad, float)
    return rad


def backbone_mask(
    atoms: _AtomArrayOrStack, backbone_definition: BackboneDefinition
) -> NDArray[np.bool_]:
    if backbone_definition == "dockq":
        mask = np.isin(atoms.atom_name, DOCKQ_BACKBONE_ATOMS)
    else:
        mask = struc.filter_peptide_backbone(atoms)
    return mask


def mask_from_res_list(
    atoms: _AtomArrayOrStack, res_list: list[int]
) -> NDArray[np.bool_]:
    mask = np.isin(atoms.res_id, res_list)
    return mask


def mask_from_resn_resi_chain_list(
    atoms: _AtomArrayOrStack, resn_resi_ch_list: List[Tuple[str, int, str]]
) -> NDArray[np.bool_]:
    struc_info = [(atm.res_name, atm.res_id, atm.chain_id) for atm in atoms]
    mask = np.array([(sinfo in resn_resi_ch_list) for sinfo in struc_info])
    return mask


def apply_mask(atoms: _AtomArrayOrStack, mask: NDArray[np.bool_]) -> _AtomArrayOrStack:
    """
    Apply a boolean mask to an AtomArray or AtomArrayStack to filter atoms.

    Parameters:
    atoms (AtomArray | AtomArrayStack): The atoms to be filtered.
    mask (NDArray[np.bool_]): The boolean mask that specifies which atoms to keep.

    Returns:
    AtomArray | AtomArrayStack: The filtered atoms.
    """
    if isinstance(atoms, AtomArray):
        return atoms[mask]
    elif isinstance(atoms, AtomArrayStack):
        return atoms[..., mask]
    else:
        raise TypeError("atoms must be an AtomArray or AtomArrayStack")


def filter_atoms(
    atoms: _AtomArrayOrStack,
    calpha_only: bool = False,
    backbone_only: bool = False,
    heavy_only: bool = True,
    backbone_definition: BackboneDefinition = BackboneDefinition.dockq,
) -> _AtomArrayOrStack:
    if calpha_only:
        atoms = apply_mask(atoms, atoms.atom_name == "CA")
    if backbone_only:
        atoms = apply_mask(atoms, backbone_mask(atoms, backbone_definition))
    if heavy_only:
        atoms = apply_mask(atoms, atoms.element != "H")
    return atoms


def atom_array_to_rdkit_mol(
    lig_atom_array: AtomArray,
    smiles: str,
    lig_resn: str = "LIG",
    chain_mapping: Optional[Dict[str, str]] = None,
    add_h: bool = True,
) -> Mol:
    # Change chain from two letters to one
    if chain_mapping is None:
        lig_atom_array.chain_id[:] = "X"
    else:
        lig_ch_id = list(set(lig_atom_array.chain_id))[0]
        lig_atom_array.chain_id[:] = chain_mapping.get(lig_ch_id, "X")

    temp1 = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    file = struc.io.pdbx.PDBxFile()
    struc.io.pdbx.set_structure(file, lig_atom_array, data_block=lig_resn)
    struc.io.save_structure(temp1.name, lig_atom_array)
    pdb_block = temp1.read()
    os.unlink(temp1.name)
    ligand_pose = Chem.MolFromPDBBlock(pdb_block)
    try:
        if smiles != "":
            # If we have smiles, use it to assign bond orders
            # Load original molecule from smiles
            template = Chem.MolFromSmiles(smiles)
            new_mol = AllChem.AssignBondOrdersFromTemplate(template, ligand_pose)
        else:
            new_mol = Chem.Mol(ligand_pose)
    except (ValueError, TypeError) as e:
        log.warning(f"Bond assignment to ligand failed...: {e}")
        # If bond assignment or addition of H fails
        new_mol = Chem.Mol(ligand_pose)
    try:
        if add_h:
            log.info("Adding Hs ...")
            new_mol_h = Chem.AddHs(new_mol, addCoords=True)
            return new_mol_h
        else:
            new_mol_h = Chem.Mol(new_mol)
            return new_mol_h
    except (ValueError, TypeError) as e:
        log.warning(f"Addition of H to ligand failed...: {e}")
        new_mol_h = Chem.Mol(new_mol)
        return new_mol_h


def protonate_protein(
    pdb_block: str,
    timeout: int = 300,
    flags: List[str] = [
        "-OH",
        "-HIS",
        "-NOHETh",
    ],
) -> Tuple[str, str]:
    """
    Protonates pdb_block with reduce

    Parameters
    ----------
    pdb_block : str
        PBD-formatted string
    timeout : int, optional
        Number of seconds to wait for completion
    flags : list of str, optional
        Flags to pass to reduce

    Returns
    -------
    Tuple[str, str]
        The hydrogenated pdb block and the standard error output of reduce

    """
    args = ["reduce"] + list(flags) + ["-"]

    proc = subprocess.Popen(
        args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    try:
        out, err = proc.communicate(pdb_block.encode(), timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()

    return out.decode(), err.decode()


def get_molecular_properties(mol: Mol) -> List[Any]:
    """
    Get molecular properites.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol:
        rdkit molecule

    Returns
    -------
    List[Any]
        List of properties

    """
    chosen_descriptors = [
        "FractionCSP3",
        "HeavyAtomCount",
        "MolLogP",
        "MolWt",
        "NumAliphaticRings",
        "NumAromaticRings",
        "NumHAcceptors",
        "NumHDonors",
        "NumRotatableBonds",
        "RingCount",
    ]
    mol_descriptor_calculator = MolecularDescriptorCalculator(chosen_descriptors)
    # Use molecular descriptor calculator on RDKit mol object
    return list(mol_descriptor_calculator.CalcDescriptors(mol))


def get_morgan_fingerprint(mol: Mol, radius: int = 3, n_bits: int = 2048) -> Any:
    """
    Get ligand morgan fingerprint

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol:
        rdkit molecule
    radius : int
        fingerprint radius
    n_bits : 2048
        number of bits

    Returns
    -------
    str
        Base64 string
    """
    f = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=radius, useFeatures=True, nBits=n_bits
    )
    return f.ToBase64()


def extract_bsite_water(
    water: AtomArray, binding_site_entity: AtomArray, thresh_distance: float = 3.0
) -> Set[Any]:
    """Get binding site water"""

    cell_list = struc.CellList(binding_site_entity, cell_size=thresh_distance)

    contacts = cell_list.get_atoms(water.coord, radius=thresh_distance)
    if len(contacts) != 0:
        contact_indices = np.where((contacts != -1).any(axis=1))[0]
        # Get protein atoms in contact
        contact_atoms = water[contact_indices]
        contact_string_list = list(
            {(atom.res_name, atom.res_id, atom.chain_id) for atom in contact_atoms}
        )

        contact_string = {
            ",".join([i[0], str(i[1]), i[2]]) for i in contact_string_list
        }
    else:
        contact_string = set({})

    return contact_string


def extract_bsite_metal(
    metal: AtomArray, binding_site_entity: AtomArray, thresh_distance: float = 3.0
) -> Set[Any]:
    """Get binding site water"""
    cell_list = struc.CellList(binding_site_entity, cell_size=thresh_distance)

    contacts = cell_list.get_atoms(metal.coord, radius=thresh_distance)
    if len(contacts) != 0:
        contact_indices = np.where((contacts != -1).any(axis=1))[0]
        # Get protein atoms in contact
        contact_atoms = metal[contact_indices]
        contact_string_list = list(
            {(atom.res_name, atom.res_id, atom.chain_id) for atom in contact_atoms}
        )

        contact_string = {
            ",".join([i[0], str(i[1]), i[2]]) for i in contact_string_list
        }
    else:
        contact_string = set({})

    return contact_string


def extract_specific_covalent_connections(
    list_of_lig_ids_res_nums: list[Any],
    componenets_graphs: dict[Any, Any],
    node_lookup: dict[Any, Any],
) -> Dict[Any, Any]:
    specific_components = {
        k: [t for l in j for t in node_lookup[l]]
        for res in list_of_lig_ids_res_nums
        for k, v in componenets_graphs.items()
        for j in v
        for p in j
        if f"{res[0]},{res[1]}" in p
    }
    specific_components = {k: v for k, v in specific_components.items() for j in v}

    return specific_components


def extract_nucleic_acid_base_pairing(
    component_graphs: Dict[Any, Any],
) -> Dict[Any, Any]:
    combined_components = {
        k: list({"-".join(list(j)) for j in v})
        for k, v in component_graphs.items()
        for j in v
    }
    return combined_components


def reformat_resn_resi(
    resi_resn: Tuple[np.ndarray[int, Any], np.ndarray[int, Any]],
) -> List[Tuple[str, int]]:
    return [(resi_resn[1][idx], resi) for idx, resi in enumerate(resi_resn[0])]


def transfer_category(
    reference_block: AtomArray,
    reference_category_list: List[Any],
    target_structure: AtomArray,
) -> AtomArray:
    target_block = target_structure.make_mmcif_block()
    for cat in reference_category_list:
        comp_atoms = reference_block.get_mmcif_category(cat)
        if cat in [
            "_pdbx_nonpoly_scheme",
            "_pdbx_poly_seq_scheme",
            "_chem_comp_atom",
            "_chem_comp_bond",
        ]:
            target_block.set_mmcif_category(cat, comp_atoms)
        else:
            target_block.set_mmcif_category(f"_added{cat}", comp_atoms)

    target_structure.update_mmcif_block(target_block)
    return target_structure


def canonicalize(smiles: str) -> Union[str, Mol]:
    try:
        return Chem.CanonSmiles(smiles)
    except:
        return smiles


def create_receptor_entity_from_struc(
    structure: AtomArray, resn_resi: List[Any], chain_id: str
) -> np.ndarray[int, Any]:
    binder_info = [(rinfo[0], rinfo[1], chain_id) for rinfo in resn_resi]
    struc_info = [(atm.res_name, atm.res_id, atm.chain_id) for atm in structure]
    return np.array([(sinfo in binder_info) for sinfo in struc_info])


def get_original_chain_mapping(
    old_chain_ids: List[str], new_chain_ids: List[str]
) -> Dict[str, str]:
    """
    Using the extreme case where original chains
    are also numbers
    ----
    Tests:
    old_ch = ["1", "2", "3","11", "22", "111"],
    new_ch = ["11", "21", "31","111", "223", "1111"]

    old_ch = ["A", "B", "AA","AB", "AC", "AD"],
    new_ch =["A2", "B3", "AA3","AB5", "AC1", "AD7"]
    ----
    """
    mapping = {}
    for new_ch in new_chain_ids:
        old_ch_candidates = []
        for old_ch in old_chain_ids:
            if new_ch.startswith(old_ch) & (new_ch[len(old_ch) :] != ""):
                old_ch_candidates.append(old_ch)
        # Deal with cases like: A -> A1; A, AA ->AA1, select the largest
        stripped = [new_ch[len(i) :] for i in old_ch_candidates]

        min_stripped = min([new_ch[len(i) :] for i in old_ch_candidates], key=len)

        old_ch_final = old_ch_candidates[stripped.index(min_stripped)]
        mapping[new_ch] = old_ch_final
    return mapping


def get_binding_site_info(
    structure: AtomArray,
    binder_chain_id: str,
    resn_resi: List[Any],
    binder_entity: AtomArray,
    binder_entity_type: str,
    binder_type: str,
    all_covalent_dict: Dict[Any, Any],
    covalent_atom_lookup: Dict[Any, Any],
    thresh_distance: float = 4.0,
) -> Any:
    binder_filter = create_receptor_entity_from_struc(
        structure, resn_resi, binder_chain_id
    )
    receptor_entity_filter = ~(binder_filter) & ~(struc.filter_solvent(structure))
    receptor_entity = structure[receptor_entity_filter]
    if (receptor_entity[struc.filter_amino_acids(receptor_entity)]) == 0:
        return ("", "", "", "", "", "", "", "", "", "", "", "", "")
    all_water = structure[structure.res_name == "HOH"]
    all_metals = structure[np.isin(structure.res_name, BINDING_SITE_METALS)]

    set_of_connected_binders = extract_specific_covalent_connections(
        resn_resi, all_covalent_dict, covalent_atom_lookup
    )

    cell_list = struc.CellList(binder_entity, cell_size=thresh_distance)

    contacts = cell_list.get_atoms(receptor_entity.coord, radius=thresh_distance)

    # Only retain atoms in the protein with contact to at least one atom of the ligand
    contact_indices = np.where((contacts != -1).any(axis=1))[0]
    # Get protein atoms in contact
    contact_atoms = receptor_entity[contact_indices]
    contact_string_list = list(
        {(atom.res_name, atom.res_id, atom.chain_id) for atom in contact_atoms}
    )
    pocket_chains = {atom.chain_id for atom in contact_atoms}

    number_of_chains = len(pocket_chains)
    contact_string = ";".join(
        [",".join([i[0], str(i[1]), i[2]]) for i in contact_string_list]
    )

    # Get bsite water
    binding_site = receptor_entity[contact_indices]
    if (len(all_water) == 0) or (len(binding_site) == 0) or (len(binder_entity) == 0):
        binder_receptor_bridging_water = ""
    else:
        bsite_water_receptor = extract_bsite_water(all_water, binding_site)
        # Get interaction with ligand

        bsite_water_binder = extract_bsite_water(all_water, binder_entity)

        binder_receptor_bridging_water = ";".join(
            list(bsite_water_receptor.intersection(bsite_water_binder))
        )

    # Extract binding site metals

    if (len(all_metals) == 0) or (len(binding_site) == 0):
        bsite_metals = ""
    else:
        bsite_metals = ";".join(list(extract_bsite_metal(all_metals, binding_site)))

    covalent_interactions = set_of_connected_binders.get("covale")
    metalic_interactions = set_of_connected_binders.get("metalc")
    hydrog_interactions = set_of_connected_binders.get("hydrog")
    if covalent_interactions is not None:
        is_cov = (
            len([cov.split(",")[1] in THREE_LETTER_AA for cov in covalent_interactions])
            > 0
        )
        cov_string = ";".join(list(covalent_interactions))
        ptm = (is_cov) & (binder_entity_type == "Polymer")
    else:
        is_cov = False
        cov_string = ""
        ptm = False
    if metalic_interactions is not None:
        met_string = ";".join(list(metalic_interactions))
    else:
        met_string = ""
    if hydrog_interactions is not None:
        hydro_string = ";".join(list(hydrog_interactions))
    else:
        hydro_string = ""
    is_hetero = len(resn_resi) == 1

    return (
        contact_string,
        cov_string,
        met_string,
        hydro_string,
        is_cov,
        number_of_chains,
        is_hetero,
        binder_type,
        binder_entity_type,
        ptm,
        binder_receptor_bridging_water,
        bsite_metals,
        pocket_chains,
    )


def reduce_atom_array(
    protein_atom_array: AtomArray, chain_mapping: Dict[Any, Any]
) -> Mol:
    if len(protein_atom_array) == 0:
        return None
    # Map chains
    else:
        for old_ch, new_chain in chain_mapping.items():
            ch_mask = protein_atom_array.chain_id == old_ch
            protein_atom_array.chain_id[ch_mask] = new_chain

        temp1 = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        file = PDBxFile()
        set_structure(file, protein_atom_array, data_block="model1")
        struc.io.save_structure(temp1.name, protein_atom_array)

        reduced_pdb_block = protonate_protein(temp1.read().decode())[0]
        os.unlink(temp1.name)
        return Chem.MolFromPDBBlock(reduced_pdb_block, removeHs=False)


def add_hydrogen_to_lig_array(
    lig_atom_array: AtomArray,
    smiles: str,
    lig_resn: str = "LIG",
    chain_mapping: Optional[Dict[Any, Any]] = None,
    add_h: bool = True,
) -> Mol:
    if len(lig_atom_array) == 0:
        return None
    else:
        # Change chain from two letters to one
        if chain_mapping is None:
            lig_atom_array.chain_id[:] = "X"
        else:
            lig_ch_id = list(set(lig_atom_array.chain_id))[0]
            lig_atom_array.chain_id[:] = chain_mapping.get(lig_ch_id, "X")

        temp1 = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        file = PDBxFile()
        set_structure(file, lig_atom_array, data_block=lig_resn)
        struc.io.save_structure(temp1.name, lig_atom_array)
        pdb_block = temp1.read()
        os.unlink(temp1.name)
        ligand_pose = Chem.MolFromPDBBlock(pdb_block)
        try:
            if smiles != "":
                # If we have smiles, use it to assign bond orders load original molecule from smiles
                template = Chem.MolFromSmiles(smiles)
                new_mol = AllChem.AssignBondOrdersFromTemplate(template, ligand_pose)
            else:
                new_mol = Chem.Mol(ligand_pose)

        except (ValueError, TypeError):
            log.warning("Bond assignment to ligand failed...: {e}")
            # If bond assignment or addition of H fails
            new_mol = Chem.Mol(ligand_pose)

        try:
            if add_h:
                log.info("Adding Hs ...")
                new_mol_h = Chem.AddHs(new_mol, addCoords=True)
                return new_mol_h
            else:
                new_mol_h = Chem.Mol(new_mol)
                return new_mol_h
        except (ValueError, TypeError) as e:
            log.warning(f"Addition of H failed to ligand...: {e}")
            new_mol_h = Chem.Mol(new_mol)
            return new_mol_h


def add_bond_order_to_protein_atom_array(atom_array: AtomArray) -> AtomArray:
    bonds = connect_via_residue_names(atom_array)
    atom_array.bonds = bonds
    return atom_array


def get_fingerprint_similarity(
    list_of_ids: List[int], fps: List[Any], threshold: float = 0.6
) -> Generator[Any, Any, Any]:
    # Compare all fp pairwise without duplicates
    for n in range(len(fps) - 1):  # -1 so the last fp will not be used
        s = cDataStructs.DataStructs.BulkTanimotoSimilarity(
            fps[n], fps[n + 1 :]
        )  # +1 compare with the next to the last fp

        for m in range(len(s)):
            if s[m] >= threshold:
                yield (list_of_ids[n], list_of_ids[n + 1 :][m], s[m])
