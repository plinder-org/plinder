# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
import random
import re
import shutil
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import biotite.structure as struc
import gemmi
import numpy as np
import pandas as pd
import requests
from biotite.structure.io.pdbx import PDBxFile, get_structure, set_structure
from pdbeccdutils.core.boundmolecule import infer_bound_molecules
from pdbeccdutils.core.clc_reader import _parse_pdb_mmcif
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol, RWMol
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

from plinder.core.utils.log import setup_logger
from plinder.data.structure.atoms import get_original_chain_mapping

log = setup_logger(__name__)

# Set global seed
random.seed(42)
RDLogger.DisableLog("rdApp.*")


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


def str_to_int(i: str) -> int:
    """
    Converts a string into integer.
    Returns 0 if a string cannot be
    converted.

    Parameters
    ----------
    i : str
        The first parameter.

    Returns
    -------
    int
       conversion of the string
    """
    try:
        return int(i)
    except ValueError:
        return 0


def get_all_indices(lst: list[Any], item: Any) -> list[Any]:
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


def bond_pdb_order(value_order: str) -> Union[Chem.rdchem.BondType, None]:
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


def remove_quote_from_string(value: str) -> str:
    """
    remove quote from string

    Parameters
    ----------
    value : str

    Returns
    -------
    str
    """
    return value.replace('"', "")


def get_all_bound_molecules(
    mmcif_file: Path, to_discard: List[str] = ["HOH"]
) -> Dict[Tuple[Any, ...], Any]:
    """
    Protonates pdb_block with reduce

    Parameters
    ----------
    mmcif_file : Path
        mmcif path
    to_discard : List[str]
        List of molecules to discard

    Returns
    -------
    Dict[Tuple[str], List[Any]]
        The hydrogenated pdb block and the standard error output of reduce

    """
    cif_block = gemmi.cif.read_file(str(mmcif_file)).sole_block()
    bms = infer_bound_molecules(str(mmcif_file), to_discard)

    # Note: Don't use PDBBlock, for some reason it doesn't work
    molecules = {}
    for bm in bms:
        bm_dict = bm.to_dict()
        key = tuple(sorted([i["id"] for i in bm_dict["residues"]]))
        bm_dict["rdkit_mol"] = _parse_pdb_mmcif(cif_block, bm.graph)
        molecules[key] = bm_dict

    return molecules


def get_specific_bound_molecules(
    mmcif_file: Path,
    residue_list: List[Tuple[str, str]],
    to_discard: List[str] = ["HOH"],
) -> Mol:
    """
    Protonates pdb_block with reduce

    Parameters
    ----------
    mmcif_file : Path
        mmcif path
    residue_list : List[Tuple[str, str]]
        list of residues to select
    to_discard : List[str]
        List of molecules to discard

    Returns
    -------
    Mol
        RDKit Mol

    """
    cif_block = gemmi.cif.read_file(str(mmcif_file)).sole_block()
    bms = infer_bound_molecules(str(mmcif_file), to_discard)
    # Note: Don't use PDBBlock, for some reason it doesn't work
    residue_info = {f"{ch}{resi}" for ch, resi in residue_list}
    for bm in bms:
        bm_dict = bm.to_dict()
        key = tuple(sorted([i["id"] for i in bm_dict["residues"]]))
        if set(key) == residue_info:
            bm_dict["rdkit_mol"] = _parse_pdb_mmcif(cif_block, bm.graph)
            return bm_dict["rdkit_mol"][0]


def get_selected_residues_pdb_block(
    mmcif_file: Path, residue_criteria: List[Tuple[str, str, str]]
) -> Any:
    """
    Selects specified residues from a structure based
    on residue number, name, and chain ID, and returns
    a string in PDB format representing these residues.

    Bond order assined with "_chem_comp_atom" and
    "_chem_comp_bond" category of nextgen rcsb mmcif
    -enriched.cif.gz file

    Parameters
    ----------
    structure : gemmi.Structure:
        The input protein structure.
    residue_criteria : List[Tuple[str, str, str]
        A list of tuples,
        each specifying a residue by
        (chain_id, residue number, residue name).

    Returns
    -------
    str
        String in PDB format representing the selected residues.

    """
    cif_block = gemmi.cif.read_file(str(mmcif_file)).sole_block()

    if (
        "_chem_comp_atom." not in cif_block.get_mmcif_category_names()
        or "_chem_comp_bond." not in cif_block.get_mmcif_category_names()
    ):
        return None
    atoms = cif_block.find("_chem_comp_atom.", ["comp_id", "atom_id"])
    atoms_dict = defaultdict(list)
    for atom in atoms:
        atoms_dict[atom["comp_id"]].append(atom["atom_id"])

    res_names = [i[0] for i in residue_criteria]
    bonds = [
        bond
        for bond in cif_block.find(
            "_chem_comp_bond.", ["comp_id", "atom_id_1", "atom_id_2", "value_order"]
        )
        if bond["comp_id"] in res_names
    ]
    bonds_dict = defaultdict(list)
    for bond in bonds:
        bonds_dict[bond["comp_id"]].append(
            (
                bond["atom_id_1"],
                bond["atom_id_2"],
                bond["value_order"],
            )
        )

    structure = gemmi.read_structure(str(mmcif_file))
    # Create a new structure to hold the selected residues
    new_structure = gemmi.Structure()

    # Iterate over the list of residue criteria
    for residue_name, residue_num, chain_id in residue_criteria:
        # Try to find the chain and residue in the original structure
        chain = structure[0][chain_id]
        if chain:
            residue = chain[str(residue_num)][residue_name]
            # Add the chain and residue to the new structure if not already added
            new_model = (
                new_structure[0]
                if new_structure
                else new_structure.add_model(gemmi.Model("model"), 0)
            )
            new_chain = new_model.add_chain(gemmi.Chain(chain_id))
            new_chain.add_residue(residue)

    # Generate a PDB-formatted string
    pdb_block = new_structure.make_pdb_string()
    mol = AllChem.MolFromPDBBlock(pdb_block, proximityBonding=False)
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
                # Extra atom
                pass

            if (
                (comp_id not in res_names)
                | atom_1.startswith("H")
                | atom_2.startswith("H")
            ):
                # Skip hydrogens
                pass
            else:
                try:
                    atom_1_ids = get_all_indices(atoms_ids, f"{comp_id}:{atom_1}")
                    atom_2_ids = get_all_indices(atoms_ids, f"{comp_id}:{atom_2}")
                    for atom_1_id, atom_2_id, order in zip(
                        atom_1_ids, atom_2_ids, np.repeat(row[2], len(atom_1_ids))
                    ):
                        bond_order = bond_pdb_order(order)
                        rw_mol.AddBond(int(atom_1_id), int(atom_2_id), bond_order)
                except ValueError:
                    log.warning(
                        f"Error perceiving {atom_1} - {atom_2} bond in _chem_comp_bond"
                    )
                except RuntimeError:
                    log.warning(f"Duplicit bond {atom_1} - {atom_2}")

    bonded_mol = rw_mol.GetMol()
    return Chem.MolToPDBBlock(bonded_mol)


def to_alpha(string: str) -> str:
    """
    Digit to alphabet
    """
    if string.isdigit():
        # If the original chain is a digit
        return string[:-1]
    else:
        return "".join(filter(lambda x: x.isalpha(), string))


def download_from_url(url: str, filename: Path) -> Path:
    """
    Download url

    Parameters
    ----------
    url : str
        rdkit molecule
    filename : Path
        output file path

    Returns
    -------
    None
    """
    if not filename.exists():
        r = requests.get(url, allow_redirects=True)
        with open(filename, "wb") as out:
            out.write(r.content)
    return filename


def download_alphafold_model(uniprot_id: str, output_path: Path) -> Path:
    url = "https://alphafold.ebi.ac.uk/files/" + "AF-{uniprot_id}-F1-model_v4.cif"
    return download_from_url(url, output_path / f"AF-{uniprot_id}.cif")


def download_extract_tar_bindingdb_affinity_data() -> Path:
    from io import BytesIO
    from urllib.request import urlopen
    from zipfile import ZipFile

    zipurl = "https://www.bindingdb.org/bind/downloads/BindingDB_All_202401_tsv.zip"
    with urlopen(zipurl) as zipresp:
        with ZipFile(BytesIO(zipresp.read())) as zfile:
            zfile.extractall()
    return Path.cwd() / "BindingDB_All_202401_tsv"


def download_extract_tar_pdbbind_affinity_data() -> Path:
    if not Path("index").exists():
        url = (
            "http://www.pdbbind.org.cn/"
            + "download/PDBbind_v2020_plain_text_index.tar.gz"
        )
        with requests.get(url, stream=True) as rx, tarfile.open(
            fileobj=rx.raw, mode="r:gz"
        ) as tarobj:
            tarobj.extractall()
        shutil.rmtree("readme")
    return Path.cwd() / "index/INDEX_general_PL_data.2020"


def sort_ccd_codes(code_list: List[str]) -> List[str]:
    """Pick long first, then alphabetical letters followed by numbers
    Args:
        code_list (Set[str]): set of CCD strings

    Returns:
        List[str]: list of sorted CCD string set
    """
    code_list = sorted(sorted(code_list), key=len, reverse=True)
    final_list = [code for code in code_list if not re.findall("([0-9])", code[0])] + [
        code for code in code_list if re.findall("([0-9])", code[0])
    ]
    return final_list


def get_ccd_synonym_frame(smidict: Dict[str, str]) -> pd.DataFrame:
    """Get Synonym dictonary for CCD SMILES
    Args:
        smidict (Dict[str, str]): CCD smiles dict from download_components_cif
        and get_ccd_smiles_dict()

    Returns:
        Dict[str, str]: dictonary mapping synonymous CCD code to preferred one
    """
    ccd_df = pd.DataFrame.from_dict(smidict, orient="index").reset_index()
    # note: SMILES assumed to be CANONICALIZED by OE read by get_ccd_smiles_dict()
    ccd_df.columns = ["ccd_code", "SMILES"]
    ccd_sets = ccd_df.groupby("SMILES").aggregate(list).reset_index()
    # ccd_sets
    ccd_sets["ccd_synonym_count"] = ccd_sets["ccd_code"].apply(lambda x: len(x))
    ccd_dups = ccd_sets[ccd_sets.ccd_synonym_count > 1].copy()
    # keep unique_ccd as sorted first entry
    ccd_dups["unique_ccd"] = ccd_dups["ccd_code"].apply(lambda x: sort_ccd_codes(x)[0])
    return ccd_dups
    # ccd_dups = ccd_dups.explode('ccd_code')
    # ccd_dups.index = ccd_dups['ccd_code']
    # return ccd_dups['unique_ccd'].to_dict()


def generate_bio_assembly(
    mmcif_filename: Path,
    output_path: Path,
    list_of_extra_categories: List[str] = [
        "refine",
        "entity",
        "entity_poly",
        "pdbx_nonpoly_scheme",
        "pdbx_entity_branch_link",
        "pdbx_poly_seq_scheme",
        "pdbx_branch_scheme",
        "chem_comp_atom",
        "chem_comp_bond",
    ],
) -> Tuple[Path, Any]:
    """
    Generate biological assemblies for the given mmCIF file.
    """
    output_path.mkdir(parents=True, exist_ok=True)
    # Read CIF file and the structure
    block = gemmi.cif.read(str(mmcif_filename))[0]
    structure = gemmi.make_structure_from_block(block)
    old_chain_id = [st.name for st in structure[0]]

    structure.transform_to_assembly(
        assembly_name="1", how=gemmi.HowToNameCopiedChain.AddNumber
    )

    fp = tempfile.NamedTemporaryFile(delete=False)
    structure.make_mmcif_document().write_file(str(fp.name))
    new_chain_id = [st.name for st in structure[0]]
    chain_map = get_original_chain_mapping(old_chain_id, new_chain_id)
    chain_map = {v: k for k, v in chain_map.items()}
    # Write out the assembly in CIF format
    entry_id = str(mmcif_filename.parent).lower()[-4:]
    new_fname = output_path / f"{entry_id}-assembly.cif"

    add_extra_loop_to_bioassembly(
        Path(mmcif_filename),
        Path(fp.name),
        chain_map,
        Path(new_fname),
        list_of_extra_categories,
    )
    return new_fname, chain_map


def add_extra_loop_to_bioassembly(
    raw_cif_path: Path,
    bio_assembly_cif_path: Path,
    chain_map: Dict[str, str],
    output_cif_path: Path,
    list_of_extra_categories: List[str] = [
        "refine",
        "entity",
        "entity_poly",
        "pdbx_nonpoly_scheme",
        "pdbx_entity_branch_link",
        "pdbx_poly_seq_scheme",
        "pdbx_branch_scheme",
        "chem_comp_atom",
        "chem_comp_bond",
    ],
) -> Path:
    cif_handle = read_mmcif_file(raw_cif_path)
    bio_cif_handle = PDBxFile.read(str(bio_assembly_cif_path))

    for aff in list_of_extra_categories:
        cat = cif_handle.get_category(aff, expect_looped=True)
        final_cat = {}
        if aff in ["pdbx_poly_seq_scheme", "pdbx_nonpoly_scheme"]:
            inner_dict = {}
            if "pdb_strand_id" in cat.keys():
                relevant_indices = [
                    idx
                    for idx, i in enumerate(cat["pdb_strand_id"])
                    if i in list(chain_map.keys())
                ]
            for key, val in cat.items():
                if key in ["pdb_strand_id"]:
                    inner_dict[key] = np.array(
                        [chain_map[i] for i in val if i in list(chain_map.keys())]
                    )
                else:
                    inner_dict[key] = val[relevant_indices]
            final_cat[aff] = inner_dict
            bio_cif_handle.set_category(aff, inner_dict)
        else:
            if cat is not None:
                bio_cif_handle.set_category(aff, cat)
    bio_cif_handle.write(output_cif_path)
    return output_cif_path


def convert_category(category: Dict[str, np.ndarray[int, Any]]) -> Dict[Any, Any]:
    """
    Convert a PDBx/mmCIF category to a dictionary indexed by sequential ids.
    with keys and values taken from the original value arrays.
    """
    category_dict: Dict[Any, Any] = {}
    if category is not None:
        for i in range(len(category[list(category.keys())[0]])):
            category_dict[i] = {}
            for key, value in category.items():
                category_dict[i][key] = value[i]
    return category_dict


def replace_with_nan(value: Any) -> Any:
    if (
        (value == "?")
        or (value == "nan")
        or (value == "NaN")
        or (value == ".")
        or (value == "")
    ):
        return np.nan
    return value


def replace_with_empty_str(value: Any) -> Any:
    if (value == "?") or (value == "NaN") or (value == "."):
        return ""
    return value


def extract_rcsb_info(pdbx_file: PDBxFile) -> Tuple[str, float, str]:
    """ "
    Use full mmcif.gz
    """
    refine_category = convert_category(
        pdbx_file.get_category("refine", expect_looped=True)
    )
    try:
        refine = {k: replace_with_nan(v) for k, v in refine_category[0].items()}
    except AttributeError:
        refine = {}

    exp_method = replace_with_empty_str(str(refine["pdbx_refine_id"]))

    try:
        resolution = float(replace_with_nan(str(refine["ls_d_res_high"])))
    except:
        resolution = np.nan

    date_dict_category = convert_category(
        pdbx_file.get_category("pdbx_database_status", expect_looped=True)
    )

    try:
        date_dict = {k: replace_with_nan(v) for k, v in date_dict_category[0].items()}
    except AttributeError:
        date_dict = {}

    date = replace_with_empty_str(date_dict["recvd_initial_deposition_date"])

    return date, resolution, exp_method


def get_ec(pdbx_file: PDBxFile, pocket_chain_ids: List[str]) -> str:
    """ "
    Use full mmcif.gz
    """
    poly_scheme = convert_category(
        pdbx_file.get_category("pdbx_poly_seq_scheme", expect_looped=True)
    )
    entity = convert_category(pdbx_file.get_category("entity", expect_looped=True))

    entity = {
        value["id"]: value["pdbx_ec"]
        for key, value in entity.items()
        if value.get("pdbx_ec") not in ["?", None]
    }
    poly_scheme_dict = dict()
    for key, value in poly_scheme.items():
        poly_scheme_dict[value["pdb_strand_id"]] = value["entity_id"]

    pocket_chain_ids = [to_alpha(i) for i in pocket_chain_ids]

    auth_chain_ec_number_dict = {
        ch: replace_with_empty_str(entity.get(poly_scheme_dict.get(ch, ""), ""))
        for ch in pocket_chain_ids
    }

    return ";".join(list(auth_chain_ec_number_dict.values()))


def convert_chain(chains: List[str]) -> Dict[str, str]:
    chains = sorted(chains)
    new_chains = {ch: chr(65 + i) for i, ch in enumerate(chains)}
    return new_chains


def get_chain_mapping(cif_path: Path) -> Dict[str, str]:
    cif_file_handle = PDBxFile.read(str(cif_path))
    structure = get_structure(cif_file_handle, model=1)
    return convert_chain(list(set(structure.chain_id)))


def get_ccd_smiles_dict(ciffile: Path) -> dict[str, str]:
    df = pd.read_parquet(ciffile.parent / "components.parquet")
    return dict(zip(df["binder_id"], df["canonical_smiles"]))


def get_entity_type(cif_path: Path, chain_id: str, res_id: int, res_name: str) -> Any:
    model = gemmi.read_structure(str(cif_path))[0]
    return model[chain_id][str(res_id)][res_name].entity_type.name


def extract_rdk_mol_from_cif(
    mmcif_file: Path, res_name: str, res_id: int, chain_id: str, new_chain_id: str = "X"
) -> Mol:
    f = PDBxFile.read(str(mmcif_file))
    atom_arr = get_structure(f, model=1)
    lig_arr = atom_arr[
        (atom_arr.res_name == res_name)
        & (atom_arr.res_id == int(res_id))
        & (atom_arr.chain_id == chain_id)
    ]
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as f:
        file = PDBxFile()
        lig_arr.chain_id = np.array([new_chain_id] * len(lig_arr))
        set_structure(file, lig_arr, data_block="LIG")
        struc.io.save_structure(f.name, lig_arr)
        mol = AllChem.MolFromPDBBlock(f.read(), proximityBonding=False, sanitize=False)
    return mol


def save_new_cif_file(
    cif_path: Path, metadata_df: pd.DataFrame, output_path: Path
) -> None:
    """Save new cif file with metadata
    pocket_chain_binder_residues: List[Tuple[List[str], str, str, str]]
    """
    output_path.mkdir(parents=True, exist_ok=True)
    pdbid = cif_path.stem[:4]
    cif_file_handle = PDBxFile.read(str(cif_path))
    structure = get_structure(cif_file_handle, model=1)
    prot = structure[struc.filter_amino_acids(structure)]

    prot_list = []
    lig_list = []
    pocket_list = metadata_df["binding_site_residues"].to_list()
    ligand_info_list = metadata_df["binder_info"].to_list()

    pocket_chain = [
        list({res.split(",")[-1] for res in pock_str.split(";")})
        for pock_str in pocket_list
    ]
    lig_chain_str = [lig_str.split(",")[-1] for lig_str in ligand_info_list]
    lig_resn_str = [lig_str.split(",")[0] for lig_str in ligand_info_list]
    lig_resi_str = [lig_str.split(",")[1] for lig_str in ligand_info_list]

    pocket_chain_binder_residues = list(
        zip(pocket_chain, lig_resn_str, lig_resi_str, lig_chain_str)
    )

    for chain_ids, lig_resn, lig_resi, lig_chain in pocket_chain_binder_residues:
        prot = prot[(np.isin(prot.chain_id, chain_ids))]
        prot_list.append(prot)

        lig = structure[
            (structure.chain_id == lig_chain)
            & (structure.res_name == lig_resn)
            & (structure.res_id == lig_resi)
        ]
        lig_list.append(lig)

    combined_list = prot_list + lig_list
    combined_list = [atm for protlig in combined_list for atm in protlig]
    combined_struc = struc.array(combined_list)
    metadata_df.fillna("?", inplace=True)
    categories = metadata_df.groupby("pdbid").agg(list).T.to_dict()[pdbid]
    categories = {k: np.array(v) for k, v in categories.items()}
    set_structure(cif_file_handle, combined_struc)
    # Add fingerprints
    cif_file_handle.set_category("plinder_metadata", categories)
    cif_file_handle.write(output_path / f"{pdbid}-processed.cif")


def load_fingerprint_from_base64(base64_fp: str, size: int) -> Any:
    fp_from_base64 = ExplicitBitVect(size)
    fp_from_base64.FromBase64(base64_fp)
    return fp_from_base64


def load_fingerprints_from_cif(
    mmcif_file: Path, lig_fp_size: int = 2024, pli_fp_size: int = 16384
) -> Tuple[str, Any, Any]:
    cif_block = gemmi.cif.read_file(str(mmcif_file)).sole_block()
    plinder_fing = cif_block.find(
        "_plinder_metadata.",
        ["pdb_binder_id", "fingerprint_id", "lig_fingerprint", "pli_fingerprint"],
    )[0]
    pdb_lig_id, finger_id, lig_base64, pli_base64 = (
        plinder_fing["pdb_binder_id"],
        plinder_fing["fingerprint_id"],
        plinder_fing["lig_fingerprint"],
        plinder_fing["pli_fingerprint"],
    )

    lig_fingerprint = load_fingerprint_from_base64(lig_base64, lig_fp_size)
    pli_fingerprint = load_fingerprint_from_base64(pli_base64, pli_fp_size)

    return (f"{pdb_lig_id[:4]}_{finger_id}", lig_fingerprint, pli_fingerprint)


def load_spliting_parameters(
    mmcif_file: Path, lig_fp_size: int = 2024, pli_fp_size: int = 16384
) -> Tuple[str, Any, Any, Any]:
    fp_id, lig_fp, pli_fp = load_fingerprints_from_cif(
        mmcif_file, lig_fp_size, pli_fp_size
    )
    cif_block = gemmi.cif.read_file(str(mmcif_file)).sole_block()
    return (
        fp_id,
        lig_fp,
        pli_fp,
        cif_block.find_value("_plinder_metadata.ecod_t_name"),
    )
