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
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol, RWMol

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