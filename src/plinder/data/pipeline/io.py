# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
Wrap all network requests in a retry decorator
and use a convention of looking for a file in a
pre-determined location before fetching it from
the network.
"""

import gzip
import json
import os
import shutil
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from subprocess import check_output
from typing import Any, Literal, Optional, TypeVar

import requests
from tqdm import tqdm

from plinder.core.utils.io import download_alphafold_cif_file, retry
from plinder.core.utils.log import setup_logger
from plinder.data.pipeline import transform

# RCSB retired its NextGen rsync service in June 2026; PDBj remains an
# official wwPDB rsync mirror for this archive.
CIF_PATH = "rsync-nextgen.pdbj.org::ftp_nextgen/data/entries/divided"
CIF_GLOB = "*-enrich.cif.gz"
VAL_PATH = "rsync.rcsb.org::ftp/validation_reports"
VAL_GLOB = "*_validation.xml.gz"
CIF_PORT = "873"
VAL_PORT = "33444"
KINDS = ["cif", "val"]
KIND_TYPES = Literal["cif", "val"]
LOG = setup_logger(__name__)
T = TypeVar("T")


@retry
def download_cofactors(
    *,
    data_dir: Path,
    url: str = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors",
    force_update: bool = False,
) -> dict[str, Any]:
    """
    Download ligand cofactor data.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    url : str
        URL to fetch data from
    force_update : bool, default=False
        if True, re-download data

    Returns
    -------
    cofactors : dict[str, Any]
        cofactor data
    """
    cofactor_path = data_dir / "dbs" / "cofactors" / "cofactors.json"
    cofactor_path.parent.mkdir(exist_ok=True, parents=True)
    if not cofactor_path.is_file() or force_update:
        LOG.info(f"download_cofactors: {url}")
        resp = requests.get(url)
        resp.raise_for_status()
        obj: dict[str, Any] = resp.json()
        with cofactor_path.open("w") as f:
            json.dump(obj, f, indent=4)
    else:
        with cofactor_path.open() as f:
            obj = json.load(f)
    return obj


@retry
def download_affinity_data(
    *,
    data_dir: Path,
    bindingdb_url: str = "https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202604_tsv.zip",
    force_update: bool = False,
) -> Any:
    """
    Download binding affinity data.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    bindinddb_url : str
        bindingdb : url
    force_update : bool, default=False
        if True, re-download data

    Returns
    -------
    binding_affinity : dict[str, Any]
        binding affinity data
    """
    from io import BytesIO
    from urllib.request import urlopen
    from zipfile import ZipFile

    affinity_path = data_dir / "dbs" / "affinity" / "affinity.json"
    bindingdb_raw_affinity_path = data_dir / "dbs" / "affinity" / "BindingDB_All.tsv"

    # Make sub directories
    bindingdb_raw_affinity_path.parent.mkdir(parents=True, exist_ok=True)
    if not affinity_path.is_file() or force_update:
        # Download BindingDB
        if (
            not bindingdb_raw_affinity_path.is_file()
            or bindingdb_raw_affinity_path.stat().st_size == 0
            or force_update
        ):
            LOG.info(f"download_bindingdb_affinity_data: {bindingdb_url}...")
            with urlopen(bindingdb_url) as zipresp:
                with ZipFile(BytesIO(zipresp.read())) as zfile:
                    zfile.extractall(path=bindingdb_raw_affinity_path.parent)

        LOG.info("transforming_affinity_data: extracting median affinity")
        binding_db_affinity_df = transform.transform_bindingdb_affinity_data(
            raw_affinity_path=bindingdb_raw_affinity_path
        )
        binding_db_affinity_df["preference"] = 1

        all_affinity_df = binding_db_affinity_df.drop_duplicates()
        all_affinity_df = all_affinity_df[all_affinity_df["pchembl"].notna()]

        all_affinity_df = all_affinity_df.loc[
            all_affinity_df.groupby("pdbid_ligid")["preference"].idxmin()
        ]
        all_affinity_df = all_affinity_df.set_index("pdbid_ligid")
        obj = {
            "pchembl": json.loads(all_affinity_df[["pchembl"]].to_json())["pchembl"],
            "target_sequence": json.loads(
                all_affinity_df[["target_sequence"]].to_json()
            )["target_sequence"],
        }
        with affinity_path.open("w") as f:
            json.dump(obj, f, indent=4)
    else:
        with affinity_path.open() as f:
            obj = json.load(f)
    return obj


@retry
def download_components_cif(
    *,
    data_dir: Path,
    url: str = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz",
    force_update: bool = False,
) -> Path:
    """
    Download components cif. Additionally aggregate
    the cif to a dataframe and store as parquet.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    url : str
        URL to fetch data from
    force_update : bool, default=False
        if True, re-download data

    Returns
    -------
    components_path : Path
        path to downloaded components data
    """
    components_path = data_dir / "dbs" / "components" / "components.cif"
    components_path.parent.mkdir(parents=True, exist_ok=True)
    if not components_path.is_file() or force_update:
        LOG.info(f"download_components_cif: {url}")
        resp = requests.get(url)
        resp.raise_for_status()
        gz = components_path.parent / "components.cif.gz"
        gz.write_bytes(resp.content)
        with gzip.open(gz, "rb") as arch:
            with components_path.open("wb") as file:
                shutil.copyfileobj(arch, file)
    components_pqt = data_dir / "dbs" / "components" / "components.parquet"
    if not components_pqt.is_file() or force_update:
        LOG.info(f"download_components_cif: transforming {components_path}")
        df = transform.transform_components_data(raw_components_path=components_path)
        df.to_parquet(components_pqt, index=False)
    return components_path


@retry
def download_seqres_data(
    *,
    data_dir: Path,
    url: str = "https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz",
    force_update: bool = False,
) -> Path:
    """
    Download input for mmseqs database.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    url : str
        URL to fetch data from
    force_update : bool, default=False
        if True, re-download data

    Returns
    -------
    seqres_path : Path
        location of downloaded seqres data
    """
    seqres_path = data_dir / "dbs" / "seqres" / "pdb_seqres.txt.gz"
    seqres_path.parent.mkdir(exist_ok=True, parents=True)
    if not seqres_path.is_file() or force_update:
        LOG.info(f"download_seqres_data: {url}")
        resp = requests.get(url)
        resp.raise_for_status()
        with seqres_path.open("wb") as f:
            f.write(resp.content)
    return seqres_path


@retry
def rsync_rcsb(
    *,
    kind: KIND_TYPES,
    two_char_code: str,
    data_dir: Path,
    pdb_id: Optional[str] = None,
) -> None:
    """
    Download PDB source files from the archive's supported rsync mirrors.

    Parameters
    ----------
    kind : str
        the kind of files to rsync ("cif" or "val")
    two_character_code : str
        two character code to sync
    data_dir : Path
        root directory for local filesystem
    """
    if kind not in KINDS:
        raise ValueError(f"kind={kind} not in {KINDS}")
    suffix = None
    if kind == "cif":
        server = CIF_PATH
        contents = CIF_GLOB
        port = CIF_PORT
        if pdb_id is not None:
            suffix = f"pdb_0000{pdb_id}"
    else:
        server = VAL_PATH
        contents = VAL_GLOB
        port = VAL_PORT
        if pdb_id is not None:
            suffix = pdb_id

    server = f"{server}/{two_char_code}/"
    if suffix is not None:
        server = f"{server}{suffix}"
    dest = f"{data_dir}/{two_char_code}/"
    Path(dest).mkdir(exist_ok=True, parents=True)

    cmd = (
        f"rsync -rlpt -z --delete --port={port} --no-perms "
        f'--include "*/" --include "{contents}" --exclude="*" '
        f"{server} {dest}"
    )
    LOG.info(f"running: {cmd}")
    check_output(cmd, shell=True)


@retry
def list_rcsb(
    *,
    kind: KIND_TYPES,
    two_char_code: Optional[str] = None,
    pdb_id: Optional[str] = None,
) -> list[str]:
    if kind not in KINDS:
        raise ValueError(f"kind={kind} not in {KINDS}")
    if kind == "cif":
        server = CIF_PATH
        port = CIF_PORT
    else:
        server = VAL_PATH
        port = VAL_PORT
    if two_char_code is not None:
        server = f"{server}/{two_char_code}/"
        if pdb_id is not None:
            server = f"{server}{pdb_id}"
    else:
        server = f"{server}/"
    cmd = f"rsync --port={port} --list-only {server}"
    LOG.info(f"running: {cmd}")
    output = check_output(cmd, shell=True, text=True).splitlines()
    return [
        line.split()[-1]
        for line in output
        if line.startswith("d") and not line.endswith(".")
    ]


def get_missing_two_char_codes(
    *,
    kind: KIND_TYPES,
    data_dir: Path,
    two_char_codes: list[str],
) -> list[str]:
    missing = []
    glob = CIF_GLOB if kind == "cif" else VAL_GLOB
    for two_char_code in tqdm(two_char_codes):
        two_char_dir = data_dir / two_char_code
        if not two_char_dir.is_dir():
            missing.append(two_char_code)
            continue
        two_char_cifs = list(two_char_dir.glob(f"*/{glob}"))
        two_char_entries = list_rcsb(kind=kind, two_char_code=two_char_code)
        if len(two_char_cifs) != len(two_char_entries):
            delta = len(two_char_entries) - len(two_char_cifs)
            LOG.info(f"two_char_code={two_char_code} missing {delta} files!")
            missing.append(two_char_code)
    return missing


def get_missing_pdb_ids(
    *,
    kind: KIND_TYPES,
    data_dir: Path,
    two_char_code: str,
) -> list[str]:
    missing = []
    glob = CIF_GLOB if kind == "cif" else VAL_GLOB
    for pdb_id in list_rcsb(kind=kind, two_char_code=two_char_code):
        pdb_dir = data_dir / two_char_code / pdb_id
        if not pdb_dir.is_dir():
            missing.append(pdb_id)
            continue
        if not len(list(pdb_dir.glob(glob))):
            missing.append(pdb_id)
    return missing


@retry
def download_alphafold_cif_files(
    *,
    data_dir: Path,
    url: str = "https://alphafold.ebi.ac.uk/files",
    uniprot_url: str = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=list&query=database:pdb",
    force_update: bool = False,
    threads: int = 10,
    timeout: int = 7200,
) -> Path:
    """
    Fetch AlphaFold CIF file

    Parameters
    ----------
    data_dir : Path
        root directory for local filesystem
    uniprot_id: str
        UniProt ID
    """
    cif_dir = data_dir / "dbs" / "alphafold"
    cif_dir.mkdir(exist_ok=True, parents=True)
    if len(os.listdir(cif_dir)) > 50_000:
        return cif_dir
    uniprot_ids_path = data_dir / "dbs" / "uniprot" / "pdb_uniprot.txt"
    (data_dir / "dbs" / "uniprot").mkdir(exist_ok=True, parents=True)
    if not uniprot_ids_path.exists() or force_update:
        resp = requests.get(uniprot_url)
        resp.raise_for_status()
        with open(uniprot_ids_path, "w") as f:
            f.write(resp.text)
    with open(uniprot_ids_path) as f:
        uniprot_ids = [uniprot_id.strip() for uniprot_id in f]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                download_alphafold_cif_file,
                uniprot_id,
                cif_dir,
                url,
                force_update,
            )
            for uniprot_id in uniprot_ids
        ]
        wait(futures, timeout=timeout, return_when=ALL_COMPLETED)
        for future in futures:
            exc = future.exception()
            if exc is not None:
                raise exc
    return cif_dir


@retry
def download_uniprot_fasta_data(
    *,
    data_dir: Path,
    url: str = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=database:pdb",
    force_update: bool = False,
) -> Path:
    """
    Download input for mmseqs database.

    Parameters
    ----------
    data_dir : Path
        the root plinder dir
    url : str
        URL to fetch data from
    force_update : bool, default=False
        if True, re-download data

    Returns
    -------
    uniprot_fasta_path : Path
        location of downloaded sequence data for UniProts matching
    """
    uniprot_fasta_path = data_dir / "dbs" / "uniprot" / "pdb_uniprot.txt.gz"
    uniprot_fasta_path.parent.mkdir(exist_ok=True, parents=True)
    if not uniprot_fasta_path.is_file() or force_update:
        resp = requests.get(url)
        resp.raise_for_status()
        with uniprot_fasta_path.open("wb") as f:
            f.write(resp.content)
    return uniprot_fasta_path
