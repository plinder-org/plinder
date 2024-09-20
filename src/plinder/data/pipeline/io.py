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

import pandas as pd
import requests
from tqdm import tqdm

from plinder.core.utils.io import download_alphafold_cif_file, retry
from plinder.core.utils.log import setup_logger
from plinder.data.pipeline import transform

CIF_PATH = "rsync-nextgen.wwpdb.org::rsync/data/entries/divided"
CIF_GLOB = "*-enrich.cif.gz"
VAL_PATH = "rsync.rcsb.org::ftp/validation_reports"
VAL_GLOB = "*_validation.xml.gz"
RCSB_PORT = "33444"
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
    bindingdb_url: str = "https://www.bindingdb.org/bind/downloads/BindingDB_All_202401_tsv.zip",
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
    papyrus_raw_affinity_path = (
        data_dir / "dbs" / "affinity" / "papyrus_affinity_raw.tar.gz"
    )
    bindingdb_raw_affinity_path = (
        data_dir / "dbs" / "affinity" / "BindingDB_All_202401.tsv"
    )
    moad_raw_affinity_path = data_dir / "dbs" / "affinity" / "moad_affinity.csv"

    # Make sub directories
    papyrus_raw_affinity_path.parent.mkdir(parents=True, exist_ok=True)
    bindingdb_raw_affinity_path.parent.mkdir(parents=True, exist_ok=True)
    moad_raw_affinity_path.parent.mkdir(parents=True, exist_ok=True)
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
        affinity_json = all_affinity_df[["pchembl"]].to_json()
        obj: dict[str, Any] = json.loads(affinity_json)
        with affinity_path.open("w") as f:
            json.dump(obj, f, indent=4)
    else:
        with affinity_path.open() as f:
            obj = json.load(f)
    return obj["pchembl"]


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
def download_ecod_data(
    *,
    data_dir: Path,
    url: str = "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt",
    force_update: bool = False,
) -> Path:
    """
    Download ECOD data.

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
    ecod_path : Path
        path to downloaded ECOD data
    """
    raw_ecod_path = data_dir / "dbs" / "ecod" / "ecod_raw.tsv"
    raw_ecod_path.parent.mkdir(parents=True, exist_ok=True)
    if not raw_ecod_path.is_file() or force_update:
        LOG.info(f"download_ecod_data: {url}")
        resp = requests.get(url)
        resp.raise_for_status()
        raw_ecod_path.write_text(resp.text)
    ecod_path = data_dir / "dbs" / "ecod" / "ecod.parquet"
    if not ecod_path.is_file() or force_update:
        LOG.info(f"download_ecod_data: transforming {raw_ecod_path}")
        ecod = transform.transform_ecod_data(raw_ecod_path=raw_ecod_path)
        ecod.to_parquet(ecod_path, index=False)
    return ecod_path


@retry
def download_panther_data(
    *,
    data_dir: Path,
    url: str = "http://data.pantherdb.org/ftp/generic_mapping/panther_classifications.tar.gz",
    force_update: bool = False,
) -> Path:
    """
    Download panther data. Also
    partitions panther by trailing character
    in uniprot code.

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
    panther_path : Path
        path to downloaded panther data
    """
    raw_panther_path = data_dir / "dbs" / "panther" / "panther_raw.tar.gz"
    raw_panther_path.parent.mkdir(parents=True, exist_ok=True)
    if not raw_panther_path.is_file() or force_update:
        LOG.info(f"download_panther_data: {url}")
        resp = requests.get(url)
        resp.raise_for_status()
        raw_panther_path.write_bytes(resp.content)
    panther_path = data_dir / "dbs" / "panther" / "panther.parquet"
    if not panther_path.is_file() or force_update:
        LOG.info(f"download_panther_data: transforming {raw_panther_path}")
        panther = transform.transform_panther_data(raw_panther_path=raw_panther_path)
        panther["shard"] = panther["uniprotac"].str[-1]
        panther.to_parquet(raw_panther_path.parent / "panther.parquet", index=False)
        for shard, grp in panther.groupby("shard"):
            panther_shard = data_dir / "dbs" / "panther" / f"panther_{shard}.parquet"
            grp.to_parquet(panther_shard, index=False)
    return panther_path


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
def download_kinase_data(
    *,
    data_dir: Path,
    url: str = "https://klifs.net/api/",
    force_update: bool = False,
) -> Path:
    """
    Download kinase data.

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
    kinase_uniprotac_path : Path
        location of downloaded kinase data
    """
    kinase_info_path = data_dir / "dbs" / "kinase" / "kinase_information.parquet"
    if not kinase_info_path.parent.exists() or force_update:
        LOG.info(f"download_kinase_data: data_dir={data_dir}")
        kinase_info_path.parent.mkdir(exist_ok=True, parents=True)
    if not kinase_info_path.is_file() or force_update:
        LOG.info(f"download_kinase_data: {url}/kinase_information")
        resp = requests.get(f"{url}/kinase_information")
        resp.raise_for_status()
        kinase_info = pd.DataFrame(resp.json(), dtype=str)
        kinase_info.to_parquet(kinase_info_path, index=False)
    kinase_ligand_path = data_dir / "dbs" / "kinase" / "kinase_ligand_ccd_codes.parquet"
    if not kinase_ligand_path.is_file() or force_update:
        LOG.info(f"download_kinase_data: {url}/ligands_list")
        resp = requests.get(f"{url}/ligands_list")
        resp.raise_for_status()
        kinase_ligand = pd.DataFrame(resp.json(), dtype=str)
        kinase_ligand.to_parquet(kinase_ligand_path, index=False)
    kinase_struc_path = kinase_info_path.parent / "kinase_structures.parquet"
    if not kinase_struc_path.is_file() or force_update:
        kinase_info = pd.read_parquet(kinase_info_path)
        outputs = []
        LOG.info(f"download_kinase_data: {url}/structures_list iteratively...")
        for kid in kinase_info["kinase_ID"].unique():
            resp = requests.get(f"{url}/structures_list", params={"kinase_ID": kid})
            if 200 <= resp.status_code < 400:
                outputs.extend(resp.json())
        kinase_struc = pd.DataFrame(outputs, dtype=str).drop_duplicates()
        kinase_struc.to_parquet(kinase_struc_path, index=False)
    kinase_uniprotac_path = kinase_info_path.parent / "kinase_uniprotac.parquet"
    if not kinase_uniprotac_path.is_file() or force_update:
        LOG.info("download_kinase_data: merging information and structure")
        kinase_info = pd.read_parquet(kinase_info_path)
        kinase_struc = pd.read_parquet(kinase_struc_path)
        df = pd.merge(
            kinase_struc,
            kinase_info[
                [
                    "kinase_ID",
                    "name",
                    "HGNC",
                    "family",
                    "group",
                    "kinase_class",
                    "full_name",
                    "uniprot",
                ]
            ],
            on="kinase_ID",
            how="left",
        )
        df["pdb"] = df["pdb"].astype(str).map(str.lower)
        df["pdbid_chainid"] = df["pdb"].str.cat(df["chain"], sep="_")
        for int_col in [
            "structure_ID",
            "kinase_ID",
            "missing_atoms",
            "missing_residues",
        ]:
            df[int_col] = df[int_col].astype(int)
        for float_col in [
            "rmsd1",
            "rmsd2",
            "resolution",
            "quality_score",
            "Grich_distance",
            "Grich_rotation",
            "Grich_angle",
        ]:
            df[float_col] = df[float_col].astype(float)
        df.to_parquet(kinase_uniprotac_path, index=False)
    return kinase_uniprotac_path


@retry
def rsync_rcsb(
    *,
    kind: KIND_TYPES,
    two_char_code: str,
    data_dir: Path,
    pdb_id: Optional[str] = None,
) -> None:
    """
    Run an RCSB rsync command.

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
        if pdb_id is not None:
            suffix = f"pdb_0000{pdb_id}"
    else:
        server = VAL_PATH
        contents = VAL_GLOB
        if pdb_id is not None:
            suffix = pdb_id

    server = f"{server}/{two_char_code}/"
    if suffix is not None:
        server = f"{server}{suffix}"
    dest = f"{data_dir}/{two_char_code}/"
    Path(dest).mkdir(exist_ok=True, parents=True)

    cmd = (
        f"rsync -rlpt -z --delete --port={RCSB_PORT} --no-perms "
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
    server = CIF_PATH if kind == "cif" else VAL_PATH
    if two_char_code is not None:
        server = f"{server}/{two_char_code}/"
        if pdb_id is not None:
            server = f"{server}{pdb_id}"
    else:
        server = f"{server}/"
    cmd = f"rsync --port={RCSB_PORT} --list-only {server}"
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
