# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import logging
import shutil
import subprocess
from itertools import zip_longest
from pathlib import Path
from typing import Any, Generator

import pandas as pd
from rdkit import Chem

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


def pad_integer_lig_ids_with_zeros(lig_id: str) -> str:
    """
    Fixes cases where integer ligand ids are striped\
    of their leading "0s". E.g  1 -> 001; 59 ->059

    Parameters
    ----------
    lig_id : str
        ligand ccd codes
    Returns
    -------
    str
    """
    try:
        lig_id = str(int(lig_id))
        return f"00{lig_id}"[-3:]
    except:
        return str(lig_id)


def make_mmp_index_from_smiles_df(
    smiles_df: pd.DataFrame,
    output_path: Path,
    database_name: str,
    scratch_dir: Path,
    nthreads: int = 20,
) -> Path:
    output_path.mkdir(parents=True, exist_ok=True)
    smiles_df = (
        smiles_df.groupby("ligand_ccd_code")
        .first()
        .reset_index()[["ligand_smiles", "ligand_ccd_code"]]
    )
    # remove rows with no smiles
    smiles_df = smiles_df[smiles_df.ligand_smiles.apply(lambda x: str(x) != "")]
    # Remove compounds containing metal-dative bond
    smiles_df = smiles_df[
        ~(
            smiles_df.ligand_smiles.str.contains(">")
            | smiles_df.ligand_smiles.str.contains("<")
        )
    ]
    smiles_df.dropna().to_csv(
        f"{output_path}/{database_name}_input.smi", sep=" ", index=False
    )

    # Split smiles file
    split_smile_cmd = (
        f"mmpdb smi_split {database_name}_input.smi --has-header -n {nthreads}"
    )

    try:
        subprocess.check_output(
            split_smile_cmd, shell=True, stderr=subprocess.STDOUT, cwd=output_path
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Smiles database split failed {e.output} ...")
        raise

    fragment_tasks = [
        f"mmpdb fragment -j 1 {smi.name}" for smi in output_path.glob("*.*.smi")
    ]
    with (scratch_dir / "fragment_tasks.txt").open("w") as f:
        f.write("\n".join(fragment_tasks))

    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "plinder.data.pipeline.mpqueue",
                f"{scratch_dir / 'fragment_tasks.txt'}",
                "--cwd",
                str(output_path),
                "--cores",
                str(nthreads),
            ],
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Database fragmentation failed {e.output} ...")
        raise
    (scratch_dir / "fragment_tasks.txt").unlink()

    # Partition fragment database
    fragdbs = list(output_path.glob("*.*.fragdb"))
    partition_tasks = [
        f"mmpdb fragdb_partition *.*.fragdb --task-id {i + 1} --num-tasks {nthreads + 1}"
        for i in range(len(fragdbs))
    ]
    with (scratch_dir / "partition_tasks.txt").open("w") as f:
        f.write("\n".join(partition_tasks))
    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "plinder.data.pipeline.mpqueue",
                f"{scratch_dir / 'partition_tasks.txt'}",
                "--cwd",
                str(output_path),
                "--cores",
                str(nthreads),
            ],
            stderr=subprocess.STDOUT,
            cwd=output_path,
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Fragment database partitioning failed {e.output} ...")
        raise
    (scratch_dir / "partition_tasks.txt").unlink()

    index_tasks = [
        f"mmpdb index {partition.name} -o {partition.name}.csv.gz"
        for partition in output_path.glob("*partition.*.fragdb")
    ]
    with (scratch_dir / "index_tasks.txt").open("w") as f:
        f.write("\n".join(index_tasks))

    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "plinder.data.pipeline.mpqueue",
                f"{scratch_dir / 'index_tasks.txt'}",
                "--cwd",
                str(output_path),
                "--cores",
                str(nthreads),
            ],
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Partition indexing  failed {e.output} ...")
        raise
    (scratch_dir / "index_tasks.txt").unlink()

    # Merge
    with (output_path / f"{database_name}.csv.gz").open("wb") as db:
        for partition in output_path.glob("partition.*.csv.gz"):
            LOG.info(f"catting {partition} into {database_name}")
            with partition.open("rb") as piece:
                shutil.copyfileobj(piece, db)

    return output_path / f"{database_name}.csv.gz"


def make_mmp_index_from_annotation_table(
    data_dir: Path,
    annotation_df: pd.DataFrame,
    database_name: str = "plinder_mms",
) -> Path:
    output_path = data_dir / "mmp"
    scratch_dir = data_dir / "scratch" / "mmp"
    output_path.mkdir(exist_ok=True, parents=True)
    scratch_dir.mkdir(exist_ok=True, parents=True)
    rename = {}
    if "ligand_smiles" not in annotation_df.columns:
        rename["ligand_rdkit_canonical_smiles"] = "ligand_smiles"
    if "ligand_ccd_code" not in annotation_df.columns:
        rename["ligand_unique_ccd_code"] = "ligand_ccd_code"
    smiles_df = annotation_df.rename(columns=rename)
    smiles_df = smiles_df[~smiles_df.ligand_smiles.isna()]
    logging.info(f"smiles Dataframe output {smiles_df}")
    return make_mmp_index_from_smiles_df(
        smiles_df, output_path, database_name, scratch_dir
    )


def add_mmp_clusters_to_data(
    load_mmp_df: pd.DataFrame,
    system_df: pd.DataFrame,
    cluster_folder: Path,
    protein_metric: str = "protein_fident_weighted_sum",
    protein_threshold: int = 95,
    protein_directed: bool = False,
    pocket_metric: str = "pocket_fident_qcov",
    pocket_threshold: int = 100,
    pocket_directed: bool = False,
    min_constant_size: int = 10,
) -> pd.DataFrame:
    """Add mmpdb data to pocket and protein similarity dataset.

    Parameters
    ----------
    load_mmp_df : Path
        mmpdb index dataframe with columns \
        ["SMILES1",  "SMILES2",  "id1",  "id2", "V1>>V2", "CONSTANT"]
    system_df : pd.DataFram
        Systems dataframe
    protein_similarity_path : Path
        protein sequence identity cluster file path
    pocket_similarity_path : Path
        pocket sequence identity cluster file path
    protein_similarity_tag : str
        protein similarity tag (e.g protein_fident_weighted_sum__0.95__weak__component)
    pocket_similarity_tag : str
        pocket similarity tag (e.g pocket_fident_weighted_sum__1.0__strong__component)
    min_constant_size : int
        minimum constant size

    Returns
    -------
    pd.DataFrame
    """
    protein_component_type = "strong" if protein_directed else "weak"
    protein_similarity_tag = (
        f"{protein_metric}__{protein_threshold}__{protein_component_type}__component"
    )
    pocket_component_type = "strong" if pocket_directed else "weak"
    pocket_similarity_tag = (
        f"{pocket_metric}__{pocket_threshold}__{pocket_component_type}__component"
    )

    # Load protein and pocket similarity data
    protein_similarity_df = pd.read_parquet(
        cluster_folder
        / "cluster=components"
        / f"directed={protein_directed}"
        / f"metric={protein_metric}"
        / f"threshold={protein_threshold}.parquet"
    )
    protein_similarity_dict = dict(
        zip(protein_similarity_df["system_id"], protein_similarity_df["label"])
    )
    pocket_similarity_df = pd.read_parquet(
        cluster_folder
        / "cluster=components"
        / f"directed={pocket_directed}"
        / f"metric={pocket_metric}"
        / f"threshold={pocket_threshold}.parquet"
    )
    pocket_similarity_dict = dict(
        zip(pocket_similarity_df["system_id"], pocket_similarity_df["label"])
    )

    # Merge protein similarity with system df
    system_df[protein_similarity_tag] = system_df["system_id"].map(
        protein_similarity_dict
    )

    # Merge pocket similarity with system df
    system_df[pocket_similarity_tag] = system_df["system_id"].map(
        pocket_similarity_dict
    )

    # Load mmp index
    # load_mmp_df = pd.read_csv(mmp_index, compression="gzip", sep="\t", header=None)
    # load_mmp_df.columns = ["SMILES1", "SMILES2", "id1", "id2", "V1>>V2", "CONSTANT"]

    # Pad interger lig ids with zeros
    load_mmp_df["id1"] = load_mmp_df.id1.apply(pad_integer_lig_ids_with_zeros)
    load_mmp_df["id2"] = load_mmp_df.id2.apply(pad_integer_lig_ids_with_zeros)

    # Consider only single position changes to the constant!
    load_mmp_df = load_mmp_df[
        load_mmp_df.CONSTANT.apply(lambda x: x.count("*")) == 1
    ].copy()

    # Select relevant columns from system df
    pocket_df = system_df[
        ["ligand_unique_ccd_code", protein_similarity_tag, pocket_similarity_tag]
    ].copy()

    pocket_df.loc[:, "prot_pocket_id"] = (
        pocket_df.loc[:, protein_similarity_tag]
        + "_"
        + pocket_df.loc[:, pocket_similarity_tag]
    )

    # Map ligand ccd code to prot_pocket_id sets
    ligand_code_pocket_mapping = (
        pocket_df.groupby("ligand_unique_ccd_code").agg(set)["prot_pocket_id"].to_dict()
    )

    load_mmp_df["prot_pocket_set_id1"] = load_mmp_df["id1"].map(
        ligand_code_pocket_mapping
    )
    load_mmp_df["prot_pocket_set_id2"] = load_mmp_df["id2"].map(
        ligand_code_pocket_mapping
    )

    load_mmp_df.dropna(inplace=True)

    # Get MMPs that share prot_pocket_id
    load_mmp_df["prot_pocket_set_shared"] = load_mmp_df[
        ["prot_pocket_set_id1", "prot_pocket_set_id2"]
    ].apply(lambda x: x.iloc[0].intersection(x.iloc[1]), axis=1)
    # remove pairs that do not share a prot_pocket id
    mmps_pocket_df = load_mmp_df[
        load_mmp_df["prot_pocket_set_shared"].apply(lambda x: len(x) > 0)
    ]

    # Group pocket across different mmps
    mmps_pocket_df1 = mmps_pocket_df.explode("prot_pocket_set_shared")

    # Identity congeneric series - MMS - group that shares
    # identical constant (with a single vector) and prot_pockets !
    grp_congeneric_df = mmps_pocket_df1.groupby(
        ["CONSTANT", "prot_pocket_set_shared"]
    ).agg(tuple)[["id1", "id2"]]
    # set to tuple for being hashable
    grp_congeneric_df["congeneric_series"] = grp_congeneric_df[["id1", "id2"]].apply(
        lambda x: tuple(sorted(set(list(x.iloc[0]) + list(x.iloc[1])))), axis=1
    )

    # NOTE: some ligands will appear in multiple instances!
    # TODO: be careful when dropping during BO to see if there were multiple and if a ligand is to be removed!
    grp_congeneric_df["mms_unique_count"] = grp_congeneric_df[
        "congeneric_series"
    ].apply(lambda x: len({i for i in x}))

    # Reset index to include "CONSTANT" as a column
    grp_congeneric_df = grp_congeneric_df.reset_index()

    # Get constant size
    grp_congeneric_df["const_ROMol"] = grp_congeneric_df.CONSTANT.apply(
        Chem.MolFromSmarts
    )
    grp_congeneric_df["const_size"] = grp_congeneric_df["const_ROMol"].apply(
        lambda x: x.GetNumHeavyAtoms()
    )

    # Extract data with min_constant_size threshold
    grp_congeneric_df_thres = grp_congeneric_df[
        grp_congeneric_df.const_size >= min_constant_size
    ].copy()

    # Set congeneric id that is unique to mms
    grp_congeneric_df_thres = grp_congeneric_df_thres.reset_index().rename(
        columns={"index": "congeneric_id"}
    )
    grp_congeneric_df_thres["congeneric_id"] = grp_congeneric_df_thres[
        "congeneric_id"
    ].apply(lambda x: f"c{x:.0f}")

    # get table for each member in congeneric_series
    grp_congeneric_df_thres["congeneric_ligand_ccd_code"] = grp_congeneric_df_thres[
        "congeneric_series"
    ].copy()
    grp_congeneric_df_thres = grp_congeneric_df_thres.explode(
        "congeneric_ligand_ccd_code"
    )

    # Merge mms data with system-level data
    grp_congeneric_df_thres["prot_pocket_lig_tag"] = (
        grp_congeneric_df_thres["prot_pocket_set_shared"]
        + "_"
        + grp_congeneric_df_thres["congeneric_ligand_ccd_code"]
    )
    system_df["prot_pocket_lig_tag"] = system_df[
        [protein_similarity_tag, pocket_similarity_tag, "ligand_unique_ccd_code"]
    ].apply(lambda x: f"{x.iloc[0]}_{x.iloc[1]}_{x.iloc[2]}", axis=1)

    final_df = (
        system_df[["prot_pocket_lig_tag", "system_id"]]
        .drop_duplicates()
        .merge(grp_congeneric_df_thres, on="prot_pocket_lig_tag", how="left")
        .drop_duplicates()
        .drop(columns=["const_ROMol", "id1", "id2", "prot_pocket_lig_tag"])
    )
    # remove those that are not mapped!
    final_df = final_df[~final_df["congeneric_id"].isna()]

    return final_df


def split_list_into_batches(
    lst: list[list[Any]], batch_size: int | None = None, num_batches: int | None = None
) -> Generator[Any, Any, Any]:
    """Split list of items into list of lists.

    Parameters
    ----------
    lst : list[list[Any]],
        list to be split
    batch_size : int | None = None,
        size of each batch
    num_batches : int | None = None,
        number of batches to return
    Returns
    -------
    Generator[Any, Any, Any]
    """
    if batch_size is None and num_batches is None:
        raise ValueError("Either batch_size or num_batches must be provided.")

    if num_batches is not None:
        # Ceiling division to ensure all items are included
        batch_size = -(-len(lst) // num_batches)

    args = [iter(lst)] * batch_size  # type: ignore
    for batch in zip_longest(*args, fillvalue=None):
        yield [item for item in batch if item is not None]
