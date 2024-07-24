# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def parse_pdb_range(pdb_range: str) -> tuple[str, str, str]:
    """
    Parse the PDB ranges from the raw ECOD dataset.
    Default to ".", "0", "0" if it's unparseable.

    Parameters
    ----------
    pdb_range : str
        "-"-delimited range values

    Returns
    -------
    chain, start, end : tuple[str, str, str]
    """
    if ":" not in pdb_range:
        return ".", "0", "0"
    try:
        chain, range_str = pdb_range.split(":")
    except ValueError:
        return ".", "0", "0"
    ranges = range_str.split("-")
    # single range value
    if len(ranges) == 1:
        fr = to = ranges[0]
    # two positive range values
    elif len(ranges) == 2:
        fr, to = ranges
    # one negative range value
    elif len(ranges) == 3:
        if not ranges[1]:
            fr, to = ranges[0], f"-{ranges[-1]}"
        else:
            fr, to = f"-{ranges[1]}", ranges[-1]
    # two negative range values
    elif len(ranges) == 4:
        fr, to = f"-{ranges[1]}", f"-{ranges[-1]}"
    else:
        raise Exception(f"Could not parse: {pdb_range}")
    return chain, fr, to


def transform_ecod_data(*, raw_ecod_path: Path) -> pd.DataFrame:
    """
    Convert the raw ECOD data to a dataset more
    amenable to the plinder annotation data.

    Parameters
    ----------
    df : pd.DataFrame
        the raw ECOD dataset

    Returns
    -------
    transformed : pd.DataFrame
        the transformed ECOD dataset
    """
    df = pd.read_csv(raw_ecod_path, sep="\t", header=1, skiprows=3)
    domains = []
    for pdb, chain, pdb_ranges, ecod_domain_id, t_name in zip(
        df["pdb"],
        df["chain"],
        df["pdb_range"],
        df["ecod_domain_id"],
        df["t_name"],
    ):
        for pdb_range in pdb_ranges.split(","):
            chain, fr, to = parse_pdb_range(pdb_range)
            domains.append((pdb, chain, ecod_domain_id, t_name, fr, to))
    cols = ["pdb", "chain", "domainid", "domain", "pdb_from", "pdb_to"]
    return pd.DataFrame(domains, columns=cols)


def transform_panther_data(*, raw_panther_path: Path) -> pd.DataFrame:
    """
    Unpack the tarball archive and collect the
    contained files to a single parquet file.

    Parameters
    ----------
    raw_panther_path : Path
        location of panther tarball

    Returns
    -------
    transformed : pd.DataFrame
        the collated panther dataset
    """
    panther_dir = raw_panther_path.parent
    with tarfile.open(raw_panther_path) as arch:
        arch.extractall(panther_dir)
    dfs = []
    names = ["uniprotac", "panther"]
    for path in (panther_dir / "panther_classifications").rglob("*"):
        dfs.append(pd.read_csv(path, delimiter="\t", header=None, names=names))
        dfs[-1]["panther_class"] = path.name
    return pd.concat(dfs).drop_duplicates()


def transform_bindingdb_affinity_data(*, raw_affinity_path: Path) -> pd.DataFrame:
    """
    Unpack the tarball archive and collect the
    contained files to a single parquet file.

    Parameters
    ----------
    raw_affinity_path : Path
        location of affinity data

    Returns
    -------
    transformed : pd.DataFrame
        median affinity dataset
    """

    def calc_pchembl(affinity: float) -> Any:
        affinity = affinity * 10**-9
        if affinity > 0:
            return -1.0 * np.log10(affinity)
        else:
            return np.nan

    cols = [
        "Ligand HET ID in PDB",
        "PDB ID(s) for Ligand-Target Complex",
        "Ki (nM)",
        # "IC50 (nM)",
        "Kd (nM)",
        "EC50 (nM)",
    ]
    df = pd.read_csv(raw_affinity_path, sep="\t", usecols=cols, low_memory=False)

    df["pchembl"] = (
        # df[["Ki (nM)", "IC50 (nM)", "Kd (nM)"]]
        df[["Ki (nM)", "Kd (nM)"]]
        .apply(set, axis=1)
        .apply(lambda x: [i for i in x if str(i) != "nan"])
    )
    df = df[df["pchembl"].apply(lambda x: x != [])]
    df["pchembl"] = df["pchembl"].apply(
        lambda x: calc_pchembl(float(str(x[0]).replace(">", "").replace("<", "")))
    )

    df = df[
        ["PDB ID(s) for Ligand-Target Complex", "Ligand HET ID in PDB", "pchembl"]
    ].drop_duplicates()
    df = df[
        (df["Ligand HET ID in PDB"].notna())
        & (df["PDB ID(s) for Ligand-Target Complex"].notna())
    ]
    df["pdb_id"] = df["PDB ID(s) for Ligand-Target Complex"].apply(
        lambda x: x.split(",")
    )
    df = df.explode(["pdb_id"]).drop_duplicates()
    df["pdbid_ligid"] = (
        df["pdb_id"].str.upper() + "_" + df["Ligand HET ID in PDB"].str.strip()
    )
    df = df[["pdbid_ligid", "pchembl"]].drop_duplicates()

    return df.groupby("pdbid_ligid").median().reset_index()


def transform_papyrus_affinity_data(*, raw_affinity_path: Path) -> pd.DataFrame:
    """
    Unpack the tarball archive and collect the
    contained files to a single parquet file.

    Parameters
    ----------
    raw_affinity_path : Path
        location of affinity data

    Returns
    -------
    transformed : pd.DataFrame
        median affinity dataset
    """
    df = pd.read_csv(raw_affinity_path, sep="\t", compression="zip")
    affinity_df = (
        df[
            [
                "accession",
                "Quality",
                "source",
                "pchembl_value_Median",
                "PDBID_ligand",
                "PDBID_protein",
            ]
        ]
        .copy()
        .rename(columns={"pchembl_value_Median": "pchembl"})
    )
    affinity_df["PDBID_protein"] = affinity_df["PDBID_protein"].apply(
        lambda x: x.split(";")
    )
    affinity_df = affinity_df.explode("PDBID_protein")
    affinity_df = affinity_df[affinity_df.pchembl.notna()]
    affinity_df["pdbid_ligid"] = (
        affinity_df["PDBID_protein"].str.upper() + "_" + affinity_df["PDBID_ligand"]
    )
    return (
        affinity_df[["pdbid_ligid", "pchembl"]]
        .groupby("pdbid_ligid")
        .median()
        .reset_index()
    )


def transform_moad_affinity_data(*, raw_affinity_path: Path) -> pd.DataFrame:
    """
    Unpack the tarball archive and collect the
    contained files to a single parquet file.

    Parameters
    ----------
    raw_affinity_path : Path
        location of affinity data

    Returns
    -------
    transformed : pd.DataFrame
        median affinity dataset
    """

    def calc_pchembl(affinity: float, unit: str) -> Any:
        if unit == "fM":
            affinity = affinity * 10**-15

            if affinity > 0:
                return -1.0 * np.log10(affinity)
            else:
                return np.nan
        elif unit == "pM":
            affinity = affinity * 10**-12
            if affinity > 0:
                return -1.0 * np.log10(affinity)
            else:
                return np.nan
        elif unit == "nM":
            affinity = affinity * 10**-9
            if affinity > 0:
                return -1.0 * np.log10(affinity)
            else:
                return np.nan
        elif unit == "uM":
            affinity = affinity * 10**-6
            if affinity > 0:
                return -1.0 * np.log10(affinity)
            else:
                return np.nan
        elif unit == "mM":
            affinity = affinity * 10**-3
            if affinity > 0:
                return -1.0 * np.log10(affinity)
            else:
                return np.nan
        elif unit == "M":
            return affinity

    with open(raw_affinity_path) as f:
        combined_list = []
        for line in f.readlines():
            line_split = line.split(",")
            tmp_enzyme_class = line_split[0]
            tmp_pdbid = line_split[2]
            if len(tmp_enzyme_class.split(".")) == 4:
                new_enzyme_class = tmp_enzyme_class
            if len(tmp_pdbid) > 0:
                if "Family" in line_split[1]:
                    family_representative = True
                else:
                    family_representative = False
                new_pdbid = tmp_pdbid
            if (line_split[3] != "") & (line_split[5] != "Ka"):
                combined_list.append(
                    [
                        new_enzyme_class,
                        family_representative,
                        new_pdbid,
                        line_split[3],
                        line_split[4],
                        line_split[7],
                        line_split[8],
                        line_split[9],
                    ]
                )
    moad_df = pd.DataFrame(
        combined_list,
        columns=[
            "ec_no.",
            "ec_family_rep",
            "pdbid",
            "binder_and_chain",
            "valid_ligand",
            "affinity",
            "unit",
            "smiles",
        ],
    )
    moad_df["pdbid"] = moad_df["pdbid"].str.lower()

    moad_df["binder_id"] = moad_df["binder_and_chain"].apply(lambda x: x.split(":")[0])
    moad_df["binder_id"] = moad_df["binder_id"].apply(lambda x: x.split())
    moad_df = moad_df.explode("binder_id")
    moad_df["pdbid_ligid"] = moad_df["pdbid"].str.upper() + "_" + moad_df["binder_id"]
    # This will set instances with undefined affinity to nan
    moad_df["pchembl"] = moad_df[["affinity", "unit"]].apply(
        lambda x: calc_pchembl(float(x[0]), x[1]) if x[0] != "" else np.nan, axis=1
    )
    return (
        moad_df[["pdbid_ligid", "pchembl"]]
        .groupby("pdbid_ligid")
        .median()
        .reset_index()
    )


def transform_components_data(*, raw_components_path: Path) -> pd.DataFrame:
    import gemmi

    data = gemmi.cif.read_file(raw_components_path.as_posix())
    rows = []
    for block in data:
        (
            binder_id,
            chemical_name,
            molecular_weight,
        ) = block.find("_chem_comp.", ["id", "name", "formula_weight"])[0]
        canonical_smiles, isomeric_smiles, inchikey = None, None, None
        for desc_row in block.find(
            "_pdbx_chem_comp_descriptor.",
            ["comp_id", "type", "program", "descriptor"],
        ):
            if (desc_row[1].strip() == "SMILES_CANONICAL") and (
                desc_row[2].strip() == '"OpenEye OEToolkits"'
            ):
                canonical_smiles = desc_row[3].strip('"').strip(";")

            if (desc_row[1].strip() == "SMILES") and (
                desc_row[2].strip() == '"OpenEye OEToolkits"'
            ):
                isomeric_smiles = desc_row[3].replace('"', "")

            if desc_row[1].strip() == "InChIKey":
                inchikey = desc_row[3]
        if any((i is None for i in (canonical_smiles, isomeric_smiles, inchikey))):
            continue
        rows.append(
            (
                binder_id,
                chemical_name,
                molecular_weight,
                canonical_smiles,
                isomeric_smiles,
                inchikey,
            )
        )
    return pd.DataFrame(
        rows,
        columns=[
            "binder_id",
            "chemical_name",
            "molecular_weight",
            "canonical_smiles",
            "isomeric_smiles",
            "inchikey",
        ],
    )
