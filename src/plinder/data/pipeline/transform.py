# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def transform_bindingdb_affinity_data(*, raw_affinity_path: Path) -> pd.DataFrame:
    """Parse BindingDB TSV into a per-(PDB, ligand) affinity table.

    Each row maps a ``pdbid_ligid`` key (e.g. ``"1ABC_ATP"``) to a
    median pChEMBL value derived from Ki/Kd measurements, plus the
    BindingDB target sequence for downstream validation.

    The BindingDB field ``PDB ID(s) for Ligand-Target Complex`` lists
    PDB IDs matched at 85% sequence identity, which can assign affinity
    values to the wrong complex (see `#94`_).  To guard against this,
    the target sequence is preserved so that callers can verify it
    against the actual PDB chain sequence before accepting the value.

    .. _#94: https://github.com/plinder-org/plinder/issues/94

    Parameters
    ----------
    raw_affinity_path : Path
        Path to the BindingDB TSV file.

    Returns
    -------
    pd.DataFrame
        Columns: ``pdbid_ligid``, ``pchembl``, ``target_sequence``.
    """

    def calc_pchembl(affinity: float) -> Any:
        # pchembl = -log10(affinity_M); naming follows the ChEMBL convention
        # for the log-transformed value, NOT the ChEMBL database itself.
        affinity = affinity * 10**-9
        if affinity > 0:
            return -1.0 * np.log10(affinity)
        else:
            return np.nan

    # BindingDB renamed the target sequence column across releases:
    #   old (<=2024): "BindingDB Target Chain  Sequence" (double space)
    #   new (>=2025): "BindingDB Target Chain Sequence 1" (numbered)
    _SEQ_COL_NEW = "BindingDB Target Chain Sequence 1"
    _SEQ_COL_OLD = "BindingDB Target Chain  Sequence"
    header = set(pd.read_csv(raw_affinity_path, sep="\t", nrows=0).columns)
    if _SEQ_COL_NEW in header:
        seq_col = _SEQ_COL_NEW
    elif _SEQ_COL_OLD in header:
        seq_col = _SEQ_COL_OLD
    else:
        raise ValueError(
            "BindingDB TSV is missing target sequence column. "
            f"Expected '{_SEQ_COL_NEW}' or '{_SEQ_COL_OLD}'. "
            "Required for target sequence validation (#94)."
        )
    cols = [
        "Ligand HET ID in PDB",
        "PDB ID(s) for Ligand-Target Complex",
        "Ki (nM)",
        "Kd (nM)",
        "EC50 (nM)",
        seq_col,
    ]
    df = pd.read_csv(raw_affinity_path, sep="\t", usecols=cols, low_memory=False)

    df["pchembl"] = (
        df[["Ki (nM)", "Kd (nM)"]]
        .apply(set, axis=1)
        .apply(lambda x: [i for i in x if str(i) != "nan"])
    )
    df = df[df["pchembl"].apply(lambda x: x != [])]
    df["pchembl"] = df["pchembl"].apply(
        lambda x: calc_pchembl(float(str(x[0]).replace(">", "").replace("<", "")))
    )

    df.rename(columns={seq_col: "target_sequence"}, inplace=True)
    df = df[
        [
            "PDB ID(s) for Ligand-Target Complex",
            "Ligand HET ID in PDB",
            "target_sequence",
            "pchembl",
        ]
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
    df = df[["pdbid_ligid", "pchembl", "target_sequence"]].drop_duplicates()

    # Per pdbid_ligid: take median pchembl, keep first non-null target sequence
    grouped = df.groupby("pdbid_ligid").agg(
        pchembl=("pchembl", "median"),
        target_sequence=("target_sequence", "first"),
    )
    return grouped.reset_index()


def transform_components_data(*, raw_components_path: Path) -> pd.DataFrame:
    import biotite.structure.io.pdbx as pdbx

    data = pdbx.CIFFile.read(str(raw_components_path))
    rows = []
    for block in data.values():
        if "chem_comp" not in block:
            continue
        chem_comp = block["chem_comp"]
        binder_id = chem_comp["id"].as_array()[0]
        chemical_name = chem_comp["name"].as_array()[0]
        molecular_weight = chem_comp["formula_weight"].as_array()[0]

        canonical_smiles, isomeric_smiles, inchikey = None, None, None
        if "pdbx_chem_comp_descriptor" in block:
            desc = block["pdbx_chem_comp_descriptor"]
            types = desc["type"].as_array()
            programs = desc["program"].as_array()
            descriptors = desc["descriptor"].as_array()
            for dtype, prog, val in zip(types, programs, descriptors):
                dtype_s = dtype.strip()
                prog_s = prog.strip().strip('"')
                if dtype_s == "SMILES_CANONICAL" and prog_s == "OpenEye OEToolkits":
                    canonical_smiles = val.strip('"').strip(";")
                if dtype_s == "SMILES" and prog_s == "OpenEye OEToolkits":
                    isomeric_smiles = val.replace('"', "")
                if dtype_s == "InChIKey":
                    inchikey = val
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
