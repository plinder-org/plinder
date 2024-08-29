from typing import Optional

import pandas as pd
from omegaconf import DictConfig

from plinder.core.index.utils import get_plindex
from plinder.core.scores.clusters import query_clusters
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

_SPLIT = None


@timeit
def get_split(
    *,
    cfg: Optional[DictConfig] = None,
) -> pd.DataFrame:
    """
    Fetch the plinder split and cache it

    Parameters
    ----------
    cfg : DictConfig, default=None
        the plinder-core config

    Returns
    -------
    pd.DataFrame
        the plinder split
    """
    global _SPLIT
    if _SPLIT is not None:
        return _SPLIT
    cfg = cfg or get_config()
    suffix = f"{cfg.data.splits}/{cfg.data.split_file}"
    split = cpl.get_plinder_path(rel=suffix)
    LOG.info(f"reading {split}")
    _SPLIT = pd.read_parquet(split)
    return _SPLIT


# TODO: this can be removed after cleanup merge
def reset_lipinski_and_other(df: pd.DataFrame) -> pd.DataFrame:
    # Reset artifacts/ions
    mask_reset = df["ligand_is_ion"] | df["ligand_is_artifact"]
    df.loc[
        mask_reset,
        [
            "ligand_is_fragment",
            "ligand_is_lipinski",
            "ligand_is_cofactor",
            "ligand_is_covalent",
            "ligand_is_oligo",
        ],
    ] = False

    # Lipinski like Ro3
    mask_ro3 = (
        (df["ligand_molecular_weight"] < 300)
        & (df["ligand_crippen_clogp"] < 3)
        & (df["ligand_num_hbd"] <= 3)
        & (df["ligand_num_hba"] <= 3)
        & ~mask_reset
    )
    df.loc[mask_ro3, ["ligand_is_fragment", "ligand_is_lipinski"]] = True

    # Lipinski like Ro5
    mask_ro5 = (
        (df["ligand_molecular_weight"] < 500)
        & (df["ligand_crippen_clogp"] < 5)
        & (df["ligand_num_hbd"] <= 5)
        & (df["ligand_num_hba"] <= 10)
        & ~mask_reset
        & ~mask_ro3
    )
    df.loc[mask_ro5, "ligand_is_lipinski"] = True

    # ligand_is_other
    df["ligand_is_other"] = ~(
        df["ligand_is_invalid"]
        | df["ligand_is_ion"]
        | df["ligand_is_oligo"]
        | df["ligand_is_artifact"]
        | df["ligand_is_cofactor"]
        | df["ligand_is_lipinski"]
        | df["ligand_is_fragment"]
        | df["ligand_is_covalent"]
    )

    return df


def get_extended_plindex(
    *,
    cfg: Optional[DictConfig] = None,
    plindex: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    if plindex is None:
        plindex = get_plindex(cfg=cfg)
    plindex = reset_lipinski_and_other(plindex)
    plindex["system_ligand_max_qed"] = plindex.groupby("system_id")[
        "ligand_qed"
    ].transform("max")
    plindex["system_pass_validation_criteria"] = plindex[
        "system_pass_validation_criteria"
    ].fillna(False)
    for n in [
        "lipinski",
        "cofactor",
        "fragment",
        "oligo",
        "artifact",
        "other",
        "covalent",
        "invalid",
        "ion",
    ]:
        plindex[f"system_ligand_has_{n}"] = plindex.groupby("system_id")[
            f"ligand_is_{n}"
        ].transform("any")
    plindex["system_interacting_protein_chains_total_length"] = plindex[
        "system_interacting_protein_chains_length"
    ].apply(lambda x: sum(int(i) for i in x.split(";")))
    plindex["ligand_is_proper"] = (
        ~plindex["ligand_is_ion"] & ~plindex["ligand_is_artifact"]
    )
    plindex["system_proper_num_ligand_chains"] = plindex.groupby("system_id")[
        "ligand_is_proper"
    ].transform("sum")
    plindex["system_proper_num_interactions"] = (
        plindex["ligand_num_interactions"]
        .where(plindex["ligand_is_proper"], other=0)
        .groupby(plindex["system_id"])
        .transform("sum")
    )
    plindex["system_proper_pocket_num_residues"] = (
        plindex["ligand_num_neighboring_residues"]
        .where(plindex["ligand_is_proper"], other=0)
        .groupby(plindex["system_id"])
        .transform("sum")
    )
    plindex["system_proper_ligand_max_molecular_weight"] = (
        plindex["ligand_molecular_weight"]
        .where(plindex["ligand_is_proper"], other=float("-inf"))
        .groupby(plindex["system_id"])
        .transform("max")
    )
    plindex["uniqueness"] = (
        plindex["system_id_no_biounit"]
        + "_"
        + plindex["pli_qcov__100__strong__component"]
    )
    plindex["system_num_ligands_in_biounit"] = plindex.groupby(
        ["entry_pdb_id", "system_biounit_id"]
    )["system_id"].transform("count")
    plindex["system_num_unique_ligands_in_biounit"] = plindex.groupby(
        ["entry_pdb_id", "system_biounit_id"]
    )["ligand_unique_ccd_code"].transform("nunique")
    plindex["system_proper_num_ligands_in_biounit"] = plindex.groupby(
        ["entry_pdb_id", "system_biounit_id"]
    )["ligand_is_proper"].transform("sum")
    ccd_dict = (
        plindex[plindex["ligand_is_proper"]]
        .groupby("system_id")["ligand_unique_ccd_code"]
        .agg(lambda x: "-".join(sorted(set(x))))
        .to_dict()
    )
    plindex["system_proper_ligand_unique_ccd_codes"] = plindex["system_id"].map(
        ccd_dict
    )
    clusters = query_clusters(
        columns=["system_id", "label", "threshold", "metric", "cluster"],
        filters=[("directed", "==", "False")],
    )
    if clusters is None:
        LOG.error("No clusters found")
        return pd.DataFrame()
    clusters = clusters.pivot_table(
        values="label",
        index="system_id",
        columns=["metric", "cluster", "threshold"],
        aggfunc="first",
    )
    clusters.columns = [
        f"{metric}__{threshold.replace('.parquet', '')}__{cluster}"
        for metric, cluster, threshold in clusters.columns
    ]
    clusters.reset_index(inplace=True)
    plindex = pd.merge(plindex, clusters, on="system_id", how="left")
    return plindex
