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


def get_extended_plindex(
    *,
    cfg: Optional[DictConfig] = None,
    plindex: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    The extended PLINDER index includes all of the cluster
    labels and additional columns that are used in the
    scoring and splitting.

    Parameters
    ----------
    cfg : DictConfig, default=None
        the plinder-core config
    plindex : pd.DataFrame, default=None
        the full PLINDER index

    Returns
    -------
    pd.DataFrame
        the extended PLINDER index
    """
    if plindex is None:
        plindex = get_plindex(cfg=cfg)
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
    plindex["system_protein_chains_total_length"] = plindex[
        "system_protein_chains_length"
    ].apply(lambda x: sum(x))

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
        f"{metric}__{threshold}__{cluster}"
        for metric, cluster, threshold in clusters.columns
    ]
    clusters.reset_index(inplace=True)
    plindex = pd.merge(plindex, clusters, on="system_id", how="left")
    return plindex
