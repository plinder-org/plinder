# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from json import load
from pathlib import Path
from typing import Any, Optional
from zipfile import ZipFile

import pandas as pd
from omegaconf import DictConfig

from plinder.core.utils import cpl, unpack
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

_PLINDEX = None
_MANIFEST = None


@timeit
def get_plindex(
    *,
    cfg: Optional[DictConfig] = None,
) -> pd.DataFrame:
    """
    Fetch the plindex and cache it

    Parameters
    ----------
    cfg : DictConfig, default=None
        the plinder-core config

    Returns
    -------
    pd.DataFrame
        the plindex
    """
    global _PLINDEX

    if _PLINDEX is not None:
        return _PLINDEX
    cfg = cfg or get_config()
    suffix = f"{cfg.data.index}/{cfg.data.index_file}"
    index = cpl.get_plinder_path(rel=suffix)
    LOG.info(f"reading {index}")
    _PLINDEX = pd.read_parquet(index)
    return _PLINDEX


def get_manifest(
    *,
    cfg: Optional[DictConfig] = None,
) -> pd.DataFrame:
    """
    Fetch the manifest and cache it

    Parameters
    ----------
    cfg : DictConfig, default=None
        the plinder-core config
    plindex : pd.DataFrame, default=None
        the plindex

    Returns
    -------
    pd.DataFrame
        the manifest
    """
    global _MANIFEST

    if _MANIFEST is not None:
        return _MANIFEST
    cfg = cfg or get_config()
    suffix = f"{cfg.data.manifest}/{cfg.data.manifest_file}"
    manifest = Path(f"{cfg.data.plinder_dir}/{suffix}")
    if not manifest.exists() or cfg.data.force_update:
        manifest.parent.mkdir(exist_ok=True, parents=True)
        plindex = get_plindex(cfg=cfg)
        plindex[["system_id", "entry_pdb_id"]].to_parquet(manifest, index=False)
    _MANIFEST = pd.read_parquet(manifest)
    return _MANIFEST


def _prune_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """
    Prune the entry as in Entry.prune

    Parameters
    ----------
    entry : dict[str, Any]
        the entry

    Returns
    -------
    dict[str, Any]
        the pruned entry
    """
    entry["systems"] = {
        id_: s
        for id_, s in entry["systems"].items()
        if any(not l["is_ion"] and not l["is_artifact"] for l in s["ligands"])
        and len(
            set(
                chain
                for ligand in s["ligands"]
                for chain in sorted(ligand["interacting_residues"].keys())
            )
        )
        <= 5
        and len(s["ligands"]) <= 5
    }
    return entry


@timeit
def load_entries(
    *,
    cfg: Optional[DictConfig] = None,
    two_char_codes: list[str] | None = None,
    pdb_ids: list[str] | None = None,
    prune: bool = True,
) -> dict[str, Any]:
    """
    Load the entries from a list of pdb IDs or two character codes.
    If no filters are provided, all entries are loaded.

    Parameters
    ----------
    cfg : DictConfig
        the plinder-core config
    two_char_codes : list[str] | None, default=None
        only consider particular two character codes
    pdb_ids : list[str] | None, default=None
        only consider particular pdb IDs

    Returns
    -------
    dict[str, Any]
        the entries
    """
    cfg = cfg or get_config()
    entry_dir = Path(cfg.data.plinder_dir) / cfg.data.entries
    entry_dir.mkdir(exist_ok=True, parents=True)

    zips = unpack.get_zips_to_unpack(
        kind=cfg.data.entries,
        cfg=cfg,
        two_char_codes=two_char_codes,
        pdb_ids=pdb_ids,
    )
    reduced: dict[str, Any] = {}
    LOG.info(f"loading entries from {len(zips)} zips")
    for zip_path, pdb_ids in zips.items():
        with ZipFile(zip_path) as archive:
            if len(pdb_ids):
                names = [f"{pdb_id}.json" for pdb_id in pdb_ids]
            else:
                names = archive.namelist()
            for name in names:
                with archive.open(name) as obj:
                    pdb_id = name.replace(".json", "")
                    # TODO: port Entry to plinder.core for model validation
                    if prune:
                        reduced[pdb_id] = _prune_entry(load(obj))
                    else:
                        reduced[pdb_id] = load(obj)
    LOG.info(f"loaded {len(reduced)} entries")
    return reduced


def download_plinder_cmd() -> None:
    """
    Download the full plinder dataset for the current configuration.
    Note that even though this is wrapped in a progress bar, the estimated
    completion time can vary wildly as it iterates over larger files vs.
    smaller ones.
    """
    cfg = get_config()
    LOG.info(
        f"""
downloading {cfg.data.plinder_remote} -> {cfg.data.plinder_dir}
if this is the first time you are running this command, it will take a while!
the estimated time on the progress bar may vary wildly based on file size
if you need to cancel this and come back to it, it will pick up where it left off"""
    )
    cpl.get_plinder_path()
