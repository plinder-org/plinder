# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from json import load
from pathlib import Path
from typing import Any, Optional
from zipfile import ZipFile

import pandas as pd
from omegaconf import DictConfig

from plinder.core.utils import cpl, gcs
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
    local = Path(f"{cfg.data.plinder_dir}/{suffix}")
    if not local.exists() or cfg.data.force_update:
        remote = f"{cfg.data.plinder_remote}/{suffix}"
        LOG.info(f"downloading {remote}")
        df = pd.read_parquet(remote)
        local.parent.mkdir(exist_ok=True, parents=True)
        df.to_parquet(local, index=False)
        return df
    LOG.info(f"reading {local}")
    _PLINDEX = pd.read_parquet(local)
    return _PLINDEX


def get_manifest(
    *,
    cfg: Optional[DictConfig] = None,
    plindex: Optional[pd.DataFrame] = None,
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
    suffix = f"{cfg.data.index}/{cfg.data.manifest_file}"
    manifest = Path(f"{cfg.data.plinder_dir}/{suffix}")
    if not manifest.exists() or cfg.data.force_update:
        if plindex is None:
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
def _load_entries_from_zips(
    *,
    cfg: DictConfig,
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

    per_zip: dict[str, list[str]] | None = None
    entry_msg = "all"
    if pdb_ids is not None:
        zip_paths_set = set()
        per_zip = {}
        for pdb_id in pdb_ids:
            code = pdb_id[-3:-1]
            zip_paths_set.add(entry_dir / f"{code}.zip")
            per_zip.setdefault(code, [])
            per_zip[code].append(f"{pdb_id}.json")
        entry_msg = str(sum((len(pz) for pz in per_zip.values())))
        zip_paths = list(zip_paths_set)
    elif two_char_codes is not None:
        zip_paths = [entry_dir / f"{code}.zip" for code in two_char_codes]
    else:
        # start with remote paths in case we have not downloaded them yet
        zip_paths = [
            entry_dir / Path(blob).name
            for blob in gcs.list_dir(
                gcs_path=(Path(cfg.data.plinder_remote) / cfg.data.entries).as_posix(),
                cfg=cfg,
            )
        ]
    gcs.download_if_not_exists(
        local_paths=zip_paths,
        remote_root=Path(cfg.data.plinder_remote) / cfg.data.entries,
        cfg=cfg,
    )
    reduced: dict[str, Any] = {}
    LOG.info(f"attempting to load {entry_msg} entries from {len(zip_paths)} zips")
    for zip_path in zip_paths:
        with ZipFile(zip_path) as archive:
            names = archive.namelist()
            if per_zip is not None:
                names = per_zip[zip_path.stem]
            for name in names:
                with archive.open(name) as obj:
                    pdb_id = name.replace(".json", "")
                    # TODO: port Entry to plinder-core for better validation
                    if prune:
                        reduced[pdb_id] = _prune_entry(load(obj))
                    else:
                        reduced[pdb_id] = load(obj)
    LOG.info(f"loaded {len(reduced)} entries")
    return reduced


def load_entries(
    *,
    pdb_ids: str | list[str] | None = None,
    two_char_codes: str | list[str] | None = None,
    cfg: Optional[DictConfig] = None,
    prune: bool = True,
) -> dict[str, Any]:
    if isinstance(pdb_ids, str):
        pdb_ids = [pdb_ids]
    if isinstance(two_char_codes, str):
        two_char_codes = [two_char_codes]
    result: dict[str, Any] = _load_entries_from_zips(
        cfg=cfg,
        two_char_codes=two_char_codes,
        pdb_ids=pdb_ids,
        prune=prune,
    )
    return result


def download_plinder_cmd() -> None:
    """
    Download the full plinder dataset for the current configuration.
    Note that even though this is wrapped in a progress bar, the estimated
    completion time can vary wildly as it iterates over larger files vs.
    smaller ones.
    """
    cfg = get_config()
    LOG.info(f"downloading {cfg.data.plinder_remote} -> {cfg.data.plinder_dir}")
    LOG.info("if this is the first time you are running this command, it will take a while!")
    LOG.info("the estimated time on the progress bar will vary wildly based on file size")
    LOG.info("if you need to cancel this and come back to it, it will pick up where it left off")
    cpl.download_many(rel="")
