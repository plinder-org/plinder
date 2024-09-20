# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from argparse import ArgumentParser
from json import load
from pathlib import Path
from shutil import rmtree
from textwrap import dedent
from time import time
from typing import Any, Optional
from zipfile import ZipFile

import pandas as pd
from omegaconf import DictConfig

from plinder.core.utils import cpl, unpack
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import get_zips_to_unpack

LOG = setup_logger(__name__)

_PLINDEX = None
_MANIFEST = None


@timeit
def get_plindex() -> pd.DataFrame:
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
    from plinder.core.scores import query_index

    global _PLINDEX

    if _PLINDEX is not None:
        return _PLINDEX
    _PLINDEX = query_index(columns=["*"], splits=["*"])
    return _PLINDEX


def get_manifest() -> pd.DataFrame:
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
    from plinder.core.scores import query_index

    global _MANIFEST

    if _MANIFEST is not None:
        return _MANIFEST
    _MANIFEST = query_index(columns=["system_id", "entry_pdb_id"], splits=["*"])
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


def _remove_old_linked_structures(data_dir: Path) -> None:
    zips = (data_dir / "linked_structures").glob("*zip")
    for zip in zips:
        if zip.stem in ["apo", "pred"]:
            continue
        zip.unlink()
        done = zip.parent / (zip.stem + "_done")
        if done.is_file():
            done.unlink()
    if (data_dir / "linked_structures" / "apo").is_dir():
        LOG.info("found old apo linked structures, removing")
        rmtree(data_dir / "linked_structures" / "apo")
    if (data_dir / "linked_structures" / "pred").is_dir():
        LOG.info("found old pred linked structures, removing")
        rmtree(data_dir / "linked_structures" / "pred")


def download_plinder_cmd(args: list[str] | None = None) -> None:
    """
    Download the full plinder dataset for the current configuration.
    Note that even though this is wrapped in a progress bar, the estimated
    completion time can vary wildly as it iterates over larger files vs.
    smaller ones.
    """
    t0 = time()
    parser = ArgumentParser(usage=download_plinder_cmd.__doc__)
    parser.add_argument("--release", default=None, help="plinder release")
    parser.add_argument("--iteration", default=None, help="plinder iteration")
    parser.add_argument("-y", "--yes", action="store_true", help="skip confirmation")
    ns, args = parser.parse_known_args(args=args)
    autodo = ns.yes
    if len(args):
        LOG.warning(f"ignoring arguments {args}")
    kwargs = None
    if ns.release is not None:
        kwargs = dict(data=dict(plinder_release=ns.release))
    if ns.iteration is not None:
        if kwargs is None:
            kwargs = dict(data=dict(plinder_iteration=ns.iteration))
        else:
            kwargs["data"]["plinder_iteration"] = ns.iteration
    cfg = get_config(config=kwargs)
    LOG.info(
        dedent(
            f"""
            Syncing {cfg.data.plinder_remote} -> {cfg.data.plinder_dir}.
            If this is the first time you are running this command, it will take a while!

            The estimated time on the progress bar may vary wildly based on varied file sizes.
            If you need to cancel this and come back to it, it will pick up where it left off.
            """
        )
    )
    LOG.debug("cleaning up old linked structures")
    _remove_old_linked_structures(Path(cfg.data.plinder_dir))
    for attr in cfg.data:
        if (
            attr.startswith("plinder_")
            or attr.endswith("_file")
            or attr in ["ingest", "validation", "force_update"]
        ):
            continue
        path = None
        if attr == "scores":
            do = (
                input("Download the full scores dataset? [Y/n] ").lower()
                in ["", "y", "yes"]
                if not autodo
                else True
            )
            if do:
                for subdb in ["apo", "pred", "holo"]:
                    msg = f"Syncing {getattr(cfg.data, attr)}/search_db={subdb}"
                    if subdb == "holo":
                        msg += ", this may take a while!"
                    LOG.info(msg)
                    if subdb == "holo":
                        LOG.info(
                            "Note that the tqdm progress bar for holo is not very useful, please be patient!"
                        )
                    cpl.get_plinder_path(
                        rel=f"{getattr(cfg.data, attr)}/search_db={subdb}",
                        force_progress=True,
                    )
            else:
                LOG.info(
                    "skipping scores download, plinder.core.scores will download it lazily on request!"
                )
        else:
            msg = f"Syncing {getattr(cfg.data, attr)}"
            do = True
            if attr in ["linked_structures", "systems"]:
                if not autodo:
                    do = input(f"Download the {attr} dataset? [Y/n] ").lower() in [
                        "",
                        "y",
                        "yes",
                    ]
                else:
                    do = True
                msg += ", this may take a while!"
            if do:
                LOG.info(msg)
                path = cpl.get_plinder_path(
                    rel=getattr(cfg.data, attr),
                    force_progress=True,
                )
            else:
                LOG.info(
                    f"skipping {attr} download, plinder.core.PlinderSystem will download lazily as needed on request!"
                )
        if path is not None and attr in ["linked_structures", "systems"]:
            LOG.info(
                f"extracting {getattr(cfg.data, attr)} archives, you may want to stretch your legs."
            )
            codes: list[str] | None = [p.stem for p in path.glob("*zip")]
            if attr == "linked_structures":
                codes = None
            get_zips_to_unpack(kind=attr, two_char_codes=codes)

    t1 = time()
    total = t1 - t0
    timing = f"{total:.2f}s"
    if total > 60:
        timing = f"{total / 60:.2f}m"
    elif total > 3600:
        timing = f"{total / 3600:.2f}h"

    LOG.info(
        dedent(
            f"""
            Sync complete in {timing}!

            If you downloaded all of the data, you can run:

                export PLINDER_OFFLINE=true

            This will avoid checking that files are still in sync when using plinder.core.
            If you didn't download all of the data, plinder.core will download it lazily when
            it's needed. By default, plinder.core will check that files are still in sync
            in case any of the files for an existing release need to be patched.
            """
        )
    )
