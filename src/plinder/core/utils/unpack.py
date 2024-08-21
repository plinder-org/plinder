# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

from omegaconf import DictConfig

from plinder.core.utils import cpl
from plinder.core.utils.config import get_config

ZIP_KINDS = Literal["entries", "linked_structures", "systems"]
ID_KINDS = Literal["system_ids", "pdb_ids", "two_char_code"]


def expand_config_context(
    *,
    system_ids: Optional[str | list[str]] = None,
    pdb_ids: Optional[str | list[str]] = None,
    two_char_codes: Optional[str | list[str]] = None,
    cfg: Optional[DictConfig] = None,
) -> tuple[str, list[str]]:
    conf = cfg or get_config()
    if system_ids is None:
        systems = conf.context.system_ids
    elif isinstance(system_ids, str):
        systems = [system_ids]
    else:
        systems = system_ids
    if pdb_ids is None:
        pdbs = conf.context.pdb_ids
    elif isinstance(pdb_ids, str):
        pdbs = [pdb_ids]
    else:
        pdbs = pdb_ids
    if two_char_codes is None:
        two_chars = conf.context.two_char_codes
    elif isinstance(two_char_codes, str):
        two_chars = [two_char_codes]
    else:
        two_chars = two_char_codes
    if systems:
        return "system_ids", systems
    if pdbs:
        return "pdb_ids", pdbs
    return "two_char_codes", two_chars


def get_zips_to_unpack(
    *,
    kind: ZIP_KINDS,
    system_ids: Optional[str | list[str]] = None,
    pdb_ids: Optional[str | list[str]] = None,
    two_char_codes: Optional[str | list[str]] = None,
    cfg: Optional[DictConfig] = None,
) -> dict[Path, list[str]]:
    conf = cfg or get_config()
    id_kind, ids = expand_config_context(
        system_ids=system_ids,
        pdb_ids=pdb_ids,
        two_char_codes=two_char_codes,
        cfg=cfg,
    )

    root = cpl.get_plinder_path(rel=getattr(conf.data, kind), download=False)
    zips: dict[Path, list[str]] = {}
    if id_kind == "system_ids":
        for system_id in ids:
            two_char_code = system_id[1:3]
            key = root / f"{two_char_code}.zip"
            zips.setdefault(key, [])
            zips[key].append(system_id)
    elif id_kind == "pdb_ids":
        for pdb_id in ids:
            code = pdb_id[-3:-1]
            key = root / f"{code}.zip"
            zips.setdefault(key, [])
            zips[key].append(pdb_id)
    elif id_kind == "two_char_codes":
        for two_char_code in ids:
            key = root / f"{two_char_code}.zip"
            zips.setdefault(key, [])

    paths = list(zips.keys())
    cpl.download_paths(paths=cpl.get_plinder_paths(paths=paths))
    return zips
