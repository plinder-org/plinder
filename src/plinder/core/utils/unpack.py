# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from typing import Optional

from omegaconf import DictConfig

from plinder.core.utils.config import get_config


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
