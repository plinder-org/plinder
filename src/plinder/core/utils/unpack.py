# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

from omegaconf import DictConfig

from plinder.core.utils import gcs
from plinder.core.utils.config import get_config

ZIP_KINDS = Literal["entries", "systems"]
ID_KINDS = Literal["system_id", "pdb_id", "two_char_code"]


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
        return "system_id", systems
    if pdbs:
        return "pdb_id", pdbs
    return "two_char_code", two_chars


def get_zips_to_unpack(
    *,
    kind: ID_KINDS,
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

    root = Path(conf.data.plinder_dir) / getattr(conf.data, kind)
    zips: dict[Path, list[str]] = {}
    if id_kind == "system_id":
        for system_id in ids:
            two_char_code = system_id[1:3]
            key = root / f"{two_char_code}.zip"
            zips.setdefault(key, [])
            zips[key].append(system_id)
    elif id_kind == "pdb_id":
        for pdb_id in ids:
            code = pdb_id[-3:-1]
            key = root / f"{code}.zip"
            zips.setdefault(key, [])
            zips[key].append(pdb_id)
    elif id_kind == "two_char_code":
        for two_char_code in ids:
            key = root / f"{two_char_code}.zip"
            zips.setdefault(key, [])

    gcs.download_if_not_exists(
        local_paths=list(zips.keys()),
        remote_root=Path(conf.data.plinder_remote) / getattr(conf.data, kind),
    )
    return zips


def unpack_zips(
    *,
    zips: dict[Path, list[str]],
    kind: str,
) -> None:
    pass


#     entry_dir = Path(conf.data.plinder_dir) / conf.data.entries
#     entry_dir.mkdir(exist_ok=True, parents=True)
#
#     per_zip: dict[str, list[str]] | None = None
#     entry_msg = "all"
#     if pdb_ids is not None:
#         zip_paths_set = set()
#         per_zip = {}
#         for pdb_id in pdb_ids:
#             code = pdb_id[-3:-1]
#             zip_paths_set.add(entry_dir / f"{code}.zip")
#             per_zip.setdefault(code, [])
#             per_zip[code].append(f"{pdb_id}.json")
#         entry_msg = str(sum((len(pz) for pz in per_zip.values())))
#         zip_paths = list(zip_paths_set)
#     elif two_char_codes is not None:
#         zip_paths = [entry_dir / f"{code}.zip" for code in two_char_codes]
#     else:
#         # start with remote paths in case we have not downloaded them yet
#         zip_paths = [
#             entry_dir / Path(blob).name
#             for blob in gcs.list_dir(
#                 gcs_path=(Path(cfg.data.plinder_remote) / cfg.data.entries).as_posix(),
#                 cfg=cfg,
#             )
#         ]
