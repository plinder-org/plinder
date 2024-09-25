# type: ignore
# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from typing import Any, Optional, Union

from omegaconf import DictConfig

from plinder.core.index.system import PlinderSystem
from plinder.core.index.utils import get_manifest
from plinder.core.utils.unpack import expand_config_context


# TODO: not used currently
def load_systems(
    *,
    system_ids: Optional[Union[str, list[str]]] = None,
    pdb_ids: Optional[Union[str, list[str]]] = None,
    two_char_codes: Optional[Union[str, list[str]]] = None,
    cfg: Optional[DictConfig] = None,
) -> dict[str, Any]:
    kind, items = expand_config_context(
        system_ids=system_ids, pdb_ids=pdb_ids, two_char_codes=two_char_codes, cfg=cfg
    )
    if kind == "system_ids":
        return {system_id: PlinderSystem(system_id=system_id) for system_id in items}
    manifest = get_manifest()
    systems = {}
    if kind == "pdb_ids":
        for pdb_id in items:
            ids = manifest[manifest["pdb_id"] == pdb_id]["system_id"].to_list()
            for system_id in ids:
                systems[system_id] = PlinderSystem(system_id=system_id)
        return systems
    manifest["two_char_code"] = manifest["entry_pdb_id"].str[1:3]
    for two_char_code in items:
        ids = manifest[manifest["two_char_code"] == two_char_code][
            "system_id"
        ].to_list()
        for system_id in ids:
            systems[system_id] = PlinderSystem(system_id=system_id)
    return systems
