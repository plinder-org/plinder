# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
from pathlib import Path
from typing import Any, Optional
from zipfile import ZipFile

import pandas as pd
from omegaconf import DictConfig

from plinder.core.index import utils
from plinder.core.scores.links import query_links
from plinder.core.utils.unpack import get_zips_to_unpack
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


class PlinderSystem:
    """
    Core class for interacting with a single system and its assets.
    Relies on the entry data to populate the system which can be looked
    up automatically using the from_system_id class method.
    """

    def __init__(
        self,
        system_id: str,
        entry: dict[str, Any],
    ) -> None:
        self.system_id = system_id
        entry_pdb_id = system_id.split("__")[0]
        self._entry = entry
        if system_id not in entry[entry_pdb_id]["systems"]:
            raise ValueError(f"system_id={system_id} not in entry")
        self.system = entry[entry_pdb_id]["systems"][system_id]
        self._archive = None

    @property
    def archive(self) -> Path:
        if self._archive is None:
            zips = get_zips_to_unpack(kind="systems", system_ids=[self.system_id])
            [self._archive] = list(zips.keys())
        return self._archive

    @property
    def system_cif(self) -> str:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/system.cif") as f:
                return f.read().decode("utf8")

    @property
    def receptor_cif(self) -> str:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/receptor.cif") as f:
                return f.read().decode("utf8")

    @property
    def receptor_pdb(self) -> str:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/receptor.pdb") as f:
                return f.read().decode("utf8")

    @property
    def system_pdb(self) -> str:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/system.pdb") as f:
                return f.read().decode("utf8")

    @property
    def sequences(self) -> str:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/sequences.fasta") as f:
                return f.read().decode("utf8")

    @property
    def chain_mapping(self) -> dict[str, Any]:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/chain_mapping.json") as f:
                return json.load(f)

    @property
    def water_mapping(self) -> dict[str, Any]:
        with ZipFile(self.archive) as arch:
            with arch.open(f"{self.system_id}/water_mapping.json") as f:
                return json.load(f)

    @property
    def ligands(self) -> dict[str, str]:
        ligands = {}
        with ZipFile(self.archive) as arch:
            for name in arch.namelist():
                if name.startswith(f"{self.system_id}/ligand_files"):
                    ligands[Path(name).stem] = arch.read(name).decode("utf8")
        return ligands

    @property
    def structures(self) -> list[str]:
        structures = []
        with ZipFile(self.archive) as arch:
                for name in arch.namelist():
                    if name.startswith(f"{self.system_id}"):
                        structures.append(name)
        return structures

    @property
    def linked_systems(self) -> pd.DataFrame:
        filters=[("query_system", "==", self.system_id)]
        print(filters)
        df = query_links(
            filters=[("query_system", "==", self.system_id)],
        )
        return df

    @classmethod
    def from_system_id(
        cls,
        system_id: str,
        cfg: Optional[DictConfig] = None,
        prune: bool = True,
    ) -> "PlinderSystem":
        entry_pdb_id = system_id.split("__")[0]
        try:
            entry = utils.load_entries(pdb_ids=[entry_pdb_id], cfg=cfg, prune=prune)
            system = entry[entry_pdb_id]["systems"].get(system_id)
        except KeyError:
            raise ValueError(f"system_id={system_id} not in entry={entry_pdb_id}")
        if system is None:
            raise ValueError(f"system_id={system_id} not in entry={entry_pdb_id}")
        return cls(
            system_id=system_id,
            entry=entry,
        )
