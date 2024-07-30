# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations
import json
from pathlib import Path
from typing import Any
from zipfile import ZipFile

import pandas as pd

from plinder.core.index import utils
from plinder.core.scores.links import query_links
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import get_zips_to_unpack

LOG = setup_logger(__name__)


class PlinderSystem:
    """
    Core class for interacting with a single system and its assets.
    Relies on the entry data to populate the system which can be looked
    up automatically using the from_system_id class method.
    """

    def __repr__(self) -> str:
        return f"PlinderSystem(system_id={self.system_id})"

    def __init__(
        self,
        *,
        system_id: str,
        prune: bool = True,
    ) -> None:
        self.system_id = system_id
        self.prune = prune
        self._system = None
        self._archive = None
        self._chain_mapping = None
        self._water_mapping = None

    @property
    def system(self) -> dict[str, Any] | None:
        if self._system is None:
            entry_pdb_id = self.system_id.split("__")[0]
            try:
                entry = utils.load_entries(pdb_ids=[entry_pdb_id], prune=self.prune)
                self._system = entry[entry_pdb_id]["systems"][self.system_id]
            except KeyError:
                raise ValueError(
                    f"system_id={self.system_id} not found in entry={entry_pdb_id}"
                )
        return self._system

    @property
    def archive(self) -> Path:
        if self._archive is None:
            zips = get_zips_to_unpack(kind="systems", system_ids=[self.system_id])
            [archive] = list(zips.keys())
            self._archive = Path(archive.fspath).parent / archive.stem
            if not self._archive.is_dir():
                with ZipFile(archive) as arch:
                    arch.extractall(path=self._archive)
        return self._archive

    @property
    def system_cif(self) -> str:
        return (self.archive / f"{self.system_id}/system.cif").as_posix()

    @property
    def receptor_cif(self) -> str:
        return (self.archive / f"{self.system_id}/receptor.cif").as_posix()

    @property
    def receptor_pdb(self) -> str:
        return (self.archive / f"{self.system_id}/receptor.pdb").as_posix()

    @property
    def system_pdb(self) -> str:
        return (self.archive / f"{self.system_id}/system.pdb").as_posix()

    @property
    def sequences(self) -> str:
        return (self.archive / f"{self.system_id}/sequences.fasta").as_posix()

    @property
    def chain_mapping(self) -> dict[str, Any]:
        if self._chain_mapping is None:
            with (self.archive / f"{self.system_id}/chain_mapping.json").open() as f:
                self._chain_mapping = json.load(f)
        return self._chain_mapping

    @property
    def water_mapping(self) -> dict[str, Any]:
        if self._water_mapping is None:
            with (self.archive / f"{self.system_id}/water_mapping.json").open() as f:
                self._water_mapping = json.load(f)
        return self._water_mapping

    @property
    def ligands(self) -> dict[str, str]:
        ligands = {}
        for ligand in (self.archive / f"{self.system_id}/ligand_files/").glob("*.sdf"):
            ligands[ligand.stem] = ligand.as_posix()
        return ligands

    @property
    def structures(self) -> list[str]:
        return [
            path.as_posix() for path in (self.archive / f"{self.system_id}/").rglob("*")
        ]

    @property
    def linked_systems(self) -> pd.DataFrame:
        return query_links(filters=[("query_system", "==", self.system_id)])
