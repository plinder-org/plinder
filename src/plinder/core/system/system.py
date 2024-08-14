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
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import get_zips_to_unpack

LOG = setup_logger(__name__)


class PlinderSystem:
    """
    Core class for interacting with a single system and its assets.
    Relies on the entry data to populate the system, system zips
    to obtain structure files, and the linked_structures archive
    for linked structures.
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
        self._entry = None
        self._system = None
        self._archive = None
        self._chain_mapping = None
        self._water_mapping = None
        self._linked_structures = None
        self._linked_archive = None

    @property
    def entry(self) -> dict[str, Any]:
        if self._entry is None:
            entry_pdb_id = self.system_id.split("__")[0]
            try:
                entry = utils.load_entries(pdb_ids=[entry_pdb_id], prune=self.prune)
                self._entry = entry[entry_pdb_id]
            except KeyError:
                raise ValueError(f"pdb_id={entry_pdb_id} not found in entries")
        return self._entry

    @property
    def system(self) -> dict[str, Any] | None:
        if self._system is None:
            try:
                self._system = self.entry["systems"][self.system_id]
            except KeyError:
                raise ValueError(f"system_id={self.system_id} not found in entry")
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
            path.as_posix()
            for path in (self.archive / f"{self.system_id}/").rglob("*")
            if path.is_file()
        ]

    @property
    def linked_structures(self) -> pd.DataFrame:
        if self._linked_structures is None:
            links = query_links(filters=[("query_system", "==", self.system_id)])
            self._linked_structures = links
        return self._linked_structures

    @property
    def linked_archive(self) -> Path:
        if self._linked_archive is None:
            cfg = get_config()
            structure_archive = cpl.get_plinder_path(rel=cfg.data.linked_structures)
            self._linked_archive = Path(structure_archive.fspath)
            if not (self._linked_archive / "apo").is_dir():
                with ZipFile(structure_archive) as arch:
                    arch.extractall(path=structure_archive)
        return self._linked_archive

    def get_linked_structure(self, link_kind: str, link_id: str, ext: str = "cif") -> str:
        if ext not in ["cif", "pdb"]:
            raise ValueError(f"extension must be cif or pdb, got {ext}")
        return (self.linked_archive / link_kind / self.system_id / link_id / "receptor.{ext}").as_posix()
