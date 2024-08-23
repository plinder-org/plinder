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
        self._linked_archive: Path | None = None

    @property
    def entry(self) -> dict[str, Any] | None:
        """
        Store a reference to the entry JSON for this system

        Returns
        -------
        dict[str, Any] | None
            entry JSON
        """
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
        """
        Return the system metadata from the original entry JSON

        Returns
        -------
        dict[str, Any] | None
            system metadata
        """
        if self._system is None:
            try:
                assert self.entry is not None
                self._system = self.entry["systems"][self.system_id]
            except KeyError:
                raise ValueError(f"system_id={self.system_id} not found in entry")
        return self._system

    @property
    def archive(self) -> Path | None:
        """
        Return the path to the directory containing the plinder system

        Returns
        -------
        Path | None
            directory containing the plinder system
        """
        if self._archive is None:
            zips = get_zips_to_unpack(kind="systems", system_ids=[self.system_id])
            [archive] = list(zips.keys())
            self._archive = archive.parent / self.system_id  # type: ignore
            assert self._archive is not None
            if not self._archive.is_dir():
                with ZipFile(archive) as arch:
                    arch.extractall(path=archive.parent)
        return self._archive

    @property
    def system_cif(self) -> str:
        """
        Path to the system.cif file

        Returns
        -------
        str
            path
        """
        assert self.archive is not None
        return (self.archive / "system.cif").as_posix()

    @property
    def receptor_cif(self) -> str:
        """
        Path to the receptor.cif file

        Returns
        -------
        str
            path
        """
        assert self.archive is not None
        return (self.archive / "receptor.cif").as_posix()

    @property
    def receptor_pdb(self) -> str:
        """
        Path to the receptor.pdb file

        Returns
        -------
        str
            path
        """
        assert self.archive is not None
        return (self.archive / "receptor.pdb").as_posix()

    @property
    def sequences(self) -> str:
        """
        Path to the sequences.fasta file

        Returns
        -------
        str
            path
        """
        assert self.archive is not None
        return (self.archive / "sequences.fasta").as_posix()

    @property
    def chain_mapping(self) -> dict[str, Any] | None:
        """
        Chain mapping metadata

        Returns
        -------
        dict[str, Any] | None
            chain mapping
        """
        if self._chain_mapping is None:
            assert self.archive is not None
            with (self.archive / "chain_mapping.json").open() as f:
                self._chain_mapping = json.load(f)
        return self._chain_mapping

    @property
    def water_mapping(self) -> dict[str, Any] | None:
        """
        Water mapping metadata

        Returns
        -------
        dict[str, Any] | None
            water mapping
        """
        if self._water_mapping is None:
            assert self.archive is not None
            with (self.archive / "water_mapping.json").open() as f:
                self._water_mapping = json.load(f)
        return self._water_mapping

    @property
    def ligands(self) -> dict[str, str]:
        """
        Return a dictionary of ligand names to paths to ligand sdf files

        Returns
        -------
        dict[str, str]
            dictionary of ligand names to paths to ligand sdf files
        """
        assert self.archive is not None
        ligands = {}
        for ligand in (self.archive / "ligand_files/").glob("*.sdf"):
            ligands[ligand.stem] = ligand.as_posix()
        return ligands

    @property
    def structures(self) -> list[str]:
        """
        Return a list of paths to all structures in the plinder system

        Returns
        -------
        list[str]
            list of paths to structures
        """
        assert self.archive is not None
        return [path.as_posix() for path in self.archive.rglob("*") if path.is_file()]

    @property
    def linked_structures(self) -> pd.DataFrame | None:
        """
        Return a dataframe of linked structures for this system. Note
        that the dataframe will include all of the scores for the linked
        structures as well, so that particular alternatives can be chosen
        accordingly.

        Returns
        -------
        pd.DataFrame | None
            dataframe of linked structures if present in plinder
        """
        if self._linked_structures is None:
            links = query_links(filters=[("reference_system_id", "==", self.system_id)])
            self._linked_structures = links
            if (
                self._linked_structures is not None
                and not self._linked_structures.empty
            ):
                zips = get_zips_to_unpack(
                    kind="linked_structures", system_ids=[self.system_id]
                )
                [archive] = list(zips.keys())
                self._linked_archive = archive.parent
                if not (self._linked_archive / "apo" / self.system_id).is_dir() or not (
                    self._linked_archive / "pred" / self.system_id
                ).is_dir():
                    with ZipFile(archive) as arch:
                        arch.extractall(path=archive.parent)
        return self._linked_structures

    @property
    def linked_archive(self) -> Path | None:
        """
        Path to linked structures archive if it exists

        Returns
        -------
        Path | None
            path to linked structures archive
        """
        self.linked_structures
        return self._linked_archive

    def get_linked_structure(self, link_kind: str, link_id: str) -> str:
        """
        Get the path to the requested linked structure

        Parameters
        ----------
        link_kind : str
            kind of linked structure ('apo' or 'pred')
        link_id : str
            id of linked structure

        Returns
        -------
        str
            path to linked structure
        """
        if self.linked_archive is None:
            raise ValueError("linked_archive is None!")
        structure = (
            self.linked_archive
            / link_kind
            / self.system_id
            / link_id
            / "superposed.cif"
        )
        if not structure.is_file():
            raise ValueError(f"structure={structure} does not exist!")
        return structure.as_posix()
