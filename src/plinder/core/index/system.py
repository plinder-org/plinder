# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    from ost import mol

from biotite.sequence.io.fasta import FastaFile

from plinder.core.index import utils
from plinder.core.scores.links import query_links
from plinder.core.structure.structure import Structure
from plinder.core.utils.cpl import get_plinder_path
from plinder.core.utils.io import (
    download_alphafold_cif_file,
    download_pdb_chain_cif_file,
)
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
        self.system_id: str = system_id
        self.prune: bool = prune
        self._entry: dict[str, Any] | None = None
        self._system: dict[str, Any] | None = None
        self._archive: Path | None = None
        self._chain_mapping: dict[str, Any] | None = None
        self._water_mapping: dict[str, Any] | None = None
        self._linked_structures: pd.DataFrame | None = None
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
            self._archive = archive.parent / self.system_id
            if not self._archive.is_dir():
                raise ValueError(f"system_id={self.system_id} not found in systems")
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
    def sequences_fasta(self) -> str:
        """
        Path to the sequences.fasta file

        Returns
        -------
        str
            path
        """
        assert self.archive is not None
        return (self.archive / "sequences.fasta").as_posix()

    @cached_property
    def sequences(self) -> dict[str, str]:
        """
        Path to the sequences.fasta file

        Returns
        -------
        str
            path
        """
        assert self.archive is not None
        return {k: v for k, v in FastaFile.read_iter(self.sequences_fasta)}

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
            try:
                with (self.archive / "water_mapping.json").open() as f:
                    self._water_mapping = json.load(f)
            except FileNotFoundError:
                return None
        return self._water_mapping

    @property
    def ligand_sdfs(self) -> dict[str, str]:
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
        if self._linked_archive is None:
            zips = get_zips_to_unpack(kind="linked_structures")
            if not len(zips):
                LOG.info("no linked_structures found, downloading now, stand by")
                get_plinder_path(rel="linked_structures")
                zips = get_zips_to_unpack(kind="linked_structures")
            archive = list(zips.keys())[0]
            self._linked_archive = archive.parent

        return self._linked_archive

    def get_linked_structure(self, link_kind: str, link_id: str) -> str:
        """
        Get the path to the requested linked structure

        Parameters
        ----------
        link_kind : str
            kind of linked structure ('apo', 'pred', 'holo')
        link_id : str
            id of linked structure

        Returns
        -------
        str
            path to linked structure
        """
        if self.linked_archive is None:
            raise ValueError("linked_archive is None!")
        allowed = ["apo", "pred", "holo"]
        assert link_kind in allowed, f"link_kind={link_kind} not in {allowed}"
        structure = self.linked_archive / f"{link_id}.cif"
        if not structure.is_file():
            if link_kind == "apo":
                pdb_id, chain_id = link_id.split("_")
                try:
                    download_pdb_chain_cif_file(pdb_id, chain_id, structure)
                except Exception as e:
                    raise ValueError(f"Unable to download {link_id}! {str(e)}")
            elif link_kind == "pred":
                uniprot_id = link_id.split("_")[0]
                cif_file_path = download_alphafold_cif_file(
                    uniprot_id, self.linked_archive
                )
                if cif_file_path is None:
                    raise ValueError(f"Unable to download {link_id}")
                cif_file_path.rename(structure)
            elif link_kind == "holo":
                structure = Path(PlinderSystem(system_id=link_id).receptor_cif)
            if structure is None or not structure.is_file():
                raise ValueError(f"structure={structure} does not exist!")
        return structure.as_posix()

    @cached_property
    def receptor_entity(self) -> "mol.EntityHandle":
        """
        Return the receptor entity handle

        Returns
        -------
        mol.EntityHandle
            receptor entity handle
        """
        try:
            from ost import io
        except ImportError:
            raise ImportError("Please install openstructure to use this property")
        return io.LoadMMCIF(self.receptor_cif)

    @cached_property
    def ligand_views(self) -> dict[str, "mol.ResidueView"]:
        """
        Return the ligand views

        Returns
        -------
        dict[str, mol.ResidueView]
        """
        try:
            from ost import io
        except ImportError:
            raise ImportError("Please install openstructure to use this property")

        ligand_views = {}
        for chain in self.ligand_sdfs:
            ligand_views[chain] = io.LoadEntity(
                self.ligand_sdfs[chain], format="sdf"
            ).Select("ele != H")
        return ligand_views

    @property
    def num_ligands(self) -> int:
        """
        Return the number of ligands in the system

        Returns
        -------
        int
        """
        return len(self.ligand_sdfs)

    @property
    def num_proteins(self) -> int:
        """
        Return the number of proteins in the system

        Returns
        -------
        int
        """
        return len(self.system_id.split("__")[2].split("_"))

    @property
    def alternate_structures(self) -> dict[str, Structure]:
        """
        load all alternate structures
        """
        # TODO: do we want to keep this as assertion?
        # better if then raise?
        assert self.linked_structures is not None

        structures = {}
        for id, kind in self.linked_structures[["id", "kind"]].values:
            protein_path = self.get_linked_structure(
                kind,
                id,
            )
            structures[id] = Structure(
                id=id,
                protein_path=Path(protein_path),
                protein_sequence=self.sequences,
                structure_type=kind,
            )
        return structures

    @property
    def holo_structure(self) -> Structure:
        """
        Load holo structure
        """
        return Structure(
            id=self.system_id,
            protein_path=Path(self.receptor_cif),
            protein_sequence=self.sequences,
            ligand_sdfs=self.ligand_sdfs,
            ligand_smiles=self.smiles,
            structure_type="holo",
        )

    @cached_property
    def smiles(self) -> dict[str, str] | None:
        if self.system is None:
            return None
        smiles_dict = {}
        try_smiles = ["rdkit_canonical_smiles", "smiles", "resolved_smiles"]
        for ligand in self.system["ligands"]:
            found = False
            for s in try_smiles:
                if s in ligand and ligand[s] is not None and len(ligand[s]):
                    smiles_dict[f"{ligand['instance']}.{ligand['asym_id']}"] = ligand[s]
                    found = True
                    break
            if not found:
                smiles_dict[f"{ligand['instance']}.{ligand['asym_id']}"] = ""
        return smiles_dict
