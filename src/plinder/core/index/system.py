# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    import biotite.structure as struc
    from rdkit import Chem

from biotite.sequence.io.fasta import FastaFile

from plinder.core.index.query import query_table
from plinder.core.release import PlinderRelease
from plinder.core.structure.structure import Structure
from plinder.core.utils.config import get_config
from plinder.core.utils.io import get_pdb_mmcif
from plinder.data.annotations.save_utils import (
    ReconstructedSystem,
    SystemReconstructionOptions,
    SystemReconstructionOutputs,
    reconstruct_system,
    save_ligands,
    save_reconstructed_chain,
    save_reconstructed_system,
)


def _ligand_sdf_groups(system_df: pd.DataFrame) -> dict[str, list[str]]:
    """Map each ligand's primary instance-chain to all its member instance-chains.

    A multi-chain covalent ligand spans several chains; its reconstructed SDF must
    contain all of them (matching the ingest writer and the stored SMILES), so we
    read the serialized ``ligand_instance_chains``. Falls back to the primary chain
    alone when that list is absent (single-chain ligands / older indexes).
    """
    groups: dict[str, list[str]] = {}
    for _, row in system_df.iterrows():
        primary = row.get("ligand_instance_chain")
        if not isinstance(primary, str) or not primary:
            primary = f"{row['ligand_instance']}.{row['ligand_asym_id']}"
        try:
            members = [str(chain) for chain in row.get("ligand_instance_chains")]
        except TypeError:  # missing column / NaN -> not iterable
            members = []
        groups[primary] = sorted(members) if members else [primary]
    return groups


def _extract_packed_ligand_sdfs(
    *, archive: Path, pdb_id: str, asym_ids: set[str]
) -> Path:
    """Extract only one system's canonical SDFs from its shard parquet."""
    folder = archive.parent / pdb_id / "ligand_files"
    missing = {
        asym_id for asym_id in asym_ids if not (folder / f"{asym_id}.sdf").is_file()
    }
    if not missing:
        return folder
    packed = pd.read_parquet(
        archive,
        columns=["ligand_asym_id", "sdf"],
        filters=[("pdb_id", "==", pdb_id), ("ligand_asym_id", "in", sorted(missing))],
    )
    observed = set(packed["ligand_asym_id"].astype(str))
    absent = sorted(missing.difference(observed))
    if absent:
        raise ValueError(f"canonical ligand SDFs are missing for {pdb_id}: {absent}")
    folder.mkdir(exist_ok=True, parents=True)
    for row in packed.itertuples(index=False):
        target = folder / f"{row.ligand_asym_id}.sdf"
        temporary = target.with_suffix(".tmp.sdf")
        temporary.write_bytes(row.sdf)
        temporary.replace(target)
    return folder


class PlinderSystem:
    """
    Core class for interacting with a single system and its assets.
    Annotation data is queried lazily for one entry or system. Release structure
    views are reconstructed from a deposited PDB mmCIF cached under the
    configured PLINDER directory; canonical ASU ligand SDFs are loaded from the
    ligand archive independently.  An explicit source mmCIF can override the
    managed cache.

    """

    def __repr__(self) -> str:
        return f"PlinderSystem(system_id={self.system_id})"

    def __init__(
        self,
        *,
        system_id: str,
        release: PlinderRelease | None = None,
        prune: bool = True,
        skip_3d_confgen: bool = False,
        complete_missing_atoms: bool = False,
        source_mmcif: Path | str | None = None,
        reconstruction_dir: Path | str | None = None,
        canonical_ligand_dir: Path | str | None = None,
        biounit_chains: pd.DataFrame | None = None,
        reconstruction_options: SystemReconstructionOptions = (
            SystemReconstructionOptions()
        ),
    ) -> None:
        self.system_id: str = system_id
        self.release = release or PlinderRelease()
        self.prune: bool = prune
        self.skip_3d_confgen: bool = skip_3d_confgen
        self.complete_missing_atoms: bool = complete_missing_atoms
        self.source_mmcif = Path(source_mmcif) if source_mmcif is not None else None
        default_reconstruction_root = (
            Path(self.release.data_dir)
            if self.release.data_dir is not None
            else Path(get_config().data.plinder_dir)
        )
        self.reconstruction_dir = (
            Path(reconstruction_dir)
            if reconstruction_dir is not None
            else default_reconstruction_root / "reconstructed_systems" / system_id
        )
        self.canonical_ligand_dir = (
            Path(canonical_ligand_dir) if canonical_ligand_dir is not None else None
        )
        self.reconstruction_options = reconstruction_options
        self._entry: pd.DataFrame | None = None
        self._system: pd.DataFrame | None = None
        self._entry_chains: pd.DataFrame | None = None
        self._linked_apo_structures: pd.DataFrame | None = None
        self._biounit_chains = (
            biounit_chains.copy() if biounit_chains is not None else None
        )
        self._reconstructed: ReconstructedSystem | None = None
        self._canonical_ligand_folder: Path | None = None

    @property
    def entry(self) -> pd.DataFrame:
        """
        Store the annotation rows for this PDB entry.

        Returns
        -------
        pd.DataFrame
            Ligand-level annotation rows for the entry.
        """
        if self._entry is None:
            entry_pdb_id = self.system_id.split("__")[0]
            self._entry = query_table(
                "annotation",
                columns=["*"],
                filters=[("entry_pdb_id", "==", entry_pdb_id)],
                joins=["entry_metadata"],
                release=self.release,
            )
            if self._entry.empty:
                raise ValueError(
                    f"pdb_id={entry_pdb_id} not found in the annotation index"
                )
        return self._entry

    @property
    def system(self) -> pd.DataFrame:
        """
        Return ligand-level annotation rows for this system.

        Returns
        -------
        pd.DataFrame
            Annotation rows for the system.
        """
        if self._system is None:
            self._system = query_table(
                "annotation",
                columns=["*"],
                filters=[("system_id", "==", self.system_id)],
                joins=["entry_metadata"],
                release=self.release,
            )
            if self._system.empty:
                raise ValueError(f"system_id={self.system_id} not found in the index")
        return self._system

    @property
    def receptor_type(self) -> str:
        """Return the receptor polymer composition recorded during ingestion."""
        column = "system_receptor_type"
        if column not in self.system.columns:
            raise ValueError(
                f"{column} is unavailable for {self.system_id}; "
                "this release predates receptor-type annotation"
            )
        values = self.system[column].dropna().astype(str).unique()
        if len(values) != 1 or not values[0]:
            raise ValueError(
                f"Expected one receptor type for {self.system_id}, "
                f"got {values.tolist()}"
            )
        return str(values[0])

    @property
    def receptor_chain_types(self) -> dict[str, str]:
        """Map each system receptor instance chain to its polymer type."""
        if self._entry_chains is None:
            path = self.release.fetch("entry_chains")
            pdb_id = self.system_id.split("__", maxsplit=1)[0]
            self._entry_chains = pd.read_parquet(
                path,
                filters=[("entry_pdb_id", "==", pdb_id)],
            )
        required = {"chain_asym_id", "chain_receptor_type"}
        missing_columns = required.difference(self._entry_chains.columns)
        if missing_columns:
            raise ValueError(
                "Entry-chain receptor types are unavailable; missing columns "
                f"{sorted(missing_columns)}"
            )
        asym_to_type = dict(
            zip(
                self._entry_chains["chain_asym_id"].astype(str),
                self._entry_chains["chain_receptor_type"].astype(str),
            )
        )
        instance_chains = list(self.system.iloc[0]["system_protein_chains_asym_id"])
        result = {}
        for instance_chain in instance_chains:
            asym_id = str(instance_chain).split(".", maxsplit=1)[-1]
            if asym_id not in asym_to_type:
                raise ValueError(
                    f"No receptor-chain type found for {instance_chain} in "
                    f"{self.system_id}"
                )
            result[str(instance_chain)] = asym_to_type[asym_id]
        return result

    @property
    def biounit_chains(self) -> pd.DataFrame:
        """Return ingested chain membership for this biological assembly."""
        if self._biounit_chains is None:
            path = self.release.fetch("entry_biounit_chains")
            row = self.system.iloc[0]
            self._biounit_chains = pd.read_parquet(
                path,
                filters=[
                    ("entry_pdb_id", "==", str(row["entry_pdb_id"])),
                    ("biounit_id", "==", str(row["system_biounit_id"])),
                ],
            )
        if self._biounit_chains.empty:
            raise ValueError(
                f"No biological-assembly chain metadata for {self.system_id}"
            )
        return self._biounit_chains

    @property
    def source_mmcif_path(self) -> Path:
        """Return the explicit or release-cached deposited PDB mmCIF."""
        if self.source_mmcif is None:
            manifest = self.release.fetch("entry_sources")
            cache_dir = (
                Path(self.release.data_dir) / get_config().data.source_mmcifs
                if self.release.data_dir is not None
                else None
            )
            self.source_mmcif = get_pdb_mmcif(
                self.system_id,
                cache_dir=cache_dir,
                manifest_path=manifest,
            )
        if not self.source_mmcif.is_file():
            raise FileNotFoundError(self.source_mmcif)
        return self.source_mmcif

    @property
    def reconstructed(self) -> ReconstructedSystem:
        """In-memory biological-assembly views reconstructed from source mmCIF."""
        if self._reconstructed is None:
            self._reconstructed = reconstruct_system(
                self.source_mmcif_path,
                self.system.iloc[0],
                biounit_chains=self.biounit_chains,
                options=self.reconstruction_options,
            )
        return self._reconstructed

    def reconstruct(
        self,
        *,
        outputs: SystemReconstructionOutputs,
        options: SystemReconstructionOptions | None = None,
        overwrite: bool = False,
    ) -> dict[str, Path]:
        """Write exactly the requested reconstructed assets."""
        selected_options = options or self.reconstruction_options
        cached = (
            self.reconstructed
            if selected_options == self.reconstruction_options
            else None
        )
        return save_reconstructed_system(
            self.source_mmcif_path,
            self.system.iloc[0],
            outputs=outputs,
            biounit_chains=self.biounit_chains,
            options=selected_options,
            overwrite=overwrite,
            reconstructed=cached,
        )

    def _ensure_standard_output(self, name: str) -> Path:
        path = self.reconstruction_dir / name
        if not path.is_file():
            output_name = {
                "system.cif": "system_cif",
                "receptor.cif": "receptor_cif",
                "sequences.fasta": "sequences_fasta",
            }[name]
            self.reconstruct(outputs=SystemReconstructionOutputs(**{output_name: path}))
        return path

    @property
    def archive(self) -> Path:
        """
        Return the path to the directory containing the plinder system

        Returns
        -------
        Path
            directory containing the plinder system
        """
        self.reconstruction_dir.mkdir(parents=True, exist_ok=True)
        return self.reconstruction_dir

    @property
    def system_cif(self) -> str:
        """
        Path to the system.cif file

        Returns
        -------
        str
            path
        """
        return self._ensure_standard_output("system.cif").as_posix()

    @property
    def receptor_cif(self) -> str:
        """
        Path to the receptor.cif file

        Returns
        -------
        str
            path
        """
        return self._ensure_standard_output("receptor.cif").as_posix()

    @property
    def sequences_fasta(self) -> str:
        """
        Path to the sequences.fasta file

        Returns
        -------
        str
            path
        """
        return self._ensure_standard_output("sequences.fasta").as_posix()

    @cached_property
    def sequences(self) -> dict[str, str]:
        """
        Parsed FASTA contents from ``sequences.fasta``.

        Returns
        -------
        dict[str, str]
            Mapping from chain ID to sequence.
        """
        return {k: v for k, v in FastaFile.read_iter(self.sequences_fasta)}

    @property
    def canonical_ligand_sdfs(self) -> dict[str, str]:
        """Map system ligand instances to the single canonical ASU SDF copy."""
        if self._canonical_ligand_folder is None:
            if self.canonical_ligand_dir is not None:
                folder = self.canonical_ligand_dir
            else:
                pdb_id = self.system_id.split("__", maxsplit=1)[0]
                asym_ids = set(self.system["ligand_asym_id"].astype(str))
                code = pdb_id[1:3]
                archive = self.release.fetch(
                    "ligand_archive",
                    shard=code,
                )
                folder = _extract_packed_ligand_sdfs(
                    archive=archive,
                    pdb_id=pdb_id,
                    asym_ids=asym_ids,
                )
            if not folder.is_dir():
                raise ValueError(f"canonical ligand directory does not exist: {folder}")
            self._canonical_ligand_folder = folder

        ligands = {}
        for _, ligand in self.system.iterrows():
            instance_chain = ligand.get("ligand_instance_chain")
            if not isinstance(instance_chain, str) or not instance_chain:
                instance_chain = (
                    f"{ligand['ligand_instance']}.{ligand['ligand_asym_id']}"
                )
            ligand_file = self._canonical_ligand_folder / (
                f"{ligand['ligand_asym_id']}.sdf"
            )
            if not ligand_file.is_file():
                raise ValueError(f"canonical ligand SDF does not exist: {ligand_file}")
            ligands[instance_chain] = ligand_file.as_posix()
        return ligands

    @property
    def ligand_sdfs(self) -> dict[str, str]:
        """
        Return a dictionary of ligand names to paths to ligand sdf files

        Returns
        -------
        dict[str, str]
            dictionary of ligand names to paths to ligand sdf files
        """
        ligand_dir = self.reconstruction_dir / "ligand_files"
        instance_chains = list(self.canonical_ligand_sdfs)
        missing = [
            chain
            for chain in instance_chains
            if not (ligand_dir / f"{chain}.sdf").is_file()
        ]
        if missing:
            # Pass the member-spanning groups (not a bare chain list) so a
            # multi-chain covalent ligand is written as one whole molecule.
            groups = _ligand_sdf_groups(self.system)
            save_ligands(
                self.reconstructed.biounit,
                {chain: groups.get(chain, [chain]) for chain in missing},
                ligand_dir,
            )
        ligands = {
            chain: (ligand_dir / f"{chain}.sdf").as_posix() for chain in instance_chains
        }
        absent = [path for path in ligands.values() if not Path(path).is_file()]
        if absent:
            raise ValueError(f"failed to reconstruct ligand SDFs: {absent}")
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
        return [path.as_posix() for path in self.archive.rglob("*") if path.is_file()]

    @property
    def linked_apo_structures(self) -> pd.DataFrame:
        """Return ranked deposited apo chains linked to this holo system."""
        if self._linked_apo_structures is None:
            path = self.release.fetch("linked_apo_structures")
            self._linked_apo_structures = pd.read_parquet(
                path,
                filters=[("reference_system_id", "==", self.system_id)],
            ).sort_values("rank", ignore_index=True)
        return self._linked_apo_structures

    def _linked_apo_row(self, linked_structure_id: str | None = None) -> pd.Series:
        links = self.linked_apo_structures
        if linked_structure_id is None:
            selected = links.head(1)
        else:
            selected = links.loc[
                links["linked_structure_id"].astype(str).eq(linked_structure_id)
            ]
        if selected.empty:
            requested = linked_structure_id or "rank 1"
            raise ValueError(
                f"No linked apo structure {requested!r} for {self.system_id}"
            )
        if len(selected) != 1:
            raise ValueError(
                f"Linked apo ID {linked_structure_id!r} is not unique for "
                f"{self.system_id}"
            )
        return selected.iloc[0]

    def _linked_apo_source_mmcif(
        self, row: pd.Series, source_mmcif: Path | str | None
    ) -> Path | str:
        if source_mmcif is not None:
            return source_mmcif
        manifest = self.release.fetch("entry_sources")
        cache_dir = (
            Path(self.release.data_dir) / get_config().data.source_mmcifs
            if self.release.data_dir is not None
            else None
        )
        return get_pdb_mmcif(
            str(row["source_entry_id"]),
            cache_dir=cache_dir,
            manifest_path=manifest,
        )

    def reconstruct_linked_apo(
        self,
        linked_structure_id: str | None = None,
        *,
        output_cif: Path | str | None = None,
        source_mmcif: Path | str | None = None,
        overwrite: bool = False,
    ) -> Path:
        """Reconstruct one linked apo chain from its deposited source mmCIF.

        When ``linked_structure_id`` is omitted, the highest-ranked link is
        selected. The result contains the exact biological-assembly chain that
        was scored, with source sequence and chemical-component metadata.
        """
        row = self._linked_apo_row(linked_structure_id)
        link_id = str(row["linked_structure_id"])
        if output_cif is None:
            output_cif = self.reconstruction_dir / "linked_apo" / f"{link_id}.cif"
        return save_reconstructed_chain(
            self._linked_apo_source_mmcif(row, source_mmcif),
            assembly_id=str(row["source_biounit_id"]),
            chain_instance=str(row["source_chain_instance"]),
            source_asym_id=str(row["source_chain_asym_id"]),
            output_cif=output_cif,
            structure_id=link_id,
            overwrite=overwrite,
        )

    def superpose_linked_apo(
        self,
        linked_structure_id: str | None = None,
        *,
        reference_chain: str | None = None,
        output_cif: Path | str | None = None,
        source_mmcif: Path | str | None = None,
        overwrite: bool = False,
    ) -> Path:
        """Reconstruct and fit one linked apo chain to a holo receptor chain.

        The reference chain is inferred when the holo receptor contains one
        protein chain. Multichain receptors require ``reference_chain`` so the
        fit is never chosen arbitrarily.
        """
        import biotite.structure as struc

        receptor = self.receptor_structure
        ca_mask = struc.filter_amino_acids(receptor) & (
            receptor.atom_name.astype(str) == "CA"
        )
        protein_chains = sorted(set(receptor.chain_id[ca_mask].astype(str)))
        if reference_chain is None:
            if len(protein_chains) != 1:
                raise ValueError(
                    "reference_chain is required when the holo receptor has "
                    f"{len(protein_chains)} protein chains: {protein_chains}"
                )
            reference_chain = protein_chains[0]
        if reference_chain not in protein_chains:
            raise ValueError(
                f"Reference chain {reference_chain!r} is not a protein chain in "
                f"{self.system_id}; available chains are {protein_chains}"
            )
        reference_atoms = receptor[
            receptor.chain_id.astype(str) == str(reference_chain)
        ]

        row = self._linked_apo_row(linked_structure_id)
        link_id = str(row["linked_structure_id"])
        if output_cif is None:
            output_cif = (
                self.reconstruction_dir / "linked_apo" / f"{link_id}_superposed.cif"
            )
        return save_reconstructed_chain(
            self._linked_apo_source_mmcif(row, source_mmcif),
            assembly_id=str(row["source_biounit_id"]),
            chain_instance=str(row["source_chain_instance"]),
            source_asym_id=str(row["source_chain_asym_id"]),
            output_cif=output_cif,
            structure_id=link_id,
            superpose_to=reference_atoms,
            overwrite=overwrite,
        )

    @cached_property
    def receptor_structure(self) -> "struc.AtomArray":
        """
        Return the receptor structure as biotite AtomArray.
        """
        import biotite.structure.io.pdbx as pdbx
        from biotite.structure import filter_heavy

        from plinder.data.annotations.cif_utils import read_mmcif_file

        cif_file = read_mmcif_file(self.receptor_cif)
        atoms = pdbx.get_structure(
            cif_file, model=1, use_author_fields=False, include_bonds=True
        )
        return atoms[filter_heavy(atoms)]

    @cached_property
    def ligand_structures(self) -> dict[str, "struc.AtomArray"]:
        """Return heavy-atom ligand structures as Biotite atom arrays."""
        from biotite.structure import filter_heavy
        from biotite.structure.io.mol import SDFile

        structures = {}
        for chain, ligand_sdf in self.ligand_sdfs.items():
            atoms = SDFile.read(ligand_sdf).record.get_structure()
            structures[chain] = atoms[filter_heavy(atoms)]
        return structures

    @cached_property
    def ligand_mols(self) -> dict[str, "Chem.Mol"]:
        """
        Return the ligand molecules as RDKit Mol objects.
        """
        # TODO(peppr): temporary local sanitize carrying boron/main-group
        # over-valence fixes not yet in a released peppr. Revert to
        # `from peppr import sanitize as peppr_sanitize` once upstream.
        # See plinder.core.utils.sanitize.
        from rdkit import Chem

        from plinder.core.utils.sanitize import sanitize as peppr_sanitize

        mols = {}
        for chain in self.ligand_sdfs:
            supplier = Chem.SDMolSupplier(self.ligand_sdfs[chain], sanitize=False)
            mol = next(supplier, None)
            if mol is not None:
                peppr_sanitize(mol)
                mol = Chem.RemoveAllHs(mol, sanitize=False)
                mols[chain] = mol
        return mols

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
            skip_3d_confgen=self.skip_3d_confgen,
            structure_type="holo",
            complete_missing_atoms=self.complete_missing_atoms,
        )

    @cached_property
    def smiles(self) -> dict[str, str] | None:
        smiles_dict = {}
        try_smiles = [
            "ligand_smiles",
            "ligand_resolved_smiles",
        ]
        for _, ligand in self.system.iterrows():
            ligand_key = f"{ligand['ligand_instance']}.{ligand['ligand_asym_id']}"
            found = False
            for s in try_smiles:
                value = ligand.get(s)
                if value is not None and isinstance(value, str) and value:
                    smiles_dict[ligand_key] = value
                    found = True
                    break
            if not found:
                smiles_dict[ligand_key] = ""
        return smiles_dict
