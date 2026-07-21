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

from plinder.core.scores import query_index
from plinder.core.scores.links import query_links
from plinder.core.scores.query import FILTER
from plinder.core.structure.structure import Structure
from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.cpl import get_plinder_path
from plinder.core.utils.io import (
    download_alphafold_cif_file,
    download_pdb_chain_cif_file,
    get_pdb_mmcif,
)
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import get_zips_to_unpack
from plinder.data.utils.annotations.save_utils import (
    ReconstructedSystem,
    SystemReconstructionOptions,
    SystemReconstructionOutputs,
    reconstruct_system,
    save_ligands,
    save_reconstructed_system,
)

LOG = setup_logger(__name__)


def _materialize_packed_ligand_sdfs(
    *, archive: Path, pdb_id: str, asym_ids: set[str]
) -> Path:
    """Materialize only one system's canonical SDFs from its shard Parquet."""
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
    Annotation data is queried lazily for one entry or system. V3 structure
    views are reconstructed from a deposited PDB mmCIF cached under the
    configured PLINDER directory; canonical ASU ligand SDFs are loaded from the
    ligand archive independently.  An explicit source mmCIF can override the
    managed cache.

    Existing local V2 system archives remain readable as a transitional
    fallback, but are never downloaded by this class.
    """

    def __repr__(self) -> str:
        return f"PlinderSystem(system_id={self.system_id})"

    def __init__(
        self,
        *,
        system_id: str,
        prune: bool = True,
        skip_3d_confgen: bool = False,
        source_mmcif: Path | str | None = None,
        reconstruction_dir: Path | str | None = None,
        canonical_ligand_dir: Path | str | None = None,
        biounit_chains: pd.DataFrame | None = None,
        reconstruction_options: SystemReconstructionOptions = (
            SystemReconstructionOptions()
        ),
    ) -> None:
        self.system_id: str = system_id
        self.prune: bool = prune
        self.skip_3d_confgen: bool = skip_3d_confgen
        self.source_mmcif = Path(source_mmcif) if source_mmcif is not None else None
        cfg = get_config()
        self.reconstruction_dir = (
            Path(reconstruction_dir)
            if reconstruction_dir is not None
            else Path(cfg.data.plinder_dir) / "reconstructed_systems" / system_id
        )
        self.canonical_ligand_dir = (
            Path(canonical_ligand_dir) if canonical_ligand_dir is not None else None
        )
        self.reconstruction_options = reconstruction_options
        self._entry: pd.DataFrame | None = None
        self._system: pd.DataFrame | None = None
        self._entry_chains: pd.DataFrame | None = None
        self._biounit_chains = (
            biounit_chains.copy() if biounit_chains is not None else None
        )
        self._archive: Path | None = None
        self._reconstructed: ReconstructedSystem | None = None
        self._canonical_ligand_folder: Path | None = None
        self._linked_structures: pd.DataFrame | None = None
        self._linked_archive: Path | None = None

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
            self._entry = query_index(
                columns=["*"],
                splits=["*"],
                filters=[FILTER(("entry_pdb_id", "==", entry_pdb_id))],
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
            self._system = query_index(
                columns=["*"],
                splits=["*"],
                filters=[FILTER(("system_id", "==", self.system_id))],
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
                f"Expected one receptor type for {self.system_id}, got {values.tolist()}"
            )
        return str(values[0])

    @property
    def receptor_chain_types(self) -> dict[str, str]:
        """Map each system receptor instance chain to its polymer type."""
        if self._entry_chains is None:
            cfg = get_config()
            path = cpl.get_plinder_path(
                rel=f"{cfg.data.index}/{cfg.data.entry_chain_file}"
            )
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
        """Return normalized chain membership for this biological assembly."""
        if self._biounit_chains is None:
            cfg = get_config()
            path = cpl.get_plinder_path(
                rel=f"{cfg.data.index}/{cfg.data.entry_biounit_chain_file}"
            )
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

    def _reconstruction_biounit_chains(self) -> pd.DataFrame | None:
        """Use normalized membership for V3 while leaving V2 assets unchanged."""
        if self._biounit_chains is not None:
            return self._biounit_chains
        if str(get_config().data.plinder_iteration).lower() != "v3":
            return None
        return self.biounit_chains

    def _legacy_archive(self) -> Path | None:
        """Return an already-local V2 archive without fetching one."""
        cfg = get_config()
        root = Path(cpl.get_plinder_path(rel=cfg.data.systems, download=False))
        extracted = root / self.system_id
        if (extracted / "receptor.cif").is_file():
            return extracted
        zip_path = root / f"{self.system_id[1:3]}.zip"
        if not zip_path.is_file():
            return None
        get_zips_to_unpack(kind="systems", system_ids=[self.system_id])
        return extracted if (extracted / "receptor.cif").is_file() else None

    def _require_source_mmcif(self) -> Path:
        if self.source_mmcif is None:
            self.source_mmcif = get_pdb_mmcif(self.system_id)
        if not self.source_mmcif.is_file():
            raise FileNotFoundError(self.source_mmcif)
        return self.source_mmcif

    @property
    def source_mmcif_path(self) -> Path:
        """Return the explicit or release-cached deposited PDB mmCIF."""
        return self._require_source_mmcif()

    @property
    def _uses_source_reconstruction(self) -> bool:
        """Use source reconstruction for explicit inputs and V3 releases."""
        return (
            self.source_mmcif is not None
            or str(get_config().data.plinder_iteration).lower() == "v3"
        )

    @property
    def reconstructed(self) -> ReconstructedSystem:
        """In-memory biological-assembly views reconstructed from source mmCIF."""
        if self._reconstructed is None:
            self._reconstructed = reconstruct_system(
                self._require_source_mmcif(),
                self.system.iloc[0],
                biounit_chains=self._reconstruction_biounit_chains(),
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
            self._require_source_mmcif(),
            self.system.iloc[0],
            outputs=outputs,
            biounit_chains=self._reconstruction_biounit_chains(),
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
    def archive(self) -> Path | None:
        """
        Return the path to the directory containing the plinder system

        Returns
        -------
        Path | None
            directory containing the plinder system
        """
        if self._archive is None:
            if self._uses_source_reconstruction:
                self.reconstruction_dir.mkdir(parents=True, exist_ok=True)
                self._archive = self.reconstruction_dir
            else:
                self._archive = self._legacy_archive()
            if self._archive is None:
                raise FileNotFoundError(
                    f"No local V2 system archive found for {self.system_id}. "
                    "Source-mmCIF reconstruction is enabled for V3 releases."
                )
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
        if self._uses_source_reconstruction:
            return self._ensure_standard_output("system.cif").as_posix()
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
        if self._uses_source_reconstruction:
            return self._ensure_standard_output("receptor.cif").as_posix()
        assert self.archive is not None
        return (self.archive / "receptor.cif").as_posix()

    @property
    def sequences_fasta(self) -> str:
        """
        Path to the sequences.fasta file

        Returns
        -------
        str
            path
        """
        if self._uses_source_reconstruction:
            return self._ensure_standard_output("sequences.fasta").as_posix()
        assert self.archive is not None
        return (self.archive / "sequences.fasta").as_posix()

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
                cfg = get_config()
                archive_root = Path(cfg.data.plinder_dir) / cfg.data.ligand_archives
                code = pdb_id[1:3]
                archive = cpl.get_plinder_path(
                    rel=f"{cfg.data.ligand_archives}/{code}.parquet"
                )
                folder = _materialize_packed_ligand_sdfs(
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
        if not self._uses_source_reconstruction:
            assert self.archive is not None
            return {
                ligand.stem: ligand.as_posix()
                for ligand in (self.archive / "ligand_files/").glob("*.sdf")
            }

        ligand_dir = self.reconstruction_dir / "ligand_files"
        instance_chains = list(self.canonical_ligand_sdfs)
        missing = [
            chain
            for chain in instance_chains
            if not (ligand_dir / f"{chain}.sdf").is_file()
        ]
        if missing:
            save_ligands(self.reconstructed.biounit, missing, ligand_dir)
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
    def receptor_structure(self) -> "struc.AtomArray":
        """
        Return the receptor structure as biotite AtomArray.
        """
        import biotite.structure.io.pdbx as pdbx

        from plinder.core.structure.atoms import is_hydrogen_isotope
        from plinder.data.utils.annotations.cif_utils import read_mmcif_file

        cif_file = read_mmcif_file(self.receptor_cif)
        atoms = pdbx.get_structure(
            cif_file, model=1, use_author_fields=False, include_bonds=True
        )
        return atoms[~is_hydrogen_isotope(atoms.element)]

    @cached_property
    def ligand_structures(self) -> dict[str, "struc.AtomArray"]:
        """Return heavy-atom ligand structures as Biotite atom arrays."""
        from biotite.structure.io.mol import SDFile

        from plinder.core.structure.atoms import is_hydrogen_isotope

        structures = {}
        for chain, ligand_sdf in self.ligand_sdfs.items():
            atoms = SDFile.read(ligand_sdf).record.get_structure()
            structures[chain] = atoms[~is_hydrogen_isotope(atoms.element)]
        return structures

    @cached_property
    def ligand_mols(self) -> dict[str, "Chem.Mol"]:
        """
        Return the ligand molecules as RDKit Mol objects.
        """
        from peppr import sanitize as peppr_sanitize
        from rdkit import Chem

        mols = {}
        for chain in self.ligand_sdfs:
            supplier = Chem.SDMolSupplier(self.ligand_sdfs[chain], sanitize=False)
            mol = next(supplier, None)
            if mol is not None:
                peppr_sanitize(mol)
                mol = Chem.RemoveAllHs(mol)
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
            skip_3d_confgen=self.skip_3d_confgen,
            structure_type="holo",
        )

    @cached_property
    def smiles(self) -> dict[str, str] | None:
        smiles_dict = {}
        try_smiles = [
            "ligand_rdkit_canonical_smiles",
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
