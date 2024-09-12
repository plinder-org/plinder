from __future__ import annotations

from functools import lru_cache
from pydantic import BaseModel

from pathlib import Path
from typing import Iterable, Optional, List, Union, Any, TYPE_CHECKING
import gzip
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolTransforms

import biotite.structure as struc
from biotite.structure.atoms import AtomArray
from biotite.structure.io.pdbx import get_structure
import biotite.structure.io as strucio
from biotite import TextFile
from numpy.typing import NDArray
from pydantic import ConfigDict
from pydantic.dataclasses import dataclass

# TODO: Decide whether to lift these from pinder or import them
from pinder.core.utils import setup_logger, constants as pc
from pinder.core.structure.superimpose import superimpose_chain
from pinder.core.utils.dataclass import stringify_dataclass
from pinder.core.structure.atoms import (
    get_seq_aligned_structures,
    get_per_chain_seq_alignments,
    invert_chain_seq_map,
    resn2seq,
    write_pdb,
)
from pinder.core.structure.contacts import get_atom_neighbors
from pinder.core.structure import surgery

from plinder.core.utils.config import get_config

if TYPE_CHECKING:
    import torch


log = setup_logger(__name__)
cfg = get_config()


def reverse_dict(mapping: dict[int, int]) -> dict[int, int]:
    return {v: k for k, v in mapping.items()}


def get_rdkit_mol(ligand: Path | str):
    if isinstance(ligand, Path):
        return Chem.MolFromSmiles(ligand)
    elif isinstance(ligand, Path):
        return next(Chem.SDMolSupplier(ligand))


def apply_tranform_to_rdkit_mol(mol, transform_matrix):
    rdMolTransforms.TransformConformer(mol.GetConformer(0), transform_matrix)
    return mol


def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import PDBxFile

    return PDBxFile


def atom_array_from_cif_file(
    structure: Path | AtomArray, use_author_fields: bool = True
) -> AtomArray | None:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_ciffile()
        try:
            if structure.suffix == ".gz":
                with gzip.open(str(structure), "rt", encoding="utf-8") as f:
                    mod = reader.read(f)
            else:
                mod = reader.read(structure)
            arr = get_structure(
                mod, model=1, use_author_fields=use_author_fields
            )  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


class Structure(BaseModel):
    protein_path: Path
    id: str
    list_ligand_sdf_or_smiles: Optional[List[Union[Path, str]]] = None
    protein_atom_array: Optional[AtomArray] = None
    ligand_mols: Optional[dict[str, Chem.rdchem.Mol]] = None
    add_ligand_hydrogen: bool = False
    structure_type: str = "holo"

    """Initialize structure.
    This dataclass provides abstraction over plinder systems holo, apo and predicted structures.
    It loads and processes receptor proteins and ligands in a way that allows us to capture the
    pecularitties of these molecule types.
    It's also intentionally design to be robust to capture structures, with or without ligand

    Parameters
    ----------
    protein_path : Path
        Protein structure cif path, could be holo, apo or pred path
    id : str
        PLINDER system or apo/pred id. If the structure is holo, this is used to load the holo structure
    list_ligand_sdf_or_smiles : list[Path | str] | None
        SDF path or ligand smiles. This is optional for apo and pred structures
    protein_atom_array : AtomArray | None = None
        Protein Biotite atom array
    ligand_mols : list[Chem.rdchem.Mol] | None = None
        Ligand rdkit mol
    add_ligand_hydrogen : bool = False
        Whether to add hydrogen to ligand or not
    structure_type : str = "holo"
        Structure type, "holo", "apo" or "pred"
    """

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def read_structure_files(
        protein_path: Path,
        list_ligand_sdf_or_smiles: list[Path | str] | None,
        add_ligand_hydrogen: bool = False,
    ) -> tuple[AtomArray | Chem.rdchem.Mol]:
        try:
            protein_arr = atom_array_from_cif_file(
                protein_path, use_author_fields=False
            )
        except:
            # Sometimes b-factor field parsing fails, try without it.
            protein_arr = atom_array_from_cif_file(
                protein_path, use_author_fields=False
            )
            annotation_arr: NDArray[np.double | np.str_] = np.repeat(
                0.0, protein_arr.shape[0]
            )
            protein_arr.set_annotation("b_factor", annotation_arr)
        if "" in set(protein_arr.element):
            try:
                protein_arr.element = struc.infer_elements(protein_arr)
            except Exception:
                log.debug(
                    f"Found missing elements in {protein_path} but failed to infer"
                )
        ligand_mols = {}
        if not (list_ligand_sdf_or_smiles is None):
            for ligand_sdf_or_smiles in list_ligand_sdf_or_smiles:
                if isinstance(ligand_sdf_or_smiles, Path) and (
                    ligand_sdf_or_smiles.suffix == ".sdf"
                ):
                    ligand_mol = next(Chem.SDMolSupplier(ligand_sdf_or_smiles))
                    if add_ligand_hydrogen:
                        ligand_mol = Chem.AddHs(ligand_mol, addCoords=True)
                    ligand_mols[ligand_sdf_or_smiles.stem] = ligand_mol
                elif isinstance(ligand_sdf_or_smiles, str):
                    ligand_mol = Chem.MolFromSmiles(ligand_sdf_or_smiles)
                    if add_ligand_hydrogen:
                        ligand_mol = Chem.AddHs(ligand_mol)
                    rdDistGeom.EmbedMolecule(ligand_mol)
                    ligand_mols[ligand_sdf_or_smiles] = ligand_mol
        return protein_arr, ligand_mols

    def __repr__(self) -> str:
        class_str: str = stringify_dataclass(self, 4)
        return class_str

    def __add__(self, other: Structure) -> Structure:
        combined_arr = self.protein_atom_array + other.protein_atom_array
        combined_id = f"{self.id}--{other.id}"
        combined_path = self.protein_path.parent / combined_id
        structure_type = "_".join(
            sorted(list({self.structure_type, other.structure_type}))
        )

        return Structure(
            protein_path=combined_path,
            id=combined_id,
            list_ligand_sdf_or_smiles=self.list_ligand_sdf_or_smiles,
            protein_atom_array=combined_arr,
            ligand_mols=self.ligand_mols,
            add_ligand_hydrogen=self.add_ligand_hydrogen,
            structure_type=structure_type,
        )

    def save_to_disk(self, filepath: Path | None = None) -> None:
        """Write Structure Atomarray to a PDB file.

        Parameters
        ----------
        filepath : Path | None
            Filepath to output PDB.
            If not provided, will write to self.protein_path,
            potentially overwriting if the file already exists!

        Returns
        -------
        None

        """
        if not filepath:
            filepath = self.protein_path
        write_pdb(self.protein_atom_array, filepath)

    def filter(
        self,
        property: str,
        mask: Iterable[bool | int | str],
        copy: bool = True,
        negate: bool = False,
    ) -> Structure | None:
        atom_mask = np.isin(getattr(self.protein_atom_array, property), list(mask))
        if negate:
            atom_mask = ~atom_mask
        if copy:
            arr = self.protein_atom_array[atom_mask].copy()
            return Structure(
                protein_path=self.protein_path,
                id=self.id,
                list_ligand_sdf_or_smiles=self.list_ligand_sdf_or_smiles,
                protein_atom_array=arr,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogen=self.add_ligand_hydrogen,
                structure_type=self.structure_type,
            )

        self.protein_atom_array = self.protein_atom_array[atom_mask]
        return None

    def get_per_chain_seq_alignments(
        self,
        other: Structure,
    ) -> dict[str, dict[int, int]]:
        self2other_seq: dict[str, dict[int, int]] = get_per_chain_seq_alignments(
            other.protein_atom_array, self.protein_atom_array
        )
        return self2other_seq

    def align_common_sequence(
        self,
        other: Structure,
        copy: bool = True,
        remove_differing_atoms: bool = True,
        renumber_residues: bool = False,
        remove_differing_annotations: bool = False,
    ) -> tuple[Structure, Structure]:
        ref_at = other.protein_atom_array.copy()
        target_at = self.protein_atom_array.copy()
        target2ref_seq = get_per_chain_seq_alignments(ref_at, target_at)
        ref2target_seq = invert_chain_seq_map(target2ref_seq)
        ref_at, target_at = get_seq_aligned_structures(ref_at, target_at)

        if remove_differing_atoms:
            # Even if atom counts are identical, annotation categories must be the same
            # First modify annotation arrays to use struc.filter_intersection,
            # then filter original structure with annotations to match res_id, res_name, atom_name
            # of intersecting structure
            ref_at_mod = ref_at.copy()
            target_at_mod = target_at.copy()
            ref_at_mod, target_at_mod = surgery.fix_annotation_mismatch(
                ref_at_mod, target_at_mod, ["element", "ins_code", "b_factor"]
            )
            ref_target_mask = struc.filter_intersection(ref_at_mod, target_at_mod)
            target_ref_mask = struc.filter_intersection(target_at_mod, ref_at_mod)
            if remove_differing_annotations:
                ref_at = ref_at_mod[ref_target_mask].copy()
                target_at = target_at_mod[target_ref_mask].copy()
            else:
                ref_at = ref_at[ref_target_mask].copy()
                target_at = target_at[target_ref_mask].copy()

        if not renumber_residues:
            target_at.res_id = np.array(
                [ref2target_seq[at.chain_id][at.res_id] for at in target_at]
            )

        if copy:
            self_struct = Structure(
                protein_path=self.protein_path,
                id=self.id,
                list_ligand_sdf_or_smiles=self.list_ligand_sdf_or_smiles,
                protein_atom_array=target_at,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogen=self.add_ligand_hydrogen,
                structure_type=self.structure_type,
            )

            other_struct = Structure(
                protein_path=other.protein_path,
                id=other.id,
                list_ligand_sdf_or_smiles=other.list_ligand_sdf_or_smiles,
                protein_atom_array=ref_at,
                ligand_mols=other.ligand_mols,
                add_ligand_hydrogen=other.add_ligand_hydrogen,
                structure_type=other.structure_type,
            )

            return self_struct, other_struct
        other.protein_atom_array = ref_at
        self.protein_atom_array = target_at
        return self, other

    def set_chain(self, chain_id: str) -> None:
        self.protein_atom_array.chain_id = np.repeat(
            chain_id, self.protein_atom_array.shape[0]
        )

    def superimpose(
        self,
        other: Structure,
    ) -> tuple[Structure, float, float]:
        # max_iterations=1 -> no outlier removal
        superimposed, _, other_anchors, self_anchors = _superimpose_common_atoms(
            other.protein_atom_array, self.protein_atom_array, max_iterations=1
        )
        raw_rmsd = struc.rmsd(
            other.protein_atom_array.coord[other_anchors],
            superimposed.coord[self_anchors],
        )

        superimposed, _, other_anchors, self_anchors = _superimpose_common_atoms(
            other.protein_atom_array, self.protein_atom_array
        )
        refined_rmsd = struc.rmsd(
            other.protein_atom_array.coord[other_anchors],
            superimposed.coord[self_anchors],
        )
        return (
            Structure(
                protein_path=self.protein_path,
                id=self.id,
                list_ligand_sdf_or_smiles=self.list_ligand_sdf_or_smiles,
                protein_atom_array=superimposed,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogen=self.add_ligand_hydrogen,
                structure_type=self.structure_type,
            ),
            raw_rmsd,
            refined_rmsd,
        )

    @property
    def protein_coords(self) -> NDArray[np.double]:
        """ndarray[np.double]: The coordinates of the protein atoms in the structure."""
        protein_coords: NDArray[np.double] = self.protein_atom_array.coord
        return protein_coords

    @property
    def ligand_coords(self) -> dict[str, NDArray[np.double]]:
        """ndarray[np.double]: The coordinates of the ligand atoms in the structure."""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.ligand_mols.items()
        }
        return ligand_coords

    @property
    def ligand_bonds(self) -> dict[str, list[tuple[int, int, str]]]:
        return {
            tag: [
                (at.GetBeginAtomIdx(), at.GetEndAtomIdx(), at.GetBondType())
                for at in mol.GetBonds()
            ]
            for tag, mol in self.ligand_mols.items()
        }

    @property
    def protein_dataframe(self) -> pd.DataFrame:
        """pd.DataFrame: The dataframe representation of the structure."""
        three_to_one = pc.three_to_one_noncanonical_mapping
        return pd.DataFrame(
            [
                {
                    "chain_id": at.chain_id,
                    "res_name": at.res_name,
                    "res_code": three_to_one.get(at.res_name, "X"),
                    "res_id": at.res_id,
                    "atom_name": at.atom_name,
                    "b_factor": at.b_factor,
                    "ins_code": at.ins_code,
                    "hetero": at.hetero,
                    "element": at.element,
                    "x": at.coord[0],
                    "y": at.coord[1],
                    "z": at.coord[2],
                }
                for at in self.protein_atom_array
            ]
        )

    @property
    def protein_backbone_mask(self) -> NDArray[np.bool_]:
        """ndarray[np.bool\_]: a logical mask for backbone atoms."""
        mask: NDArray[np.bool_] = struc.filter_peptide_backbone(self.protein_atom_array)
        return mask

    @property
    def protein_calpha_mask(self) -> NDArray[np.bool_]:
        """ndarray[np.bool\_]: a logical mask for alpha carbon atoms."""
        mask: NDArray[np.bool_] = self.protein_atom_array.atom_name == "CA"
        return mask

    @property
    def protein_n_atoms(self) -> int:
        """int: The number of atoms in the structure."""
        n: int = self.protein_atom_array.shape[0]
        return n

    @property
    def protein_chains(self) -> list[str]:
        """list[str]: The list of chain IDs in the structure."""
        ch_list = self._attr_from_atom_array(
            self.protein_atom_array, "chain_id", distinct=True, sort=True
        )
        return [str(ch) for ch in ch_list]

    @property
    def protein_chain_sequence(self) -> dict[str, list[str]]:
        """dict[str, list[str]]: The chain sequence dictionary, where keys are chain IDs and values are lists of residue codes."""
        ch_seq: dict[str, list[str]] = (
            self.protein_dataframe[["chain_id", "res_code", "res_name", "res_id"]]
            .drop_duplicates()
            .groupby("chain_id")["res_code"]
            .apply(list)
            .to_dict()
        )
        return ch_seq

    @property
    def protein_sequence(self) -> str:
        """str: The amino acid sequence of the structure."""
        numbering, resn = struc.get_residues(self.protein_atom_array)
        seq: str = resn2seq(resn)
        return seq

    @property
    def protein_fasta(self) -> str:
        """str: The fasta representation of the structure sequence."""
        fasta_str: str = "\n".join(
            [f">{self.protein_path.stem}", self.protein_sequence]
        )
        return fasta_str

    @property
    def protein_tokenized_sequence(self) -> "torch.Tensor":
        """torch.Tensor: The tokenized sequence representation of the structure sequence."""
        import torch

        seq_encoding = torch.tensor([pc.AA_TO_INDEX[x] for x in self.protein_sequence])
        tokenized: torch.Tensor = seq_encoding.long()
        return tokenized

    @property
    def protein_residue_names(self) -> list[str]:
        """list[str]: The list of distinct residue names in the structure."""
        res_list = self._attr_from_atom_array(
            self.protein_atom_array, "res_name", distinct=True, sort=True
        )
        return [str(r) for r in res_list]

    @property
    def protein_residues(self) -> list[int]:
        """list[int]: The list of distinct residue IDs in the structure."""
        res_list = self._attr_from_atom_array(
            self.protein_atom_array, "res_id", distinct=True, sort=True
        )
        return [int(r) for r in res_list]

    @property
    def protein_atom_names(self) -> list[str]:
        """list[str]: The list of distinct atom names in the structure."""
        at_list = self._attr_from_atom_array(
            self.protein_atom_array, "atom_name", distinct=True, sort=True
        )
        return [str(a) for a in at_list]

    @property
    def protein_b_factor(self) -> list[float]:
        """list[float]: A list of B-factor values for each atom in the structure."""
        b_factor = self._attr_from_atom_array(
            self.protein_atom_array, "b_factor", distinct=False, sort=False
        )
        return [float(b) for b in b_factor]

    @staticmethod
    def _attr_from_atom_array(
        array: AtomArray, attr: str, distinct: bool = False, sort: bool = False
    ) -> list[str] | list[int] | list[float]:
        prop = getattr(array, attr)
        if distinct:
            prop = set(prop)
        if sort:
            prop = sorted(prop)
        return list(prop)

    def model_post_init(self, __context: Any) -> None:
        # pydantic v2 renames this to dataclass post_init
        if (self.structure_type) == "holo":
            self.protein_path = (
                Path(cfg.data.plinder_dir) / "systems" / f"{self.id}" / "receptor.cif"
            )
            ligand_dir = (
                Path(cfg.data.plinder_dir) / "systems" / f"{self.id}" / "ligand_files"
            )
            self.list_ligand_sdf_or_smiles = list(ligand_dir.glob("*sdf"))

        if (self.protein_atom_array is None) & (self.ligand_mols is None):
            self.protein_atom_array, self.ligand_mols = self.read_structure_files(
                protein_path=self.protein_path,
                list_ligand_sdf_or_smiles=self.list_ligand_sdf_or_smiles,
                add_ligand_hydrogen=self.add_ligand_hydrogen,
            )
        try:
            getattr(self.protein_atom_array, "b-factor")
        except AttributeError:
            b_factors: NDArray[np.double | np.str_] = np.repeat(
                0.0, self.protein_atom_array.shape[0]
            )
            self.protein_atom_array.set_annotation("b_factor", b_factors)


def find_potential_interchain_bonded_atoms(
    structure: Structure,
    interface_res: dict[str, list[int]] | None = None,
    radius: float = 2.3,
) -> AtomArray:
    if interface_res is None:
        interface_res = structure.get_interface_residues(calpha_mask=False)
    interface_mask = structure.get_interface_mask(interface_res, calpha_only=False)
    interface = structure.atom_array[interface_mask].copy()
    assert set(interface_res.keys()) == {"R", "L"}
    interface_R = interface[(interface.chain_id == "R") & (interface.element != "H")]
    interface_L = interface[(interface.chain_id == "L") & (interface.element != "H")]
    L_neigh = get_atom_neighbors(interface_L, interface_R, radius=radius)
    R_neigh = get_atom_neighbors(interface_R, interface_L, radius=radius)
    interchain_at = R_neigh + L_neigh
    return interchain_at


def mask_common_uniprot(
    mono_A: Structure, mono_B: Structure
) -> tuple[Structure, Structure]:
    # Ensure mapping only contains those residues that are actually
    # present in the PDB file
    map_A = mono_A.resolved_mapping
    map_B = mono_B.resolved_mapping
    mask_A: list[int] | None = None
    mask_B: list[int] | None = None
    if all(isinstance(df, pd.DataFrame) for df in [map_A, map_B]):
        # Case when apo and holo uniprot mapping exists
        assert isinstance(map_A, pd.DataFrame)
        assert isinstance(map_B, pd.DataFrame)
        map_A.loc[:, "uniprot_uuid"] = (
            map_A.resi_uniprot.astype(str) + "-" + map_A.uniprot_acc.astype(str)
        )
        map_B.loc[:, "uniprot_uuid"] = (
            map_B.resi_uniprot.astype(str) + "-" + map_B.uniprot_acc.astype(str)
        )
        uniprot_mask = set(map_A.uniprot_uuid).intersection(set(map_B.uniprot_uuid))
        map_A = map_A[map_A.uniprot_uuid.isin(uniprot_mask)]
        map_B = map_B[map_B.uniprot_uuid.isin(uniprot_mask)]
        mask_A = sorted(list(set(map_A.resi.astype(int))))
        mask_B = sorted(list(set(map_B.resi.astype(int))))
    elif isinstance(map_A, pd.DataFrame):
        # Case where mono_B is predicted and already in uniprot numbering
        # Ensure uniprot used in holo mapping corresponds to predicted (could be chimeric)
        if isinstance(mono_B.pinder_id, str):
            mono_B_uniprot = mono_B.pinder_id.split("__")[1].split("-")[0]
            map_A = map_A.query(f"uniprot_acc == '{mono_B_uniprot}'").reset_index(
                drop=True
            )
        mask_B_resolved = set(map_A.resi_uniprot.astype(int)).intersection(
            set(mono_B.atom_array.res_id)
        )
        mask_B = sorted(list(mask_B_resolved))

        mask_A_resolved = set(map_A.resi.astype(int)).intersection(
            set(mono_A.atom_array.res_id)
        )
        mask_A = sorted(list(mask_A_resolved))

    elif isinstance(map_B, pd.DataFrame):
        # Case where mono_A is predicted and already in uniprot numbering
        # Ensure uniprot used in holo mapping corresponds to predicted (could be chimeric)
        if isinstance(mono_A.pinder_id, str):
            mono_A_uniprot = mono_A.pinder_id.split("__")[1].split("-")[0]
            map_B = map_B.query(f"uniprot_acc == '{mono_A_uniprot}'").reset_index(
                drop=True
            )
        mask_A_resolved = set(map_B.resi_uniprot.astype(int)).intersection(
            set(mono_A.atom_array.res_id)
        )
        mask_A = sorted(list(mask_A_resolved))
        map_B = map_B[map_B["resi_uniprot"].isin(mask_A)].reset_index(drop=True)
        mask_B_resolved = set(mono_B.atom_array.res_id).intersection(
            set(map_B.resi.astype(int))
        )
        mask_B = sorted(list(mask_B_resolved))

    if not (mask_A and mask_B):
        # Could happen if different domains are crystallized
        # Or if our mapping is incorrect
        log.error(
            "no common residues found! " f"{mono_A.pinder_id}--{mono_B.pinder_id}"
        )
        return mono_A, mono_B

    assert len(mask_A) == len(mask_B)

    mono_A_common = mono_A.filter("res_id", mask_A)
    mono_B_common = mono_B.filter("res_id", mask_B)
    assert isinstance(mono_A_common, Structure)
    assert isinstance(mono_B_common, Structure)
    return mono_A_common, mono_B_common


@lru_cache(maxsize=1)
def canonical_atom_type_mask(atom_tys: tuple[str]) -> "torch.Tensor":
    """Canonical atom masks for 21 (20 standard + missing) residue types"""
    import torch

    msk = torch.zeros(len(pc.INDEX_TO_AA_THREE), len(atom_tys))
    for aa_idx in range(len(pc.INDEX_TO_AA_THREE)):
        aa = pc.INDEX_TO_AA_THREE[aa_idx]
        for atom_idx, atom in enumerate(atom_tys):
            if atom in pc.BB_ATOMS:
                msk[aa_idx, atom_idx] = 1
            elif atom in pc.AA_TO_SC_ATOMS[aa]:
                msk[aa_idx, atom_idx] = 1
    atom_mask: torch.Tensor = msk.bool()
    return atom_mask


@lru_cache(maxsize=1)
def backbone_atom_tensor(atom_tys: tuple[str]) -> "torch.Tensor":
    import torch

    atoms: torch.Tensor = torch.tensor(
        [i for i, at in enumerate(atom_tys) if at in pc.BB_ATOMS]
    ).long()
    return atoms


def _superimpose_common_atoms(
    fixed: AtomArray, mobile: AtomArray, max_iterations: int = 10
) -> tuple[
    AtomArray,
    struc.AffineTransformation,
    NDArray[np.int_],
    NDArray[np.int_],
]:
    """Try to superimpose two structures based on homology.
    If this fails due to a lack of common anchors (e.g. in case of very short peptides),
    fall back to superimposing corresponding atoms with common atom annotations.
    """
    try:
        return_value: tuple[
            AtomArray,
            struc.AffineTransformation,
            NDArray[np.int_],
            NDArray[np.int_],
        ] = superimpose_chain(fixed, mobile, max_iterations=max_iterations)
        return return_value
    except ValueError as error:
        # Only run fallback if number of anchors is the issue
        if "anchor" not in str(error):
            raise
        fixed_common_mask = struc.filter_intersection(fixed, mobile)
        fixed_coord = fixed.coord[fixed_common_mask]
        mobile_common_mask = struc.filter_intersection(mobile, fixed)
        mobile_coord = mobile.coord[mobile_common_mask]
        _, transformation = struc.superimpose(fixed_coord, mobile_coord)
        mobile_superimposed = transformation.apply(mobile)
        return (
            mobile_superimposed,
            transformation,
            np.where(fixed_common_mask)[0],
            np.where(mobile_common_mask)[0],
        )
