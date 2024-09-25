from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Iterable

import biotite.structure as struc
import numpy as np
from biotite.structure.atoms import AtomArray
from biotite.structure.io.mol import SDFile
from numpy.typing import NDArray
from pydantic import BaseModel, model_validator
from rdkit import Chem

from plinder.core.structure import surgery
from plinder.core.structure.atoms import (
    _stack_atom_array_features,
    atom_array_from_cif_file,
    generate_input_conformer,
    get_residue_index_mapping_mask,
    make_atom_mask,
    match_ligands,
)
from plinder.core.structure.superimpose import superimpose_chain
from plinder.core.structure.vendored import (
    get_per_chain_seq_alignments,
    get_seq_aligned_structures,
    invert_chain_seq_map,
    resn2seq,
    write_pdb,
)
from plinder.core.utils import constants as pc
from plinder.core.utils.dataclass import stringify_dataclass
from plinder.core.utils.log import setup_logger

if TYPE_CHECKING:
    import torch


log = setup_logger(__name__)


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


def reverse_dict(mapping: dict[int, int]) -> dict[int, int]:
    return {v: k for k, v in mapping.items()}


def get_rdkit_mol(ligand: Path | str) -> Chem.rdchem.Mol:
    if isinstance(ligand, Path):
        return Chem.MolFromSmiles(ligand)
    elif isinstance(ligand, Path):
        return next(Chem.SDMolSupplier(ligand))


class Structure(BaseModel):
    id: str
    protein_path: Path
    protein_sequence: dict[str, str] | None = None
    ligand_sdfs: dict[str, str] | None = None
    ligand_smiles: dict[str, str] | None = None
    protein_atom_array: AtomArray | None = None
    ligand_mols: (
        dict[
            str,
            tuple[
                Chem.Mol,
                Chem.Mol,
                # tuple[NDArray, NDArray],
                Chem.Mol,
                # tuple[NDArray, NDArray],
                tuple[NDArray, NDArray],
            ],
        ]
        | None
    ) = None
    add_ligand_hydrogens: bool = False
    structure_type: str = "holo"

    """Initialize structure.
    This dataclass provides abstraction over plinder systems holo, apo and predicted structures.
    It loads and processes receptor proteins and ligands in a way that allows us to capture the
    pecularitties of these molecule types.
    It's also intentionally design to be robust to capture structures, with or without ligand

    Parameters
    ----------
    id : str
        PLINDER system or apo/pred id. If the structure is holo, this is used to load the holo structure
    protein_path : Path
        Protein structure cif path, could be holo, apo or pred path
    protein_sequence : dict[str, str] | None = None
        Protein sequence dictionary with chain id as key and sequence as value
        Set from protein_path if not provided
    ligand_sdfs : dict[str, str] | None = None
        Dictionary of ligand sdf file paths with chain id as key and sdf file path as value
    ligand_smiles : dict[str, str] | None = None
        Dictionary of ligand smiles with chain id as key and smiles as value
    protein_atom_array : AtomArray | None = None
        Protein Biotite atom array
    ligand_mols : ligand_mols: Optional[
        dict[
            str,
            tuple[Chem.Mol, Chem.Mol, Chem.Mol, tuple[NDArray, NDArray]]
        ]
    ]
        Dictionary of ligand molecule id to
            molecule loaded from smiles (2D),
            template (random) conformer generated from 2D,
            resolved (holo) mol conformer,
            paired stacked arrays (template vs holo) mapping atom order by index
    add_ligand_hydrogens : bool = False
        Whether to add hydrogen to ligand or not
    structure_type : str = "holo"
        Structure type, "holo", "apo" or "pred"
    """

    class Config:
        arbitrary_types_allowed = True

    @model_validator(mode="after")
    def initialize(self) -> Structure:
        if self.protein_atom_array is None:
            self.load_protein()
        if self.protein_sequence is None:
            self.load_sequence()
        if self.ligand_sdfs is not None or self.ligand_smiles is not None:
            self.load_ligands()
        return self

    def load_protein(self) -> None:
        protein_arr = atom_array_from_cif_file(
            self.protein_path, use_author_fields=False
        )
        if protein_arr is None:
            raise ValueError("Protein atom array could not be loaded")
        annotation_arr: NDArray[np.double | np.str_] = np.repeat(
            0.0, protein_arr.shape[0]
        )
        protein_arr.set_annotation("b_factor", annotation_arr)
        if "" in set(protein_arr.element):
            try:
                protein_arr.element = struc.infer_elements(protein_arr)
            except Exception:
                log.debug(
                    f"Found missing elements in {self.protein_path} but failed to infer"
                )

        # Remove water
        self.protein_atom_array = protein_arr[protein_arr.res_name != "HOH"]

        # Set b-factor again
        try:
            getattr(self.protein_atom_array, "b-factor")
        except AttributeError:
            b_factors: NDArray[np.double | np.str_] = np.repeat(
                0.0, self.protein_atom_array.shape[0]
            )
            self.protein_atom_array.set_annotation("b_factor", b_factors)

    def load_sequence(self) -> None:
        if self.protein_atom_array is None:
            raise ValueError("Protein atom array not loaded")
        self.protein_sequence = {}
        for chain in self.protein_chain_ordered:
            self.protein_sequence[chain] = struc.to_sequence(
                self.protein_atom_array[self.protein_atom_array.chain_id == chain]
            )
        if not len(self.protein_sequence):
            raise ValueError("Protein sequence could not be loaded")

    def load_ligands(self) -> None:
        if self.ligand_sdfs is None:
            raise ValueError("Ligand SDFs not provided")
        if self.ligand_smiles is None:
            raise ValueError("Ligand SMILES not provided")
        if set(self.ligand_smiles.keys()) != set(self.ligand_sdfs.keys()):
            raise ValueError("Ligand SMILES and SDFs keys do not match")
        self.ligand_mols = {}
        for name, sdf in self.ligand_sdfs.items():
            # Match molecules
            (
                template_mol,
                resolved_mol,
                atoms_order_stacks,
            ) = match_ligands(self.ligand_smiles[name], sdf)

            # get input_conformer with matches
            (template_mol_conformer) = generate_input_conformer(
                template_mol, addHs=self.add_ligand_hydrogens
            )

            self.ligand_mols[name] = (
                template_mol,
                template_mol_conformer,
                resolved_mol,
                atoms_order_stacks,
            )

    def __repr__(self) -> str:
        class_str: str = stringify_dataclass(self, 4)
        return class_str

    def __add__(self, other: Structure) -> Structure | None:
        if (self.protein_atom_array is not None) and (
            other.protein_atom_array is not None
        ):
            combined_arr = self.protein_atom_array + other.protein_atom_array
            combined_id = f"{self.id}--{other.id}"
            combined_path = self.protein_path.parent / combined_id
            structure_type = "_".join(
                sorted(list({self.structure_type, other.structure_type}))
            )
            return Structure(
                id=combined_id,
                protein_path=combined_path,
                protein_sequence=self.protein_sequence,
                ligand_sdfs=self.ligand_sdfs,
                ligand_smiles=self.ligand_smiles,
                protein_atom_array=combined_arr,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
                structure_type=structure_type,
            )
        else:
            return None

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
            if self.protein_atom_array is not None:
                arr = self.protein_atom_array[atom_mask].copy()

            return Structure(
                id=self.id,
                protein_path=self.protein_path,
                protein_sequence=self.protein_sequence,
                ligand_sdfs=self.ligand_sdfs,
                ligand_smiles=self.ligand_smiles,
                protein_atom_array=arr,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
                structure_type=self.structure_type,
            )
        assert self.protein_atom_array is not None
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
        assert (other.protein_atom_array is not None) and (
            self.protein_atom_array is not None
        )
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
                id=self.id,
                protein_path=self.protein_path,
                protein_sequence=self.protein_sequence,
                ligand_sdfs=self.ligand_sdfs,
                ligand_smiles=self.ligand_smiles,
                protein_atom_array=target_at,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
                structure_type=self.structure_type,
            )

            other_struct = Structure(
                id=other.id,
                protein_path=other.protein_path,
                protein_sequence=self.protein_sequence,
                ligand_sdfs=other.ligand_sdfs,
                ligand_smiles=other.ligand_smiles,
                protein_atom_array=ref_at,
                ligand_mols=other.ligand_mols,
                add_ligand_hydrogens=other.add_ligand_hydrogens,
                structure_type=other.structure_type,
            )

            return self_struct, other_struct
        other.protein_atom_array = ref_at
        self.protein_atom_array = target_at
        return self, other

    def set_chain(self, chain_id: str) -> None:
        if self.protein_atom_array is not None:
            self.protein_atom_array.chain_id = np.repeat(
                chain_id, self.protein_atom_array.shape[0]
            )

    def superimpose(
        self,
        other: Structure,
        strip_ligands: bool = True,
    ) -> tuple[Structure, float, float]:
        """transforms self to superimpose on other"""

        if self.protein_atom_array is None:
            raise ValueError("both structures must have protein atom arrays")
        if other.protein_atom_array is None:
            raise ValueError("both structures must have protein atom arrays")

        # max_iterations=1 -> no outlier removal
        # first structure - fixed (other), second - mobile (self)
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
                id=self.id,
                protein_path=self.protein_path,
                protein_sequence=self.protein_sequence,
                ligand_sdfs=None if strip_ligands else self.ligand_sdfs,
                ligand_smiles=None if strip_ligands else self.ligand_smiles,
                protein_atom_array=superimposed,
                ligand_mols=None if strip_ligands else self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
                structure_type=self.structure_type,
            ),
            raw_rmsd,
            refined_rmsd,
        )

    @property
    def input_sequence_residue_mask_stacked(self) -> list[list[int]]:
        """Input sequence stacked by chain"""
        # TODO: do we want to keep this as assertion?
        # better if then raise?
        assert self.protein_atom_array is not None
        assert self.protein_sequence is not None

        seqres_masks = get_residue_index_mapping_mask(
            self.protein_sequence, self.protein_atom_array
        )
        return [seqres_masks[ch] for ch in self.protein_chain_ordered]

    @property
    def input_sequence_list_ordered_by_chain(self) -> list[str] | None:
        """List of protein chains ordered the way it is in structure."""
        if self.protein_chain_ordered and self.protein_sequence:
            return [self.protein_sequence[ch] for ch in self.protein_chain_ordered]
        else:
            return []

    @property
    def protein_chain_ordered(self) -> list[str]:
        """List of protein chains ordered the way it is in structure"""
        assert self.protein_atom_array is not None
        chain_order = []
        for ch_id in self.protein_atom_array.chain_id:
            if ch_id not in chain_order:
                chain_order.append(ch_id)
        return chain_order

    @property
    def input_sequence_full_atom_feat(self) -> dict[str, list[str]]:
        """Resolved sequence full atom features."""
        # TODO: do we want to keep this as assertion?
        # better if then raise?
        assert self.protein_sequence is not None
        assert self.protein_chain_ordered is not None

        seq_res_atom_dict = {
            ch: [
                atm
                for res in self.protein_sequence[ch]
                for atm in pc.ORDERED_AA_FULL_ATOM[pc.ONE_TO_THREE[res]]
            ]
            for ch in self.protein_chain_ordered
        }
        return seq_res_atom_dict

    @property
    def sequence_atom_mask(self) -> list[list[int]]:
        """Sequence mask indicating which residues are resolved"""
        # TODO: do we want to keep this as assertion?
        # better if then raise?
        assert self.protein_atom_array is not None
        assert self.protein_sequence is not None

        seqres_masks = get_residue_index_mapping_mask(
            self.protein_sequence, self.protein_atom_array
        )
        return [
            make_atom_mask(
                self.protein_atom_array[self.protein_atom_array.chain_id == ch],
                self.protein_sequence[ch],
                seqres_masks[ch],
            )
            for ch in self.protein_chain_ordered
        ]

    @property
    def ligand_chain_ordered(self) -> list[str]:
        """List of Ligand chains sorted in order."""
        if self.ligand_mols is not None:
            return sorted(self.ligand_mols.keys())
        else:
            return []

    @property
    def protein_coords(self) -> list[NDArray]:
        """list[NDArray]: The coordinates of the protein atoms in the structure."""
        assert self.protein_atom_array is not None

        protein_coords: list[NDArray] = [
            coord
            for coord in _stack_atom_array_features(
                self.protein_atom_array, "coord", self.protein_chain_ordered
            )
        ]

        return protein_coords

    @property
    def input_ligand_conformer_atom_array(self) -> dict[str, AtomArray]:
        """dict[str, AtomArray]: The coordinates of the input 3D conformer generated from input SMILES"""
        ligands = {}
        for c in self.input_ligand_conformers:
            with tempfile.NamedTemporaryFile(suffix=".sdf") as tmp_file:
                Chem.SDWriter(tmp_file.name).write(self.input_ligand_conformers[c])
                ligands[c] = SDFile.read(tmp_file.name).record.get_structure()
        return ligands

    @property
    def ligand_atom_array(self) -> dict[str, AtomArray]:
        """dict[str, AtomArray]: The coordinates of the input 3D conformer generated from input SMILES"""
        if self.ligand_sdfs is None:
            return {}
        ligands = {}
        for c in self.ligand_sdfs:
            ligands[c] = SDFile.read(self.ligand_sdfs[c]).record.get_structure()
        return ligands

    @property
    def input_ligand_templates(self) -> dict[str, Chem.Mol]:
        """Ligand 2D mol objects from input SMILES"""
        assert self.ligand_mols is not None
        return {tag: mol_tuple[0] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def protein_calpha_coords(self) -> NDArray[np.double]:
        assert self.protein_atom_array is not None
        """list[NDArray]: The coordinates of the protein clapha atoms in the structure."""
        protein_calpha_coords: list[NDArray] = [
            coord
            for coord in _stack_atom_array_features(
                self.protein_atom_array[self.protein_atom_array.atom_name == "CA"],
                "coord",
                self.protein_chain_ordered,
            )
        ]
        return protein_calpha_coords

    @property
    def input_ligand_conformers(self) -> dict[str, Chem.Mol]:
        """Ligand 3D mol conformer objects from input SMILES"""
        assert self.ligand_mols is not None
        return {tag: mol_tuple[1] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def input_ligand_conformers_coords(self) -> dict[str, NDArray[np.double]]:
        """dict[NDArray]: The coordinates of the input 3D conformer generated from input SMILES"""
        ligand_coords: dict[str, NDArray[np.double]] = {}
        for tag, mol in self.input_ligand_conformers.items():
            try:
                ligand_coords[tag] = mol.GetConformer().GetPositions()
            except ValueError:
                ligand_coords[tag] = np.array(len(mol.GetAtoms()) * [[0.0, 0.0, 0.0]])
        return ligand_coords

    @property
    def resolved_ligand_mols(self) -> dict[str, Chem.Mol]:
        """Resolved holo ligand mol objects"""
        assert self.ligand_mols is not None
        return {tag: mol_tuple[2] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def ligand_template2resolved_atom_order_stacks(
        self,
    ) -> dict[str, tuple[NDArray, NDArray]]:
        """for every ligand this gets a pair of atom order array stacks providing index sort to match template atoms to holo conformer atoms"""
        return (
            {tag: mol_tuple[3] for tag, mol_tuple in self.ligand_mols.items()}
            if self.ligand_mols
            else {}
        )

    @property
    def resolved_ligand_mols_coords(self) -> dict[str, NDArray[np.double]]:
        """Dictionary of ligands and their resolved atom coordinates in the holo structure."""
        ligand_coords: dict[str, NDArray[np.double]] = {}
        if self.ligand_mols is not None:
            ligand_coords = {
                tag: mol.GetConformer().GetPositions()
                for tag, mol in self.resolved_ligand_mols.items()
            }
        return ligand_coords

    @property
    def protein_backbone_mask(self) -> NDArray[np.bool_]:
        """ndarray[np.bool\_]: a logical mask for backbone atoms."""
        assert self.protein_atom_array is not None
        mask: NDArray[np.bool_] = struc.filter_peptide_backbone(self.protein_atom_array)
        return mask

    @property
    def protein_calpha_mask(self) -> NDArray[np.bool_]:
        """ndarray[np.bool\_]: a logical mask for alpha carbon atoms."""
        assert self.protein_atom_array is not None
        mask: NDArray[np.bool_] = self.protein_atom_array.atom_name == "CA"
        return mask

    @property
    def protein_n_atoms(self) -> int:
        """int: The number of atoms in the structure."""
        assert self.protein_atom_array is not None
        n: int = self.protein_atom_array.shape[0]
        return n

    @property
    def protein_chains(self) -> list[str]:
        """list[str]: The list of chain IDs in the structure."""
        if self.protein_atom_array is not None:
            ch_list = self._attr_from_atom_array(
                self.protein_atom_array, "chain_id", distinct=True, sort=True
            )
            return [str(ch) for ch in ch_list]
        else:
            return []

    @property
    def protein_sequence_from_structure(self) -> str:
        """str: The amino acid sequence of the structure."""
        assert self.protein_atom_array is not None
        numbering, resn = struc.get_residues(self.protein_atom_array)
        seq: str = resn2seq(resn)
        return seq

    @property
    def protein_structure_tokenized_sequence(self) -> "torch.Tensor":
        """torch.Tensor: The tokenized sequence representation of the structure sequence."""
        import torch

        seq_encoding = torch.tensor(
            [pc.AA_TO_INDEX[x] for x in self.protein_sequence_from_structure]
        )
        tokenized: torch.Tensor = seq_encoding.long()
        return tokenized

    @property
    def protein_unique_residue_names(self) -> list[str]:
        """list[str]: The list of distinct residue names in the structure."""
        assert self.protein_atom_array is not None
        res_list = self._attr_from_atom_array(
            self.protein_atom_array, "res_name", distinct=True, sort=True
        )
        return [str(r) for r in res_list]

    @property
    def protein_unique_residue_ids(self) -> list[int]:
        """list[int]: The list of distinct residue IDs in the structure."""
        assert self.protein_atom_array is not None
        res_list = self._attr_from_atom_array(
            self.protein_atom_array, "res_id", distinct=True, sort=True
        )
        return [int(r) for r in res_list]

    @property
    def protein_unique_atom_names(self) -> list[str]:
        """list[str]: The list of distinct atom names in the structure."""
        assert self.protein_atom_array is not None
        at_list = self._attr_from_atom_array(
            self.protein_atom_array, "atom_name", distinct=True, sort=True
        )
        return [str(a) for a in at_list]

    @property
    def protein_structure_b_factor(self) -> list[float]:
        """list[float]: A list of B-factor values for each atom in the structure."""
        assert self.protein_atom_array is not None
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

    @classmethod
    def get_properties(cls) -> list[str]:
        return [name for name in dir(cls) if isinstance(getattr(cls, name), property)]
