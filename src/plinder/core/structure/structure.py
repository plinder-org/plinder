from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Optional

import biotite.structure as struc
import numpy as np
import pandas as pd
from biotite import TextFile
from biotite.structure.atoms import AtomArray
from numpy.typing import NDArray
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from plinder.core.structure import surgery
from plinder.core.structure.atoms import (
    atom_array_from_cif_file,
    generate_conformer,
    get_per_chain_seq_alignments,
    get_per_chain_seq_to_structure_alignments,
    get_seq_aligned_structures,
    invert_chain_seq_map,
    match_ligands,
    resn2seq,
    write_pdb,
)
from plinder.core.structure.superimpose import superimpose_chain
from plinder.core.utils import constants as pc
from plinder.core.utils.config import get_config
from plinder.core.utils.dataclass import stringify_dataclass

# TODO: Decide whether to lift these from pinder or import them
from plinder.core.utils.log import setup_logger

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


class Structure(BaseModel):
    id: str
    protein_path: Path
    protein_sequence: Path
    list_ligand_sdf_and_resolved_smiles: Optional[list[tuple[Path, str]]] = None
    protein_atom_array: Optional[AtomArray] = None
    renumbered_protein_array_and_seqs: Optional[
        tuple[AtomArray, dict[str, str], dict[str, list[int]]]
    ] = None
    ligand_mols: Optional[
        dict[
            str,
            tuple[Chem.rdchem.Mol, Chem.rdchem.Mol, Chem.rdchem.Mol, tuple[int, ...]],
        ]
    ] = None
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
    list_ligand_sdf_and_resolved_smiles : Optional[list[tuple[Path, str]]]
        SDF path and ligand smiles. This is optional for apo and pred structures
    protein_atom_array : AtomArray | None = None
        Protein Biotite atom array
    renumbered_protein_array_and_seqs: Optional[
        tuple[AtomArray, dict[str, str], dict[str, list[int]]]] = None
        Protein array, dictionary of aligned sequences and dictionary of mapped indices
    ligand_mols : ligand_mols: Optional[
        dict[
            str,
            tuple[Chem.rdchem.Mol, Chem.rdchem.Mol, Chem.rdchem.Mol, tuple[int, ...]],
        ]
    ]
        Dictionary of tuple of unresolved ligand mol,
        resolved aligned mol (via rdkit ConstrainedEmbed), resolved mol randon conformer and tuple of matched sbstructure indices
    add_ligand_hydrogen : bool = False
        Whether to add hydrogen to ligand or not
    structure_type : str = "holo"
        Structure type, "holo", "apo" or "pred"
    """

    class Config:
        arbitrary_types_allowed = True

    @property
    def resolved_sequences(self):
        resolved_sequences = {}
        with open(self.protein_sequence) as f_seq:
            lines = f_seq.readlines()
            for idx, line in enumerate(lines):
                if ">" in line:
                    chain = line[1:].strip()
                    seq = lines[idx + 1].strip().replace("\n", "")
                    # canonicalize
                    seq = seq.translate(str.maketrans(pc.non_canonical_aa))

                    resolved_sequences[chain] = seq

        return resolved_sequences

    @classmethod
    def load_structure(
        cls,
        id: str,
        protein_path: Path,
        protein_sequence: Path,
        list_ligand_sdf_and_resolved_smiles: Optional[list[tuple[Path, str]]] = None,
        protein_atom_array: Optional[AtomArray] = None,
        renumbered_protein_array_and_seqs: Optional[
            tuple[AtomArray, dict[str, str], dict[str, list[int]]]
        ] = None,
        ligand_mols: Optional[
            dict[
                str,
                tuple[
                    Chem.rdchem.Mol, Chem.rdchem.Mol, Chem.rdchem.Mol, tuple[int, ...]
                ],
            ]
        ] = None,
        add_ligand_hydrogen: bool = False,
        structure_type: str = "holo",
    ) -> Structure | None:
        structure = cls(
            id=id,
            protein_path=protein_path,
            protein_sequence=protein_sequence,
            list_ligand_sdf_and_resolved_smiles=list_ligand_sdf_and_resolved_smiles,
            protein_atom_array=protein_atom_array,
            renumbered_protein_array_and_seqs=renumbered_protein_array_and_seqs,
            ligand_mols=ligand_mols,
            add_ligand_hydrogen=add_ligand_hydrogen,
            structure_type=structure_type,
        )

        if (structure.protein_atom_array is None) & (structure.ligand_mols is None):
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

            if structure_type == "holo":
                structure.renumbered_protein_array_and_seqs = (
                    get_per_chain_seq_to_structure_alignments(
                        structure.resolved_sequences, protein_arr
                    )
                )
            else:
                structure.renumbered_protein_array_and_seqs = (protein_arr, {}, {})

            structure.ligand_mols = {}
            if not (list_ligand_sdf_and_resolved_smiles is None):
                for ligand_sdf, resolved_smiles in list_ligand_sdf_and_resolved_smiles:
                    # Match molecules
                    (
                        matches,
                        resolved_ligand_mol,
                        original_unresolved_mol,
                    ) = match_ligands(
                        resolved_smiles,
                        ligand_sdf,
                        add_hydrogen=add_ligand_hydrogen,
                    )

                    resolved_ligand_mol_conformer = generate_conformer(
                        resolved_ligand_mol
                    )
                    structure.ligand_mols[ligand_sdf.stem] = (
                        original_unresolved_mol,
                        resolved_ligand_mol,
                        resolved_ligand_mol_conformer,
                        matches,
                    )
        elif structure.protein_atom_array is not None:
            if structure_type == "holo":
                structure.renumbered_protein_array_and_seqs = (
                    get_per_chain_seq_to_structure_alignments(
                        structure.resolved_sequences, protein_arr
                    )
                )
            else:
                structure.renumbered_protein_array_and_seqs = (protein_arr, {}, {})

        structure.protein_atom_array = structure.renumbered_protein_array_and_seqs[0]
        try:
            getattr(structure.protein_atom_array, "b-factor")
        except AttributeError:
            b_factors: NDArray[np.double | np.str_] = np.repeat(
                0.0, structure.protein_atom_array.shape[0]
            )
            structure.protein_atom_array.set_annotation("b_factor", b_factors)

        return structure

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
            id=combined_id,
            protein_path=combined_path,
            protein_sequence=self.protein_sequence,
            list_ligand_sdf_and_resolved_smiles=self.list_ligand_sdf_and_resolved_smiles,
            protein_atom_array=combined_arr,
            renumbered_protein_array_and_seqs=self.renumbered_protein_array_and_seqs,
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
                id=self.id,
                protein_path=self.protein_path,
                protein_sequence=self.protein_sequence,
                list_ligand_sdf_and_resolved_smiles=self.list_ligand_sdf_and_resolved_smiles,
                protein_atom_array=arr,
                renumbered_protein_array_and_seqs=self.renumbered_protein_array_and_seqs,
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
                id=self.id,
                protein_path=self.protein_path,
                protein_sequence=self.protein_sequence,
                list_ligand_sdf_and_resolved_smiles=self.list_ligand_sdf_and_resolved_smiles,
                protein_atom_array=target_at,
                renumbered_protein_array_and_seqs=self.renumbered_protein_array_and_seqs,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogen=self.add_ligand_hydrogen,
                structure_type=self.structure_type,
            )

            other_struct = Structure(
                id=other.id,
                protein_path=other.protein_path,
                protein_sequence=self.protein_sequence,
                list_ligand_sdf_and_resolved_smiles=other.list_ligand_sdf_and_resolved_smiles,
                protein_atom_array=ref_at,
                renumbered_protein_array_and_seqs=self.renumbered_protein_array_and_seqs,
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
                id=self.id,
                protein_path=self.protein_path,
                protein_sequence=self.protein_sequence,
                list_ligand_sdf_and_resolved_smiles=self.list_ligand_sdf_and_resolved_smiles,
                protein_atom_array=superimposed,
                renumbered_protein_array_and_seqs=self.renumbered_protein_array_and_seqs,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogen=self.add_ligand_hydrogen,
                structure_type=self.structure_type,
            ),
            raw_rmsd,
            refined_rmsd,
        )

    @property
    def aligned_unresolved_seqs(self):
        return self.renumbered_protein_array_and_seqs[1]

    @property
    def unresolved_aligned_indices(self):
        return self.renumbered_protein_array_and_seqs[2]

    @property
    def protein_coords(self) -> NDArray[np.double]:
        """ndarray[np.double]: The coordinates of the protein atoms in the structure."""
        protein_coords: NDArray[np.double] = self.protein_atom_array.coord
        return protein_coords

    @property
    def original_unresolved_mols(self) -> dict[str, Chem.rdchem.Mol]:
        """Original ligand mol objects."""
        return {tag: mol_tuple[0] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def resolved_ligand_mols(self) -> dict[str, Chem.rdchem.Mol]:
        """Original ligand mol objects."""
        return {tag: mol_tuple[1] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def resolved_ligand_conformers(self) -> dict[str, Chem.rdchem.Mol]:
        """Original ligand mol objects."""
        return {tag: mol_tuple[2] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def resolved_ligand_conformers_coords(self) -> dict[str, NDArray[np.double]]:
        """ndarray[np.double]: The coordinates of the reinitialized resolved ligand atoms in the structure."""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.resolved_ligand_conformers.items()
        }
        return ligand_coords

    @property
    def resolved_ligand_mols_coords(self) -> dict[str, NDArray[np.double]]:
        """ndarray[np.double]: The coordinates of the aligned resolved ligand atoms in the structure."""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.resolved_ligand_mols.items()
        }
        return ligand_coords

    @property
    def original_unresolved_mols_coords(self) -> dict[str, NDArray[np.double]]:
        """ """
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.original_unresolved_mols.items()
        }
        return ligand_coords

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
    def protein_chain_unresolved_sequence(self) -> dict[str, list[str]]:
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
    def unresolved_protein_sequence(self) -> str:
        """str: The amino acid sequence of the structure."""
        numbering, resn = struc.get_residues(self.protein_atom_array)
        seq: str = resn2seq(resn)
        return seq

    @property
    def unresolved_protein_fasta(self) -> str:
        """str: The fasta representation of the structure sequence."""
        fasta_str: str = "\n".join(
            [f">{self.protein_path.stem}", self.unresolved_protein_sequence]
        )
        return fasta_str

    @property
    def protein_tokenized_sequence(self) -> "torch.Tensor":
        """torch.Tensor: The tokenized sequence representation of the structure sequence."""
        import torch

        seq_encoding = torch.tensor(
            [pc.AA_TO_INDEX[x] for x in self.unresolved_protein_sequence]
        )
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

    @classmethod
    def get_properties(cls):
        return [name for name in dir(cls) if isinstance(getattr(cls, name), property)]


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


def _align_monomers_with_mask(
    monomer1: Structure,
    monomer2: Structure,
    remove_differing_atoms: bool = True,
    renumber_residues: bool = False,
    remove_differing_annotations: bool = False,
) -> tuple[Structure, Structure]:
    monomer2, monomer1 = monomer2.align_common_sequence(
        monomer1,
        remove_differing_atoms=remove_differing_atoms,
        renumber_residues=renumber_residues,
        remove_differing_annotations=remove_differing_annotations,
    )
    monomer2, _, _ = monomer2.superimpose(monomer1)
    return monomer1, monomer2


def _align_structures_with_mask(
    multi_chain_structure: Structure,
    map_of_alt_monomer_structures: dict[str, Structure],
    remove_differing_atoms: bool = True,
    renumber_residues: bool = False,
    remove_differing_annotations: bool = False,
) -> tuple[Structure, Structure]:
    cropped_and_superposed_ref_dict = {}
    cropped_and_superposed_target_dict = {}
    for ref_chain, target_monomer_structure in map_of_alt_monomer_structures.items():
        ref_monomer_structure = multi_chain_structure.filter(
            property="chain_id", mask=[ref_chain]
        )
        target_monomer_structure.set_chain(ref_chain)

        ref_monomer_structure, target_monomer_structure = _align_monomers_with_mask(
            ref_monomer_structure,
            target_monomer_structure,
            remove_differing_atoms=remove_differing_atoms,
            renumber_residues=renumber_residues,
            remove_differing_annotations=remove_differing_annotations,
        )
        cropped_and_superposed_ref_dict[ref_chain] = ref_monomer_structure
        cropped_and_superposed_target_dict[ref_chain] = target_monomer_structure
    target_values = list(cropped_and_superposed_target_dict.values())
    ref_values = list(cropped_and_superposed_ref_dict.values())
    new_merged_target = target_values[0]
    new_merged_ref = ref_values[0]
    for target_val in target_values[1:]:
        new_merged_target += target_val
    for ref_val in ref_values[1:]:
        new_merged_ref += ref_val
    # Transfer ligand
    new_merged_target.ligand_mols = multi_chain_structure.ligand_mols
    new_merged_ref.ligand_mols = multi_chain_structure.ligand_mols

    # merge monomers (and/or part of multi_chain_structure) into single new structure
    reference_chains = set(multi_chain_structure.protein_atom_array.chain_id)
    chains_not_mapped = reference_chains - set(map_of_alt_monomer_structures.keys())
    if len(chains_not_mapped) == 0:
        return new_merged_ref, new_merged_target
    else:
        for chain in chains_not_mapped:
            unmapped_structure = multi_chain_structure.filter(
                property="chain_id", mask=chain
            )
            new_merged_target += unmapped_structure
            new_merged_ref += unmapped_structure

    return cropped_and_superposed_ref_dict, cropped_and_superposed_target_dict
