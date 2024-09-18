from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Optional

import biotite.structure as struc
import numpy as np
from biotite import TextFile
from biotite.structure.atoms import AtomArray
from numpy.typing import NDArray
from pinder.core.loader.structure import _superimpose_common_atoms
from pinder.core.structure import surgery
from pinder.core.structure.atoms import (
    get_per_chain_seq_alignments,
    get_seq_aligned_structures,
    invert_chain_seq_map,
    resn2seq,
    write_pdb,
)
from pinder.core.utils.dataclass import stringify_dataclass
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from plinder.core.structure.atoms import (
    atom_array_from_cif_file,
    generate_input_conformer,
    get_ligand_atom_index_mapping_mask,
    get_residue_index_mapping_mask,
    make_atom_mask,
    make_one_hot_atom_features,
    match_ligands,
)
from plinder.core.utils import constants as pc
from plinder.core.utils.config import get_config

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
    list_ligand_sdf_and_input_smiles: Optional[list[tuple[Path, str]]] = None
    protein_atom_array: Optional[AtomArray] = None
    ligand_mols: Optional[
        dict[
            str,
            tuple[
                Chem.rdchem.Mol,
                Chem.rdchem.Mol,
                tuple[int, ...],
                Chem.rdchem.Mol,
                tuple[int, ...],
            ],
        ]
    ] = None
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
    list_ligand_sdf_and_input_smiles : Optional[list[tuple[Path, str]]]
        SDF path and ligand smiles. This is optional for apo and pred structures
    protein_atom_array : AtomArray | None = None
        Protein Biotite atom array
    ligand_mols : ligand_mols: Optional[
        dict[
            str,
            tuple[Chem.rdchem.Mol, Chem.rdchem.Mol, Chem.rdchem.Mol, tuple[int, ...]],
        ]
    ]
    # TODO: this is a WRONG description!!
        Dictionary of template molecule loaded from smiles (2D), template conformer,
        array of conformer-template matching atom ids, resolved mol, resolved atom matches
          template_mol,
                        template_mol_conformer,
                        _conformer_atom_matches,
                        resolved_mol,
                        _resolved_atom_matches,
    add_ligand_hydrogens : bool = False
        Whether to add hydrogen to ligand or not
    structure_type : str = "holo"
        Structure type, "holo", "apo" or "pred"
    """

    class Config:
        arbitrary_types_allowed = True

    @property
    def input_sequences(self) -> dict[str, str]:
        _resolvd_sequences = {}
        with open(self.protein_sequence) as f_seq:
            lines = f_seq.readlines()
            for idx, line in enumerate(lines):
                if ">" in line:
                    chain = line[1:].strip()
                    seq = lines[idx + 1].strip().replace("\n", "")
                    # canonicalize
                    seq = seq.translate(str.maketrans(pc.non_canonical_aa))

                    _resolvd_sequences[chain] = seq
        return _resolvd_sequences

    @classmethod
    def load_structure(
        cls,
        id: str,
        protein_path: Path,
        protein_sequence: Path,
        list_ligand_sdf_and_input_smiles: Optional[list[tuple[Path, str]]] = None,
        protein_atom_array: Optional[AtomArray] = None,
        ligand_mols: Optional[
            dict[
                str,
                tuple[
                    Chem.rdchem.Mol,
                    Chem.rdchem.Mol,
                    tuple[int, ...],
                    Chem.rdchem.Mol,
                    tuple[int, ...],
                ],
            ]
        ] = None,
        add_ligand_hydrogens: bool = False,
        structure_type: str = "holo",
    ) -> Structure | None:
        structure = cls(
            id=id,
            protein_path=protein_path,
            protein_sequence=protein_sequence,
            list_ligand_sdf_and_input_smiles=list_ligand_sdf_and_input_smiles,
            protein_atom_array=protein_atom_array,
            ligand_mols=ligand_mols,
            add_ligand_hydrogens=add_ligand_hydrogens,
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
            # Remove water
            structure.protein_atom_array = protein_arr[protein_arr.res_name != "HOH"]
            structure.ligand_mols = {}
            if not (list_ligand_sdf_and_input_smiles is None):
                for ligand_sdf, input_smiles in list_ligand_sdf_and_input_smiles:
                    # Match molecules
                    (
                        template_mol,
                        resolved_mol,
                        resolved_atom_matches,
                    ) = match_ligands(input_smiles, ligand_sdf)

                    # get input_conformer with matches
                    (
                        template_mol_conformer,
                        conformer_atom_matches,
                    ) = generate_input_conformer(
                        template_mol, addHs=add_ligand_hydrogens
                    )
                    # TODO: Double check with VO
                    # Extract the template indices that are matches
                    # _conformer_atom_matches = [i[1] for i in conformer_atom_matches[0]]
                    # _resolved_atom_matches = [i[1] for i in resolved_atom_matches[0]]
                    print("A", conformer_atom_matches)
                    _conformer_atom_matches = tuple(conformer_atom_matches[0].values())
                    _resolved_atom_matches = tuple(resolved_atom_matches[0].values())

                    structure.ligand_mols[ligand_sdf.stem] = (
                        template_mol,
                        template_mol_conformer,
                        _conformer_atom_matches,
                        resolved_mol,
                        _resolved_atom_matches,
                    )
                print(structure.ligand_mols)

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
            list_ligand_sdf_and_input_smiles=self.list_ligand_sdf_and_input_smiles,
            protein_atom_array=combined_arr,
            ligand_mols=self.ligand_mols,
            add_ligand_hydrogens=self.add_ligand_hydrogens,
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
                list_ligand_sdf_and_input_smiles=self.list_ligand_sdf_and_input_smiles,
                protein_atom_array=arr,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
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
                list_ligand_sdf_and_input_smiles=self.list_ligand_sdf_and_input_smiles,
                protein_atom_array=target_at,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
                structure_type=self.structure_type,
            )

            other_struct = Structure(
                id=other.id,
                protein_path=other.protein_path,
                protein_sequence=self.protein_sequence,
                list_ligand_sdf_and_input_smiles=other.list_ligand_sdf_and_input_smiles,
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
                list_ligand_sdf_and_input_smiles=self.list_ligand_sdf_and_input_smiles,
                protein_atom_array=superimposed,
                ligand_mols=self.ligand_mols,
                add_ligand_hydrogens=self.add_ligand_hydrogens,
                structure_type=self.structure_type,
            ),
            raw_rmsd,
            refined_rmsd,
        )

    @property
    def input_sequence_stacked_mask(self) -> NDArray[np.int]:
        """Input sequence stacked by chain"""
        seqres_masks = get_residue_index_mapping_mask(
            self.input_sequences, self.protein_atom_array
        )
        return [seqres_masks[ch] for ch in self.protein_chain_ordered]

    # TODO: review this!!!
    # @property
    # def resolved_ligand_mask(self):
    #     masks = []
    #     for ch in self.ligand_chain_ordered:
    #         mol = self.ligand_mols[ch][1]
    #         matching_idxs = self.ligand_mols[ch][-1]
    #         masks.append(get_ligand_atom_index_mapping_mask(mol, matching_idxs))
    #     return masks
    @property
    def resolved_smiles_ligand_mask(self) -> list[list[int]]:
        """List of protein chains ordered the way it is in structure."""
        masks = []
        for ch in self.ligand_chain_ordered:
            mol = self.ligand_mols[ch][1]
            matching_idxs = self.ligand_mols[ch][-1]
            masks.append(get_ligand_atom_index_mapping_mask(mol, matching_idxs))
        return masks

    @property
    def input_sequence_list_ordered_by_chain(self) -> list[str]:
        """List of protein chains ordered the way it is in structure."""
        return [self.input_sequences[ch] for ch in self.protein_chain_ordered]

    @property
    def protein_chain_ordered(self) -> list[list[int]]:
        """List of protein chains ordered the way it is in structure"""
        chain_order = []
        for ch_id in self.protein_atom_array.chain_id:
            if ch_id not in chain_order:
                chain_order.append(ch_id)
        return chain_order

    @property
    def input_sequence_full_atom_feat(self) -> list[list[int]]:
        """Resolved sequence full atom features."""
        seq_res_atom_list = [
            [
                pc.ORDERED_AA_FULL_ATOM[pc.ONE_TO_THREE[res]]
                for res in self.input_sequences[ch]
            ]
            for ch in self.protein_chain_ordered
        ]
        return [
            [make_one_hot_atom_features(res) for res in ch] for ch in seq_res_atom_list
        ]

    @property
    def sequence_atom_mask(self) -> list[list[int]]:
        """Sequence mask indicating which residues are resolved"""
        seqres_masks = get_residue_index_mapping_mask(
            self.input_sequences, self.protein_atom_array
        )
        return [
            make_atom_mask(
                self.protein_atom_array[self.protein_atom_array.chain_id == ch],
                self.input_sequences[ch],
                seqres_masks[ch],
            )
            for ch in self.protein_chain_ordered
        ]

    @property
    def ligand_chain_ordered(self) -> list[str]:
        """List of Ligand chains sorted in order."""
        return [
            path_and_smiles[0].stem
            for path_and_smiles in self.list_ligand_sdf_and_input_smiles
        ]

    @property
    def protein_coords(self) -> NDArray[np.double]:
        """ndarray[np.double]: The coordinates of the protein atoms in the structure."""
        protein_coords: NDArray[np.double] = self.protein_atom_array.coord
        return protein_coords

    @property
    def input_ligand_templates(self) -> dict[str, Chem.rdchem.Mol]:
        """Ligand 2D mol objects from input SMILES"""
        return {tag: mol_tuple[0] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def protein_calpha_coords(self) -> NDArray[np.double]:
        """ndarray[np.double]: The coordinates of the protein clapha atoms in the structure."""
        protein_calpha_coords: NDArray[np.double] = self.protein_atom_array[
            self.protein_atom_array.atom_name == "CA"
        ].coord
        return protein_calpha_coords

    @property
    def input_ligand_conformers(self) -> dict[str, Chem.rdchem.Mol]:
        """Ligand 3D mol conformer objects from input SMILES"""
        return {tag: mol_tuple[1] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def input_ligand_conformer_atom_index_maps(self) -> dict[str, list[dict[int, int]]]:
        """List of ligand conformer atom index maps to input SMILES"""
        return {tag: mol_tuple[2] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def input_ligand_conformer_masks(self) -> dict[str, list[tuple[int]]]:
        """
        Map matched template ligands to input ligands conformer
        and convert the indices to binary mask."""
        masks = []
        for ch in self.ligand_chain_ordered:
            mol = self.ligand_mols[ch][3]

            matching_idxs = self.ligand_mols[ch][2]
            masks.append(get_ligand_atom_index_mapping_mask(mol, matching_idxs))
        return masks

    @property
    def input_ligand_template_masks(self) -> dict[str, list[tuple[int]]]:
        """
        Map matched template ligands to input resolved ligands
        and convert the indices to binary mask."""
        masks = []
        for ch in self.ligand_chain_ordered:
            mol = self.ligand_mols[ch][0]

            matching_idxs = self.ligand_mols[ch][4]
            masks.append(get_ligand_atom_index_mapping_mask(mol, matching_idxs))
        return masks

    @property
    def input_ligand_template_coords(self) -> dict[str, NDArray[np.double]]:
        """ndarray[np.double]: The coordinates of the input 2D conformer generated from input SMILES"""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.input_ligand_templates.items()
        }
        return ligand_coords

    @property
    def input_ligand_conformer_coords(self) -> dict[str, NDArray[np.double]]:
        """ndarray[np.double]: The coordinates of the input 3D conformer generated from input SMILES"""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.input_ligand_conformers.items()
        }
        return ligand_coords

    @property
    def resolved_ligand_mols(self) -> dict[str, Chem.rdchem.Mol]:
        """Resolved holo ligand mol objects"""
        return {tag: mol_tuple[3] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def resolved_ligand_structure_atom_index_maps(
        self,
    ) -> dict[str, list[dict[int, int]]]:
        """List of resolved holo structure atom index maps to input SMILES"""
        return {tag: mol_tuple[4] for tag, mol_tuple in self.ligand_mols.items()}

    @property
    def resolved_ligand_structure_coords(self) -> dict[str, NDArray[np.double]]:
        """ndarray[np.double]: The coordinates of the holo resolved ligand atoms in the holo structure."""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.resolved_ligand_mols.items()
        }
        return ligand_coords

    @property
    def conformer2resolved_mapping(self) -> dict[str, list[dict[int, int]]]:
        """for every ligand it generates a list of dictionaries providing pairwise mapping between conformer and holo"""
        from itertools import product

        conf2holo_mappings = {}
        for (
            tag,
            smi2holo_maps,
        ) in self.resolved_ligand_structure_atom_index_maps.items():
            smi2conf_maps = self.input_ligand_conformer_atom_index_maps[tag]
            conf2holo_mappings[tag] = [
                {v: smi2holo[k] for k, v in smi2conf.items() if k in smi2holo}
                for smi2conf, smi2holo in product(smi2conf_maps, smi2holo_maps)
            ]
        return conf2holo_mappings

    # TODO: complete
    # @property
    # def resolved_to_conformer_masks(self) -> dict[str, list[tuple[int]]]:
    #     masks = []
    #     # map the two
    #     resolved_ligand_structure_masks
    #     input_ligand_conformer_masks
    #     return masks
    @property
    def resolved_ligand_mols_coords(self) -> dict[str, NDArray[np.double]]:
        """Holo PDB ligand coordinate"""
        ligand_coords: NDArray[np.double] = {
            tag: mol.GetConformer().GetPositions()
            for tag, mol in self.resolved_ligand_mols.items()
        }
        return ligand_coords

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
    def protein_chain_uninput_sequence(self) -> dict[str, list[str]]:
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
    def protein_sequence_from_structure(self) -> str:
        """str: The amino acid sequence of the structure."""
        numbering, resn = struc.get_residues(self.protein_atom_array)
        seq: str = resn2seq(resn)
        return seq

    @property
    def protein_structure_sequence_fasta(self) -> str:
        """str: The fasta representation of the structure sequence."""
        fasta_str: str = "\n".join(
            [f">{self.protein_path.stem}", self.protein_sequence_from_structure]
        )
        return fasta_str

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
    def protein_structure_residue_names(self) -> list[str]:
        """list[str]: The list of distinct residue names in the structure."""
        res_list = self._attr_from_atom_array(
            self.protein_atom_array, "res_name", distinct=True, sort=True
        )
        return [str(r) for r in res_list]

    @property
    def protein_structure_residues(self) -> list[int]:
        """list[int]: The list of distinct residue IDs in the structure."""
        res_list = self._attr_from_atom_array(
            self.protein_atom_array, "res_id", distinct=True, sort=True
        )
        return [int(r) for r in res_list]

    @property
    def protein_structure_atom_names(self) -> list[str]:
        """list[str]: The list of distinct atom names in the structure."""
        at_list = self._attr_from_atom_array(
            self.protein_atom_array, "atom_name", distinct=True, sort=True
        )
        return [str(a) for a in at_list]

    @property
    def protein_structure_b_factor(self) -> list[float]:
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


def _align_monomers_with_mask(
    monomer1: Structure,
    monomer2: Structure,
    remove_differing_atoms: bool = True,
    renumber_residues: bool = False,
    remove_differing_annotations: bool = False,
) -> tuple[Structure, Structure]:
    """Align and superpose monomers.

    Parameters
    ----------
    monomer1: Structure
        Monomer structure 1
    monomer2: Structure
        Monomer structure 2
    remove_differing_atoms: bool = True
        Remove differing atoms between monomers
    renumber_residues: bool = False
        renumber of residues

    Returns
    -------
    tuple[Structure, Structure]
        Align and superposed structures

    """
    if (monomer1 is None) | (monomer2 is None):
        return monomer1, monomer2
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
) -> tuple[dict[str, Structure], dict[str, Structure]]:
    """Align and crop two structures.

    Parameters
    ----------
    multi_chain_structure : Structure
        Reference Structure to align to.
        This could be multi-chain structure (could be holo structure)
    map_of_alt_monomer_structures : dict[str, Structure]
        Dictionary of chains of target structures to align (could be chains of apo/pred)
    remove_differing_atoms : bool = True
        Remove atoms different between reference and target
    renumber_residues : bool = False
        Renumber residue
    remove_differing_annotations : bool = False
        Remove different annotations from AtomArray

    Returns
    -------
    tuple[dict[str, Structure], dict[str, Structure]]
        dictionary of chains of aligned and cropped reference and target chains

    """
    cropped_and_superposed_ref_dict = {}
    cropped_and_superposed_target_dict = {}
    for ref_chain, target_monomer_structure in map_of_alt_monomer_structures.items():
        ref_monomer_structure = multi_chain_structure.filter(
            property="chain_id", mask=[ref_chain]
        )
        if target_monomer_structure is not None:
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
    if new_merged_target is not None:
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
            if new_merged_target is not None:
                new_merged_target += unmapped_structure
            new_merged_ref += unmapped_structure

    return cropped_and_superposed_ref_dict, cropped_and_superposed_target_dict
