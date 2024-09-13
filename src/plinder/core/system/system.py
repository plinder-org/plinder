# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np
from numpy.typing import NDArray
import pandas as pd
from rdkit import Chem
import torch
from torch import Tensor

if TYPE_CHECKING:
    from ost import mol


from plinder.core.index import utils
from plinder.core.scores.links import query_links
from plinder.core.utils.cpl import get_plinder_path
from plinder.core.utils.io import (
    download_alphafold_cif_file,
    download_pdb_chain_cif_file,
)
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import get_zips_to_unpack
from plinder.core.structure.structure import Structure
from plinder.core.utils import constants as pc
from plinder.core.loader.ligand_featurizer import (
    lig_atom_featurizer,
    generate_conformer,
)

LOG = setup_logger(__name__)

RES_IDX_PAD_VALUE = -99
COORDS_PAD_VALUE = -100
ATOM_TYPE_PAD_VALUE = -1


def structure2tensor(
    protein_atom_coordinates: NDArray[np.double] | None = None,
    protein_atom_types: NDArray[np.str_] | None = None,
    protein_residue_coordinates: NDArray[np.double] | None = None,
    protein_residue_ids: NDArray[np.int_] | None = None,
    protein_residue_types: NDArray[np.str_] | None = None,
    ligand_mols: Chem.rdchem.Mol = None,
    dtype: torch.dtype = torch.float32,
    use_ligand_conformer: bool = True,
) -> dict[str, torch.Tensor]:
    property_dict = {}
    if protein_atom_types is not None:
        types_array_ele = np.zeros(
            (len(protein_atom_types), len(set(list(pc.ELE2NUM.values()))))
        )
        for i, name in enumerate(protein_atom_types):
            types_array_ele[i, pc.ELE2NUM.get(name, "C")] = 1.0

        property_dict["protein_atom_types"] = torch.tensor(types_array_ele).type(dtype)

    if protein_residue_types is not None:
        unknown_name_idx = max(pc.AA_TO_INDEX.values()) + 1
        types_array_res = np.zeros((len(protein_residue_types), 1))
        for i, name in enumerate(protein_residue_types):
            types_array_res[i] = pc.AA_TO_INDEX.get(name, unknown_name_idx)
        property_dict["protein_residue_types"] = torch.tensor(types_array_res).type(
            dtype
        )
    if protein_atom_coordinates is not None:
        property_dict["protein_atom_coordinates"] = torch.tensor(
            protein_atom_coordinates, dtype=dtype
        )
    if protein_residue_coordinates is not None:
        property_dict["protein_residue_coordinates"] = torch.tensor(
            protein_residue_coordinates, dtype=dtype
        )
    if protein_residue_ids is not None:
        property_dict["protein_residue_ids"] = torch.tensor(
            protein_residue_ids, dtype=dtype
        )
    if ligand_mols is not None:
        property_dict["ligand_features"] = torch.tensor(
            [lig_atom_featurizer(ligand_mol) for ligand_mol in ligand_mols]
        )
        if use_ligand_conformer:
            ligand_mols = [generate_conformer(ligand_mol) for ligand_mol in ligand_mols]
        lig_coords = torch.tensor(
            [ligand_mol.GetConformer().GetPositions() for ligand_mol in ligand_mols]
        ).float()

        property_dict["ligand_atom_coordinates"] = lig_coords

    return property_dict


def pad_to_max_length(
    mat: Tensor,
    max_length: int | Sequence[int] | Tensor,
    dims: Sequence[int],
    value: int | float | None = None,
) -> Tensor:
    """Takes a tensor and pads it to maximum length with right padding on the specified dimensions.

    Parameters:
        mat (Tensor): The tensor to pad. Can be of any shape
        max_length (int | Sequence[int] | Tensor): The size of the tensor along specified dimensions after padding.
        dims (Sequence[int]): The dimensions to pad. Must have the same number of elements as `max_length`.
        value (int, optional): The value to pad with, by default None

    Returns:
        Tensor : The padded tensor. Below are examples of input and output shapes
            Example 1:
                input: (2, 3, 4), max_length: 5, dims: [0, 2]
                output: (5, 3, 5)
            Example 2:
                input: (2, 3, 4), max_length: 5, dims: [0]
                output: (5, 3, 4)
            Example 3:
                input: (2, 3, 4), max_length: [5, 7], dims: [0, 2]
                output: (5, 3, 7)

    """
    if not isinstance(max_length, int):
        assert len(dims) == len(max_length)

    num_dims = len(mat.shape)
    pad_idxs = [(num_dims - i) * 2 - 1 for i in dims]

    if isinstance(max_length, int):
        pad_sizes = [
            max_length - mat.shape[int(-(i + 1) / 2 + num_dims)] if i in pad_idxs else 0
            for i in range(num_dims * 2)
        ]
    else:
        max_length_list = (
            list(max_length) if not isinstance(max_length, list) else max_length
        )
        pad_sizes = [
            (
                max_length_list[int(-(i + 1) / 2 + num_dims)]
                - mat.shape[int(-(i + 1) / 2 + num_dims)]
                if i in pad_idxs
                else 0
            )
            for i in range(num_dims * 2)
        ]

    return torch.nn.functional.pad(input=mat, pad=tuple(pad_sizes), value=value)


def pad_and_stack(
    tensors: list[Tensor],
    dim: int = 0,
    dims_to_pad: list[int] | None = None,
    value: int | float | None = None,
) -> Tensor:
    """Pads a list of tensors to the maximum length observed along each dimension and then stacks them along a new dimension (given by `dim`).

    Parameters:
        tensors (list[Tensor]): A list of tensors to pad and stack
        dim (int): The new dimension to stack along.
        dims_to_pad (list[int] | None): The dimensions to pad
        value (int | float | None, optional): The value to pad with, by default None

    Returns:
        Tensor: The padded and stacked tensor. Below are examples of input and output shapes
            Example 1: Sequence features (although redundant with torch.rnn.utils.pad_sequence)
                input: [(2,), (7,)], dim: 0
                output: (2, 7)
            Example 2: Pair features (e.g., pairwise coordinates)
                input: [(4, 4, 3), (7, 7, 3)], dim: 0
                output: (2, 7, 7, 3)

    """
    assert (
        len({t.ndim for t in tensors}) == 1
    ), f"All `tensors` must have the same number of dimensions."

    # Pad all dims if none are specified
    if dims_to_pad is None:
        dims_to_pad = list(range(tensors[0].ndim))

    # Find the max length of the dims_to_pad
    shapes = torch.tensor([t.shape for t in tensors])
    envelope = shapes.max(dim=0).values
    max_length = envelope[dims_to_pad]

    padded_matrices = [
        pad_to_max_length(t, max_length, dims_to_pad, value) for t in tensors
    ]
    return torch.stack(padded_matrices, dim=dim)


def collate_complex(
    structures: list[dict[str, Tensor]],
    coords_pad_value: int = COORDS_PAD_VALUE,
    atom_type_pad_value: int = ATOM_TYPE_PAD_VALUE,
    residue_id_pad_value: int = RES_IDX_PAD_VALUE,
) -> dict[str, Tensor]:
    protein_atom_types = []
    protein_residue_types = []
    protein_atom_coordinates = []
    protein_residue_coordinates = []
    protein_residue_ids = []
    ligand_features = []
    ligand_atom_coordinates = []

    for x in structures:
        protein_atom_types.append(x["protein_atom_types"])
        protein_residue_types.append(x["rprotein_esidue_types"])
        protein_atom_coordinates.append(x["protein_atom_coordinates"])
        protein_residue_coordinates.append(x["protein_residue_coordinates"])
        protein_residue_ids.append(x["protein_residue_ids"])
        ligand_features.append(x["ligand_features"])
        ligand_atom_coordinates.append(x["ligand_atom_coordinates"])
    return {
        "protein_atom_types": pad_and_stack(
            protein_atom_types, dim=0, value=atom_type_pad_value
        ),
        "protein_residue_types": pad_and_stack(
            protein_residue_types, dim=0, value=atom_type_pad_value
        ),
        "protein_atom_coordinates": pad_and_stack(
            protein_atom_coordinates, dim=0, value=coords_pad_value
        ),
        "protein_residue_coordinates": pad_and_stack(
            protein_residue_coordinates, dim=0, value=coords_pad_value
        ),
        "protein_residue_ids": pad_and_stack(
            protein_residue_ids, dim=0, value=residue_id_pad_value
        ),
        "ligand_features": pad_and_stack(
            ligand_features, dim=0, value=atom_type_pad_value
        ),
        "ligand_atom_coordinates": pad_and_stack(
            ligand_atom_coordinates, dim=0, value=coords_pad_value
        ),
    }


def collate_batch(
    batch: list[dict[str, dict[str, Tensor] | str]],
) -> dict[str, dict[str, Tensor] | list[str]]:
    """Collate a batch of PlinderDataset items into a merged mini-batch of Tensors.

    Used as the default collate_fn for the torch DataLoader consuming PlinderDataset.

    Parameters:
        batch (list[dict[str, dict[str, Tensor] | str]]): A list of dictionaries
        containing the data for each item in the batch.

    Returns:
        dict[str, dict[str, Tensor] | list[str]]: A dictionary containing
        the merged Tensors for the batch.

    """
    sample_ids: list[str] = []
    target_ids: list[str] = []
    target_structures: list[dict[str, Tensor]] = []
    feature_structures: list[dict[str, Tensor]] = []
    for x in batch:
        assert isinstance(x["id"], str)
        assert isinstance(x["input_id"], str)
        assert isinstance(x["target_id"], str)
        assert isinstance(x["target"], dict)
        assert isinstance(x["input_features"], dict)

        sample_ids.append(x["input_id"])
        target_ids.append(x["target_id"])
        target_structures.append(x["target"])
        feature_structures.append(x["input_feature"])

    collated_batch: dict[str, dict[str, Tensor] | list[str]] = {
        "target": collate_complex(target_structures),
        "input_features": collate_complex(feature_structures),
        "input_id": sample_ids,
        "target_id": target_ids,
    }
    return collated_batch


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
            property="chain_id", mask={ref_chain}
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
        self._best_linked_structures: pd.DataFrame | None = None
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
        for chain in self.ligands:
            ligand_views[chain] = io.LoadEntity(
                self.ligands[chain], format="sdf"
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
        return len(self.ligands)

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
    def best_linked_structures_paths(self) -> dict[str, str] | None:
        """
        Return single best apo and pred by sort score.

        Returns
        -------
        pd.DataFrame | None
            dataframe of linked structures if present in plinder
        """
        # TODO: This assumes single protein chain holo system.
        # Extend this to make it more general
        if self._best_linked_structures is None:
            links = query_links(filters=[("reference_system_id", "==", self.system_id)])
            best_apo = (
                links[(links.kind == "apo")]
                .sort_values(by="sort_score")
                .id.to_list()[0]
            )
            best_pred = (
                links[(links.kind == "pred")]
                .sort_values(by="sort_score", ascending=False)
                .id.to_list()[0]
            )
            chain_id = self.system_id.split("__")[2]
            self._best_linked_structures = {
                "apo": {chain_id: self.get_linked_structure("apo", best_apo)},
                "pred": {chain_id: self.get_linked_structure("pred", best_pred)},
            }
        return self._best_linked_structures

    def create_masked_bound_unbound_complexes(
        self,
        remove_differing_atoms: bool = True,
        renumber_residues: bool = False,
        remove_differing_annotations: bool = False,
    ) -> tuple[Structure, Structure, Structure]:
        """Create complexes for apo and predicted cropped to common holo substructures.

        The method applies a pairwise masking procedure which crops both unbound and bound structures
        such that they have equal numbers of residues and atom counts.

        Note: this method may result in very distorted holo (ground-truth) structures if
        the unbound monomer structures have little to no sequence and atoms in common.
        Unless you need all monomer types to be equal shapes, the `PinderSystem.create_complex` method
        or pure-superposition without masking (Structure.superimpose) is more appropriate.

        Parameters:
            monomer_types (Sequence[str]): The unbound monomer types to consider (apo, predicted, or both).
            remove_differing_atoms (bool):
                Whether to remove non-overlappings atoms that may still be present even after sequence-based alignment.
            renumber_residues (bool):
                Whether to renumber the residues in the receptor and ligand `Structure`'s to match numbering of the holo counterparts.
            remove_differing_annotations (bool):
                Whether to remove annotation categories (set to empty str or default value for the category type).
                This is useful if you need to perform biotite.structure.filter_intersection on the resulting structures.
                Note: this can have unintended side-effects like stripping the `element` attribute on structures.
                By default, the annotation categories are removed if they don't match in order to define the intersecting atom mask,
                after which the original structure annotations are preserved by applying the intersecting mask to the original AtomArray.
                Default is False.

        Returns:
            tuple[Structure, Structure, Structure]: A tuple of the cropped holo, apo, and predicted Structure objects, respectively.

        """
        holo_structure = self.holo_structure
        alt_structure = self.alt_structures
        apo_structure_map = alt_structure["apo"]
        pred_structure_map = alt_structure["pred"]

        if apo_structure_map is not None:
            # Ensure same number of chains in apo and holo
            holo_structure, apo_structure = _align_structures_with_mask(
                multi_chain_structure=holo_structure,
                map_of_alt_monomer_structures=apo_structure_map,
                remove_differing_atoms=remove_differing_atoms,
                renumber_residues=renumber_residues,
                remove_differing_annotations=remove_differing_annotations,
            )

        if pred_structure_map is not None:
            # Ensure same number of chains in pred and holo
            holo_structure, pred_structure = _align_structures_with_mask(
                multi_chain_structure=holo_structure,
                map_of_alt_monomer_structures=pred_structure_map,
                remove_differing_atoms=remove_differing_atoms,
                renumber_residues=renumber_residues,
                remove_differing_annotations=remove_differing_annotations,
            )

        return holo_structure, apo_structure, pred_structure

    @property
    def holo_structure(self) -> Structure | None:
        """
        Load holo structure
        """
        return Structure(protein_path=self.receptor_cif, id=self.system_id)

    @property
    def alt_structures(self) -> Structure | None:
        """
        Load apo structure
        """
        best_structures = self.best_linked_structures_paths
        return {
            kind: {
                chain: Structure(
                    protein_path=alt_path,
                    id=Path(alt_path).parent.name,
                    structure_type=kind,
                )
            }
            for kind, alts in best_structures.items()
            for chain, alt_path in alts.items()
        }
