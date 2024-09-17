# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np
import pandas as pd
import torch
from biotite.structure.atoms import AtomArray
from numpy.typing import NDArray
from pinder.core.loader.dataset import pad_and_stack
from rdkit import Chem
from torch import Tensor

if TYPE_CHECKING:
    from ost import mol


from plinder.core.index import utils
from plinder.core.scores.links import query_links
from plinder.core.structure.diffdock_utils import lig_atom_featurizer
from plinder.core.structure.structure import Structure, _align_structures_with_mask
from plinder.core.utils import constants as pc
from plinder.core.utils.cpl import get_plinder_path
from plinder.core.utils.io import (
    download_alphafold_cif_file,
    download_pdb_chain_cif_file,
)
from plinder.core.utils.log import setup_logger
from plinder.core.utils.unpack import get_zips_to_unpack

LOG = setup_logger(__name__)

RES_IDX_PAD_VALUE = -99
COORDS_PAD_VALUE = -100
ATOM_TYPE_PAD_VALUE = -1


def stack_atom_array_features(atom_arr, atom_arr_feat, chain_order_list):
    return [
        getattr(atom_arr[atom_arr.chain_id == chain], atom_arr_feat)
        for chain in chain_order_list
    ]


def stack_ligand_feat(feat_dict, chain_order_list):
    return [feat_dict[chain] for chain in chain_order_list]


def structure2tensor(
    protein_atom_array: AtomArray,
    resolved_sequence_residue_types: NDArray[np.str_] | None = None,
    resolved_sequence_mask: NDArray[np.int_] | None = None,
    resolved_ligand_mols: dict[str, Chem.rdchem.Mol] | None = None,
    resolved_ligand_mols_coords: dict[str, NDArray[np.float_]] | None = None,
    resolved_ligand_conformers_coords: dict[str, NDArray[np.float_]] = None,
    resolved_ligand_mask: list[NDArray[np.int_]] | None = None,
    protein_chain_in_order: list[str] = None,
    ligand_chain_in_order: list[str] = None,
    coords_pad_value: int = COORDS_PAD_VALUE,
    atom_type_pad_value: int = ATOM_TYPE_PAD_VALUE,
    residue_id_pad_value: int = RES_IDX_PAD_VALUE,
    dtype: torch.dtype = torch.float32,
) -> dict[str, torch.Tensor]:
    protein_atom_coordinates = [
        torch.tensor(coord)
        for coord in stack_atom_array_features(
            protein_atom_array, "coord", protein_chain_in_order
        )
    ]
    protein_atom_types = stack_atom_array_features(
        protein_atom_array, "element", protein_chain_in_order
    )
    protein_residue_coordinates = [
        torch.tensor(coord)
        for coord in stack_atom_array_features(
            protein_atom_array, "coord", protein_chain_in_order
        )
    ]
    protein_residue_ids = stack_atom_array_features(
        protein_atom_array, "res_id", protein_chain_in_order
    )
    protein_residue_types = stack_atom_array_features(
        protein_atom_array, "res_name", protein_chain_in_order
    )
    property_dict = {}

    if protein_atom_types is not None:
        types_array_ele = []
        for chain_atom_types in protein_atom_types:
            types_array_ele_by_chain = np.zeros(
                (len(chain_atom_types), len(set(list(pc.ELE2NUM.values()))))
            )
            for res_index, name in enumerate(chain_atom_types):
                types_array_ele_by_chain[res_index, pc.ELE2NUM.get(name, "C")] = 1.0
            types_array_ele.append(torch.tensor(types_array_ele_by_chain).type(dtype))
        property_dict["protein_atom_types"] = pad_and_stack(
            types_array_ele, dim=0, value=atom_type_pad_value
        )

    if protein_residue_types is not None:
        types_array_res = []
        unknown_name_idx = max(pc.AA_TO_INDEX.values()) + 1
        for chain_residue_types in protein_residue_types:
            types_array_res_by_chain = np.zeros((len(chain_residue_types), 1))

            for res_index, name in enumerate(chain_atom_types):
                types_array_res_by_chain[res_index] = pc.AA_TO_INDEX.get(
                    name, unknown_name_idx
                )
            types_array_res.append(torch.tensor(types_array_res_by_chain).type(dtype))
        property_dict["protein_residue_types"] = pad_and_stack(
            types_array_res, dim=0, value=atom_type_pad_value
        )
    if resolved_sequence_residue_types is not None:
        types_array_res_resolved = []
        unknown_name_idx = max(pc.AA_TO_INDEX.values()) + 1
        for chain_residue_types_resolved in resolved_sequence_residue_types:
            types_array_res_by_chain_resolved = np.zeros(
                (len(chain_residue_types_resolved), 1)
            )
            chain_residue_types_resolved = [
                pc.ONE_TO_THREE.get(x) for x in chain_residue_types_resolved
            ]
            for res_index, name in enumerate(chain_residue_types_resolved):
                types_array_res_by_chain_resolved[res_index] = pc.AA_TO_INDEX.get(
                    name, unknown_name_idx
                )
            types_array_res_resolved.append(
                torch.tensor(types_array_res_by_chain_resolved).type(dtype)
            )

        property_dict["resolved_protein_residue_types"] = pad_and_stack(
            types_array_res_resolved, dim=0, value=atom_type_pad_value
        )

    if protein_atom_coordinates is not None:
        property_dict["protein_atom_coordinates"] = pad_and_stack(
            protein_atom_coordinates, dim=0, value=coords_pad_value
        )
    if protein_residue_coordinates is not None:
        property_dict["protein_residue_coordinates"] = pad_and_stack(
            protein_atom_coordinates, dim=0, value=coords_pad_value
        )

    if protein_residue_ids is not None:
        property_dict["protein_residue_ids"] = pad_and_stack(
            [torch.tensor(res_id) for res_id in protein_residue_ids],
            dim=0,
            value=residue_id_pad_value,
        )

    if resolved_sequence_mask is not None:
        property_dict["resolved_sequence_mask"] = pad_and_stack(
            [torch.tensor(mask) for mask in resolved_sequence_mask],
            dim=0,
            value=residue_id_pad_value,
        )

    if resolved_ligand_mols_coords is not None:
        lig_coords = pad_and_stack(
            [
                torch.tensor(coord)
                for coord in stack_ligand_feat(
                    resolved_ligand_mols_coords, ligand_chain_in_order
                )
            ],
            dim=0,
            value=coords_pad_value,
        )
        property_dict["resolved_ligand_mols_coords"] = lig_coords
    if resolved_ligand_conformers_coords is not None:
        lig_coords = pad_and_stack(
            [
                torch.tensor(coord)
                for coord in stack_ligand_feat(
                    resolved_ligand_conformers_coords, ligand_chain_in_order
                )
            ],
            dim=0,
            value=coords_pad_value,
        )
        property_dict["ligand_conformer_atom_coordinates"] = lig_coords
    if resolved_ligand_mols is not None:
        ligand_feat = {
            ch: lig_atom_featurizer(ligand_mol)
            for ch, ligand_mol in resolved_ligand_mols.items()
        }

        property_dict["ligand_features"] = pad_and_stack(
            [
                torch.tensor(feat)
                for feat in stack_ligand_feat(ligand_feat, ligand_chain_in_order)
            ],
            dim=0,
            value=coords_pad_value,
        )

    if resolved_ligand_mask is not None:
        lig_masks = pad_and_stack(
            [torch.tensor(mask) for mask in resolved_ligand_mask],
            dim=0,
            value=coords_pad_value,
        )
        property_dict["resolved_ligand_mask"] = lig_masks

    return property_dict


def collate_complex(
    structures: list[dict[str, Tensor]],
    coords_pad_value: int = COORDS_PAD_VALUE,
    atom_type_pad_value: int = ATOM_TYPE_PAD_VALUE,
    residue_id_pad_value: int = RES_IDX_PAD_VALUE,
) -> dict[str, Tensor]:
    protein_atom_types = []
    protein_residue_types = []
    resolved_protein_residue_types = []
    protein_atom_coordinates = []
    protein_residue_coordinates = []
    protein_residue_ids = []
    resolved_sequence_mask = []
    ligand_features = []
    ligand_conformer_atom_coordinates = []
    resolved_ligand_mols_coords = []
    resolved_ligand_mask = []

    for x in structures:
        protein_atom_types.append(x["protein_atom_types"])
        protein_residue_types.append(x["protein_residue_types"])
        resolved_protein_residue_types.append(x["resolved_protein_residue_types"])
        protein_atom_coordinates.append(x["protein_atom_coordinates"])
        protein_residue_coordinates.append(x["protein_residue_coordinates"])
        protein_residue_ids.append(x["protein_residue_ids"])
        resolved_sequence_mask.append(x["resolved_sequence_mask"])
        ligand_features.append(x["ligand_features"])
        ligand_conformer_atom_coordinates.append(x["ligand_conformer_atom_coordinates"])
        resolved_ligand_mols_coords.append(x["resolved_ligand_mols_coords"])
        resolved_ligand_mask.append(x["resolved_ligand_mask"])

    return {
        "protein_atom_types": pad_and_stack(
            protein_atom_types, dim=0, value=atom_type_pad_value
        ),
        "protein_residue_types": pad_and_stack(
            protein_residue_types, dim=0, value=atom_type_pad_value
        ),
        "resolved_protein_residue_type": pad_and_stack(
            resolved_protein_residue_types, dim=0, value=coords_pad_value
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
        "resolved_sequence_mask": pad_and_stack(
            resolved_sequence_mask, dim=0, value=residue_id_pad_value
        ),
        "ligand_features": pad_and_stack(
            ligand_features, dim=0, value=atom_type_pad_value
        ),
        "ligand_conformer_atom_coordinates": pad_and_stack(
            ligand_conformer_atom_coordinates, dim=0, value=coords_pad_value
        ),
        "resolved_ligand_mols_coords": pad_and_stack(
            resolved_ligand_mols_coords, dim=0, value=coords_pad_value
        ),
        "resolved_ligand_mask": pad_and_stack(
            resolved_ligand_mask, dim=0, value=coords_pad_value
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
        assert isinstance(x["input_id"], str)
        assert isinstance(x["target_id"], str)
        assert isinstance(x["target"], dict)
        assert isinstance(x["input_features"], dict)

        sample_ids.append(x["input_id"])
        target_ids.append(x["target_id"])
        target_structures.append(x["target"])
        feature_structures.append(x["input_features"])

    collated_batch: dict[str, dict[str, Tensor] | list[str]] = {
        "target": collate_complex(target_structures),
        "input_features": collate_complex(feature_structures),
        "input_id": sample_ids,
        "target_id": target_ids,
    }
    return collated_batch


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
        resolved_smiles_dict: dict[str, str],
        prune: bool = True,
    ) -> None:
        self.system_id: str = system_id
        self.resolved_smiles_dict: dict[str, str] = resolved_smiles_dict
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
                links[(links.kind == "apo")].sort_values(by="sort_score").id.to_list()
            )
            if len(best_apo) > 0:
                best_apo = best_apo[0]
            else:
                best_apo = None
            best_pred = (
                links[(links.kind == "pred")]
                .sort_values(by="sort_score", ascending=False)
                .id.to_list()
            )
            if len(best_pred) > 0:
                best_pred = best_pred[0]
            else:
                best_pred = None
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
        """ """
        holo_structure = self.holo_structure
        alt_structure = self.alt_structures
        apo_structure_map = alt_structure.get("apo")
        pred_structure_map = alt_structure.get("pred")
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
        list_ligand_sdf_and_resolved_smiles = [
            (Path(sdf_path), self.resolved_smiles_dict[chain])
            for chain, sdf_path in self.ligands.items()
        ]
        return Structure.load_structure(
            id=self.system_id,
            protein_path=self.receptor_cif,
            protein_sequence=Path(self.sequences),
            list_ligand_sdf_and_resolved_smiles=list_ligand_sdf_and_resolved_smiles,
        )

    @property
    def alt_structures(self) -> Structure | None:
        """
        Load apo structure
        """
        best_structures = self.best_linked_structures_paths
        alt_structure_dict = {}
        for kind, alts in best_structures.items():
            alt_chain = {}
            for chain, alt_path in alts.items():
                try:
                    alt_chain[chain] = Structure.load_structure(
                        id=Path(alt_path).parent.name,
                        protein_path=alt_path,
                        protein_sequence=Path(self.sequences),
                        structure_type=kind,
                    )
                except Exception:
                    alt_chain[chain] = None
            alt_structure_dict[kind] = alt_chain
        return alt_structure_dict
