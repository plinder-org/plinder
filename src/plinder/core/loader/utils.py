from typing import Any, Sequence

import numpy as np
import pandas as pd
import torch
from numpy.typing import NDArray
from rdkit import Chem
from torch import Tensor

from plinder.core.structure.atoms import _AtomArrayOrStack
from plinder.core.structure.diffdock_utils import lig_atom_featurizer
from plinder.core.structure.structure import Structure
from plinder.core.utils import constants as pc

RES_IDX_PAD_VALUE = -99
COORDS_PAD_VALUE = -100
ATOM_TYPE_PAD_VALUE = -1

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
            max_length_list[int(-(i + 1) / 2 + num_dims)]
            - mat.shape[int(-(i + 1) / 2 + num_dims)]
            if i in pad_idxs
            else 0
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


def stack_atom_array_features(
    atom_arr: _AtomArrayOrStack,
    atom_arr_feat: str,
    chain_order_list: list[str] | None,
) -> list[NDArray[np.int_ | np.str_ | np.float_]]:
    assert chain_order_list is not None
    return [
        getattr(atom_arr[atom_arr.chain_id == chain], atom_arr_feat)
        for chain in chain_order_list
    ]


def stack_ligand_feat(
    feat_dict: dict[str, Any], chain_order_list: list[str] | None
) -> list[list[list[int]]]:
    assert chain_order_list is not None
    return [feat_dict[chain] for chain in chain_order_list]


def structure2tensor(
    protein_atom_array: _AtomArrayOrStack,
    input_sequence_residue_types: NDArray[np.str_] | None = None,
    input_sequence_mask: NDArray[np.int_] | None = None,
    input_sequence_full_atom_feat: NDArray[np.int_] | None = None,
    resolved_ligand_mols_coords: dict[str, NDArray[np.float_]] | None = None,
    input_ligand_conformer_coords: dict[str, NDArray[np.float_]] | None = None,
    input_ligand_conformers: dict[str, Chem.rdchem.Mol] | None = None,
    ligand_conformer2resolved_mask: dict[str, NDArray[np.int_]] | None = None,
    protein_chain_ordered: list[str] | None = None,
    ligand_chain_ordered: list[str] | None = None,
    atom_type_pad_value: int = ATOM_TYPE_PAD_VALUE,
    residue_id_pad_value: int = RES_IDX_PAD_VALUE,
    coords_pad_value: int = COORDS_PAD_VALUE,
    dtype: torch.dtype = torch.float32,
) -> dict[str, torch.Tensor]:
    protein_atom_coordinates = [
        torch.tensor(coord)
        for coord in stack_atom_array_features(
            protein_atom_array, "coord", protein_chain_ordered
        )
    ]
    protein_atom_types = stack_atom_array_features(
        protein_atom_array, "element", protein_chain_ordered
    )
    protein_residue_coordinates = [
        torch.tensor(coord)
        for coord in stack_atom_array_features(
            protein_atom_array, "coord", protein_chain_ordered
        )
    ]
    protein_calpha_coordinates = [
        torch.tensor(coord)
        for coord in stack_atom_array_features(
            protein_atom_array[protein_atom_array.atom_name == "CA"],
            "coord",
            protein_chain_ordered,
        )
    ]
    protein_residue_ids = stack_atom_array_features(
        protein_atom_array, "res_id", protein_chain_ordered
    )
    protein_residue_types = stack_atom_array_features(
        protein_atom_array, "res_name", protein_chain_ordered
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

            for res_index, name in enumerate(chain_residue_types):
                types_array_res_by_chain[res_index] = pc.AA_TO_INDEX.get(
                    name, unknown_name_idx
                )
            types_array_res.append(torch.tensor(types_array_res_by_chain).type(dtype))
        property_dict["protein_residue_types"] = pad_and_stack(
            types_array_res, dim=0, value=atom_type_pad_value
        )
    if input_sequence_residue_types is not None:
        types_array_res_resolved = []
        unknown_name_idx = max(pc.AA_TO_INDEX.values()) + 1
        for chain_residue_types_resolved in input_sequence_residue_types:
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
    if input_sequence_full_atom_feat is not None:
        property_dict["input_sequence_full_atom_feat"] = pad_and_stack(
            [torch.tensor(feat) for feat in input_sequence_full_atom_feat],
            dim=0,
            value=residue_id_pad_value,
        )

    if protein_atom_coordinates is not None:
        property_dict["protein_atom_coordinates"] = pad_and_stack(
            protein_atom_coordinates, dim=0, value=coords_pad_value
        )
    if protein_residue_coordinates is not None:
        property_dict["protein_residue_coordinates"] = pad_and_stack(
            protein_atom_coordinates, dim=0, value=coords_pad_value
        )

    if protein_calpha_coordinates is not None:
        property_dict["protein_calpha_coordinates"] = pad_and_stack(
            protein_calpha_coordinates, dim=0, value=coords_pad_value
        )

    if protein_residue_ids is not None:
        property_dict["protein_residue_ids"] = pad_and_stack(
            [torch.tensor(res_id) for res_id in protein_residue_ids],
            dim=0,
            value=residue_id_pad_value,
        )

    if input_sequence_mask is not None:
        property_dict["input_sequence_mask"] = pad_and_stack(
            [torch.tensor(mask) for mask in input_sequence_mask],
            dim=0,
            value=residue_id_pad_value,
        )

    if resolved_ligand_mols_coords is not None:
        lig_coords = pad_and_stack(
            [
                torch.tensor(coord)
                for coord in stack_ligand_feat(
                    resolved_ligand_mols_coords, ligand_chain_ordered
                )
            ],
            dim=0,
            value=coords_pad_value,
        )
        property_dict["resolved_ligand_mols_coords"] = lig_coords

    if input_ligand_conformer_coords is not None:
        lig_coords = pad_and_stack(
            [
                torch.tensor(coord)
                for coord in stack_ligand_feat(
                    input_ligand_conformer_coords, ligand_chain_ordered
                )
            ],
            dim=0,
            value=coords_pad_value,
        )
        property_dict["ligand_conformer_atom_coordinates"] = lig_coords

    if input_ligand_conformers is not None:
        ligand_feat = {
            ch: lig_atom_featurizer(ligand_mol)
            for ch, ligand_mol in input_ligand_conformers.items()
        }

        property_dict["ligand_features"] = pad_and_stack(
            [
                torch.tensor(feat)
                for feat in stack_ligand_feat(ligand_feat, ligand_chain_ordered)
            ],
            dim=0,
            value=coords_pad_value,
        )
    if ligand_conformer2resolved_mask is not None:
        ligand_mask = {ch: mask for ch, mask in ligand_conformer2resolved_mask.items()}

        property_dict["ligand_conformer2resolved_mask"] = pad_and_stack(
            [
                torch.tensor(m)
                for m in stack_ligand_feat(ligand_mask, ligand_chain_ordered)
            ],
            dim=0,
            value=coords_pad_value,
        )

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
    input_sequence_mask = []

    ligand_features = []
    ligand_conformer_atom_coordinates = []
    resolved_ligand_mols_coords = []
    input_sequence_full_atom_feat = []
    protein_calpha_coordinates = []
    ligand_conformer2resolved_mask = []

    for x in structures:
        protein_atom_types.append(x["protein_atom_types"])
        protein_residue_types.append(x["protein_residue_types"])
        resolved_protein_residue_types.append(x["resolved_protein_residue_types"])
        protein_atom_coordinates.append(x["protein_atom_coordinates"])
        protein_residue_coordinates.append(x["protein_residue_coordinates"])
        protein_residue_ids.append(x["protein_residue_ids"])
        input_sequence_mask.append(x["input_sequence_mask"])
        ligand_features.append(x["ligand_features"])
        ligand_conformer_atom_coordinates.append(x["ligand_conformer_atom_coordinates"])
        resolved_ligand_mols_coords.append(x["resolved_ligand_mols_coords"])
        ligand_conformer2resolved_mask.append(x["ligand_conformer2resolved_mask"])
        input_sequence_full_atom_feat.append(x["input_sequence_full_atom_feat"])
        protein_calpha_coordinates.append(x["protein_calpha_coordinates"])

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
        "input_sequence_mask": pad_and_stack(
            input_sequence_mask, dim=0, value=residue_id_pad_value
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
        "input_sequence_full_atom_feat": pad_and_stack(
            input_sequence_full_atom_feat, dim=0, value=coords_pad_value
        ),
        "protein_calpha_coordinates": pad_and_stack(
            protein_calpha_coordinates, dim=0, value=coords_pad_value
        ),
        "ligand_conformer2resolved_mask": pad_and_stack(
            ligand_conformer2resolved_mask, dim=0, value=coords_pad_value
        ),
    }


def collate_batch(
    batch: list[
        dict[
            str,
            list[str]
            | list[Structure]
            | list[dict[str, dict[str, Structure]]]
            | list[Path]
            | list[dict[str, Tensor]],
        ]
    ],
) -> dict[
    str,
    list[str]
    | list[Structure]
    | list[dict[str, dict[str, Structure]]]
    | list[Path]
    | list[dict[str, Tensor]],
]:
    """Collate a batch of PlinderDataset items into a merged mini-batch of Tensors.

    Used as the default collate_fn for the torch DataLoader consuming PlinderDataset.

    Parameters:
        batch (list[dict[str, dict[str, Tensor] | str]]): A list of dictionaries
        containing the data for each item in the batch.

    Returns:
        dict[str, dict[str, Tensor] | list[str]]: A dictionary containing
        the merged Tensors for the batch.

    """
    system_ids: list[str] = []
    holo_structures: list[Structure] = []
    alternate_structures: list[dict[str, dict[str, Structure]]] = []
    feature_and_coords: list[dict[str, Tensor]] = []
    paths: list[Path] = []
    for x in batch:
        assert isinstance(x["system_id"], str)
        assert isinstance(x["holo_structure"], Structure)
        assert isinstance(x["alternate_structures"], dict)
        assert isinstance(x["features_and_coords"], dict)
        assert isinstance(x["path"], str)

        system_ids.append(x["system_id"])
        holo_structures.append(x["holo_structure"])
        alternate_structures.append(x["alternate_structures"])
        feature_and_coords.append(x["features_and_coords"])
        paths.append(x["path"])

    collated_batch: dict[
        str,
        list[str]
        | list[Structure]
        | list[dict[str, dict[str, Structure]]]
        | list[Path]
        | list[dict[str, Tensor]],
    ] = {
        "system_ids": system_ids,
        "holo_structures": holo_structures,
        "alternate_structures": alternate_structures,
        "paths": paths,
        "features_and_coords": collate_complex(feature_and_coords),  # type: ignore
    }
    return collated_batch
