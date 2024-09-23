from __future__ import annotations

from pathlib import Path
from typing import Sequence

import torch
from torch import Tensor

from plinder.core.structure.structure import Structure

PAD_VALUE = -100


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
    ), "All `tensors` must have the same number of dimensions."

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
    batch_features: list[dict[str, Tensor]],
    pad_value: int = PAD_VALUE,
) -> dict[str, Tensor]:
    collated_and_padded_properties = {}
    batch_size = len(batch_features)

    feature_names = batch_features[0].keys()
    for feat_name in feature_names:
        collated_properties = []
        for idx in range(batch_size):
            feat = batch_features[idx][feat_name]
            collated_properties.append(feat)
        collated_and_padded_properties[feat_name] = pad_and_stack(
            collated_properties, dim=0, value=pad_value
        )
    return collated_and_padded_properties


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
        batch list[
        dict[
            str,
            list[str]
            | list[Structure]
            | list[dict[str, dict[str, Structure]]]
            | list[Path]
            | list[dict[str, Tensor]],
        ]
    ]: A list of dictionaries
        containing the data for each item in the batch.

    Returns:
        dict[
    str,
    list[str]
    | list[Structure]
    | list[dict[str, dict[str, Structure]]]
    | list[Path]
    | list[dict[str, Tensor]]
    A dictionary containing the merged items in the batch.

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
