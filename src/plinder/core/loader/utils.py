from __future__ import annotations

from typing import Sequence

import torch
from torch import Tensor

from plinder.core.index.system import PlinderSystem

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
    all_collated_and_padded_properties = {}
    batch_size = len(batch_features)

    feature_groups_and_names = {k: list(v.keys()) for k, v in batch_features[0].items()}
    for feat_group, feature_names in feature_groups_and_names.items():
        collated_and_padded_properties = {}
        for feat_name in feature_names:
            collated_properties = []
            for idx in range(batch_size):
                feat = batch_features[idx][feat_group][feat_name]
                collated_properties.append(torch.tensor(feat))
            collated_and_padded_properties[feat_name] = pad_and_stack(
                collated_properties, dim=0, value=pad_value
            )
        all_collated_and_padded_properties[feat_group] = collated_and_padded_properties
    return collated_and_padded_properties


def collate_batch(
    batch: list[dict[str, str | PlinderSystem | dict[str, dict[str, Tensor]]]]
) -> dict[
    str,
    list[str] | list[PlinderSystem] | list[dict[str, dict[str, Tensor]]],
]:
    """Collate a batch of PlinderDataset items into a merged mini-batch of Tensors.

    Used as the default collate_fn for the torch DataLoader consuming PlinderDataset.

    Parameters:
       batch: list[
        dict[
            str,
            str
            | PlinderSystem
            | dict[str, dict[str, Tensor]]
            | Path,
        ],
    ]
        A list of dictionaries
            containing the data for each item in the batch.

    Returns:
        dict[
    str,
    | list[PlinderSystem]
    | list[dict[str, dict[str, Tensor]]]
    | list[Path]]
    A dictionary containing the merged items in the batch.

    """
    system_ids: list[str] = []
    plinder_system: list[PlinderSystem] = []
    feature_and_coords: list[dict[str, Tensor]] = []
    paths: list[str] = []
    for x in batch:
        assert isinstance(x["system_id"], str)
        assert isinstance(x["plinder_system"], PlinderSystem)
        assert isinstance(x["features_and_coords"], dict)
        assert isinstance(x["path"], str)

        system_ids.append(x["system_id"])
        plinder_system.append(x["plinder_system"])
        feature_and_coords.append(x["features_and_coords"])
        paths.append(x["path"])

    collated_batch: dict[
        str,
        list[str] | list[PlinderSystem] | list[dict[str, dict[str, Tensor]]],
    ] = {
        "system_ids": system_ids,
        "plinder_system": plinder_system,
        "paths": paths,
        "features_and_coords": collate_complex(feature_and_coords),  # type: ignore
    }
    return collated_batch
