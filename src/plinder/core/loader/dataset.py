# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from typing import Any, Callable

import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset

from plinder.core.index.system import PlinderSystem
from plinder.core.loader.utils import collate_batch, structure2tensor
from plinder.core.scores.index import query_index
from plinder.core.scores.query import FILTERS
from plinder.core.structure.structure import Structure


def structure2tensor_transform(structure: Structure) -> dict[str, torch.Tensor]:
    props: dict[str, torch.Tensor] = structure2tensor(
        protein_atom_array=structure.protein_atom_array,
        input_sequence_residue_types=structure.input_sequence_list_ordered_by_chain,
        input_sequence_mask=structure.input_sequence_stacked_mask,
        input_sequence_full_atom_feat=structure.input_sequence_full_atom_feat,
        resolved_ligand_mols_coords=structure.resolved_ligand_mols_coords,
        input_ligand_conformers=structure.input_ligand_conformers,
        ligand_conformer2resolved_mask=structure.ligand_conformer2resolved_mask,
        input_ligand_conformer_coords=structure.input_ligand_conformer_coords,
        protein_chain_ordered=structure.protein_chain_ordered,
        ligand_chain_ordered=structure.ligand_chain_ordered,
        dtype=torch.float32,
    )
    return props


class PlinderDataset(Dataset):  # type: ignore
    """
    Creates a dataset from plinder systems

    Parameters
    ----------
    df : pd.DataFrame | None
        the split to use
    split : str
        the split to sample from
    input_structure_priority : str, default="apo"
        Which alternate structure to proritize
    transform: Callable[
        [Structure], torch.Tensor | dict[str, torch.Tensor]
    ] = structure2tensor_transform,
        Transformation to turn structure to input tensors
    """

    def __init__(
        self,
        split: str,
        filters: FILTERS = None,
        transform: Callable[
            [Structure], torch.Tensor | dict[str, torch.Tensor]
        ] = structure2tensor_transform,
        **kwargs: Any,
    ):
        self._system_ids = list(
            set(query_index(splits=[split], filters=filters)["system_id"])
        )
        self._num_examples = len(self._system_ids)
        self._transform = transform

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(
        self, index: int
    ) -> dict[str, int | str | pd.DataFrame | dict[str, str | pd.DataFrame]]:
        if not 0 <= index < self._num_examples:
            raise IndexError(index)
        s = PlinderSystem(system_id=self._system_ids[index])

        holo_structure = s.holo_structure
        if self._transform is not None:
            features_and_coords = self._transform(holo_structure)

        item: dict[str, Any] = {
            "system_id": holo_structure.id,
            "holo_structure": holo_structure,
            "alternate_structures": s.alternate_structures,
            "features_and_coords": features_and_coords,
            "path": s.system_cif,
        }
        return item


def get_torch_loader(
    dataset: PlinderDataset,
    batch_size: int = 2,
    shuffle: bool = True,
    sampler: None = None,  # None: 'Sampler[PlinderDataset]' | None = None,
    num_workers: int = 1,
    collate_fn: Callable[[list[dict[str, Any]]], dict[str, Any]] = collate_batch,
    **kwargs: Any,
) -> DataLoader[PlinderDataset]:
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        sampler=sampler,
        collate_fn=collate_fn,
        **kwargs,
    )
