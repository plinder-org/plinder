# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from typing import Any, Callable

import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset

from plinder.core.index.system import PlinderSystem
from plinder.core.loader.featurizer import structure_featurizer
from plinder.core.loader.utils import collate_batch
from plinder.core.scores import query_index
from plinder.core.scores.query import FILTERS
from plinder.core.structure.structure import Structure
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)


class PlinderDataset(Dataset):  # type: ignore
    """
    Creates a dataset from plinder systems

    Parameters
    ----------
    split : str
        the split to sample from
    filters: FILTERS, default=None
        Index filter to select specific system ids
    use_alternate_structures: bool, default=True
        Whether to load alternate structures
    featurizer: Callable[
            [Structure, int], dict[str, torch.Tensor]
    ] = structure_featurizer,
        Transformation to turn structure to input tensors
    """

    def __init__(
        self,
        split: str,
        filters: FILTERS = None,
        use_alternate_structures: bool = True,
        featurizer: Callable[
            [Structure], torch.Tensor | dict[str, torch.Tensor]
        ] = structure_featurizer,
        **kwargs: Any,
    ):
        index = query_index(splits=[split], filters=filters)
        LOG.info(f"Loading {index.system_id.nunique()} systems")
        self._system_ids = list(set(index["system_id"]))
        self._num_examples = len(self._system_ids)

        self._featurizer = featurizer
        self._use_alternate_structures = use_alternate_structures

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(
        self, index: int
    ) -> dict[str, int | str | pd.DataFrame | dict[str, str | pd.DataFrame]]:
        if not 0 <= index < self._num_examples:
            raise IndexError(index)
        s = PlinderSystem(system_id=self._system_ids[index])

        holo_structure = s.holo_structure
        features_and_coords = None
        if self._featurizer is not None:
            features_and_coords = self._featurizer(holo_structure)

        item: dict[str, Any] = {
            "system_id": holo_structure.id,
            "holo_structure": holo_structure,
            "alternate_structures": s.alternate_structures
            if self._use_alternate_structures
            else {},
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
