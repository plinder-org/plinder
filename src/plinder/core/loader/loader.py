# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path
from typing import Callable

import atom3d.util.formats as fo
import pandas as pd
from torch.utils.data import Dataset

from plinder.core.index import utils
from plinder.core.system import system


class PlinderDataset(Dataset):
    """
    Creates a dataset from plinder systems


    Parameters
    ----------
    transform : Callable, default=None
        transformation function for data augmentation
    """

    def __init__(
        self,
        file_with_system_ids: str | Path,
        transform: Callable | None = None,
        store_file_path: bool = True,
        load_alternative_structures: bool = False,
    ):
        if file_with_system_ids is None:
            self._system_ids: list[str] = utils.get_manifest()["system_id"].to_list()
        else:
            self._system_ids = pd.read_csv(file_with_system_ids)["system_id"].to_list()
        self._num_examples = len(self._system_ids)
        self._transform = transform
        self._store_file_path = store_file_path
        self.load_alternative_structures = load_alternative_structures

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(self, index: int):
        if not 0 <= index < self._num_examples:
            raise IndexError(index)

        s = system.PlinderSystem(system_id=self._system_ids[index])

        item = {
            "id": index,
            "system_id": s.system_id,
            "df": fo.bp_to_df(fo.read_any(s.system_cif)),
            "alternative_structures": {},
        }

        if self.load_alternative_structures:
            links = s.linked_systems
            raise NotImplementedError
            # item["alternative_structures"] = {
            #     tag: {
            #         chain: fo.bp_to_df(fo.read_any(apo_pred_path / tag / apo_pred_name))
            #         for chain, apo_pred_name in chain_map.items
            #     }
            #     for tag, chain_map in APO_PRED_JSON[id].items
            # }
        if self._transform:
            item = self._transform(item)
        return item
