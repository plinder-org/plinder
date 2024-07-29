# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations
from pathlib import Path
from typing import Callable

from torch.utils.data import Dataset
import atom3d.util.formats as fo
from biotite.structure.io.pdbx import PDBxFile

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
        transform: Callable | None = None,
        store_file_path: bool = True,
        load_alternative_structures: bool = False,
    ):
        self._system_ids = utils.get_manifest()["system_id"].to_list()
        self._num_examples = len(self._system_ids)
        self._transform = transform
        self._store_file_path = store_file_path
        self.load_alternative_structures = load_alternative_structures

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(self, index: int):
        if not 0 <= index < self._num_examples:
            raise IndexError(index)

        s = system.PlinderSystem.from_system_id(self._system_ids[index])

        buf = s.system_cif
        print("buf", buf)
        r = PDBxFile.read(buf)
        #s.system_cif)
        print("read", r)


        item = {
            "id": s.system_id,
            "system": s,
            "struc": r, # r.read(s.system_cif),
            "atoms": None, # s.
            "alternative_structures": {},
        }


        # file_path = self._file_list[index]
        # apo_pred_path = file_path.parent
        #
        # id = file_path.name
        # item = {
        #     "atoms": fo.bp_to_df(fo.read_any(file_path)),
        #     "id": id,
        #     "alternative_structures": {},
        # }
        # if self.load_alternative_structures:
        #     item["alternative_structures"] = {
        #         tag: {
        #             chain: fo.bp_to_df(fo.read_any(apo_pred_path / tag / apo_pred_name))
        #             for chain, apo_pred_name in chain_map.items
        #         }
        #         for tag, chain_map in APO_PRED_JSON[id].items
        #     }
        # if self._store_file_path:
        #     item["file_path"] = str(file_path)
        # if self._transform:
        #     item = self._transform(item)
        return item
