# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import atom3d.util.formats as fo
import pandas as pd
from torch.utils.data import Dataset

from plinder.core.split.utils import get_split
from plinder.core.system import system


class PlinderDataset(Dataset):  # type: ignore
    """
    Creates a dataset from plinder systems

    Parameters
    ----------
    df : pd.DataFrame | None
        the split to use
    split : str
        the split to sample from
    file_with_system_ids : str | Path
        path to a file containing a list of system ids (default: full index)
    store_file_path : bool, default=True
        if True, include the file path of the source structures in the dataset
    load_alternative_structures : bool, default=False
        if True, include alternative structures in the dataset
    num_alternative_structures : int, default=1
        number of alternative structures (apo and pred) to include
    """

    def __init__(
        self,
        df: pd.DataFrame | None = None,
        split: str = "train",
        split_parquet_path: str | Path | None = None,
        store_file_path: bool = True,
        load_alternative_structures: bool = False,
        num_alternative_structures: int = 1,
    ):
        if df is None:
            if split_parquet_path is None:
                df = get_split()
            else:
                df = pd.read_parquet(split_parquet_path)
        self._system_ids = df.loc[df["split"] == split, "system_id"].to_list()
        self._num_examples = len(self._system_ids)
        self._store_file_path = store_file_path
        self.load_alternative_structures = load_alternative_structures
        self.num_alternative_structures = num_alternative_structures

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(
        self, index: int
    ) -> dict[str, int | str | pd.DataFrame | dict[str, str | pd.DataFrame]]:
        if not 0 <= index < self._num_examples:
            raise IndexError(index)

        s = system.PlinderSystem(system_id=self._system_ids[index])

        item = {
            "id": index,
            "system_id": s.system_id,
            "df": fo.bp_to_df(fo.read_any(s.system_cif)),
            "alternative_structures": {},
        }
        if self._store_file_path:
            item["path"] = s.system_cif

        if self.load_alternative_structures:
            if s.linked_structures is not None and not s.linked_structures.empty:
                links = s.linked_structures.groupby("kind")
                for kind, group in links:
                    for link_id in group["id"].values[
                        : self.num_alternative_structures
                    ]:
                        structure = s.get_linked_structure(
                            link_kind=str(kind), link_id=link_id
                        )
                        item["alternative_structures"][
                            f"{kind}_{link_id}_df"
                        ] = fo.bp_to_df(fo.read_any(structure))
                        if self._store_file_path:
                            item["alternative_structures"][
                                f"{kind}_{link_id}_path"
                            ] = structure
        return item
