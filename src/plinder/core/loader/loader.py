# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
from pathlib import Path

from torch.utils.data import Dataset
import atom3d.util.formats as fo


@cache
def load_apo_and_prep_data(
    representative_apo_pred_json: Path
) -> dict[str, dict[str, dict[str, str]]]:
    """
    Sample
    {
    "sys_id1": {
        "apo": {"A": "apo_id1",
                "B": "apo_id2"},
        "pred": {"A": "pred_id1",
            "B": "pred_id2"}
        },
    "sys_id2": {
        "apo": {"A": "apo_id1",
                "B": "apo_id2"},
        "pred": {"A": "pred_id1",
            "B": "pred_id2"}
        },
    }
    """
    with open(representative_apo_pred_json, "rb") as f:
        data = json.load(f)
    return data


APO_PRED_PATH = Path("")
APO_PRED_JSON = load_apo_and_prep_data(APO_PRED_PATH)


class PlinderDataset(Dataset):
    """
    Creates a dataset from a list of PDB files.
    :param file_list: path to LMDB file containing dataset
    :type file_list: list[Union[str, Path]]
    :param transform: transformation function for data augmentation, defaults to None
    :type transform: function, optional
    """

    def __init__(
        self,
        file_list,
        transform=None,
        store_file_path=True,
        load_alternative_structures=False,
    ):
        """constructor"""
        self._file_list = [Path(x).absolute() for x in file_list]
        self._num_examples = len(self._file_list)
        self._transform = transform
        self._store_file_path = store_file_path
        self.load_alternative_structures = load_alternative_structures

    def __len__(self) -> int:
        return self._num_examples

    def __getitem__(self, index: int):
        if not 0 <= index < self._num_examples:
            raise IndexError(index)

        file_path = self._file_list[index]
        apo_pred_path = file_path.parent

        id = file_path.name
        item = {
            "atoms": fo.bp_to_df(fo.read_any(file_path)),
            "id": id,
            "alternative_structures": {},
        }
        if self.load_alternative_structures:
            item["alternative_structures"] = {
                tag: {
                    chain: fo.bp_to_df(fo.read_any(apo_pred_path / tag / apo_pred_name))
                    for chain, apo_pred_name in chain_map.items
                }
                for tag, chain_map in APO_PRED_JSON[id].items
            }
        if self._store_file_path:
            item["file_path"] = str(file_path)
        if self._transform:
            item = self._transform(item)
        return item
