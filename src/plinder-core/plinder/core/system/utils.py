from typing import Any, Optional, Union

from omegaconf import DictConfig

from plinder.core.index.utils import _load_entries_from_zips


def load_systems(
    *,
    system_ids: Optional[Union[str, list[str]]] = None,
    pdb_ids: Optional[Union[str, list[str]]] = None,
    two_char_codes: Optional[Union[str, list[str]]] = None,
    cfg: Optional[DictConfig] = None,
) -> dict[str, Any]:
    if isinstance(pdb_ids, str):
        pdb_ids = [pdb_ids]
    if isinstance(two_char_codes, str):
        two_char_codes = [two_char_codes]
    result: dict[str, Any] = _load_entries_from_zips(
        cfg=cfg,
        two_char_codes=two_char_codes,
        pdb_ids=pdb_ids,
    )
    return result
