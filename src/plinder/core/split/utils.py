from typing import Optional

import pandas as pd
from omegaconf import DictConfig

from plinder.core.utils import cpl
from plinder.core.utils.config import get_config
from plinder.core.utils.dec import timeit
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

_SPLIT = None


@timeit
def get_split(
    *,
    cfg: Optional[DictConfig] = None,
) -> pd.DataFrame:
    """
    Fetch the plinder split and cache it

    Parameters
    ----------
    cfg : DictConfig, default=None
        the plinder-core config

    Returns
    -------
    pd.DataFrame
        the plinder split
    """
    global _SPLIT
    if _SPLIT is not None:
        return _SPLIT
    cfg = cfg or get_config()
    suffix = f"{cfg.data.splits}/{cfg.data.split_file}"
    split = cpl.get_plinder_path(rel=suffix)
    LOG.info(f"reading {split}")
    _SPLIT = pd.read_parquet(split)
    return _SPLIT
