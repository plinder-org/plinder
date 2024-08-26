# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from textwrap import dedent

_root = Path(__file__).parent.parent

from plinder.core.index.utils import get_manifest, get_plindex
from plinder.core.split.utils import get_split
from plinder.core.system.system import PlinderSystem
from plinder.core.utils.config import get_config

try:
    from plinder.core.loader.loader import PlinderDataset
except (ImportError, ModuleNotFoundError):
    print(
        dedent(
            """\
            plinder.core.PlinderDataset requires pytorch and atom3d.

            please run:

                pip install plinder[loader]

            to enable the data loader
            """
        )
    )
    PlinderDataset = None  # type: ignore

__all__ = [
    "get_config",
    "get_plindex",
    "get_manifest",
    "get_split",
    "PlinderSystem",
    "PlinderDataset",
]
