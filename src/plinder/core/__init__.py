# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
The plinder.core package collects useful functions and classes for interacting
with the PLINDER dataset. It manages app configuration and will automatically
download (and / or sync) the dataset to a local cache in a lazy manner, when
particular assets are requested. Downloads are checked against the release
manifest for size and MD5 before replacing a local file. Cached files with the
expected size are reused without rehashing; force a refresh to repair same-size
local corruption.

Note
----
Set the environment variable `PLINDER_OFFLINE=true` to use local files without
network access.
"""

from plinder.core.index.interface import PlinderInterface
from plinder.core.index.query import query_table
from plinder.core.index.system import PlinderSystem
from plinder.core.release import (
    RELEASE_PATHS,
    RELEASE_TABLES,
    PlinderRelease,
)
from plinder.core.utils.config import get_config
from plinder.core.utils.io import download_pdb_mmcifs, get_pdb_mmcif

__all__ = [
    "get_config",
    "get_pdb_mmcif",
    "download_pdb_mmcifs",
    "PlinderInterface",
    "PlinderSystem",
    "query_table",
    "PlinderRelease",
    "RELEASE_PATHS",
    "RELEASE_TABLES",
]
