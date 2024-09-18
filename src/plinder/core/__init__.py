# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
The plinder.core package collects useful functions and classes for interacting
with the PLINDER dataset. It manages app configuration and will automatically
download (and / or sync) the dataset to a local cache in a lazy manner, when
particular assets are requested. One side effect of this is that plinder.core
will (by default) compare the MD5 checksums of files on disk and files in cloud
storage when they are accessed.

Note
----
You can disable the MD5 checksum comparison between local files and remote files
by setting the environment variable `PLINDER_OFFLINE=true`.
"""
from plinder.core.index.system import PlinderSystem
from plinder.core.index.utils import get_manifest, get_plindex
from plinder.core.split.utils import get_split
from plinder.core.utils.config import get_config

__all__ = [
    "get_config",
    "get_plindex",
    "get_manifest",
    "get_split",
    "PlinderSystem",
]
