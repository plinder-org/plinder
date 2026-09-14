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
