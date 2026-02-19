# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path

from ._version import _get_version

_root = Path(__file__).parent
__version__ = _get_version()
