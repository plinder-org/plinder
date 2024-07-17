# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Namespace package root for plinder-data."""
from pathlib import Path

from plinder.data._version import _get_version

__version__ = _get_version()

_root = Path(__file__).parent
