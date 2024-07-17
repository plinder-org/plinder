"""Namespace package root for plinder-data."""
from pathlib import Path

from plinder.data._version import _get_version

__version__ = _get_version()

_root = Path(__file__).parent
