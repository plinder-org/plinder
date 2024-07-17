# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
Module for supplying version information.

This module provides the function `_get_version()`, which gets the version if
either setuptools_scm or our package is installed, and returns "unknown"
otherwise. Getting the version from setuptools_scm is primarily useful for
Docker image building (where it doesn't make sense to make the host install our
package just to obtain the version information to pass to the image build
process) and for editable installs (where having setuptools_scm installed is
the only way to get accurate version information).

When our package is installed, its version (the result of `get_version()`) can
be accessed as `plinder.__version__`.
"""

import warnings


def _get_version() -> str:
    try:
        # Our first choice would be to get the version from setuptools_scm if it
        # is installed (only way that works with editable installs)
        from setuptools_scm import get_version as scm_get_version

        version: str = scm_get_version(root="../..", relative_to=__file__)
    except (ImportError, LookupError):
        from importlib.metadata import PackageNotFoundError
        from importlib.metadata import version as importlib_version

        try:
            # Our second choice is to try to get the version from importlib
            version = importlib_version("plinder")
        except PackageNotFoundError:
            # We will land here if our package isn't actually installed
            warnings.warn("Neither our package nor setuptools_scm are installed")
            version = "unknown"
            version = "unknown"
    return version
