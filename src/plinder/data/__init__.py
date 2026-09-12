# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""PLINDER ingest and dataset-generation modules.

Optional dependencies are imported by the modules that use them.  Keeping the
package initializer dependency-free also lets :mod:`plinder.core` reuse the
mmCIF reconstruction utilities in a standard pip installation.
"""
