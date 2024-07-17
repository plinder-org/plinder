# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest

try:
    import geotite
except ImportError:
    pytest.skip(allow_module_level=True)


def test_plinder_core():
    from plinder.core import geotite

    assert geotite is not None
