import pytest

try:
    import geotite
except ImportError:
    pytest.skip(allow_module_level=True)


def test_plinder_core():
    from plinder.core import geotite

    assert geotite is not None
