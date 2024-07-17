# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
def test_ost():
    try:
        import openstructure as ost
    except (ImportError, ModuleNotFoundError):
        import ost

    assert ost is not None
