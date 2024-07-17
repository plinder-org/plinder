def test_ost():
    try:
        import openstructure as ost
    except (ImportError, ModuleNotFoundError):
        import ost

    assert ost is not None
