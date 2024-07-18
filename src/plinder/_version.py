def _get_version() -> str:
    try:
        from setuptools_scm import get_version as scm_get_version

        version: str = scm_get_version(
            root="../..",
            relative_to=__file__,
        )
    except (ImportError, LookupError):
        try:
            from importlib.metadata import PackageNotFoundError
            from importlib.metadata import version as importlib_version

            try:
                version = importlib_version(__name__.split(".")[0])
            except PackageNotFoundError:
                version = "unknown"
        except (ImportError, ModuleNotFoundError):
            version = "unknown"

    return version
