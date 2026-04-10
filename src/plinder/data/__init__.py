# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from textwrap import dedent

try:
    import networkit  # noqa
except (ImportError, ModuleNotFoundError):
    raise ImportError(
        dedent(
            """\
            plinder.data requires networkit >= 11.0 to be installed.
            Please refer to the documentation for installation instructions.
            See details here:

                https://plinder-org.github.io/plinder/contribution/development.html#creating-the-conda-environment
            """
        )
    )
