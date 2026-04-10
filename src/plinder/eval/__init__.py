# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from textwrap import dedent

try:
    import ost  # noqa
except (ImportError, ModuleNotFoundError):
    raise ImportError(
        dedent(
            """\
            plinder.eval requires OpenStructure >= 2.8.0 (ost).
            Install with: pip install plinder[eval]

            Note: OpenStructure requires numpy<2. Data generation
            (plinder.data) does NOT require OpenStructure.

            See: https://plinder-org.github.io/plinder/contribution/development.html#creating-the-conda-environment
            """
        )
    )
