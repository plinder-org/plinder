# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from textwrap import dedent

try:
    import ost
except (ImportError, ModuleNotFoundError):
    raise ImportError(
        dedent(
            """\
            plinder.data requires the OpenStructureToolkit (ost) to be installed.
            Please refer to the documentation for installation instructions and current limitations.
            See the note here: https://github.com/plinder-org/plinder?tab=readme-ov-file#-getting-started
            """
        )
    )

try:
    import networkit
except (ImportError, ModuleNotFoundError):
    raise ImportError(
        dedent(
            """\
            plinder.data requires the networkit library to be installed.
            Please refer to the documentation for installation instructions and current limitations.
            See the note here: https://github.com/plinder-org/plinder?tab=readme-ov-file#-getting-started
            """
        )
    )
