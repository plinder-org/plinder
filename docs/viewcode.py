# The code in this file is based on the file with the same name in the Biotite project
# licensed under BSD-3-clause license.
from __future__ import annotations

import inspect
from importlib import import_module
from typing import Any

import plinder


def linkcode_resolve(domain: str, info: dict[str, Any]) -> str | None:
    """
    See https://www.sphinx-doc.org/en/master/usage/extensions/linkcode.html.
    """
    version = plinder.__version__
    base_url = f"https://github.com/plinder-org/plinder/blob/v{version}/src/"

    if domain != "py":
        return None

    package_name = info["module"]
    attr_name = info["fullname"]

    package = import_module(package_name)
    try:
        attr = getattr(package, attr_name)
    except AttributeError:
        # The attribute is not defined within PLINDER or is part of a class
        # -> do not provide a link
        return None
    attr = getattr(package, attr_name)
    module = inspect.getmodule(attr)

    try:
        source_lines, first = inspect.getsourcelines(attr)
    except TypeError:
        # The attribute is some special object, e.g. a 'partial' object
        return None
    last = first + len(source_lines) - 1

    return base_url + f"{module.__name__.replace('.', '/')}.py#L{first}-L{last}"
