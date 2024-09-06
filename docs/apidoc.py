# The code in this file is based on the file with the same name in the Biotite project
# licensed under BSD-3-clause license.

import enum
import shutil
import types
from importlib import import_module
from pathlib import Path
from textwrap import dedent
from types import ModuleType
from typing import Any

from sphinx.application import Sphinx

_INDENT = " " * 4


def generate_api_reference(package_name: str, output_dir: Path) -> None:
    """
    Create the API reference for the given package and store it in the output directory.

    This creates `.rst` files containing `autodoc` directives.

    Parameters
    ----------
    package_name : str
        The name of the package to document.
    output_dir : Path
        The directory to store the generated files in.

    Notes
    -----
    Use `.rst` files are `apidoc`and `numpydoc` produce reStructuredText formatted text.
    """
    package = import_module(package_name)
    class_names = []
    function_names = []
    for attr_name in package.__all__:
        if _is_public_class(package, attr_name):
            class_names.append(attr_name)
        elif _is_public_function(package, attr_name):
            function_names.append(attr_name)

    # Remove existing rst files
    if output_dir.exists():
        shutil.rmtree(output_dir)
    # Create rst files
    output_dir.mkdir(parents=True, exist_ok=True)
    _create_package_page(
        output_dir / "index.rst", package_name, class_names, function_names
    )
    for class_name in class_names:
        _create_class_page(
            output_dir / f"{package_name}.{class_name}.rst", package_name, class_name
        )
    for function_name in function_names:
        _create_function_page(
            output_dir / f"{package_name}.{function_name}.rst",
            package_name,
            function_name,
        )


def _create_package_page(
    output_path: Path,
    package_name: str,
    classes: list[str],
    functions: list[str],
):
    """
    Create the `.rst` file for the package overview.

    Parameters
    ----------
    output_path : Path
        The path to the output file.
    package_name : str
        The name of the package as it would be imported.
    classes, functions : list of str
        The names of the classes and functions to be documented, that are attributes of
        the package, respectively.
    """
    attributes = classes + functions
    # Enumerate classes and functions
    summary_string = dedent(
        """

        .. autosummary::
            :nosignatures:
            :toctree:

    """
    )
    summary_string += "\n".join([_INDENT + attr for attr in attributes])

    # Assemble page
    file_content = (
        dedent(
            f"""

        ``{package_name}``
        {"=" * (len(package_name) + 4)}
        .. currentmodule:: {package_name}

        .. automodule:: {package_name}

        .. currentmodule:: {package_name}

    """
        )
        + summary_string
    )
    with open(output_path, "w") as f:
        f.write(file_content)


def _create_class_page(output_path, package_name, class_name):
    """
    Create the `.rst` file for the given class.

    Parameters
    ----------
    output_path : Path
        The path to the output file.
    package_name : str
        The name of the package as it would be imported.
    class_name : str
        The name of the class to document, that is part of the package.
    """
    file_content = dedent(
        f"""
        :sd_hide_title: true

        ``{class_name}``
        {"=" * (len(class_name) + 4)}
        .. autoclass:: {package_name}.{class_name}
            :show-inheritance:
            :members:
            :member-order: bysource
            :undoc-members:
            :inherited-members:
    """
    )
    with open(output_path, "w") as f:
        f.write(file_content)


def _create_function_page(output_path, package_name, function_name):
    """
    Create the `.rst` file for the given function.

    Parameters
    ----------
    output_path : Path
        The path to the output file.
    package_name : str
        The name of the package as it would be imported.
    function_name : str
        The name of the function to document, that is part of the package.
    """
    file_content = dedent(
        f"""
        :sd_hide_title: true

        ``{function_name}``
        {"=" * (len(function_name) + 4)}
        .. autofunction:: {package_name}.{function_name}
    """
    ).strip()
    with open(output_path, "w") as f:
        f.write(file_content)


def _is_public_class(module: ModuleType, attr_name: str) -> bool:
    """
    Check if the attribute is a public class.
    """
    return _is_public(attr_name) and isinstance(getattr(module, attr_name), type)


def _is_public_function(module: ModuleType, attr_name: str) -> bool:
    """
    Check if the attribute is a public function.
    """
    return _is_public(attr_name) and callable(getattr(module, attr_name))


def _is_public(attr_name: str) -> bool:
    """
    Check if the attribute is public.
    """
    return not attr_name.startswith("_")


def skip_nonrelevant(
    app: Sphinx, what: str, name: str, obj: Any, skip: bool, options: dict
) -> bool:
    """
    Skip all class members, that are not methods, enum values or inner
    classes, since other attributes are already documented in the class
    docstring.

    Furthermore, skip all class members, that are inherited from
    non-PLINDER base classes.

    See
    https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#event-autodoc-skip-member
    for more information.
    """
    if skip:
        return True
    if not _is_relevant_type(obj):
        return True
    if obj.__module__ is None:
        # Some built-in functions have '__module__' set to None
        return True
    package_name = obj.__module__.split(".")[0]
    if package_name != "plinder":
        return True
    return False


def _is_relevant_type(obj: Any) -> bool:
    """
    Check if the given object is an attribute that is worth documenting.
    """
    if type(obj).__name__ == "method_descriptor":
        # These are some special built-in Python methods
        return False
    return (
        (
            # Functions
            type(obj)
            in [types.FunctionType, types.BuiltinFunctionType, types.MethodType]
        )
        | (
            # Enum instance
            isinstance(obj, enum.Enum)
        )
        | (
            # Inner class
            isinstance(obj, type)
        )
    )
