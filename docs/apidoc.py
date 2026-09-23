# The code in this file is based on the file with the same name in the Biotite
# project, licensed under the BSD-3-Clause license.

from __future__ import annotations

import ast
import shutil
from importlib import import_module
from pathlib import Path
from textwrap import dedent


def generate_api_reference(
    package_name: str,
    output_dir: Path,
    *,
    excluded_modules: set[str] | None = None,
) -> None:
    """Create module reference pages for an installed PLINDER package.

    Parameters
    ----------
    package_name : str
        Import name of the package to document.
    output_dir : Path
        Directory in which the generated ``.rst`` files are written.
    excluded_modules : set of str, optional
        Module names, or package prefixes, that require optional dependencies.
    """
    modules = _discover_modules(package_name, excluded_modules or set())

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    _write_package_index(output_dir / "index.rst", package_name, modules)
    for module_name, source_path in modules:
        _write_module_page(output_dir / f"{module_name}.rst", module_name, source_path)


def clear_api_reference(output_dir: Path) -> None:
    """Remove generated API directories while retaining hand-written pages."""
    if not output_dir.exists():
        return
    for path in output_dir.iterdir():
        if path.is_dir():
            shutil.rmtree(path)


def _discover_modules(
    package_name: str,
    excluded_modules: set[str],
) -> list[tuple[str, Path]]:
    """Return source modules below *package_name* in import-name order."""
    package = import_module(package_name)
    package_root = Path(package.__file__).parent
    modules: list[tuple[str, Path]] = []

    for source_path in package_root.rglob("*.py"):
        relative_path = source_path.relative_to(package_root)
        if source_path.name == "__init__.py":
            continue
        suffix = relative_path.with_suffix("").parts
        module_name = ".".join((package_name, *suffix))
        if not _is_excluded(module_name, excluded_modules):
            modules.append((module_name, source_path))

    return sorted(set(modules))


def _is_excluded(module_name: str, excluded_modules: set[str]) -> bool:
    return any(
        module_name == excluded or module_name.startswith(f"{excluded}.")
        for excluded in excluded_modules
    )


def _write_package_index(
    output_path: Path,
    package_name: str,
    modules: list[tuple[str, Path]],
) -> None:
    title = f"``{package_name}``\n" f'{"=" * (len(package_name) + 4)}\n'
    if modules:
        entries = "\n".join(f"    {module_name}" for module_name, _ in modules)
        output_path.write_text(
            title + "\n.. toctree::\n" + "    :maxdepth: 1\n\n" + f"{entries}\n"
        )
    else:
        output_path.write_text(title)


def _write_module_page(
    output_path: Path,
    module_name: str,
    source_path: Path,
) -> None:
    imported_names = _imported_names(source_path)
    excluded_members = ""
    if imported_names:
        excluded_members = f"    :exclude-members: {', '.join(imported_names)}\n"
    output_path.write_text(
        dedent(
            f"""
            :sd_hide_title: true

            ``{module_name}``
            {"=" * (len(module_name) + 4)}

            .. automodule:: {module_name}
                :members:
                :undoc-members:
                :show-inheritance:
                :member-order: bysource
            """
        ).rstrip()
        + "\n"
        + excluded_members
    )


def _imported_names(source_path: Path) -> list[str]:
    """Return names bound by top-level imports in a Python source file."""
    tree = ast.parse(source_path.read_text())
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            names.update(
                alias.asname or alias.name.split(".")[0] for alias in node.names
            )
        elif isinstance(node, ast.ImportFrom):
            names.update(alias.asname or alias.name for alias in node.names)
    return sorted(name for name in names if not name.startswith("_"))
