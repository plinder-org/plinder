from collections.abc import Iterable, Mapping
from dataclasses import fields, is_dataclass
from typing import Any

import pandas as pd
from biotite.structure.atoms import AtomArray


def atom_array_summary_markdown_repr(array: AtomArray) -> str:
    df = pd.DataFrame(
        {
            k: array.get_annotation(k)
            for k in array.get_annotation_categories()
            if k not in ["element", "atom_id", "b_factor", "atom_name"]
        }
    ).drop_duplicates()
    markdown: str = df.to_markdown(index=False)
    return markdown


def stringify_dataclass(
    obj: Any, indent: int = 4, _indents: int = 0, verbose_atom_array: bool = False
) -> str:
    """Pretty repr (or print) a (possibly deeply-nested) dataclass.
    Each new block will be indented by `indent` spaces (default is 4).

    https://stackoverflow.com/questions/66807878/pretty-print-dataclasses-prettier-with-line-breaks-and-indentation
    """
    if isinstance(obj, str):
        return f"'{obj}'"

    if not is_dataclass(obj) and not isinstance(obj, (Mapping, Iterable)):
        return str(obj)

    if hasattr(obj, "shape"):
        if isinstance(obj, AtomArray) and verbose_atom_array:
            return "\n" + atom_array_summary_markdown_repr(obj)

        return f"{type(obj)} with shape {obj.shape}"

    this_indent = indent * _indents * " "
    next_indent = indent * (_indents + 1) * " "
    # dicts, lists, and tuples will re-assign this
    start, end = f"{type(obj).__name__}(", ")"

    if is_dataclass(obj):
        body = "\n".join(
            f"{next_indent}{field.name}="
            f"{stringify_dataclass(getattr(obj, field.name), indent, _indents + 1)},"
            for field in fields(obj)
        )

    elif isinstance(obj, Mapping):
        if isinstance(obj, dict):
            start, end = "{", "}"

        body = "\n".join(
            f"{next_indent}{stringify_dataclass(key, indent, _indents + 1)}: "
            f"{stringify_dataclass(value, indent, _indents + 1)},"
            for key, value in obj.items()
        )

    else:  # is Iterable
        if isinstance(obj, list):
            start, end = "[", "]"
        elif isinstance(obj, tuple):
            start = "("

        body = "\n".join(
            f"{next_indent}{stringify_dataclass(item, indent, _indents + 1)},"
            for item in obj
        )

    return f"{start}\n{body}\n{this_indent}{end}"
