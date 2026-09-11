# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import Any, Generator

from pydantic import BaseModel

EXCLUDE_DESCRIPTION_PREFIX = "[EXCLUDE]"
CUSTOM_EXPORT_DESCRIPTION_PREFIX = "[CUSTOM_EXPORT]"


def _description_text(description: str | None) -> str:
    if not description:
        return ""
    return str(description).lstrip().replace("\n", " ").strip()


def description_excluded_from_column_docs(description: str | None) -> bool:
    """Return whether a field is intentionally absent from column docs."""
    return _description_text(description).startswith(EXCLUDE_DESCRIPTION_PREFIX)


def description_excluded_from_flat_export(description: str | None) -> bool:
    """Return whether automatic flat export must skip a field."""
    return _description_text(description).startswith((
        EXCLUDE_DESCRIPTION_PREFIX,
        CUSTOM_EXPORT_DESCRIPTION_PREFIX,
    ))


def column_description_text(description: str | None) -> str:
    """Return public prose with any export marker removed."""
    text = _description_text(description)
    for prefix in (EXCLUDE_DESCRIPTION_PREFIX, CUSTOM_EXPORT_DESCRIPTION_PREFIX):
        if text.startswith(prefix):
            return text.removeprefix(prefix).lstrip()
    return text


class DocBaseModel(BaseModel):
    @classmethod
    def get_descriptions_and_types(cls) -> dict[str, tuple[str | None, str | None]]:
        """
        Returns a dictionary mapping attribute and property names to their descriptions and types.

        Returns:
        --------
        dict[str, str | None]
            A dictionary mapping attribute and property names to their descriptions and types.
        """
        descriptions = {}
        annotations = cls.__annotations__
        for name, value in cls.model_fields.items():
            descriptions[name] = (value.description, annotations[name])

        for name, prop in cls.__dict__.items():
            if isinstance(prop, cached_property) or isinstance(prop, property):
                dtype = None
                if hasattr(prop, "func"):
                    dtype = prop.func.__annotations__.get("return", None)
                descriptions[name] = (prop.__doc__, dtype)

        return descriptions

    @classmethod
    def document_properties(
        cls, prefix: str
    ) -> Generator[tuple[str, str | None, str], Any, Any]:
        for field, field_info in cls.get_descriptions_and_types().items():
            description, dtype = field_info
            if description_excluded_from_column_docs(description):
                continue
            if field.startswith(prefix):
                name = field
            else:
                name = f"{prefix}_{field}"
            if description:
                descr = column_description_text(description)
            else:
                descr = "[DESCRIPTION MISSING]"
            if "pass_criteria" in name and "validation" not in name:
                name = name.replace("pass_criteria", "pass_validation_criteria")
            yield (name, dtype, descr.strip())

    @classmethod
    def document_properties_to_tsv(
        cls, prefix: str, filename: Path, nested: bool = False
    ) -> None:
        with open(filename, "w") as tsv:
            tsv.write("\t".join(["Name", "Type", "Description"]) + "\n")
            for name, dtype, descr in cls.document_properties(prefix):
                if nested:
                    tsv.write("\t".join([name, f"list[{str(dtype)}]", descr]) + "\n")
                else:
                    tsv.write("\t".join([name, str(dtype), descr]) + "\n")
