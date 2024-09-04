# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from functools import cached_property

from pydantic import BaseModel


class DocBaseModel(BaseModel):
    @classmethod
    def get_descriptions(cls) -> dict[str, str | None]:
        """
        Returns a dictionary mapping attribute and property names to their descriptions.

        Returns:
        --------
        dict[str, str | None]
            A dictionary mapping attribute and property names to their descriptions.
        """
        descriptions = {}
        for name, value in cls.model_fields.items():
            descriptions[name] = value.description
        for name, prop in cls.__dict__.items():
            if isinstance(prop, cached_property) or isinstance(prop, property):
                descriptions[name] = prop.__doc__
        return descriptions
