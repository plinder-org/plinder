# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import abc
from enum import Enum
from typing import Any, List

from pydantic import BaseModel, ConfigDict, PlainSerializer
from typing_extensions import Annotated

DOCKQ_BACKBONE_ATOMS = ["C", "CA", "N", "O"]


class Base(BaseModel):
    model_config = ConfigDict(
        extra="forbid", validate_assignment=True, use_enum_values=True
    )


class ChainConfig(Base):
    """Preparation configuration settings.

    Default configuration for preparation of normalized monomers,
    ready for use in PPI generation tasks.
    """

    decoy_receptor: List[str] = ["R"]
    decoy_ligand: List[str] = ["L"]
    native_receptor: List[str] = ["R"]
    native_ligand: List[str] = ["L"]


class BackboneDefinition(str, Enum):
    # Biotite (standard) C, CA, N
    biotite = "biotite"
    # DockQ definition: C, CA, N, O
    dockq = "dockq"


class DatasetName(str, Enum):
    plinder_s = "plinder_s"
    plinder_xl = "plinder_xl"
    plinder_af2 = "plinder_af2"


class MonomerName(str, Enum):
    holo = "holo"
    apo = "apo"
    predicted = "predicted"


class ExcludeComputedModelBase(BaseModel, metaclass=abc.ABCMeta):
    def model_dump(self, **kwargs) -> dict[str, Any]:  # type: ignore
        model_dict = super().model_dump(**kwargs)
        model_fields = set(self.model_fields.keys())
        return {key: value for key, value in model_dict.items() if key in model_fields}

    # @model_serializer()
    # def exclude_computed(self) -> dict[str, Any]:
    #     model_dict = super().model_dump()
    #     model_fields = set(self.model_fields.keys())
    #     return {key: value for key, value in model_dict.items() if key in model_fields}


def manual_model_dump(model_to_dump: ExcludeComputedModelBase) -> dict[str, Any]:
    return model_to_dump.model_dump()


ExcludeComputedModel = Annotated[
    ExcludeComputedModelBase, PlainSerializer(manual_model_dump)
]
