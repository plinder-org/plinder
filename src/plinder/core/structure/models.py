from __future__ import annotations

from enum import Enum
from typing import List

from pydantic import BaseModel, ConfigDict

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
