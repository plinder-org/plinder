from enum import Enum

from pydantic import BaseModel, ConfigDict


class Base(BaseModel):
    model_config = ConfigDict(
        extra="forbid", validate_assignment=True, use_enum_values=True
    )


class ChainConfig(Base):
    """Preparation configuration settings.

    Default configuration for preparation of normalized monomers,
    ready for use in PPI generation tasks.
    """

    decoy_receptor: list[str] = ["R"]
    decoy_ligand: list[str] = ["L"]
    native_receptor: list[str] = ["R"]
    native_ligand: list[str] = ["L"]


class BackboneDefinition(str, Enum):
    # Biotite (standard) C, CA, N
    biotite = "biotite"
    # DockQ definition: C, CA, N, O
    dockq = "dockq"
