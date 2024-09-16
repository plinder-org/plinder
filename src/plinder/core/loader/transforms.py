from plinder.core import PlinderSystem
from plinder.core.structure.structure import Structure


class TransformBase:
    def __init__(self) -> None:
        pass

    def __call__(self, dimer: PlinderSystem) -> PlinderSystem:
        return self.transform(dimer)

    def transform(self, dimer: PlinderSystem) -> PlinderSystem:
        raise NotImplementedError


class StructureTransform:
    def __call__(self, structure: Structure) -> Structure:
        return self.transform(structure)

    def transform(self, structure: Structure) -> Structure:
        raise NotImplementedError

    def __repr__(self) -> str:
        return self.__class__.__name__


class SelectAtomTypes(StructureTransform):
    def __init__(self, atom_types: list[str] = ["CA"]) -> None:
        self.atom_types = atom_types

    def transform(self, structure: Structure) -> Structure:
        return structure.filter("atom_name", self.atom_types)
