from plinder.core.structure.structure import Structure


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
        return structure.filter("atom_name", self.atom_types)  # type: ignore
