# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from biotite.structure import AtomArray
from plinder.data.structure.atoms import atom_array_from_cif_file


# TODO: likely not needed
def test_atom_array_from_cif_file(cif_1qz5_unzipped):
    arr = atom_array_from_cif_file(cif_1qz5_unzipped)
    assert isinstance(arr, AtomArray)
