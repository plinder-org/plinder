import pytest
from plinder.core import PlinderSystem
from plinder.core.loader.transforms import (
    SelectAtomTypes,
    StructureTransform,
)
from plinder.core.utils import constants as pc


def test_transform_abc(read_plinder_mount):
    s = PlinderSystem(system_id="1avd__1__1.A__1.C").holo_structure
    with pytest.raises(NotImplementedError):
        StructureTransform().transform(s)


@pytest.mark.parametrize(
    "system_id, atom_types",
    [
        ("1avd__1__1.A__1.C", ["CA"]),
        ("1avd__1__1.A__1.C", ["CA", "N", "C", "O"]),
        ("1avd__1__1.A__1.C", ["foo"]),
    ],
)
def test_select_atom_types_structure_transform(
    read_plinder_mount, system_id, atom_types
):
    valid_atom_names = set(pc.ALL_ATOMS)
    expected_atom_names = set(atom_types).intersection(valid_atom_names)
    s = PlinderSystem(system_id=system_id).holo_structure
    t = SelectAtomTypes(atom_types=atom_types)
    assert len(s.protein_unique_atom_names) > len(expected_atom_names)
    s = t.transform(s)
    assert set(s.protein_unique_atom_names) == set(expected_atom_names)
