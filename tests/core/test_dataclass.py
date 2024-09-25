from plinder.core.utils.dataclass import (
    atom_array_summary_markdown_repr,
    stringify_dataclass,
)


def test_stringify_dataclass(read_plinder_mount):
    from plinder.core import PlinderSystem

    system_id = "19hc__1__1.A_1.B__1.V_1.X_1.Y"
    system = PlinderSystem(system_id=system_id)
    struct = system.holo_structure
    assert isinstance(stringify_dataclass(struct), str)


def test_markdown_repr(read_plinder_mount):
    from plinder.core import PlinderSystem

    system_id = "19hc__1__1.A_1.B__1.V_1.X_1.Y"
    system = PlinderSystem(system_id=system_id)
    struct = system.holo_structure
    markdown = atom_array_summary_markdown_repr(struct.protein_atom_array)
    assert isinstance(markdown, str)
