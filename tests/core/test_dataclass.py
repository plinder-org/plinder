from plinder.core.utils.dataclass import (
    atom_array_summary_markdown_repr,
    stringify_dataclass,
)


def test_stringify_dataclass(cached_plinder_system):
    system_id = "1avd__1__1.A__1.C"
    system = cached_plinder_system(system_id)
    struct = system.holo_structure
    assert isinstance(stringify_dataclass(struct), str)


def test_markdown_repr(cached_plinder_system):
    system_id = "1avd__1__1.A__1.C"
    system = cached_plinder_system(system_id)
    struct = system.holo_structure
    markdown = atom_array_summary_markdown_repr(struct.protein_atom_array)
    assert isinstance(markdown, str)
