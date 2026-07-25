# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from plinder.data.annotations.get_ligand_validation import _select_altcode


def test_select_altcode_uses_selected_source_conformer() -> None:
    assert _select_altcode({"D", "C", "A", "B"}, "C") == "C"
    assert _select_altcode({".", "A"}, ".") == "."


def test_select_altcode_has_deterministic_fallback() -> None:
    assert _select_altcode({"D", "C", "A", "B"}, "missing") == "A"
    assert _select_altcode({".", "A"}, "missing") == "."
