from types import SimpleNamespace

import pytest

from plinder.data.annotations.get_ligand_validation import ResidueValidation
from plinder.data.annotations.ligand_validation import (
    format_selection_validation,
    validation_columns,
)


def _residue(rscc):
    return SimpleNamespace(
        validation=ResidueValidation(
            altcode=".",
            inscode=".",
            rsr=0.1,
            rsrz=0.2,
            rscc=rscc,
            average_occupancy=1.0,
            average_b_factor=20.0,
            unknown_residue=False,
            atom_count=5,
            unknown_atom_count=0,
            heavy_atom_count=5,
            num_unresolved_heavy_atoms=0,
            is_outlier=dict.fromkeys(
                ["geometry", "density", "chirality", "clashes"], False
            ),
            is_atom_count_consistent=True,
            has_clashing_partial_occupancy_atoms=False,
            alt_count=1,
        )
    )


def _chain(*residues):
    return SimpleNamespace(residues=dict(enumerate(residues, 1)), validation=None)


def test_ligand_summary_excludes_other_system_ligands():
    chains = {"F": _chain(_residue(0.852)), "E": _chain(_residue(0.984))}
    result = format_selection_validation({"1.F": [1]}, chains, prefix="ligand")
    assert result["ligand_validation_average_rscc"] == pytest.approx(0.852)
    assert result["ligand_validation_num_residues"] == 1
    assert set(result) == set(validation_columns("ligand"))


def test_multichain_summary_pools_residues_not_chain_means():
    chains = {"A": _chain(_residue(0.9), _residue(0.6)), "B": _chain(_residue(0.3))}
    result = format_selection_validation(
        {"1.A": [1, 2, 2], "1.B": [1]}, chains, prefix="ligand"
    )
    assert result["ligand_validation_average_rscc"] == pytest.approx(0.6)
    assert result["ligand_validation_num_residues"] == 3


def test_missing_validation_is_null_and_partial_coverage_is_visible():
    chains = {"A": _chain(_residue(0.9), SimpleNamespace(validation=None))}
    missing = format_selection_validation({"1.A": [2]}, chains, prefix="ligand")
    assert all(value is None for value in missing.values())
    partial = format_selection_validation(
        {"1.A": [1, 2]}, chains, prefix="ligand_pocket"
    )
    assert partial["ligand_pocket_validation_num_processed_residues"] == 1
    assert partial["ligand_pocket_validation_num_residues"] == 2
    assert partial["ligand_pocket_validation_percent_processed_residues"] == 50
    assert partial["ligand_pocket_validation_num_unresolved_heavy_atoms"] is None


def test_collation_keeps_ligand_validation_independent_of_system():
    from plinder.data.pipeline.collate import _annotation_columns_to_publish

    columns = [
        "ligand_id",
        *validation_columns("ligand"),
        *validation_columns("ligand_pocket"),
        "system_ligand_validation_average_rscc",
    ]
    assert _annotation_columns_to_publish(columns) == columns[:-1]


def test_ambiguous_validation_residue_remains_missing(monkeypatch):
    def ambiguous(*args, **kwargs):
        raise AssertionError("Duplicate alternative configurations")

    monkeypatch.setattr(
        "PDBValidation.Residue.Residue.CreateFromMmCIFPosition", ambiguous
    )
    assert ResidueValidation.from_residue("A", 1, "1", object()) is None
