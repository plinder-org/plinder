# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from types import SimpleNamespace
from unittest.mock import Mock

import pandas as pd
import pytest
from plinder.data.annotations.get_ligand_validation import (
    ResidueListValidation,
    ResidueValidation,
    ResidueValidationThresholds,
    _select_altcode,
)


def test_select_altcode_uses_selected_source_conformer() -> None:
    assert _select_altcode({"D", "C", "A", "B"}, "C") == "C"
    assert _select_altcode({".", "A"}, ".") == "."


def test_select_altcode_has_deterministic_fallback() -> None:
    assert _select_altcode({"D", "C", "A", "B"}, "missing") == "A"
    assert _select_altcode({".", "A"}, "missing") == "."


@pytest.mark.parametrize(
    ("reference_count", "polymer", "expected"),
    [
        (None, False, None),
        (None, True, None),
        (5, False, 0),
        (8, False, 3),
        (6, False, 1),
        (4, False, None),
        (0, False, None),
        (6, True, 0),
        (4, True, 0),
        (8, True, 3),
        (3, True, None),
    ],
)
def test_unresolved_atoms_require_a_reference(
    monkeypatch, reference_count, polymer, expected
):
    residue = Mock(
        **{
            "countAtomsHeavyPDBX.return_value": 5,
            "countAtomsHeavyConop.return_value": reference_count,
            "isType.return_value": polymer,
            "getAltCode.return_value": ".",
            "getInsertionCode.return_value": ".",
            "getRSR.return_value": 0.1,
            "getRSRZ.return_value": 0.2,
            "getRSCC.return_value": 0.9,
            "getOccupancies.return_value": [1.0],
            "getBFactors.return_value": [20.0],
            "getResName.return_value": "LIG",
            "countAtomsPDBX.return_value": 5,
            "countUnknownAtoms.return_value": 0,
            "isOutlier.return_value": False,
            "isAtomCountConsistent.return_value": True,
            "hasClashingPartialOccupancyAtoms.return_value": False,
        }
    )
    if reference_count is None:
        residue.isAtomCountConsistent.side_effect = TypeError("missing reference")
    monkeypatch.setattr(
        "PDBValidation.Residue.Residue.CreateFromMmCIFPosition",
        lambda *_args: SimpleNamespace(
            listAlt=lambda: {"."}, getAlt=lambda _alt: residue
        ),
    )
    validation = ResidueValidation.from_residue("A", 1, "1", object())
    assert validation is not None
    assert validation.num_unresolved_heavy_atoms == expected
    assert validation.heavy_atom_count == 5
    assert validation.rscc == 0.9
    if reference_count is None:
        assert validation.is_atom_count_consistent is None
        residue.isAtomCountConsistent.assert_not_called()
    else:
        residue.isAtomCountConsistent.assert_called_once_with()


def _validation(unresolved: int | None) -> ResidueValidation:
    return ResidueValidation(
        altcode=".",
        inscode=".",
        rsr=0.1,
        rsrz=0.2,
        rscc=0.9,
        average_occupancy=1.0,
        average_b_factor=20.0,
        unknown_residue=False,
        atom_count=5,
        unknown_atom_count=0,
        heavy_atom_count=5,
        num_unresolved_heavy_atoms=unresolved,
        is_outlier={
            name: False for name in ["geometry", "density", "chirality", "clashes"]
        },
        is_atom_count_consistent=True,
        has_clashing_partial_occupancy_atoms=False,
        alt_count=1,
    )


@pytest.mark.parametrize(
    ("counts", "expected"),
    [
        ([0, 2], 2),
        ([0, 0], 0),
        ([None, 2], None),
        ([None, None], None),
        ([-5, 2], None),
    ],
)
def test_unresolved_total_preserves_unknown_counts(counts, expected, tmp_path):
    validation = ResidueListValidation.from_residues(
        [_validation(count) for count in counts], ResidueValidationThresholds()
    )
    assert validation is not None
    assert validation.num_unresolved_heavy_atoms == expected
    assert validation.average_rscc == pytest.approx(0.9)
    assert validation.num_processed_residues == 2
    restored = ResidueListValidation.model_validate_json(validation.model_dump_json())
    assert restored.num_unresolved_heavy_atoms == expected
    path = tmp_path / "validation.parquet"
    pd.DataFrame([restored.format()]).to_parquet(path, index=False)
    value = pd.read_parquet(path).loc[0, "validation_num_unresolved_heavy_atoms"]
    assert pd.isna(value) if expected is None else value == expected


def test_unresolved_total_is_unknown_when_a_residue_was_not_validated():
    validation = ResidueListValidation.from_residues(
        [_validation(0), None], ResidueValidationThresholds()
    )
    assert validation is not None
    assert validation.num_unresolved_heavy_atoms is None
    assert validation.num_processed_residues == 1
    assert validation.percent_processed_residues == 50


@pytest.mark.parametrize(("unresolved", "passes"), [(0, True), (None, False)])
def test_unknown_pocket_atom_count_cannot_pass_system_validation(unresolved, passes):
    from plinder.data.annotations.aggregate_annotations import QualityCriteria, System
    from plinder.data.annotations.ligand_utils import Ligand

    system = System(
        pdb_id="1abc",
        biounit_id="1",
        receptor_type="protein",
        ligands=[
            Ligand(
                pdb_id="1abc",
                biounit_id="1",
                asym_id="L",
                instance=1,
                num_heavy_atoms=5,
                num_unresolved_heavy_atoms=0,
            )
        ],
        ligand_validation=ResidueListValidation.from_residues(
            [_validation(0)], ResidueValidationThresholds()
        ),
        pocket_validation=ResidueListValidation.from_residues(
            [_validation(unresolved)], ResidueValidationThresholds()
        ),
    )
    result = system.format_validation(True, QualityCriteria())
    assert result["system_pass_validation_criteria"] is passes
