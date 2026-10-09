"""Known interface sites can be located without an input binding partner."""

from types import SimpleNamespace

import pandas as pd
import pytest

from plinder.core.scores.custom import calculate_custom_interface_side_scores
from plinder.core.scores.entries import InterfaceView
from plinder.core.scores.interface import (
    INTERFACE_SIDE_COLUMNS,
    calculate_interface_side_scores,
)
from plinder.core.scores.mapping import pack_residue_identities


def interface():
    return InterfaceView(
        id="1abc__1__1.B--2.D",
        pdb_id="1abc",
        biounit_id="1",
        chain_1="1.B",
        chain_2="2.D",
        chain_1_residue_number_to_index={10: 0, 20: 1, 30: 2, 40: 3},
        chain_2_residue_number_to_index={1: 0, 2: 1},
        num_contact_residue_pairs=4,
    )


def alignment(**changes):
    return {
        "query_entry": "1abc",
        "target_entry": "monomer",
        "query_chain_mapped": "B",
        "target_chain_mapped": "A.with.dot",
        "source": "mmseqs",
        "query_selected_residue_numbers": [10, 20, 30, 99],
        "target_selected_residue_numbers": [101, 102, -1, 109],
        "selected_residue_identity_bits": pack_residue_identities(
            [True, False, True, True]
        ),
        **changes,
    }


def score(rows):
    known = interface()
    return calculate_interface_side_scores(
        pd.DataFrame(rows), interfaces={known.id: known}
    )


def test_monomer_hits_have_directional_site_coverage_and_identity():
    result = score([alignment()])
    assert len(result) == 1
    hit = result.iloc[0]
    assert hit.custom_structure_id == "monomer"
    assert hit.custom_chain == "A.with.dot"
    assert hit.plinder_interface_id == interface().id
    assert hit.plinder_chain == "1.B"
    assert hit.plinder_partner_chain == "2.D"
    assert hit.interface_side_qcov == 50
    assert hit.interface_side_fident == 25
    assert hit.plinder_residue_numbers == [10, 20]
    assert hit.custom_residue_numbers == [101, 102]


def test_alternative_alignments_are_not_unioned_and_backends_stay_separate():
    result = score(
        [
            alignment(),
            alignment(
                query_selected_residue_numbers=[30, 40],
                target_selected_residue_numbers=[103, 104],
                selected_residue_identity_bits=pack_residue_identities([True, True]),
            ),
            alignment(source="foldseek"),
        ]
    )
    assert len(result) == 2
    assert result.interface_side_qcov.tolist() == [50, 50]
    mmseqs = result.set_index("source").loc["mmseqs"]
    assert mmseqs.interface_side_fident == 50
    assert mmseqs.plinder_residue_numbers == [30, 40]


def test_either_interface_side_can_match():
    result = score(
        [
            alignment(
                query_chain_mapped="D",
                query_selected_residue_numbers=[1, 2],
                target_selected_residue_numbers=[101, 102],
            )
        ]
    )
    assert result.iloc[0].plinder_chain == "2.D"
    assert result.iloc[0].plinder_partner_chain == "1.B"
    assert result.iloc[0].interface_side_qcov == 100


@pytest.mark.parametrize(
    "rows",
    [
        [],
        [alignment(query_chain_mapped="unrelated")],
        [alignment(target_selected_residue_numbers=[-1, -1, -1, -1])],
    ],
)
def test_no_site_match_keeps_output_columns(rows):
    result = score(rows)
    assert result.empty
    assert list(result.columns) == INTERFACE_SIDE_COLUMNS


def test_mismatched_residue_maps_fail():
    with pytest.raises(ValueError):
        score([alignment(target_selected_residue_numbers=[101])])


def test_duplicate_residues_do_not_inflate_coverage():
    result = score([alignment(query_selected_residue_numbers=[10, 10, 20, 20])])
    assert result.iloc[0].interface_side_qcov == 50


def test_sequence_output_preserves_original_fasta_identifier(tmp_path, monkeypatch):
    from plinder.core.scores import custom

    known = interface()
    monkeypatch.setattr(
        custom,
        "_load_release_entry_views",
        lambda assets, *, pdb_ids: {
            "1abc": SimpleNamespace(interfaces={known.id: known})
        },
    )
    path = tmp_path / "alignments.parquet"
    pd.DataFrame([alignment()]).to_parquet(path)
    manifest = tmp_path / "manifest.parquet"
    pd.DataFrame(
        {"structure_id": ["monomer"], "sequence_id": ["my_original_sequence"]}
    ).to_parquet(manifest)
    output = tmp_path / "sites.parquet"
    calculate_custom_interface_side_scores(
        {"mmseqs": path}, assets=None, output_path=output, chain_manifest=manifest
    )
    result = pd.read_parquet(output)
    assert result.sequence_id.tolist() == ["my_original_sequence"]
    assert result.custom_residue_numbers.iloc[0].tolist() == [101, 102]
