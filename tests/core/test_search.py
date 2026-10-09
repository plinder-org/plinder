"""Public search modes select the requested comparisons."""

import importlib

import pandas as pd
import pytest


@pytest.mark.parametrize(
    "mode,ligands,interfaces",
    [
        ("auto", None, True),
        ("pockets", False, False),
        ("ligands", True, False),
        ("interfaces", False, True),
        ("all", None, True),
    ],
)
@pytest.mark.parametrize("include_monomers", [True, False])
def test_search_mode_dispatch(
    tmp_path, monkeypatch, mode, ligands, interfaces, include_monomers
):
    module = importlib.import_module("plinder.core.scores.search")
    source = tmp_path / "query.cif"
    source.write_text("data_query\n")
    table = pd.DataFrame({"input_id": ["query"], "structure_path": [str(source)]})
    monkeypatch.setattr(module, "score_custom_cif_files", lambda sources, **kw: kw)

    result = module.search(
        table,
        output_dir=tmp_path / "output",
        mode=mode,
        include_monomers=include_monomers,
    )

    assert result["include_ligands"] is ligands
    assert result["include_interfaces"] is interfaces
    assert result["include_monomers"] is include_monomers


@pytest.mark.parametrize("mode", ["auto", "pockets", "interfaces", "all"])
def test_fasta_discovers_sites(tmp_path, monkeypatch, mode):
    module = importlib.import_module("plinder.core.scores.search")
    monkeypatch.setattr(module, "score_custom_sequence_file", lambda source, **kw: kw)
    result = module.search(
        tmp_path / "query.fa", output_dir=tmp_path / "output", mode=mode
    )
    assert result["include_interfaces"] is (mode != "pockets")
    assert result["backends"] == ("mmseqs",)
    assert result["include_monomers"] is True


@pytest.mark.parametrize("mode", ["all", "ligands"])
def test_structure_file_does_not_require_a_table(tmp_path, monkeypatch, mode):
    module = importlib.import_module("plinder.core.scores.search")
    monkeypatch.setattr(module, "score_custom_cif_files", lambda sources, **kw: kw)
    source = tmp_path / "query.cif"
    source.write_text("data_query\n")
    result = module.search(source, output_dir=tmp_path / "output", mode=mode)
    assert result["include_ligands"] is (True if mode == "ligands" else None)
    assert result["include_interfaces"] is (mode != "ligands")


def test_fasta_cannot_compare_ligands(tmp_path):
    module = importlib.import_module("plinder.core.scores.search")
    with pytest.raises(ValueError, match="FASTA inputs"):
        module.search(
            tmp_path / "query.fa", output_dir=tmp_path / "output", mode="ligands"
        )


@pytest.mark.parametrize("mode", ["invalid", "both"])
def test_invalid_mode_lists_all(tmp_path, mode):
    module = importlib.import_module("plinder.core.scores.search")
    with pytest.raises(ValueError, match="interfaces or all"):
        module.search(tmp_path / "query.cif", output_dir=tmp_path / "output", mode=mode)


def test_default_runs_applicable_features(tmp_path, monkeypatch):
    module = importlib.import_module("plinder.core.scores.search")
    monkeypatch.setattr(module, "score_custom_sequence_file", lambda source, **kw: kw)
    result = module.search(tmp_path / "query.fa", output_dir=tmp_path / "output")
    assert result["include_interfaces"] is True
