# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License, Version 2.0
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from plinder import cli


def test_discover_mmcif_files_is_nonrecursive_by_default(tmp_path):
    direct = tmp_path / "direct.CIF"
    direct.touch()
    nested = tmp_path / "nested" / "model.mmcif.gz"
    nested.parent.mkdir()
    nested.touch()
    (tmp_path / "notes.txt").touch()

    assert cli._discover_mmcif_files(tmp_path, recursive=False) == [direct]
    assert cli._discover_mmcif_files(tmp_path, recursive=True) == [direct, nested]


def test_component_mapping_rejects_duplicate_component():
    with pytest.raises(ValueError, match="repeats component"):
        cli._component_mapping(
            ["LIG=CCO", "LIG=CCC"],
            option="--ligand-smiles",
        )


def test_link_cli_uses_public_assets_and_auto_mode(tmp_path, monkeypatch):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    cif_file = input_dir / "model.cif"
    cif_file.touch()
    output_dir = tmp_path / "links"
    public_data = tmp_path / "public-release"
    public_data.mkdir()
    observed = {}

    def fake_score(cif_files, **kwargs):
        observed["cif_files"] = list(cif_files)
        observed.update(kwargs)
        kwargs["work_dir"].mkdir(parents=True)
        protein = kwargs["work_dir"] / "protein_scores.parquet"
        interface = kwargs["work_dir"] / "interface_scores.parquet"
        residue_pairs = kwargs["work_dir"] / "aligned_pocket_residues.parquet"
        pd.DataFrame({"value": [1, 2]}).to_parquet(protein, index=False)
        pd.DataFrame({"value": [1]}).to_parquet(interface, index=False)
        pd.DataFrame({"value": [1, 2, 3]}).to_parquet(
            residue_pairs,
            index=False,
        )
        return SimpleNamespace(
            protein_scores=protein,
            ligand_scores=None,
            interface_scores=interface,
            aligned_pocket_residues=residue_pairs,
        )

    monkeypatch.setattr(cli, "score_custom_cif_files", fake_score)
    args = cli.build_parser().parse_args(
        [
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--data-dir",
            str(public_data),
            "--ligand-ccd",
            "LIG=ATP",
            "--save-aligned-pocket-residues",
        ]
    )
    summary = cli._run_link(args)

    assert cli.build_parser().prog == "plinder_link"
    assert observed["cif_files"] == [cif_file]
    assert observed["data_dir"] == public_data
    assert observed["include_ligands"] is None
    assert observed["include_interfaces"] is None
    assert observed["ligand_ccd_code_dict"] == {"LIG": "ATP"}
    assert observed["store_aligned_pocket_residues"] is True
    assert summary["score_modes"] == ["protein", "interface"]
    assert summary["outputs"]["protein_scores"]["rows"] == 2
    assert summary["outputs"]["aligned_pocket_residues"]["rows"] == 3
    assert Path(summary["summary"]).is_file()
