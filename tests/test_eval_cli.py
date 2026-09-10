"""The command line forwards to the same evaluation API and reports failures."""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from plinder.eval import cli


@pytest.mark.parametrize("failures", [[], [{"error": "missing reference"}]])
def test_cli_forwards_options_and_reports_status(
    monkeypatch, tmp_path, capsys, failures
):
    observed = {}

    def evaluate(predictions, **kwargs):
        observed.update(predictions=predictions, **kwargs)
        return {
            "ligands": pd.DataFrame({"status": ["unassigned"]}),
            "failures": pd.DataFrame(failures),
        }

    monkeypatch.setattr(cli, "evaluate_predictions", evaluate)
    code = cli.main(
        [
            str(tmp_path / "predictions"),
            "--output-dir",
            str(tmp_path / "results"),
            "--mode",
            "ligands",
            "--num-workers",
            "3",
            "--data-dir",
            str(tmp_path / "release"),
            "--no-posebusters",
            "--include-all-ligands",
            "--ligand-smiles",
            json.dumps({"LIG": "CCO"}),
            "--ligand-ccd-codes",
            json.dumps({"OTHER": "ATP"}),
            "--ligand-chain",
            "P",
            "--ligand-chain",
            "Q",
            "--ligand-options=--min-pep-length 4",
            "--interface-options=--custom-mapping 'A:long chain'",
            "--ost-executable",
            "/custom/bin/ost",
        ]
    )
    assert code == int(bool(failures))
    assert observed["predictions"] == tmp_path / "predictions"
    assert observed["output_dir"] == tmp_path / "results"
    assert observed["release"].data_dir == tmp_path / "release"
    assert observed["mode"] == "ligands" and observed["num_workers"] == 3
    assert observed["posebusters"] is False
    assert observed["include_all_ligands"] is True
    assert observed["ligand_smiles"] == {"LIG": "CCO"}
    assert observed["ligand_ccd_codes"] == {"OTHER": "ATP"}
    assert observed["ligand_chains"] == ["P", "Q"]
    assert observed["ligand_options"] == ["--min-pep-length", "4"]
    assert observed["interface_options"] == ["--custom-mapping", "A:long chain"]
    assert observed["ost_executable"] == "/custom/bin/ost"
    assert "ligands: 1 rows" in capsys.readouterr().out


@pytest.mark.parametrize("value", ["bad JSON", "[]", '{"LIG": 1}', '{"LIG": ""}'])
def test_cli_rejects_invalid_chemistry_mapping(tmp_path, value):
    with pytest.raises(SystemExit) as exc:
        cli.main(
            [
                str(tmp_path),
                "--output-dir",
                str(tmp_path / "out"),
                "--ligand-smiles",
                value,
            ]
        )
    assert exc.value.code == 2


def test_cli_rejects_invalid_workers_before_loading_release(tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(
            [str(tmp_path), "--output-dir", str(tmp_path / "out"), "--num-workers", "0"]
        )
    assert exc.value.code == 2
    assert not (tmp_path / "out").exists()


def test_cli_help_without_optional_evaluation_packages():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            """
import importlib.abc
import runpy
import sys

class NoEvaluationPackages(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'ost', 'posebusters', 'PDBValidation', 'networkit'}:
            raise ModuleNotFoundError(fullname)

sys.meta_path.insert(0, NoEvaluationPackages())
sys.argv = ['plinder_eval', '--help']
runpy.run_module('plinder.eval.cli', run_name='__main__')
""",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--num-workers" in result.stdout
    assert "--ligand-ccd-codes" in result.stdout
    assert "--prediction_file" not in result.stdout


def test_packaged_command_uses_new_evaluator():
    import tomllib

    config = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())
    assert config["project"]["scripts"]["plinder_eval"] == "plinder.eval.cli:main"
