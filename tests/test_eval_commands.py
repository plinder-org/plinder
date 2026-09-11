"""CLI contracts and small comparisons using the installed evaluation tools."""

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from plinder.eval.commands import run_openstructure, run_posebusters


@pytest.mark.parametrize("failure", ["exit", "absent", "status", "malformed"])
def test_failed_ost_run_preserves_previous_result(tmp_path, monkeypatch, failure):
    output = tmp_path / "result.json"
    output.write_text('{"previous": true}')

    def run(command, *, stdout, stderr, check):
        pending = Path(command[-1])
        stdout.write("diagnostic from tool")
        if failure == "status":
            pending.write_text('{"status": "FAILURE", "exception": "bad CIF"}')
        elif failure == "malformed":
            pending.write_text("not JSON")
        return subprocess.CompletedProcess(command, 1 if failure == "exit" else 0)

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises((RuntimeError, json.JSONDecodeError)):
        run_openstructure(
            "model.cif", "reference.cif", output, action="compare-structures"
        )
    assert output.read_text() == '{"previous": true}'
    assert output.with_suffix(".json.log").read_text() == "diagnostic from tool"
    assert not list(tmp_path.glob(".evaluation-*"))


def test_ost_retains_native_fields_and_arguments(tmp_path, monkeypatch):
    native = {
        "status": "SUCCESS",
        "dockq": [],
        "dockq_ave_full": 0.0,
        "future_ost_field": {"value": 3},
    }
    model = tmp_path / "model with spaces.cif"

    def run(command, **kwargs):
        assert command[command.index("--model") + 1] == str(model)
        assert "--min-pep-length" in command
        Path(command[-1]).write_text(json.dumps(native))
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(subprocess, "run", run)
    result = run_openstructure(
        model,
        "reference.cif",
        tmp_path / "result.json",
        action="compare-structures",
        options=["--min-pep-length", "4"],
    )
    assert result == native


@pytest.mark.parametrize("failure", ["missing_rows", "unloaded_receptor", "empty"])
def test_posebusters_rejects_incomplete_results(tmp_path, monkeypatch, failure):
    posebusters = pytest.importorskip("posebusters")
    from rdkit import Chem

    ligand = tmp_path / "pose.sdf"
    receptor = Chem.MolFromSmiles("C")
    receptor.GetAtomWithIdx(0).SetMonomerInfo(Chem.AtomPDBResidueInfo())
    ligand.touch()

    def bust(self, **kwargs):
        file = "another.sdf" if failure == "missing_rows" else str(ligand)
        return pd.DataFrame(
            []
            if failure == "empty"
            else [[file, 0, True, failure != "unloaded_receptor"]],
            columns=["file", "position", "mol_pred_loaded", "mol_cond_loaded"],
        )

    monkeypatch.setattr(posebusters.PoseBusters, "bust", bust)
    with pytest.raises(RuntimeError):
        run_posebusters([ligand], receptor, tmp_path / "result.csv")
    assert not (tmp_path / "result.csv").exists()


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
def test_native_ligand_comparison_retains_displaced_pose(test_dir, tmp_path):
    from rdkit import Chem

    system = test_dir / "reconstructed_systems/1avd__1__1.A__1.C"
    ligand = system / "ligand_files/1.C.sdf"
    molecule = next(Chem.SDMolSupplier(str(ligand)))
    shifted = tmp_path / "displaced.sdf"
    conformer = molecule.GetConformer()
    for index in range(molecule.GetNumAtoms()):
        position = conformer.GetAtomPosition(index)
        conformer.SetAtomPosition(index, (position.x + 100, position.y, position.z))
    with Chem.SDWriter(str(shifted)) as writer:
        writer.write(molecule)
    result = run_openstructure(
        system / "receptor.cif",
        system / "receptor.cif",
        tmp_path / "ligands.json",
        action="compare-ligand-structures",
        options=[
            "--model-ligands",
            str(ligand),
            str(shifted),
            "--reference-ligands",
            str(ligand),
        ],
    )
    assert len(result["model_ligands"]) == 2
    matched = result["rmsd"]["assigned_scores"][0]
    assert matched["score"] == pytest.approx(0, abs=1e-3)
    assert matched["bb_rmsd"] == pytest.approx(0, abs=1e-3)
    assert matched["lddt_lp"] == pytest.approx(1)
    assert str(shifted) in result["rmsd"]["model_ligand_unassigned_reason"]


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
def test_native_interface_comparison_keeps_missing_chain_penalty(test_dir, tmp_path):
    root = test_dir / "reconstructed_systems"
    reference = root / "1avd__1__1.A_2.A__1.D/receptor.cif"
    for model, expected in [
        (reference, 1.0),
        (root / "1avd__1__1.A__1.C/receptor.cif", 0.0),
    ]:
        result = run_openstructure(
            model,
            reference,
            tmp_path / "interface.json",
            action="compare-structures",
        )
        assert result["dockq_ave_full"] == pytest.approx(expected, abs=1e-3)
        assert {"ilddt", "qs_global", "qs_best"}.issubset(result)
        if model == reference:
            assert result["ilddt"] == pytest.approx(1)
            assert result["qs_global"] == pytest.approx(1)


def test_native_posebusters_checks_clashes_with_mmcif_receptor(test_dir, tmp_path):
    pytest.importorskip("posebusters")
    from biotite.interface import rdkit
    from biotite.structure.io import pdbx
    from rdkit import Chem

    system = test_dir / "reconstructed_systems/1avd__1__1.A__1.C"
    atoms = pdbx.get_structure(
        pdbx.CIFFile.read(system / "receptor.cif"), model=1, include_bonds=True
    )
    atoms.set_annotation("chain_id", ["long_chain_id"] * len(atoms))
    receptor = rdkit.to_mol(atoms)
    assert (
        receptor.GetAtomWithIdx(0).GetPDBResidueInfo().GetChainId() == "long_chain_id"
    )
    ligand = system / "ligand_files/1.C.sdf"
    molecule = next(Chem.SDMolSupplier(str(ligand)))
    conformer = molecule.GetConformer()
    coords = conformer.GetPositions()
    # Put a ligand atom exactly on a receptor atom to force a protein clash.
    coords += atoms.coord[0] - coords[0]
    for index, position in enumerate(coords):
        conformer.SetAtomPosition(index, position)
    clashing = tmp_path / "clashing.sdf"
    with Chem.SDWriter(str(clashing)) as writer:
        writer.write(molecule)
    result = run_posebusters([ligand, clashing], receptor, tmp_path / "validation.csv")
    assert len(result) == 2
    assert result["mol_cond_loaded"].all()
    assert result["mol_pred_loaded"].all()
    assert result["minimum_distance_to_protein"].notna().all()
    assert not result.set_index("file").loc[
        str(clashing), "minimum_distance_to_protein"
    ]


def test_posebusters_rejects_receptor_without_residue_metadata(tmp_path):
    pytest.importorskip("posebusters")
    from rdkit import Chem

    with pytest.raises(ValueError, match="residue metadata"):
        run_posebusters(["pose.sdf"], Chem.MolFromSmiles("CC"), tmp_path / "result.csv")
