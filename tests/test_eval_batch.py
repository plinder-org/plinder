"""Folder evaluation keeps failures, native assignments and bounded workers."""

import gzip
import json
import multiprocessing
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from types import SimpleNamespace

import pandas as pd
import pytest

from plinder.core.release import PlinderRelease
from plinder.eval import batch, evaluate_predictions
from plinder.eval.batch import _ligand_rows, _Reference


def _initialize_with_loaded_blas():
    from threadpoolctl import threadpool_info, threadpool_limits

    # pandas/NumPy were imported before the ProcessPool initializer. Ensure the
    # loaded runtime really has more than one thread, independent of node defaults.
    threadpool_limits(limits=2)
    global _pools_before_initializer
    _pools_before_initializer = threadpool_info()
    batch._worker_threads()


def _observe_worker_blas():
    import numpy as np
    from threadpoolctl import threadpool_info

    values = np.ones((16, 16))
    assert (values @ values).sum() == 4096
    return _pools_before_initializer, threadpool_info()


def test_worker_initializer_limits_already_loaded_blas():
    from threadpoolctl import threadpool_info

    parent_environment = dict(os.environ)
    parent_pools = threadpool_info()
    with ProcessPoolExecutor(
        max_workers=1,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_initialize_with_loaded_blas,
    ) as pool:
        before, after = pool.submit(_observe_worker_blas).result(timeout=60)
    before_blas = [entry for entry in before if entry["user_api"] == "blas"]
    after_blas = [entry for entry in after if entry["user_api"] == "blas"]
    assert before_blas and all(entry["num_threads"] == 2 for entry in before_blas)
    assert after_blas and all(entry["num_threads"] == 1 for entry in after_blas)
    assert all(entry["num_threads"] == 1 for entry in after)
    assert os.environ == parent_environment
    assert threadpool_info() == parent_pools


@pytest.fixture
def reference(test_dir):
    root = test_dir / "reconstructed_systems/1avd__1__1.A__1.C"
    return _Reference(
        "ligands",
        "1avd__1__1.A__1.C",
        root / "receptor.cif",
        {"reference_ligand": root / "ligand_files/1.C.sdf"},
    )


def test_ligand_summary_keeps_independent_assignments_and_missing_ligands(tmp_path):
    paths = {"one": tmp_path / "one.sdf", "missing": tmp_path / "missing.sdf"}
    ref = _Reference("ligands", "reference", tmp_path / "ref.cif", paths)

    def metric(model):
        return {
            "assigned_scores": [
                {
                    "reference_ligand": str(paths["one"]),
                    "model_ligand": model,
                    "score": 0.8,
                    "coverage": 1.0,
                    "bb_rmsd": 0.2,
                    "lddt_lp": 0.9,
                }
            ],
            "reference_ligand_unassigned_reason": {
                str(paths["missing"]): ["no_match", "No matching model ligand"]
            },
        }

    native = {
        "reference_ligands": list(map(str, paths.values())),
        "rmsd": metric("model1"),
        "lddt_pli": metric("model2"),
    }
    rows = _ligand_rows(native, ref, {"model1": "A", "model2": "B"})
    assert rows[0]["rmsd_model_ligand"] == "A"
    assert rows[0]["lddt_pli_model_ligand"] == "B"
    assert rows[1] == {
        "ligand_id": "missing",
        "status": "unassigned",
        "rmsd_unassigned": "no_match",
        "lddt_pli_unassigned": "no_match",
    }
    native["rmsd"]["assigned_scores"].append(native["rmsd"]["assigned_scores"][0])
    with pytest.raises(ValueError, match="twice"):
        _ligand_rows(native, ref, {"model1": "A", "model2": "B"})


def test_reference_lookup_uses_release_and_proper_selection(monkeypatch, tmp_path):
    release = PlinderRelease(data_dir=tmp_path)
    calls = []

    def query(table, **kwargs):
        calls.append((table, kwargs))
        return pd.DataFrame({
            "system_id": ["system", "system"]
            if table == "annotation"
            else ["interface"]
        })

    monkeypatch.setattr(batch, "query_table", query)

    def system(**kwargs):
        assert kwargs["release"] is release
        return SimpleNamespace(receptor_cif=tmp_path / "receptor.cif")

    def selected(reference, *, include_all_ligands):
        assert include_all_ligands is False
        return {"proper": tmp_path / "whole.sdf"}

    monkeypatch.setattr(batch, "PlinderSystem", system)
    monkeypatch.setattr(batch, "reference_ligands", selected)
    monkeypatch.setattr(
        batch,
        "PlinderInterface",
        lambda **kwargs: SimpleNamespace(interface_cif=tmp_path / "interface.cif"),
    )
    refs, failures = batch._references("1avd", "both", release, False)
    assert failures == []
    assert [r.system_id for r in refs] == ["system", "interface"]
    assert refs[0].ligands == {"proper": tmp_path / "whole.sdf"}
    assert all(
        kwargs["release"] is release
        and kwargs["filters"] == [("entry_pdb_id", "==", "1avd")]
        for _, kwargs in calls
    )


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
def test_native_folder_evaluation_with_two_workers(
    reference, test_dir, tmp_path, monkeypatch
):
    pytest.importorskip("posebusters")
    input_dir = tmp_path / "predictions"
    folder = input_dir / "1avd"
    folder.mkdir(parents=True)
    system = reference.receptor.parent
    shutil.copyfile(system / "system.cif", folder / "complete.CIF")
    shutil.copyfile(reference.receptor, folder / "missing.cif")
    interface = _Reference(
        "interfaces",
        "dimer",
        test_dir / "reconstructed_systems/1avd__1__1.A_2.A__1.D/receptor.cif",
    )
    monkeypatch.setattr(
        batch, "_references", lambda *args: ([reference, interface], [])
    )
    tables = evaluate_predictions(
        input_dir, output_dir=tmp_path / "evaluation", num_workers=2
    )
    assert tables["failures"].empty, tables["failures"].to_dict("records")
    ligands = tables["ligands"].set_index("prediction")
    assert len(ligands) == 2
    assert ligands.loc["1avd/complete.CIF", "rmsd"] == pytest.approx(0, abs=1e-3)
    assert ligands.loc["1avd/complete.CIF", "lddt_pli"] == pytest.approx(1)
    assert ligands.loc["1avd/missing.cif", "status"] == "unassigned"
    assert pd.isna(ligands.loc["1avd/missing.cif", "rmsd"])
    assert tables["interfaces"]["dockq_ave_full"].eq(0).all()
    assert set(tables["posebusters"]["model_ligand"]) == {"1.C"}
    assert tables["posebusters"]["mol_cond_loaded"].all()
    output = tmp_path / "evaluation"
    for key in ("ligands", "interfaces", "posebusters"):
        pd.testing.assert_frame_equal(
            pd.read_parquet(output / f"{key}.parquet"), tables[key]
        )
    assert (output / "failures.tsv").is_file()
    assert list((output / "details").rglob("*.json"))


def test_worker_failure_preserves_expected_rows(reference, tmp_path, monkeypatch):
    folder = tmp_path / "predictions" / "1avd"
    folder.mkdir(parents=True)
    (folder / "broken.cif").write_text("not a CIF")
    interface = _Reference("interfaces", "interface", reference.receptor)
    monkeypatch.setattr(
        batch, "_references", lambda *args: ([reference, interface], [])
    )
    tables = evaluate_predictions(
        folder.parent,
        output_dir=tmp_path / "evaluation",
        num_workers=1,
        ost_executable=str(tmp_path / "absent-ost"),
    )
    assert tables["ligands"]["status"].tolist() == ["error"]
    assert tables["interfaces"]["status"].tolist() == ["error"]
    assert set(tables["failures"]["stage"]) == {"prepare", "interfaces"}
    assert tables["ligands"]["rmsd"].isna().all()


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
@pytest.mark.parametrize("suffix", [".pdb", ".PDB", ".pdb.gz", ".PDB.gz"])
def test_native_pdb_interface_evaluation(test_dir, tmp_path, monkeypatch, suffix):
    from biotite.structure.io import pdb

    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )

    reference = _Reference(
        "interfaces",
        "dimer",
        test_dir / "reconstructed_systems/1avd__1__1.A_2.A__1.D/receptor.cif",
    )
    atoms = get_structure_with_altloc(read_mmcif_file(reference.receptor))
    # PDB needs single-character chains; OST must map these to the reference.
    atoms.chain_id = ["R" if chain == "1.A" else "L" for chain in atoms.chain_id]
    model = pdb.PDBFile()
    model.set_structure(atoms)
    folder = tmp_path / "predictions" / "dimer"
    folder.mkdir(parents=True)
    path = folder / f"model{suffix}"
    if suffix.endswith(".gz"):
        with gzip.open(path, "wt") as stream:
            model.write(stream)
    else:
        model.write(path)
    original = path.read_bytes()
    monkeypatch.setattr(batch, "_references", lambda *args: ([reference], []))
    tables = evaluate_predictions(
        folder.parent, output_dir=tmp_path / "evaluation", mode="interfaces"
    )
    assert tables["failures"].empty, tables["failures"].to_dict("records")
    row = tables["interfaces"].iloc[0]
    assert row["prediction"] == f"dimer/model{suffix}"
    assert row["status"] == "success"
    for metric in ("lddt", "ilddt", "qs_global", "qs_best", "dockq_ave_full"):
        assert row[metric] == pytest.approx(1)
    assert tables["ligands"].empty
    assert tables["posebusters"].empty
    assert path.read_bytes() == original
    assert not list((tmp_path / "evaluation" / "details").rglob("*.cif"))


def test_pdb_ligand_preparation_reports_format_error(reference, tmp_path, monkeypatch):
    folder = tmp_path / "predictions" / "1avd"
    folder.mkdir(parents=True)
    (folder / "model.pdb").write_text("HEADER    PREDICTION\nEND\n")
    monkeypatch.setattr(batch, "_references", lambda *args: ([reference], []))
    tables = evaluate_predictions(
        folder.parent, output_dir=tmp_path / "evaluation", mode="ligands"
    )
    assert tables["ligands"]["status"].tolist() == ["error"]
    assert tables["ligands"]["rmsd"].isna().all()
    assert tables["failures"]["stage"].tolist() == ["prepare"]
    assert "receptor with ligand SDFs" in tables["failures"].iloc[0]["error"]


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
@pytest.mark.parametrize("suffix", [".pdb", ".cif"])
def test_native_paired_ligand_evaluation(reference, tmp_path, monkeypatch, suffix):
    from biotite.structure.io import pdb
    from rdkit import Chem

    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )

    pytest.importorskip("posebusters")
    folder = tmp_path / "predictions" / "1avd"
    folder.mkdir(parents=True)
    receptor = folder / f"close{suffix}"
    if suffix == ".pdb":
        atoms = get_structure_with_altloc(read_mmcif_file(reference.receptor))
        atoms.chain_id[:] = "R"
        model = pdb.PDBFile()
        model.set_structure(atoms)
        model.write(receptor)
    else:
        shutil.copyfile(reference.receptor, receptor)
    for name in ("far", "missing", "different"):
        shutil.copyfile(receptor, folder / f"{name}{suffix}")
        (folder / f"{name}.ligands").mkdir()
    original_sdf = reference.ligands["reference_ligand"]
    shutil.copyfile(original_sdf, folder / "close.SDF")
    molecule = next(Chem.SDMolSupplier(str(original_sdf), removeHs=False))
    conformer = molecule.GetConformer()
    for index in range(molecule.GetNumAtoms()):
        position = conformer.GetAtomPosition(index)
        conformer.SetAtomPosition(index, (position.x + 100, position.y, position.z))
    with Chem.SDWriter(str(folder / "far.ligands/LIG.sdf")) as writer:
        writer.write(molecule)
    other = Chem.MolFromSmiles("CCO")
    conformer = Chem.Conformer(other.GetNumAtoms())
    for index in range(other.GetNumAtoms()):
        conformer.SetAtomPosition(index, (100 + index * 1.5, 0, 0))
    other.AddConformer(conformer)
    with Chem.SDWriter(str(folder / "different.ligands/LIG.sdf")) as writer:
        writer.write(other)
    monkeypatch.setattr(batch, "_references", lambda *args: ([reference], []))
    tables = evaluate_predictions(
        folder.parent, output_dir=tmp_path / "evaluation", mode="ligands", num_workers=2
    )
    assert tables["failures"].empty, tables["failures"].to_dict("records")
    rows = tables["ligands"].set_index("prediction")
    assert len(rows) == 4
    assert rows.loc[f"1avd/close{suffix}", "rmsd"] == pytest.approx(0, abs=1e-3)
    assert rows.loc[f"1avd/close{suffix}", "lddt_pli"] == pytest.approx(1)
    assert rows.loc[f"1avd/far{suffix}", "status"] == "success"
    assert pd.isna(rows.loc[f"1avd/far{suffix}", "rmsd"])
    assert rows.loc[f"1avd/far{suffix}", "rmsd_unassigned"] == "model_binding_site"
    assert rows.loc[f"1avd/far{suffix}", "lddt_pli"] == pytest.approx(0)
    assert rows.loc[f"1avd/far{suffix}", "lddt_pli_model_ligand"] == "LIG.sdf"
    assert rows.loc[f"1avd/missing{suffix}", "status"] == "unassigned"
    assert rows.loc[f"1avd/different{suffix}", "status"] == "unassigned"
    assert rows.loc[f"1avd/different{suffix}", "rmsd_unassigned"] == "identity"
    # Different chemistry with the same input filename stays prediction-specific.
    prepared_other = (
        tmp_path
        / "evaluation"
        / "details"
        / "1avd"
        / f"different{suffix}"
        / "model"
        / "ligand_0000.sdf"
    )
    assert Chem.MolToSmiles(next(Chem.SDMolSupplier(str(prepared_other)))) == "CCO"
    checks = tables["posebusters"].set_index("prediction")
    assert len(checks) == 3
    assert checks["mol_cond_loaded"].all()
    assert checks.loc[f"1avd/close{suffix}", "model_ligand"] == "close.SDF"
    assert checks.loc[f"1avd/far{suffix}", "model_ligand"] == "LIG.sdf"
    assert checks.loc[f"1avd/different{suffix}", "model_ligand"] == "LIG.sdf"


def test_worker_count_is_bounded_and_results_do_not_reuse_old_files(
    reference, tmp_path, monkeypatch
):
    """Use real spawned workers and a tiny tool double to observe concurrency."""
    folder = tmp_path / "predictions" / "1avd"
    folder.mkdir(parents=True)
    for index in range(5):
        (folder / f"model{index}.cif").write_text("input handled by test tool")
    events = tmp_path / "events"
    events.mkdir()
    tool = tmp_path / "ost"
    tool.write_text(
        "#!/usr/bin/env python\n"
        + f"""
import json, os, pathlib, sys, time
args = sys.argv
events = pathlib.Path({str(events)!r})
start = time.monotonic_ns()
time.sleep(0.1)
end = time.monotonic_ns()
(events / str(os.getpid())).write_text(json.dumps([start, end, os.environ['OMP_NUM_THREADS']]))
pathlib.Path(args[args.index('--output') + 1]).write_text(json.dumps({{'status': 'SUCCESS', 'lddt': 1, 'ilddt': 1, 'qs_global': 1, 'qs_best': 1, 'dockq': [1], 'dockq_ave_full': 1, 'dockq_wave_full': 1}}))
"""
    )
    tool.chmod(0o755)
    monkeypatch.setattr(
        batch,
        "_references",
        lambda *args: ([_Reference("interfaces", "interface", reference.receptor)], []),
    )
    output = tmp_path / "evaluation"
    output.mkdir()
    pd.DataFrame({"old": [True]}).to_parquet(output / "interfaces.parquet")
    tables = evaluate_predictions(
        folder.parent,
        output_dir=output,
        mode="interfaces",
        num_workers=2,
        ost_executable=tool,
    )
    assert tables["failures"].empty
    assert len(tables["interfaces"]) == 5
    intervals = [json.loads(path.read_text()) for path in events.iterdir()]
    assert len(intervals) == 5 and all(threads == "1" for _, _, threads in intervals)
    timeline = sorted(
        [(start, 1) for start, _, _ in intervals]
        + [(end, -1) for _, end, _ in intervals]
    )
    active = maximum = 0
    for _, delta in timeline:
        active += delta
        maximum = max(maximum, active)
    assert 1 <= maximum <= 2
    assert "old" not in pd.read_parquet(output / "interfaces.parquet")


@pytest.mark.parametrize("workers", [0, -1, True, 1.5])
def test_invalid_worker_count(tmp_path, workers):
    with pytest.raises(ValueError, match="positive integer"):
        evaluate_predictions(tmp_path, output_dir=tmp_path / "out", num_workers=workers)


def test_failed_reference_does_not_hide_other_references(tmp_path, monkeypatch):
    def query(table, **kwargs):
        if table == "interface_annotations":
            raise FileNotFoundError("interface table")
        return pd.DataFrame({"system_id": ["broken", "good"]})

    def system(*, system_id, release):
        if system_id == "broken":
            raise ValueError("missing source")
        return SimpleNamespace(receptor_cif=tmp_path / "receptor.cif")

    monkeypatch.setattr(batch, "query_table", query)
    monkeypatch.setattr(batch, "PlinderSystem", system)
    monkeypatch.setattr(
        batch,
        "reference_ligands",
        lambda *args, **kwargs: {"ligand": tmp_path / "pose.sdf"},
    )
    references, failures = batch._references(
        "1avd", "both", PlinderRelease(data_dir=tmp_path), False
    )
    assert [ref.system_id for ref in references] == ["good"]
    assert {failure["system_id"] for failure in failures} == {"broken", "1avd"}


def test_unavailable_references_produce_failure_table(tmp_path, monkeypatch):
    folder = tmp_path / "predictions" / "1avd"
    folder.mkdir(parents=True)
    (folder / "prediction.cif").touch()

    def unavailable(*args, **kwargs):
        raise FileNotFoundError("missing release")

    monkeypatch.setattr(batch, "query_table", unavailable)
    output = tmp_path / "evaluation"
    tables = evaluate_predictions(folder.parent, output_dir=output, mode="ligands")
    assert tables["ligands"].empty
    assert tables["failures"]["stage"].tolist() == ["reference"]
    assert pd.read_csv(output / "failures.tsv", sep="\t")["prediction"].tolist() == [
        "1avd/prediction.cif"
    ]


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
def test_native_evaluation_from_shared_table(reference, tmp_path, monkeypatch):
    folder = tmp_path / "inputs"
    folder.mkdir()
    shutil.copyfile(reference.receptor, folder / "receptor.cif")
    shutil.copyfile(reference.ligands["reference_ligand"], folder / "pose.sdf")
    table = folder / "inputs.tsv"
    pd.DataFrame({
        "input_id": ["prediction_1"],
        "reference_id": ["1avd"],
        "structure_path": ["receptor.cif"],
        "ligand_path": ["pose.sdf"],
    }).to_csv(table, sep="\t", index=False)
    monkeypatch.setattr(batch, "_references", lambda *args: ([reference], []))
    results = evaluate_predictions(
        table, output_dir=tmp_path / "results", mode="ligands", posebusters=False
    )
    assert results["failures"].empty, results["failures"].to_dict("records")
    row = results["ligands"].iloc[0]
    assert row["prediction"] == "prediction_1"
    assert row["rmsd_model_ligand"] == "pose.sdf"
    assert row["rmsd"] == pytest.approx(0, abs=1e-3)
    assert row["lddt_pli"] == pytest.approx(1)
