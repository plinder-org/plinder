"""Evaluate folders of predictions against PLINDER reference structures."""

from __future__ import annotations

import multiprocessing
import os
import sys
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Literal
from urllib.parse import quote

import pandas as pd

from plinder.core.index.interface import PlinderInterface
from plinder.core.index.query import query_table
from plinder.core.index.system import PlinderSystem
from plinder.core.release import PlinderRelease
from plinder.eval.commands import run_openstructure, run_posebusters
from plinder.eval.inputs import prepare_prediction, reference_ligands

_INTERFACE_METRICS = (
    "lddt",
    "ilddt",
    "qs_global",
    "qs_best",
    "dockq",
    "dockq_ave_full",
    "dockq_wave_full",
)
_LIGAND_COLUMNS = (
    "ligand_id",
    "rmsd",
    "bb_rmsd",
    "lddt_lp",
    "rmsd_coverage",
    "rmsd_model_ligand",
    "rmsd_unassigned",
    "lddt_pli",
    "lddt_pli_coverage",
    "lddt_pli_model_ligand",
    "lddt_pli_unassigned",
)


@dataclass(frozen=True)
class _Reference:
    kind: str
    system_id: str
    receptor: Path
    ligands: dict[str, Path] = field(default_factory=dict)


def _references(
    name: str, mode: str, release: PlinderRelease, include_all: bool
) -> tuple[list[_Reference], list[dict[str, str]]]:
    """Resolve and reconstruct references in the parent, avoiding cache races."""
    references, failures = [], []
    column = "entry_pdb_id" if len(name) == 4 and name.isalnum() else "system_id"
    for kind, table in (
        ("ligands", "annotation"),
        ("interfaces", "interface_annotations"),
    ):
        if mode not in (kind, "both"):
            continue
        try:
            rows = query_table(
                table,
                columns=["system_id"],
                filters=[(column, "==", name)],
                release=release,
            )
        except Exception as exc:
            failures.append(
                {"system_id": name, "error": f"{table}: {type(exc).__name__}: {exc}"}
            )
            continue
        for system_id in sorted(set(rows["system_id"])):
            try:
                if kind == "ligands":
                    system = PlinderSystem(system_id=system_id, release=release)
                    reference = _Reference(
                        kind,
                        system_id,
                        Path(system.receptor_cif).resolve(),
                        reference_ligands(system, include_all_ligands=include_all),
                    )
                else:
                    interface = PlinderInterface(system_id=system_id, release=release)
                    reference = _Reference(
                        kind, system_id, Path(interface.interface_cif).resolve()
                    )
                references.append(reference)
            except Exception as exc:
                failures.append(
                    {"system_id": system_id, "error": f"{type(exc).__name__}: {exc}"}
                )
    if not references and not failures:
        failures.append(
            {"system_id": name, "error": f"No {mode} references found for {name}"}
        )
    return references, failures


def _ligand_rows(
    result: dict[str, Any], reference: _Reference, model_ids: dict[str, str]
) -> list[dict[str, Any]]:
    """Keep each reference ligand and each metric's independent assignment."""
    reference_ids = {
        str(path.resolve()): ligand_id for ligand_id, path in reference.ligands.items()
    }
    if set(result["reference_ligands"]) != set(reference_ids):
        raise ValueError("OST reference ligands differ from the selected reference")
    rows = {
        key: {"ligand_id": value, "status": "success"}
        for key, value in reference_ids.items()
    }
    for metric in ("rmsd", "lddt_pli"):
        assigned = set()
        for score in result[metric]["assigned_scores"]:
            key = score["reference_ligand"]
            if key in assigned:
                raise ValueError(f"OST assigned {key} twice for {metric}")
            assigned.add(key)
            row = rows[key]
            row[metric] = score["score"]
            row[f"{metric}_coverage"] = score["coverage"]
            row[f"{metric}_model_ligand"] = model_ids[score["model_ligand"]]
            if metric == "rmsd":
                row.update({key: score[key] for key in ("bb_rmsd", "lddt_lp")})
        unassigned = result[metric]["reference_ligand_unassigned_reason"]
        if assigned & unassigned.keys() or assigned | unassigned.keys() != rows.keys():
            raise ValueError(f"OST returned incomplete {metric} assignments")
        for key, reason in unassigned.items():
            rows[key][f"{metric}_unassigned"] = reason[0]
    for row in rows.values():
        if "rmsd" not in row and "lddt_pli" not in row:
            row["status"] = "unassigned"
    return list(rows.values())


def _worker_threads() -> None:
    from threadpoolctl import threadpool_limits

    # OST subprocesses inherit these limits; PoseBusters' own pool is disabled.
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = "1"
    # Spawn imports pandas/NumPy before the initializer. Limit their already
    # loaded BLAS/OpenMP runtimes, not just libraries in future subprocesses.
    threadpool_limits(limits=1)
    if "numexpr" in sys.modules:
        sys.modules["numexpr"].set_num_threads(1)


def _evaluate_one(
    task: tuple[Path, str, list[_Reference]],
    *,
    output_dir: Path,
    posebusters: bool,
    ligand_smiles: dict[str, str],
    ligand_ccd_codes: dict[str, str],
    ligand_chains: tuple[str, ...],
    ligand_options: tuple[str, ...],
    interface_options: tuple[str, ...],
    ost_executable: str,
) -> dict[str, list[dict[str, Any]]]:
    model, prediction, references = task
    details = output_dir / "details" / prediction
    results: dict[str, list[dict[str, Any]]] = {
        key: [] for key in ("ligands", "interfaces", "posebusters", "failures")
    }

    def failure(stage: str, exc: Exception, system_id: str = "") -> None:
        results["failures"].append(
            {
                "prediction": prediction,
                "system_id": system_id,
                "stage": stage,
                "error": f"{type(exc).__name__}: {exc}",
            }
        )

    # Interface comparisons do not depend on ligand chemistry or PoseBusters.
    for ref in references:
        if ref.kind != "interfaces":
            continue
        row = {"prediction": prediction, "system_id": ref.system_id, "status": "error"}
        try:
            native = run_openstructure(
                model,
                ref.receptor,
                details / f"interface-{quote(ref.system_id, safe='._-')}.json",
                action="compare-structures",
                options=interface_options,
                executable=ost_executable,
            )
            row.update({key: native[key] for key in _INTERFACE_METRICS})
            row["status"] = "success"
        except Exception as exc:
            failure("interfaces", exc, ref.system_id)
        results["interfaces"].append(row)

    ligand_refs = [ref for ref in references if ref.kind == "ligands"]
    if not ligand_refs:
        return results
    prepared = None
    try:
        prepared = prepare_prediction(
            model,
            details / "model",
            ligand_smiles=ligand_smiles,
            ligand_ccd_codes=ligand_ccd_codes,
            ligand_chains=ligand_chains,
        )
    except Exception as exc:
        failure("prepare", exc)
    for ref in ligand_refs:
        rows = [
            {"ligand_id": ligand_id, "status": "error"} for ligand_id in ref.ligands
        ]
        if prepared is not None:
            try:
                native = run_openstructure(
                    prepared.receptor,
                    ref.receptor,
                    details / f"ligands-{quote(ref.system_id, safe='._-')}.json",
                    action="compare-ligand-structures",
                    options=[
                        *ligand_options,
                        "--model-ligands",
                        *map(str, prepared.ligands.values()),
                        "--reference-ligands",
                        *(str(path.resolve()) for path in ref.ligands.values()),
                    ],
                    executable=ost_executable,
                )
                rows = _ligand_rows(
                    native,
                    ref,
                    {str(path): key for key, path in prepared.ligands.items()},
                )
            except Exception as exc:
                failure("ligands", exc, ref.system_id)
        results["ligands"].extend(
            {"prediction": prediction, "system_id": ref.system_id, **row}
            for row in rows
        )
    if posebusters and prepared is not None and prepared.ligands:
        try:
            checks = run_posebusters(
                list(prepared.ligands.values()),
                prepared.receptor_molecule,
                details / "posebusters.csv",
            )
            ids = {str(path): key for key, path in prepared.ligands.items()}
            checks["model_ligand"] = checks["file"].map(ids)
            if checks["model_ligand"].isna().any() or checks["position"].ne(0).any():
                raise ValueError("PoseBusters poses differ from the prepared ligands")
            checks["prediction"] = prediction
            results["posebusters"] = checks.drop(
                columns=["file", "molecule"], errors="ignore"
            ).to_dict("records")
        except Exception as exc:
            failure("posebusters", exc)
    return results


def evaluate_predictions(
    predictions: str | Path,
    *,
    output_dir: str | Path,
    release: PlinderRelease | None = None,
    mode: Literal["ligands", "interfaces", "both"] = "both",
    num_workers: int = 1,
    include_all_ligands: bool = False,
    posebusters: bool = True,
    ligand_smiles: Mapping[str, str] | None = None,
    ligand_ccd_codes: Mapping[str, str] | None = None,
    ligand_chains: Sequence[str] = (),
    ligand_options: Sequence[str] = (),
    interface_options: Sequence[str] = (),
    ost_executable: str | Path = "ost",
) -> dict[str, pd.DataFrame]:
    """Evaluate ``predictions/<reference ID>/<model>.cif`` folders.

    A reference ID is a PLINDER system/interface ID, or a four-character PDB ID
    to compare against all its systems/interfaces. Each mmCIF (also .mmcif or
    gzip-compressed) is a separate prediction. ``mode`` selects ligand or
    interface evaluation, or both. References come from the selected release.

    ``num_workers`` limits concurrent prediction processes. PoseBusters runs
    without its own process pool. When calling from a Python script, put the
    call inside ``if __name__ == '__main__':``.

    Ligand rows correspond to proper reference ligands by default; set
    ``include_all_ligands`` to include reference ions/artifacts too. Predictions
    are not filtered by artifact status or receptor proximity. Missing poses
    remain as unassigned reference rows. Preparation/tool errors retain error
    rows and are listed separately, not converted to scientific scores.

    CCD/SMILES mappings and ``ligand_chains`` are passed to ``prepare_prediction``
    for every model. ``ligand_options`` and ``interface_options`` pass additional
    native OST flags to the respective actions.

    Returns tables keyed by ``ligands``, ``interfaces``, ``posebusters``, and
    ``failures``; also writes the first three as Parquet and ``failures.tsv``.
    Ligand metrics keep separate ``rmsd_model_ligand``/``lddt_pli_model_ligand``
    assignments. Join either to PoseBusters on prediction and model_ligand.
    PoseBusters covers all candidate predicted ligands, including unmatched ones.
    Native JSON, CSV, logs and prepared coordinates remain under ``details/``.
    Every call reruns comparisons and replaces the summary tables; old details
    are never read as results.
    """
    if (
        isinstance(num_workers, bool)
        or not isinstance(num_workers, int)
        or num_workers < 1
    ):
        raise ValueError("num_workers must be a positive integer")
    if mode not in {"ligands", "interfaces", "both"}:
        raise ValueError("mode must be ligands, interfaces, or both")
    predictions, output_dir = Path(predictions).resolve(), Path(output_dir).resolve()
    if not predictions.is_dir():
        raise NotADirectoryError(predictions)
    tasks: list[tuple[Path, str, list[_Reference]]] = []
    collected: dict[str, list[dict[str, Any]]] = {
        key: [] for key in ("ligands", "interfaces", "posebusters", "failures")
    }
    for folder in sorted(predictions.iterdir()):
        if not folder.is_dir() or folder.name.startswith("."):
            continue
        models = [
            path
            for path in sorted(folder.iterdir())
            if path.is_file()
            and path.name.lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz"))
        ]
        if not models:
            continue
        refs, failures = _references(
            folder.name, mode, release or PlinderRelease(), include_all_ligands
        )
        collected["failures"].extend(
            {
                "prediction": str(model.relative_to(predictions)),
                "stage": "reference",
                **failure,
            }
            for model in models
            for failure in failures
        )
        tasks.extend(
            (model, str(model.relative_to(predictions)), refs)
            for model in models
            if refs
        )
    if not tasks and not collected["failures"]:
        raise ValueError("No mmCIF predictions found in reference-ID subdirectories")
    output_dir.mkdir(parents=True, exist_ok=True)
    worker = partial(
        _evaluate_one,
        output_dir=output_dir,
        posebusters=posebusters,
        ligand_smiles=dict(ligand_smiles or {}),
        ligand_ccd_codes=dict(ligand_ccd_codes or {}),
        ligand_chains=tuple(ligand_chains),
        ligand_options=tuple(ligand_options),
        interface_options=tuple(interface_options),
        ost_executable=str(ost_executable),
    )
    with ProcessPoolExecutor(
        max_workers=num_workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_worker_threads,
    ) as pool:
        futures = {pool.submit(worker, task): task for task in tasks}
        for future in as_completed(futures):
            try:
                for key, rows in future.result().items():
                    collected[key].extend(rows)
            except Exception as exc:
                _, prediction, refs = futures[future]
                collected["failures"].append(
                    {
                        "prediction": prediction,
                        "system_id": "",
                        "stage": "worker",
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                )
                for ref in refs:
                    for ligand in ref.ligands if ref.kind == "ligands" else [None]:
                        collected[ref.kind].append(
                            {
                                "prediction": prediction,
                                "system_id": ref.system_id,
                                "status": "error",
                                **({"ligand_id": ligand} if ligand is not None else {}),
                            }
                        )
    schemas = {
        "ligands": ["prediction", "system_id", "status", *_LIGAND_COLUMNS],
        "interfaces": ["prediction", "system_id", "status", *_INTERFACE_METRICS],
        "posebusters": ["prediction", "model_ligand", "position"],
        "failures": ["prediction", "system_id", "stage", "error"],
    }
    tables = {}
    with TemporaryDirectory(dir=output_dir, prefix=".tables-") as temporary:
        for key, rows in collected.items():
            columns = list(
                dict.fromkeys(
                    [*schemas[key], *(column for row in rows for column in row)]
                )
            )
            frame = pd.DataFrame(rows).reindex(columns=columns)
            frame = frame.where(frame.notna(), None)
            order = [
                column
                for column in (
                    "prediction",
                    "system_id",
                    "ligand_id",
                    "model_ligand",
                    "stage",
                )
                if column in frame
            ]
            tables[key] = frame.sort_values(order).reset_index(drop=True)
            pending = Path(temporary) / (
                "failures.tsv" if key == "failures" else f"{key}.parquet"
            )
            if key == "failures":
                tables[key].to_csv(pending, sep="\t", index=False)
            else:
                tables[key].to_parquet(pending, index=False)
        for pending in Path(temporary).iterdir():
            pending.replace(output_dir / pending.name)
    return tables
