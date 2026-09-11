"""Run OpenStructure and PoseBusters on prepared structures.

Results retain the tools' field names and assignments. Reference selection and
PLINDER's proper-ligand selection belong to the calling workflow.
"""

from __future__ import annotations

import json
import subprocess
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any, Literal

import pandas as pd

if TYPE_CHECKING:
    from rdkit.Chem import Mol


@contextmanager
def _run(command: list[str], output: Path) -> Iterator[Path]:
    """Keep a tool log and replace results only after successful parsing."""
    output.parent.mkdir(parents=True, exist_ok=True)
    log = output.with_suffix(output.suffix + ".log")
    with TemporaryDirectory(dir=output.parent, prefix=".evaluation-") as temporary:
        pending = Path(temporary) / output.name
        with log.open("w") as stream:
            completed = subprocess.run(
                [*command, "--output", str(pending)],
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode or not pending.is_file():
            raise RuntimeError(
                f"{command[0]} failed (exit {completed.returncode}); see {log}: "
                f"{log.read_text()[-2000:]}"
            )
        yield pending
        pending.replace(output)


def run_openstructure(
    model: str | Path,
    reference: str | Path,
    output: str | Path,
    *,
    action: Literal["compare-structures", "compare-ligand-structures"],
    options: Sequence[str] = (),
    executable: str | Path = "ost",
) -> dict[str, Any]:
    """Run an OST comparison and return its complete JSON result.

    ``model`` and ``reference`` are coordinate files accepted by OST. Ligand
    comparisons include lDDT-PLI, binding-site-superposed RMSD, backbone RMSD
    (``bb_rmsd``) and lDDT-LP. Structure comparisons include lDDT, iLDDT,
    QS-score and DockQ, including penalties for unmapped reference interfaces.
    ``options`` passes additional CLI arguments directly
    to OST, for example ``["--model-ligands", "pose.sdf"]``. Embedded mmCIF
    ligands use OST's own detection and chemistry rules.

    The JSON file and a sibling ``.log`` file are retained. Failed comparisons
    raise an error, including when OST reports FAILURE with a zero exit code.
    Requires the ``ost`` executable from Bioconda's OpenStructure package.
    """
    metrics = {
        "compare-structures": ["--lddt", "--ilddt", "--qs-score", "--dockq"],
        "compare-ligand-structures": ["--lddt-pli", "--rmsd"],
    }
    if action not in metrics:
        raise ValueError(f"Unsupported OST action: {action}")
    command = [
        str(executable),
        action,
        "--model",
        str(Path(model).resolve()),
        "--reference",
        str(Path(reference).resolve()),
        *metrics[action],
        *options,
    ]
    with _run(command, Path(output)) as pending:
        result = json.loads(pending.read_text())
        if not isinstance(result, dict) or result.get("status") != "SUCCESS":
            raise RuntimeError(f"OpenStructure comparison failed: {result}")
    return result


def run_posebusters(
    ligands: Sequence[str | Path],
    receptor: Mol,
    output: str | Path,
    *,
    full_report: bool = False,
) -> pd.DataFrame:
    """Check SDF poses against a receptor RDKit molecule with residue metadata.

    Create the receptor using ``biotite.interface.rdkit.to_mol(receptor_atoms)``.
    Its residue names and hetero flags let PoseBusters distinguish receptor,
    cofactor and water atoms. An SDF round trip loses that distinction.
    Ligands must have lowercase ``.sdf`` extensions. The returned table and CSV
    retain PoseBusters' columns, file names and SDF record positions.
    Requires ``pip install plinder[eval]``.

    Failed plausibility checks are returned as False. A missing/unreadable
    result or an unloaded receptor raises an error, because intermolecular
    checks would otherwise be absent.
    """
    from posebusters import PoseBusters
    from rdkit.Chem import Mol

    if not isinstance(receptor, Mol):
        raise TypeError("receptor must be an RDKit molecule with residue metadata")
    if receptor.GetNumAtoms() == 0 or any(
        atom.GetPDBResidueInfo() is None for atom in receptor.GetAtoms()
    ):
        raise ValueError("receptor atoms must have residue metadata for PoseBusters")
    if isinstance(ligands, (str, Path)) or not ligands:
        raise ValueError("Provide a nonempty list of ligand SDF paths")
    paths = [Path(path).resolve() for path in ligands]
    if any(path.suffix != ".sdf" for path in paths):
        raise ValueError("PoseBusters ligand files must have the .sdf extension")
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    result = (
        PoseBusters(config="dock", max_workers=0)
        .bust(mol_pred=paths, mol_cond=receptor, full_report=full_report)
        .reset_index()
    )
    required = {"file", "position", "mol_pred_loaded", "mol_cond_loaded"}
    if result.empty or not required.issubset(result.columns):
        raise RuntimeError("PoseBusters returned incomplete results")
    if set(result["file"].astype(str)) != {str(path) for path in paths}:
        raise RuntimeError("PoseBusters did not report all requested ligand files")
    if not result["mol_cond_loaded"].eq(True).all():
        raise RuntimeError("PoseBusters could not load the receptor")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=".evaluation-") as temporary:
        pending = Path(temporary) / output.name
        result.to_csv(pending, index=False)
        pending.replace(output)
    return result
