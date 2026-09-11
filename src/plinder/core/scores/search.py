"""Search PLINDER using an input table, protein structures, or sequences."""

from __future__ import annotations

import shutil
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Literal

import pandas as pd

from plinder.core.release import PlinderRelease
from plinder.core.scores.custom import (
    CIF_SEARCH_BACKENDS,
    CustomProteinSearchConfig,
    CustomScoringResult,
    CustomSequenceScoringResult,
    score_custom_cif_files,
    score_custom_sequence_file,
)
from plinder.core.structure.inputs import read_structure_table, write_search_structure


def search(
    inputs: str | Path | pd.DataFrame,
    *,
    output_dir: str | Path,
    release: PlinderRelease | None = None,
    mode: Literal["auto", "pockets", "ligands", "interfaces", "both"] = "auto",
    backends: Iterable[str] | None = None,
    threads: int = 1,
    search_config: CustomProteinSearchConfig | None = None,
    plinder_entry_ids: Iterable[str] | None = None,
    store_aligned_pocket_residues: bool = False,
) -> CustomScoringResult | CustomSequenceScoringResult:
    """Search with the same structure table accepted by evaluation.

    Tables contain ``input_id``, ``structure_path`` and optional ``ligand_path``.
    A receptor and its SDFs share an input ID; paths are relative to the table.
    Complete mmCIFs use an empty ligand_path. ``reference_id`` is ignored here.

    For ligand-free queries, pass a FASTA, one PDB/mmCIF, or a folder of those
    structures. FASTA searches use MMseqs; structure searches use Foldseek and
    MMseqs by default. ``mode='pockets'`` compares protein/pocket features,
    ``'ligands'`` includes ligand similarities, and ``'interfaces'`` includes
    protein interfaces. ``'both'`` includes both; ``'auto'`` uses applicable
    features (ligand comparisons require a table).

    PDB receptors contribute amino-acid residues, including modified amino acids;
    other PDB components are omitted. Supply ligand poses through ``ligand_path``.

    The returned result contains paths to score tables and alignment files.
    Input IDs become custom structure IDs in these outputs. Generated ligand
    chain IDs map to the original SDF paths in ``inputs/<input_id>.ligands.tsv``.
    """
    if mode not in {"auto", "pockets", "ligands", "interfaces", "both"}:
        raise ValueError("mode must be auto, pockets, ligands, interfaces or both")
    output_dir = Path(output_dir).resolve()
    common: dict[str, Any] = dict(
        work_dir=output_dir,
        data_dir=release.data_dir if release is not None else None,
        threads=threads,
        search_config=search_config,
        plinder_entry_ids=plinder_entry_ids,
        store_aligned_pocket_residues=store_aligned_pocket_residues,
    )
    path = None if isinstance(inputs, pd.DataFrame) else Path(inputs).resolve()
    if path is not None and path.name.lower().endswith((
        ".fasta",
        ".fa",
        ".faa",
        ".fasta.gz",
        ".fa.gz",
        ".faa.gz",
    )):
        if mode not in {"auto", "pockets"}:
            raise ValueError(
                "FASTA inputs support pocket search; use structures for ligand/interface comparisons"
            )
        return score_custom_sequence_file(
            path,
            backends=tuple(backends) if backends is not None else ("mmseqs",),
            **common,
        )
    is_table = path is None or path.suffix.lower() in {".csv", ".tsv", ".parquet"}
    if is_table:
        structures = read_structure_table(inputs)
    else:
        if mode in {"ligands", "both"}:
            raise ValueError("Use an input table for ligand searches")
        assert path is not None
        files = sorted(path.iterdir()) if path.is_dir() else [path]
        rows = []
        for file in files:
            if file.is_file() and file.name.lower().endswith((
                ".pdb",
                ".pdb.gz",
                ".cif",
                ".cif.gz",
                ".mmcif",
                ".mmcif.gz",
            )):
                name = file.stem if file.suffix.lower() == ".gz" else file.name
                rows.append({"input_id": Path(name).stem, "structure_path": str(file)})
        if not rows:
            raise ValueError(f"No PDB/mmCIF structures found in {path}")
        if len({row["input_id"] for row in rows}) != len(rows):
            raise ValueError("Structure filenames must have unique basenames")
        structures = read_structure_table(pd.DataFrame(rows))
    include_ligands = is_table and mode in {"auto", "ligands", "both"}
    sources = []
    for input_id, structure in structures.items():
        original = Path(structure.coordinates)
        destination = output_dir / "inputs" / f"{input_id}.cif"
        if structure.ligand_sdfs is None and original.suffix.lower() in {
            ".cif",
            ".mmcif",
        }:
            destination.parent.mkdir(parents=True, exist_ok=True)
            if original.resolve() != destination.resolve():
                shutil.copyfile(original, destination)
        else:
            write_search_structure(
                structure, destination, include_ligands=include_ligands
            )
        sources.append(destination)
    return score_custom_cif_files(
        sources,
        backends=tuple(backends)
        if backends is not None
        else (("mmseqs",) if plinder_entry_ids is not None else CIF_SEARCH_BACKENDS),
        include_ligands=None if mode == "auto" and include_ligands else include_ligands,
        include_interfaces=None if mode == "auto" else mode in {"interfaces", "both"},
        **common,
    )
