# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License, Version 2.0
"""Command-line access to public PLINDER workflows."""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from plinder.core.scores.custom import (
    SEARCH_BACKENDS,
    CustomProteinSearchConfig,
    score_custom_cif_files,
)

MMCIF_SUFFIXES = (".cif", ".mmcif", ".cif.gz", ".mmcif.gz")


def _is_mmcif(path: Path) -> bool:
    return path.name.lower().endswith(MMCIF_SUFFIXES)


def _discover_mmcif_files(path: Path, *, recursive: bool) -> list[Path]:
    source = path.resolve()
    if source.is_file():
        if not _is_mmcif(source):
            raise ValueError(f"input is not an mmCIF file: {source}")
        return [source]
    if not source.is_dir():
        raise FileNotFoundError(f"input file or directory does not exist: {source}")
    candidates = source.rglob("*") if recursive else source.iterdir()
    files = sorted(
        candidate.resolve()
        for candidate in candidates
        if candidate.is_file() and _is_mmcif(candidate)
    )
    if not files:
        scope = "recursively" if recursive else "directly"
        raise ValueError(f"no mmCIF files found {scope} under {source}")
    return files


def _component_mapping(values: Sequence[str], *, option: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        component, separator, replacement = value.partition("=")
        component = component.strip()
        replacement = replacement.strip()
        if not separator or not component or not replacement:
            raise ValueError(f"{option} expects COMPONENT=VALUE, got {value!r}")
        if component in result:
            raise ValueError(f"{option} repeats component {component!r}")
        result[component] = replacement
    return result


def _row_count(path: Path | None) -> int | None:
    return pq.ParquetFile(path).metadata.num_rows if path is not None else None


def _default_output_dir(input_path: Path) -> Path:
    name = input_path.resolve().name
    lowered = name.lower()
    for suffix in MMCIF_SUFFIXES:
        if lowered.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return Path.cwd() / f"{name}_plinder_links"


def _run_link(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = (args.output_dir or _default_output_dir(args.input)).resolve()
    cif_files = [
        path
        for path in _discover_mmcif_files(args.input, recursive=args.recursive)
        if not path.is_relative_to(output_dir)
    ]
    if not cif_files:
        raise ValueError("the output directory contains every discovered mmCIF")
    ligand_smiles = _component_mapping(
        args.ligand_smiles,
        option="--ligand-smiles",
    )
    ligand_ccd = _component_mapping(args.ligand_ccd, option="--ligand-ccd")
    overlap = sorted(set(ligand_smiles).intersection(ligand_ccd))
    if overlap:
        raise ValueError(
            "a component cannot have both SMILES and CCD overrides: "
            f"{overlap}"
        )
    include_ligands: bool | None
    if args.mode == "auto":
        include_ligands = None
    else:
        include_ligands = args.mode == "ligand"
    backends = tuple(args.backend or SEARCH_BACKENDS)
    result = score_custom_cif_files(
        cif_files,
        work_dir=output_dir,
        scratch_dir=args.scratch_dir,
        data_dir=args.data_dir,
        structure_mode=args.structure_mode,
        assembly_ids=args.assembly_id or None,
        ligand_smiles_dict=ligand_smiles or None,
        ligand_ccd_code_dict=ligand_ccd or None,
        include_ligands=include_ligands,
        include_interfaces=False if args.skip_interfaces else None,
        include_shape=not args.skip_shape,
        interface_annotate_prodigy=not args.skip_prodigy,
        search_config=CustomProteinSearchConfig(
            evalue=args.evalue,
            sensitivity=args.sensitivity,
            max_seqs=args.max_seqs,
            coverage=args.coverage,
            min_seq_id=args.min_seq_id,
        ),
        backends=backends,
        threads=args.threads,
        shape_score_threads=args.shape_score_threads,
        store_aligned_pocket_residues=args.save_aligned_pocket_residues,
    )
    outputs = {
        "protein_scores": result.protein_scores,
        "ligand_scores": result.ligand_scores,
        "interface_scores": result.interface_scores,
        "aligned_pocket_residues": result.aligned_pocket_residues,
    }
    modes = [
        name.removesuffix("_scores")
        for name, path in outputs.items()
        if path is not None and name.endswith("_scores")
    ]
    summary: dict[str, Any] = {
        "status": "complete",
        "input_files": [str(path) for path in cif_files],
        "output_dir": str(output_dir),
        "public_data_dir": (
            str(args.data_dir.resolve()) if args.data_dir is not None else None
        ),
        "score_modes": modes,
        "outputs": {
            name: (
                {"path": str(path), "rows": _row_count(path)}
                if path is not None
                else None
            )
            for name, path in outputs.items()
        },
    }
    summary_path = output_dir / "link_summary.json"
    summary["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plinder_link",
        description="score custom mmCIF structures against the public PLINDER release",
    )
    parser.add_argument("input", type=Path, help="an mmCIF file or directory")
    parser.add_argument("-o", "--output-dir", type=Path)
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="downloaded public PLINDER release root; otherwise use the release cache",
    )
    parser.add_argument("--scratch-dir", type=Path)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument(
        "--mode",
        choices=["auto", "ligand", "protein"],
        default="auto",
        help="auto scores ligands only when proper ligands are present",
    )
    parser.add_argument(
        "--structure-mode",
        choices=["as_is", "pdb"],
        default="as_is",
    )
    parser.add_argument("--assembly-id", action="append", default=[])
    parser.add_argument(
        "--ligand-smiles",
        action="append",
        default=[],
        metavar="COMPONENT=SMILES",
        help="assign chemistry to a custom component; repeat as needed",
    )
    parser.add_argument(
        "--ligand-ccd",
        action="append",
        default=[],
        metavar="COMPONENT=CCD",
        help="map a custom component to a CCD code; repeat as needed",
    )
    parser.add_argument(
        "--backend",
        action="append",
        choices=list(SEARCH_BACKENDS),
        help="search backend; repeat to select both (the default)",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--shape-score-threads", type=int, default=1)
    parser.add_argument("--max-seqs", type=int, default=10_000)
    parser.add_argument("--evalue", type=float, default=0.01)
    parser.add_argument("--sensitivity", type=float, default=11.0)
    parser.add_argument("--coverage", type=float, default=0.0)
    parser.add_argument("--min-seq-id", type=float, default=0.0)
    parser.add_argument("--skip-shape", action="store_true")
    parser.add_argument("--skip-interfaces", action="store_true")
    parser.add_argument("--skip-prodigy", action="store_true")
    parser.add_argument(
        "--save-aligned-pocket-residues",
        action="store_true",
        help="write the PLINDER/custom residue pairs behind pocket identity scores",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = _run_link(args)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        parser.exit(2, f"plinder_link: error: {exc}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
