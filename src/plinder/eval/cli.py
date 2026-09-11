"""Command-line entry point for folder-based prediction evaluation."""

from __future__ import annotations

import argparse
import json
import shlex
from collections.abc import Sequence
from pathlib import Path

from plinder.core.release import PlinderRelease
from plinder.eval.batch import evaluate_predictions


def _component_mapping(value: str) -> dict[str, str]:
    try:
        mapping = json.loads(value)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected a JSON component mapping: {exc}"
        ) from exc
    if not isinstance(mapping, dict) or any(
        not isinstance(key, str)
        or not key.strip()
        or not isinstance(item, str)
        or not item.strip()
        for key, item in mapping.items()
    ):
        raise argparse.ArgumentTypeError('Use a JSON object such as {"LIG": "ATP"}')
    return mapping


def main(argv: Sequence[str] | None = None) -> int:
    """Return zero on completion, one for recorded errors, two for invalid input."""
    parser = argparse.ArgumentParser(
        description="Evaluate an input table or predictions/<reference ID>/ folder with OpenStructure and PoseBusters. Tables use input_id, structure_path, reference_id and optional ligand_path. Reference IDs may be PLINDER system/interface IDs or four-character PDB IDs."
    )
    parser.add_argument(
        "predictions",
        type=Path,
        help="CSV/TSV/Parquet input table, or a folder with reference-ID subdirectories",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--mode", choices=("ligands", "interfaces", "both"), default="both"
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=1,
        help="Maximum concurrent prediction processes (default: 1)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Local PLINDER release; otherwise use the configured release",
    )
    parser.add_argument(
        "--include-all-ligands",
        action="store_true",
        help="Include reference ions/artifacts, not only proper ligands",
    )
    parser.add_argument(
        "--no-posebusters", action="store_true", help="Skip PoseBusters checks"
    )
    parser.add_argument(
        "--ligand-smiles",
        type=_component_mapping,
        help='JSON component-to-SMILES mapping, e.g. {"LIG": "CCO"}',
    )
    parser.add_argument(
        "--ligand-ccd-codes",
        type=_component_mapping,
        help='JSON component-to-CCD mapping, e.g. {"LIG": "ATP"}',
    )
    parser.add_argument(
        "--ligand-chain",
        action="append",
        default=[],
        help="Also treat this polymer label asym ID as a ligand; repeat for multiple chains",
    )
    parser.add_argument(
        "--ligand-options",
        type=shlex.split,
        default=(),
        help="Additional OST ligand flags as one quoted string; use --ligand-options='--flag value'",
    )
    parser.add_argument(
        "--interface-options",
        type=shlex.split,
        default=(),
        help="Additional OST interface flags as one quoted string; use --interface-options='--flag value'",
    )
    parser.add_argument(
        "--ost-executable",
        default="ost",
        help="Path to the OST executable (default: ost)",
    )
    args = parser.parse_args(argv)
    try:
        tables = evaluate_predictions(
            args.predictions,
            output_dir=args.output_dir,
            release=PlinderRelease(data_dir=args.data_dir)
            if args.data_dir is not None
            else None,
            mode=args.mode,
            num_workers=args.num_workers,
            include_all_ligands=args.include_all_ligands,
            posebusters=not args.no_posebusters,
            ligand_smiles=args.ligand_smiles,
            ligand_ccd_codes=args.ligand_ccd_codes,
            ligand_chains=args.ligand_chain,
            ligand_options=args.ligand_options,
            interface_options=args.interface_options,
            ost_executable=args.ost_executable,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    for name, table in tables.items():
        print(f"{name}: {len(table)} rows")
    print(f"Results: {args.output_dir.resolve()}")
    return int(not tables["failures"].empty)


if __name__ == "__main__":
    raise SystemExit(main())
