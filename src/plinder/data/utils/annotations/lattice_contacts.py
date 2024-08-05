# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import multiprocessing
import os
import re
import subprocess
from pathlib import Path

import pandas as pd

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)

# Need to have pymol env installed
# works with Python 3.8.13 (datascience-1.0__1 environment)
pymol_install_instructions = """
# NOTE: this might need sudo permissions
apt-get update
apt install -y libgl1-mesa-glx
bash
conda create -n pymol -c schrodinger pymol-bundle -y
conda activate pymol
# test import
python -c 'import pymol'
"""
PYMOL_PYTHON = os.environ["CONDA_EXE"].split("bin")[0] + "envs/pymol/bin/python"
SCRIPT_NAME = "pymol_lattice_contacts_script.py"
NEWLINE = "\n"


def pymol_lattice_contacts_process(
    system_id: str,
    ligand_sdf: Path,
    receptor_cif: Path,
    output_file: Path,
    cutoff: float = 4,
) -> None:
    # run through a subprocess to avoid memory problems with PyMOL.
    proc = subprocess.Popen(
        f"{PYMOL_PYTHON} -u {SCRIPT_NAME} {system_id} {ligand_sdf} {receptor_cif} {output_file} {cutoff}",
        shell=True,
        stderr=subprocess.PIPE,
    )
    if proc.stderr:
        stderr = proc.stderr.read().decode()
    # check if everything went ok with the processing
    if bool(re.search("ImportError", stderr)):
        raise ImportError(
            f"""Check if PYMOL_PYTHON={PYMOL_PYTHON} is correctly installed: {stderr}"""
        )
    elif bool(re.search("Error", stderr)):
        raise ValueError(f"Processing error occured:{NEWLINE}{stderr}")


def run_system_lattice_contacts(
    system_id: str,
    data_dir: Path,
    output_dir: Path,
    cutoff: float = 4,
    overwrite: bool = False,
) -> None:
    pdb_id = system_id.split("_")[0]
    system_dir = data_dir / "systems" / pdb_id[-3:-1] / system_id
    receptor_cif = system_dir / "receptor.cif"
    ligand_dir = system_dir / "ligand_files"
    for lig_sdf in ligand_dir.glob("*.sdf"):
        contact_filename = f"{system_id}_{lig_sdf.name.split('.sdf')[0]}_{cutoff}.csv"
        output_file = output_dir / "lattice_contacts_tmp" / contact_filename
        if overwrite or not output_file.exists():
            try:
                # run the pymol process for one system
                pymol_lattice_contacts_process(
                    system_id, lig_sdf, receptor_cif, output_file, cutoff=cutoff
                )
                if not output_file.exists():
                    raise ValueError(f"Cannot find output {output_file}")
            except Exception as e:
                LOG.info(f"Error running pymol process for {system_id}:{NEWLINE}{e}")


def get_lattice_contacts_annotations(
    df: pd.DataFrame,
    data_dir: Path,
    output_dir: Path,
    cutoff: float,
) -> pd.DataFrame:
    num_cont_lig_list = []
    num_cont_symm_list = []
    for system_id in df["system_id"]:
        sys_num_cont_lig = 0
        sys_num_cont_symm = 0
        for output_file in output_dir.glob(
            f"lattice_contacts_tmp/{system_id}_*_{cutoff}.csv"
        ):
            # NOTE: for multi ligand systems, we sum up the number of contacts
            num_cont_lig, num_cont_symm = pd.read_csv(
                output_file, names=["num_cont_lig", "num_cont_symm"]
            ).values[0]
            sys_num_cont_lig += num_cont_lig
            sys_num_cont_symm += num_cont_symm
        num_cont_lig_list.append(sys_num_cont_lig)
        num_cont_symm_list.append(sys_num_cont_symm)
    df["num_cont_lig"] = num_cont_lig_list
    df["num_cont_symm"] = num_cont_symm_list
    return df


def run_annotations_for_split(
    split_file: Path,
    data_dir: Path,
    output_dir: Path,
    cutoff: float = 4,
    test_label: str = "test",
    num_processes: int = 8,
    overwrite: bool = False,
) -> None:
    if split_file.name.endswith(".csv"):
        df = pd.read_csv(split_file)
    else:
        df = pd.read_parquet(split_file)
    assert all(x in df.columns for x in ["split", "system_id"])

    output_dir.mkdir(exist_ok=True)
    # adding a place to write intermediate files
    output_tmp_dir = output_dir / "lattice_contacts_tmp"
    output_tmp_dir.mkdir(exist_ok=True)
    # run the pymol process for all test systems
    df = df[df["split"] == test_label].copy()

    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(num_processes) as p:
        p.starmap(
            run_system_lattice_contacts,
            [
                (system_id, data_dir, output_dir, cutoff, overwrite)
                for system_id in df.system_id
            ],
        )

    df = get_lattice_contacts_annotations(df, data_dir, output_dir, cutoff)
    df.to_parquet(output_dir / f"lattice_contacts_{test_label}.parquet")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Annotate test set for lattice contacts"
    )
    parser.add_argument(
        "--split_file",
        type=Path,
        help="Path to split file with [system_id, split] as columns",
    )
    parser.add_argument(
        "--data_dir",
        type=Path,
        help="Path to plinder data",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Path to output folder where lattice contact data are saved",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=4.0,
        help="Cutoff for lattice contacts",
    )
    parser.add_argument(
        "--test_label",
        type=str,
        default="test",
        help="split=<test_label> is used to get test systems",
    )
    parser.add_argument(
        "--num_processes",
        type=int,
        default=8,
        help="Number of processes to use",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite crystal contacts files",
    )

    args = parser.parse_args()

    run_annotations_for_split(
        split_file=args.split_file,
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        cutoff=args.cutoff,
        test_label=args.test_label,
        num_processes=args.num_processes,
        overwrite=args.overwrite,
    )


# python lattice_contacts.py --split_file split.csv --data_dir ~/plinder_local_data/v2/ --output_dir .

if __name__ == "__main__":
    main()
