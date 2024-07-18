# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
import subprocess as sp
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional

import pandas as pd

from plinder.core.utils.log import setup_logger

if TYPE_CHECKING:
    from plinder.data.utils.annotations.aggregate_annotations import Entry

LOG = setup_logger(__name__)


def run(cmd: list[str]) -> None:
    LOG.info(" ".join(cmd))
    sp.check_call(cmd, stderr=sp.STDOUT)


def make_db(
    input_dir: Path,
    output_dir: Path,
    db: str,
    tmp_dir: Path = Path("tmp"),
) -> None:
    """
    Create full databases. input_dir is a misnomer
    for mmseqs because it wants the full path to the
    pdb_seqres.txt.gz or pdb_uniprot.fasta. output_dir is a misnomer in both
    cases because it wants the full path to the filename
    of the database, not just the directory it writes to.

    Parameters
    ----------
    input_dir : Path
        location of input files for database
        (foldseek: input_dir = adir [/ **/*-enrich.cif.gz] for apo/holo and [/AF-*-F1-model_v4.cif] for pred)
        (mmseqs: input_dir = adir / seqres / pdb_seqres.txt.gz for apo/holo and / uniprot / pdb_uniprot.txt.gz for pred)
    output_dir : Path
        location of full final database (including file name)
    db : str
        name of db in ["foldseek", "mmseqs"]
    tmp_dir : Path, default=Path("tmp")
        scratch directory
    """
    full_db = str(output_dir / db)
    output_dir.mkdir(exist_ok=True, parents=True)
    cmd = [
        db,
        "createdb",
        str(input_dir),
        full_db,
    ]
    if db == "foldseek":
        cmd.extend(["--chain-name-mode", "1"])
    run(cmd)
    run(
        [
            db,
            "createindex",
            full_db,
            str(tmp_dir) + f"_{db}",
        ]
    )


def make_sub_db(
    entry_chain_ids: set[str],
    full_db: Path,
    sub_db: Path,
    aln_type: str,
) -> list[str]:
    """
    Requires a list of already generated annotation entries.
    Create reduced indexes from the full databases based
    on the passed entry_chain_ids.

    Parameters
    ----------
    entry_chain_ids : set[str]
        output from get_db_ids
    full_db : Path
        path to full database
    sub_db : Path
        path to sub database
    aln_type : str
        database name
    """
    # Map chain ids to full database ids
    sub_db_lookup_file = sub_db / "ids.tsv"
    found = set()
    with open(full_db.parent / f"{full_db.name}.lookup", "r") as f:
        with open(sub_db_lookup_file, "w") as fw:
            for line in f:
                fields = line.split()
                if fields[1] in entry_chain_ids:
                    fw.write(line)
                    found.add(fields[1])
    # Create subdb
    if aln_type == "foldseek":
        run(
            [
                "foldseek",
                "createsubdb",
                str(sub_db_lookup_file),
                str(full_db),
                str(sub_db / sub_db.name),
                "--subdb-mode",
                "1",
            ]
        )
        run(
            [
                "foldseek",
                "createsubdb",
                str(sub_db_lookup_file),
                f"{full_db}_ss",
                f"{sub_db / sub_db.name}_ss",
                "--subdb-mode",
                "1",
            ]
        )
        run(
            [
                "foldseek",
                "createsubdb",
                str(sub_db_lookup_file),
                f"{full_db}_ca",
                f"{sub_db / sub_db.name}_ca",
                "--subdb-mode",
                "1",
            ]
        )
    elif aln_type == "mmseqs":
        run(
            [
                "mmseqs",
                "createsubdb",
                str(sub_db_lookup_file),
                str(full_db),
                str(sub_db / sub_db.name),
            ]
        )

    missing = set(entry_chain_ids) - found
    return list(missing)


def get_ids_in_db(data_dir: Path, search_db: str, aln_type: str) -> pd.DataFrame:
    sub_db_file = data_dir / "dbs" / "subdbs" / f"{search_db}_{aln_type}" / "ids.tsv"
    return pd.read_csv(sub_db_file, sep="\t", header=None)


def get_db_ids(
    entries: Dict[str, "Entry"],
    search_db: str,
    aln_type: str,
    entry_ids: Optional[list[str]] = None,
) -> set[str]:
    """
    get chain IDs for a given database from
    the provided entries.

    Parameters
    ----------
    entries : Dict[str, Entry]
        the collection of entries
    search_db : str
        search_db in ["apo", "holo", "pred"]
    aln_type : str
        aln_type in ["foldseek", "mmseqs"]

    Returns
    -------
    db_ids : set[str]
        db IDs in the given database
    """
    db_ids = set()
    if entry_ids is None:
        entry_ids = list(entries.keys())
    for entry_id in entry_ids:
        for chain in entries[entry_id].chains_for_alignment(search_db, aln_type):
            db_ids.add(chain)
    return db_ids


def make_sub_dbs(
    db_dir: Path,
    db_sources: Dict[str, Path],
    entries: Dict[str, "Entry"],
) -> None:
    """
    Create the apo/holo subdbs for score
    generation.

    Parameters
    ----------
    db_dir : Path
        directory where subdbs get written
    db_sources : Dict[str, Path]
        map of database name to path to full database
    entries : Dict[str, Entry]
        map of all the entries
    """

    db_dir.mkdir(exist_ok=True)
    report = {}
    for search_db_aln_type, full_db in db_sources.items():
        search_db, aln_type = search_db_aln_type.split("_")
        subdb = db_dir / search_db_aln_type
        subdb.mkdir(exist_ok=True)
        report[search_db_aln_type] = make_sub_db(
            get_db_ids(entries, search_db, aln_type),
            full_db,
            subdb,
            aln_type,
        )
    with (db_dir / "missing.json").open("w") as f:
        json.dump(report, f)
