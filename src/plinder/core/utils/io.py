# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
Wrap all network requests in a retry decorator
and use a convention of looking for a file in a
pre-determined location before fetching it from
the network.
"""

from functools import wraps
from pathlib import Path
from time import sleep
from typing import Any, Callable, Optional, TypeVar

import requests
from biotite.database.rcsb import fetch
from biotite.structure.io.pdbx import CIFFile, get_structure, set_structure

from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)
T = TypeVar("T")


def retry(func: Callable[..., T]) -> Callable[..., T]:
    @wraps(func)
    def inner(*args: Any, **kwargs: Any) -> T:
        name = func.__name__
        mod = func.__module__
        log = setup_logger(".".join([mod, name]))
        retries = 5
        exc = None
        for i in range(1, retries + 1):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                wait = 2**i
                log.error(f"failed: {repr(e)}, retry in: {wait}s")
                exc = e
                sleep(wait)
        raise Exception(f"Timeout error {exc}")

    return inner


@retry
def download_alphafold_cif_file(
    uniprot_id: str,
    output_folder: Path,
    url: str = "https://alphafold.ebi.ac.uk/files",
    force_update: bool = False,
) -> Optional[Path]:
    cif_file_path = output_folder / f"AF-{uniprot_id}-F1-model_v4.cif"
    if not cif_file_path.is_file() or force_update:
        resp = requests.get(f"{url}/{cif_file_path.name}")
        if resp.status_code == 404:
            LOG.info(f"UniProt ID {uniprot_id} not in AlphaFold database")
            return None
        resp.raise_for_status()
        with open(cif_file_path, "w") as f:
            f.write(resp.text)
    return cif_file_path


@retry
def download_pdb_chain_cif_file(pdb_id: str, chain_id: str, filename: Path) -> Path:
    structure = get_structure(
        CIFFile.read(
            fetch(
                pdb_ids=pdb_id,
                format="cif",
                overwrite=False,
            )
        ),
        model=1,
        use_author_fields=False,
    )
    write_file = CIFFile()
    set_structure(write_file, structure[structure.chain_id == chain_id])
    write_file.write(filename.as_posix())
    return filename
