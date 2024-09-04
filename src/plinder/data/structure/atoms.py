# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
import os
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

import biotite.structure as struc
from biotite import TextFile
from biotite.structure import connect_via_residue_names
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from rdkit import Chem
from rdkit.Chem import AllChem, Mol

from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)

# TODO: this is likely a redundant function - only used for atom_array_from_cif_file below
def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import PDBxFile

    return PDBxFile

# TODO: this is likely a redundant function - only used for interface gaps
def atom_array_from_cif_file(
    structure: Path | AtomArray, use_author_fields: bool = True
) -> AtomArray | None:
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_ciffile()
        try:
            if structure.suffix == ".gz":
                with gzip.open(str(structure), "rt", encoding="utf-8") as f:
                    mod = reader.read(f)
            else:
                mod = reader.read(structure)
            arr = get_structure(mod, model=1, use_author_fields=use_author_fields)  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure
