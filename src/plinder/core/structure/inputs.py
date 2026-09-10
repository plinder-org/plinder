"""Structure input files for search and evaluation workflows."""

from __future__ import annotations

import gzip
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    import biotite.structure as struc
    from biotite.structure.io import pdbx
    from rdkit.Chem import Mol


@dataclass(frozen=True)
class StructureInput:
    """A complete structure, or a receptor with separate ligand poses.

    ``coordinates`` is a PDB or mmCIF file, optionally gzip-compressed.
    ``ligand_sdfs=None`` identifies a complete structure. A supplied sequence
    identifies separate ligand poses, one molecule per SDF, in the receptor's
    coordinate frame. An empty sequence represents zero predicted ligand poses.

    Examples
    --------
    ``StructureInput("complex.cif")``
    ``StructureInput("receptor.pdb", ligand_sdfs=["pose.sdf"])``
    """

    coordinates: str | Path
    ligand_sdfs: Sequence[str | Path] | None = None
    reference_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "coordinates", Path(self.coordinates))
        if self.ligand_sdfs is not None:
            if isinstance(self.ligand_sdfs, (str, Path)):
                raise TypeError("ligand_sdfs must be a sequence of SDF paths")
            object.__setattr__(
                self, "ligand_sdfs", tuple(Path(path) for path in self.ligand_sdfs)
            )

    @classmethod
    def from_path(cls, coordinates: str | Path) -> StructureInput:
        """Pair coordinates with ``model.sdf`` or ``model.ligands/*.sdf``.

        With no matching SDF or ligand folder, use the complete-structure form.
        Matching uses the coordinate basename, including for compressed files.
        """
        model = Path(coordinates)
        name = model.stem if model.suffix.lower() == ".gz" else model.name
        stem = Path(name).stem
        folder = model.parent / f"{stem}.ligands"
        singles = sorted(
            path
            for path in model.parent.iterdir()
            if path.is_file() and path.stem == stem and path.suffix.lower() == ".sdf"
        )
        if len(singles) > 1 or (singles and folder.exists()):
            raise ValueError(f"Use either {stem}.sdf or {stem}.ligands/ for {model}")
        if folder.exists():
            if not folder.is_dir():
                raise NotADirectoryError(folder)
            return cls(
                model,
                ligand_sdfs=sorted(
                    path
                    for path in folder.iterdir()
                    if path.is_file() and path.suffix.lower() == ".sdf"
                ),
            )
        return cls(model, ligand_sdfs=singles or None)


def read_structure_table(
    source: str | Path | pd.DataFrame, *, require_reference: bool = False
) -> dict[str, StructureInput]:
    """Read CSV/TSV/Parquet inputs, grouped by ``input_id``.

    Required columns are ``input_id`` and ``structure_path``. ``ligand_path``
    holds one SDF or a folder of SDFs; repeat the input row for multiple ligands.
    Evaluation also requires ``reference_id``. Relative file paths are resolved
    from the table's directory (the working directory for a DataFrame).
    """
    if isinstance(source, pd.DataFrame):
        frame, base = source.copy(), Path.cwd()
    else:
        path = Path(source).resolve()
        base = path.parent
        if path.suffix.lower() == ".parquet":
            frame = pd.read_parquet(path)
        elif path.suffix.lower() in {".csv", ".tsv"}:
            frame = pd.read_csv(
                path,
                sep="\t" if path.suffix.lower() == ".tsv" else ",",
                dtype=str,
                keep_default_na=False,
            )
        else:
            raise ValueError("Input tables must be CSV, TSV or Parquet")
    required = {"input_id", "structure_path"}
    if require_reference:
        required.add("reference_id")
    if missing := required - set(frame.columns):
        raise ValueError(f"Missing input table columns: {sorted(missing)}")
    if frame.empty:
        raise ValueError("Input table is empty")
    for column in required | {"ligand_path", "reference_id"}:
        if column not in frame:
            frame[column] = ""
        frame[column] = frame[column].fillna("").astype(str).str.strip()
        if column in required and frame[column].eq("").any():
            raise ValueError(f"Input table has empty {column} values")
    inputs = {}
    for input_id, rows in frame.groupby("input_id", sort=False):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", input_id):
            raise ValueError(
                f"Invalid input_id: {input_id!r}; use letters, digits, _, . or -"
            )
        for column in ("structure_path", "reference_id"):
            if rows[column].nunique() != 1:
                raise ValueError(f"Conflicting {column} values for {input_id}")
        ligands: list[Path] = []
        supplied = rows.loc[rows["ligand_path"].ne(""), "ligand_path"]
        for name in supplied:
            path = (base / name).resolve()
            ligands.extend(
                sorted(
                    p
                    for p in path.iterdir()
                    if p.is_file() and p.suffix.lower() == ".sdf"
                )
                if path.is_dir()
                else [path]
            )
        if len(set(ligands)) != len(ligands):
            raise ValueError(f"Repeated ligand files for {input_id}")
        inputs[input_id] = StructureInput(
            (base / rows["structure_path"].iloc[0]).resolve(),
            ligand_sdfs=ligands if len(supplied) else None,
            reference_id=rows["reference_id"].iloc[0] or None,
        )
    return inputs


def read_ligand_sdf(path: str | Path) -> Mol:
    """Read one supplied pose, preserving its coordinates and bond information."""
    import numpy as np
    from rdkit.Chem import SDMolSupplier

    path = Path(path)
    if path.suffix.lower() != ".sdf":
        raise ValueError(f"Expected a ligand SDF: {path}")
    supplier = SDMolSupplier(str(path), sanitize=False, removeHs=False)
    if len(supplier) != 1 or (molecule := supplier[0]) is None:
        raise ValueError(f"Expected one readable ligand per SDF: {path}")
    if (
        molecule.GetNumAtoms() == 0
        or molecule.GetNumConformers() == 0
        or not np.isfinite(molecule.GetConformer().GetPositions()).all()
    ):
        raise ValueError(f"Ligand SDF needs atomic coordinates: {path}")
    return molecule


def read_input_structure(
    path: str | Path
) -> tuple[struc.AtomArray, pdbx.CIFBlock | None]:
    """Read the first coordinate model and its mmCIF metadata where available."""
    from biotite.structure.io import pdb

    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )

    path = Path(path)
    if path.name.lower().endswith((".pdb", ".pdb.gz")):
        with (
            gzip.open(path, "rt") if path.suffix.lower() == ".gz" else path.open("r")
        ) as stream:
            return pdb.PDBFile.read(stream).get_structure(
                model=1, include_bonds=True, extra_fields=["charge"]
            ), None
    if not path.name.lower().endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        raise ValueError(f"Expected PDB or mmCIF coordinates: {path}")
    cif = read_mmcif_file(path)
    return get_structure_with_altloc(cif, include_bonds=True), cif.block


def write_search_structure(
    structure: StructureInput, destination: Path, *, include_ligands: bool
) -> Path:
    """Write a search mmCIF with the supplied SDF chemistry and pose coordinates.

    PDB receptors contribute amino-acid residues, including modified amino acids.
    Other PDB components are omitted; ligand poses can be supplied as SDFs.
    """
    import biotite.structure as struc
    import numpy as np
    from biotite.interface.rdkit import from_mol
    from biotite.structure.io import pdbx
    from rdkit import Chem

    from plinder.data.annotations.save_utils import save_cif_file

    atoms, block = read_input_structure(structure.coordinates)
    if block is None:
        if struc.filter_nucleotides(atoms).any():
            raise ValueError(
                "Use mmCIF with polymer sequence metadata for DNA/RNA search inputs"
            )
        # PDB author chains may mix protein, water and other components. Only
        # protein residues belong in the sequence of the generated polymer entity.
        # CCD-based selection retains modified amino acids recorded as HETATM.
        atoms = atoms[struc.filter_amino_acids(atoms)]
        if len(atoms) == 0:
            raise ValueError(
                f"No protein residues found in PDB receptor: {structure.coordinates}"
            )
    author_numbers = atoms.res_id.copy()
    if block is None:
        atoms.res_id = struc.create_continuous_res_ids(atoms, restart_each_chain=True)
    if "charge" not in atoms.get_annotation_categories():
        atoms.set_annotation("charge", np.zeros(len(atoms), dtype=int))
    mapping = {}
    used_chains, used_components = set(atoms.chain_id), set(atoms.res_name)
    if include_ligands and structure.ligand_sdfs is not None:
        for index, ligand_path in enumerate(structure.ligand_sdfs, start=1):
            molecule = read_ligand_sdf(ligand_path)
            Chem.SanitizeMol(molecule)
            ligand = from_mol(molecule, conformer_id=0, add_hydrogen=False)
            counter = index
            while f"L{counter:04d}" in used_chains | used_components:
                counter += 1
            name = f"L{counter:04d}"
            used_chains.add(name)
            used_components.add(name)
            ligand.chain_id = np.full(len(ligand), name)
            ligand.res_name = np.full(len(ligand), name)
            ligand.res_id[:] = 1
            ligand.ins_code[:] = ""
            ligand.hetero[:] = True
            ligand.atom_name = struc.create_atom_names(ligand)
            author_numbers = np.concatenate((author_numbers, ligand.res_id))
            atoms += ligand
            mapping[name] = str(Path(ligand_path).resolve())
    destination.parent.mkdir(parents=True, exist_ok=True)
    save_cif_file(
        atoms,
        destination.stem,
        destination,
        source_block=block,
        source_asym_ids={chain: chain for chain in used_chains},
    )
    if block is None:
        cif = pdbx.CIFFile.read(destination)
        cif.block["atom_site"]["auth_seq_id"] = author_numbers
        cif.write(destination)
    if mapping:
        pd.DataFrame(
            [
                {"chain_id": chain, "ligand_path": value}
                for chain, value in mapping.items()
            ]
        ).to_csv(destination.with_suffix(".ligands.tsv"), sep="\t", index=False)
    return destination
