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
from plinder.core.structure.inputs import (
    find_structure_files,
    read_structure_table,
    write_search_structure,
)


def search(
    inputs: str | Path | pd.DataFrame,
    *,
    output_dir: str | Path,
    release: PlinderRelease | None = None,
    mode: Literal["auto", "all", "pockets", "ligands", "interfaces"] = "auto",
    backends: Iterable[str] | None = None,
    threads: int = 1,
    search_config: CustomProteinSearchConfig | None = None,
    plinder_entry_ids: Iterable[str] | None = None,
    store_aligned_pocket_residues: bool = False,
    include_monomers: bool = True,
) -> CustomScoringResult | CustomSequenceScoringResult:
    """Search with the same structure table accepted by evaluation.

    Tables contain ``input_id``, ``structure_path`` and optional ``ligand_path``.
    A receptor and its SDFs share an input ID; paths are relative to the table.
    Complete mmCIFs use an empty ligand_path. ``reference_id`` is ignored here.

    Pass a FASTA, one PDB/mmCIF, a structure folder, or a table pairing receptors
    with ligand SDFs. FASTA searches use MMseqs2; structures use Foldseek and
    MMseqs2 by default. Every search returns protein-chain hits.

    ``'pockets'`` maps PLINDER ligand pockets onto similar regions of your
    proteins, ignoring query ligands. ``'ligands'`` also compares query ligands
    and their interactions. ``'interfaces'`` maps PLINDER interface sites onto
    similar regions of each input chain, and also compares complete interfaces
    when the input contains interacting chains.
    ``'auto'`` (the default) runs the comparisons supported by the input;
    ``'all'`` selects all feature families, skipping unavailable comparisons.
    FASTA supports
    pocket and interface-site discovery, but not ligand or complete-interface
    comparisons.

    PDB receptors contribute amino-acid residues, including modified amino acids;
    other PDB components are omitted. Supply ligand poses through ``ligand_path``.

    The returned result contains paths to score tables and alignment files.
    ``chain_similarity_scores`` gives direct protein-chain hits, with query and
    target coverage, identity, sequence similarity, and Foldseek lDDT as
    integer percentages. FASTA results include the original ``sequence_id``.
    With ``include_monomers=True``, chain hits also cover proteins outside
    ligand pockets and protein interfaces. Pocket and interface scores still
    use their own targets.
    ``interface_side_scores`` locates potential interface sites on input chains;
    its coverage and identity percentages use the entire known interface side
    as the denominator. ``interface_scores`` compares two-chain interfaces and
    is absent for monomer inputs. Site matches do not establish that a protein
    binds the partner observed in the matched PLINDER interface.
    ``search_config`` controls search filters. Its shared ``evalue`` defaults to
    0.01; ``foldseek_evalue``, ``mmseqs_evalue`` and ``steam_evalue`` override it
    independently for their respective backends.
    Input IDs become custom structure IDs in these outputs. Generated ligand
    chain IDs map to the original SDF paths in ``inputs/<input_id>.ligands.tsv``.
    """
    if mode not in {"auto", "pockets", "ligands", "interfaces", "all"}:
        raise ValueError("mode must be auto, pockets, ligands, interfaces or all")
    output_dir = Path(output_dir).resolve()
    common: dict[str, Any] = dict(
        work_dir=output_dir,
        data_dir=release.data_dir if release is not None else None,
        threads=threads,
        search_config=search_config,
        plinder_entry_ids=plinder_entry_ids,
        store_aligned_pocket_residues=store_aligned_pocket_residues,
        include_monomers=include_monomers,
    )
    path = None if isinstance(inputs, pd.DataFrame) else Path(inputs).resolve()
    if path is not None and path.name.lower().endswith(
        (".fasta", ".fa", ".faa", ".fasta.gz", ".fa.gz", ".faa.gz")
    ):
        if mode == "ligands":
            raise ValueError(
                "FASTA inputs support pocket and interface-site discovery; use structures for ligand comparisons"
            )
        return score_custom_sequence_file(
            path,
            backends=tuple(backends) if backends is not None else ("mmseqs",),
            include_interfaces=mode in {"auto", "interfaces", "all"},
            **common,
        )
    is_table = path is None or path.suffix.lower() in {".csv", ".tsv", ".parquet"}
    if is_table:
        structures = read_structure_table(inputs)
    else:
        assert path is not None
        rows = []
        for file in find_structure_files(path):
            name = file.stem if file.suffix.lower() == ".gz" else file.name
            rows.append({"input_id": Path(name).stem, "structure_path": str(file)})
        if not rows:
            raise ValueError(f"No PDB/mmCIF structures found in {path}")
        if len({row["input_id"] for row in rows}) != len(rows):
            raise ValueError("Structure filenames must have unique basenames")
        structures = read_structure_table(pd.DataFrame(rows))
    include_ligands = mode in {"auto", "ligands", "all"}
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
        include_ligands=None if mode in {"auto", "all"} else include_ligands,
        include_interfaces=mode in {"auto", "interfaces", "all"},
        **common,
    )
