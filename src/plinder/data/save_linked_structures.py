# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import pandas as pd

from plinder.core import PlinderSystem, scores
from plinder.core.structure.atoms import is_hydrogen_isotope
from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.cif_utils import (
    _cif_scalar,
    read_mmcif_container,
    read_mmcif_file,
)
from plinder.data.utils.annotations.save_utils import save_cif_file
from plinder.eval.docking import utils

LOG = setup_logger(__name__)


def get_resolution(cif_file: Path) -> float | None:
    if not cif_file.exists():
        LOG.info(f"no such file {cif_file}")
        return None
    block = read_mmcif_container(cif_file)
    res = _cif_scalar(block, "refine", "ls_d_res_high")
    if res is None:
        res = _cif_scalar(block, "em_3d_reconstruction", "resolution")
    if res is not None:
        return float(res)
    return None


def get_plddt(cif_file: Path) -> float | None:
    # may only be present in pred_mmseqs so fail gracefully
    if not cif_file.exists():
        LOG.info(f"no such file {cif_file}")
        return None
    block = read_mmcif_container(cif_file)
    val = _cif_scalar(block, "ma_qa_metric_global", "metric_value")
    if val is not None:
        return float(val)
    return None


def superpose_to_system(
    system_atoms: struc.AtomArray,
    target_cif_file: Path,
    save_folder: Path,
    target_chain: str | None = None,
) -> None:
    """
    Superpose a target asymmetric unit and chain to a system.

    Parameters
    ----------
    system_atoms : AtomArray
        Reference system atoms (receptor).
    target_cif_file : Path
        Path to the target asymmetric unit cif file.
    save_folder : Path
        Folder to save the superposed target mmCIF file.
    target_chain : str, optional
        Chain of the target to superpose.
    """
    # Load target
    cif_file_obj = read_mmcif_file(target_cif_file)
    target_atoms = pdbx.get_structure(
        cif_file_obj, model=1, use_author_fields=False, include_bonds=True
    )
    target_atoms = target_atoms[~is_hydrogen_isotope(target_atoms.element)]

    if target_chain is not None:
        target_atoms = target_atoms[target_atoms.chain_id == target_chain]
    target_atoms = target_atoms[~struc.filter_solvent(target_atoms)]

    # Superpose target to system
    ref_ca = system_atoms[
        struc.filter_amino_acids(system_atoms) & (system_atoms.atom_name == "CA")
    ]
    target_ca = target_atoms[
        struc.filter_amino_acids(target_atoms) & (target_atoms.atom_name == "CA")
    ]

    if len(ref_ca) > 0 and len(target_ca) > 0:
        # Match by sequence alignment
        fitted, transformation = struc.superimpose(ref_ca, target_ca)
        # Apply transformation to all target atoms
        target_atoms = struc.superimpose_apply(target_atoms, transformation)
        rmsd = struc.rmsd(ref_ca, fitted)
        LOG.info(f"target_cif {target_cif_file} rmsd: {rmsd:.2f}")

    # Save superposed target
    save_cif_file(target_atoms, "superposed", save_folder / "superposed.cif")


@dataclass
class LinkedStructureConfig:
    num_per_system: int = (
        5  # Maximum number of apo/pred/cross structures to keep per system
    )
    filter_criteria: dict[str, int] = field(
        default_factory=lambda: {
            "pocket_fident": 95,
            "protein_fident_weighted_sum": 95,
            "protein_fident_qcov_weighted_sum": 80,
            "protein_lddt_weighted_sum": 20,
        }
    )  # Filter criteria for deciding whether to keep a linked structure (AND logic)


def make_linked_structures_data_file(
    data_dir: Path,
    search_db: str,
    superposed_folder: Path,
    output_file: Path,
    cfg: LinkedStructureConfig = LinkedStructureConfig(),
    num_processes: int = 8,
) -> None:
    (superposed_folder / search_db).mkdir(exist_ok=True, parents=True)
    filters = []
    for metric, threshold in cfg.filter_criteria.items():
        filters.append([("metric", "==", metric), ("similarity", ">=", threshold)])
    output_file.parent.mkdir(exist_ok=True, parents=True)
    if (output_file.parent / f"{output_file.stem}_intermediate.parquet").is_file():
        links = pd.read_parquet(
            output_file.parent / f"{output_file.stem}_intermediate.parquet"
        )
    else:
        links = scores.query_protein_similarity(
            search_db=search_db,
            columns=["query_system", "target_system", "metric", "similarity"],
            filters=filters,
        )
        assert links is not None
        links = links.iloc[
            links.groupby(["query_system", "target_system", "metric"], observed=True)[
                "similarity"
            ].idxmax()
        ]
        links = links[
            links["query_system"].str[:4] != links["target_system"].str[:4]
        ].reset_index(drop=True)
        links = links.pivot(
            index=["query_system", "target_system"],
            columns="metric",
            values="similarity",
        ).reset_index()
        query = " and ".join([f"{m} >= {t}" for (m, t) in cfg.filter_criteria.items()])
        links = links.query(query)
        if search_db == "holo":
            links["target_id"] = links["target_system"]
        else:
            links["target_id"] = links["target_system"].map(lambda x: x.split("_")[0])
        links.to_parquet(
            output_file.parent / f"{output_file.stem}_intermediate.parquet", index=False
        )

    targets = set(links["target_id"])
    if search_db == "holo":
        target_files = [get_original_pdb_mmcif(data_dir, x) for x in targets]
    else:
        target_files = [get_cif_file(data_dir, search_db, x) for x in targets]
    del links

    sort_scores = {}
    if search_db == "pred":
        func = get_plddt
        ascending = False
    else:
        func = get_resolution
        ascending = True
    if num_processes == 1:
        resolutions = list(map(func, target_files))
    else:
        with multiprocessing.get_context("spawn").Pool(num_processes) as pool:
            resolutions = pool.map(func, target_files)
    sort_scores = dict(zip(targets, resolutions))

    nonnull = sum((v for v in sort_scores.values() if v is not None))
    LOG.info(f"non null scores: {nonnull}")
    links = pd.read_parquet(
        output_file.parent / f"{output_file.stem}_intermediate.parquet"
    )
    links["sort_score"] = links["target_id"].map(sort_scores)
    links = (
        links[links["sort_score"].notna()]
        .sort_values("sort_score", ascending=ascending)
        .groupby("query_system")
        .head(cfg.num_per_system)
        .reset_index(drop=True)
    )
    links.rename(
        columns={"query_system": "reference_system_id", "target_system": "id"},
        inplace=True,
    )
    links["receptor_file"] = links.apply(
        lambda row: (
            superposed_folder
            / search_db
            / row.reference_system_id
            / row.id
            / "superposed.cif"
        ).as_posix(),
        axis=1,
    )
    links.to_parquet(output_file, index=False)


def get_original_pdb_mmcif(data_dir: Path, system_or_pdb_id: str) -> Path:
    """Locate the original PDB mmCIF used by the ingest pipeline."""
    pdb_id = system_or_pdb_id[:4]
    return (
        data_dir
        / "ingest"
        / pdb_id[1:3]
        / f"pdb_0000{pdb_id}"
        / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
    )


def get_canonical_ligand_dir(data_dir: Path, system_or_pdb_id: str) -> Path:
    """Locate the canonical, asymmetric-unit ligand SDF directory."""
    pdb_id = system_or_pdb_id[:4]
    return data_dir / "raw_entries" / pdb_id[1:3] / pdb_id / "ligand_files"


def get_cif_file(data_dir: Path, search_db: str, system: str) -> Path:
    if search_db == "apo":
        return get_original_pdb_mmcif(data_dir, system)
    elif search_db == "pred":
        return (
            data_dir
            / "dbs"
            / "alphafold"
            / f"AF-{system.split('_')[0]}-F1-model_v4.cif"
        )
    elif search_db == "holo":
        return Path(
            PlinderSystem(
                system_id=system,
                source_mmcif=get_original_pdb_mmcif(data_dir, system),
                reconstruction_dir=data_dir / "reconstructed_systems" / system,
                canonical_ligand_dir=get_canonical_ligand_dir(data_dir, system),
            ).receptor_cif
        )
    else:
        raise ValueError("search_db much be apo, holo or pred")


def save_superposition(
    *,
    data_dir: Path,
    save_folder: Path,
    search_db: str,
    link: pd.Series,
    reference_system: PlinderSystem,
    overwrite: bool = False,
) -> bool:
    if not overwrite and (save_folder / "superposed.cif").exists():
        LOG.warning(f"save_superposition: {save_folder / 'superposed.cif'} exists")
        return True
    target_cif_file = get_cif_file(data_dir, search_db, link.id)
    if not target_cif_file.exists():
        LOG.error(
            f"get_transplanted_ligand_scores_system: {link.id} cif file doesn't exist"
        )
        return False
    target_chain = None
    if search_db != "holo":
        target_chain = link.id.split("_")[-1]
    try:
        superpose_to_system(
            system_atoms=reference_system.receptor_structure,
            target_cif_file=target_cif_file,
            save_folder=save_folder,
            target_chain=target_chain,
        )
        return True
    except Exception as e:
        LOG.error(f"save_superposition: Error in superpose_to_system: {repr(e)} {e}")
        return False


def system_save_and_score_representative(
    link: pd.Series,
    reference_system: PlinderSystem,
    data_dir: Path,
    search_db: str,
    output_folder: Path,
    overwrite: bool = False,
) -> None:
    save_folder = output_folder / search_db / link.reference_system_id / link.id
    save_folder.mkdir(exist_ok=True, parents=True)
    carry_on = save_superposition(
        data_dir=data_dir,
        save_folder=save_folder,
        search_db=search_db,
        link=link,
        reference_system=reference_system,
        overwrite=overwrite,
    )
    if not carry_on:
        return

    try:
        scores = utils.ModelScores.from_model_files(
            link.id,
            save_folder / "superposed.cif",
            [Path(path) for path in reference_system.ligand_sdfs.values()],
            reference_system,
            score_protein=True,
        ).summarize_scores()
        with open(save_folder / "scores.json", "w") as f:
            json.dump(scores, f)
    except Exception as e:
        LOG.error(
            f"system_save_and_score_representative: Error for {link.reference_system_id} against {link.id}: {e}"
        )


def system_save_and_score_representatives(
    system: str,
    saver_input: pd.DataFrame,
    data_dir: Path,
    search_db: str,
    output_folder: Path,
    overwrite: bool = False,
) -> None:
    try:
        reference_system = PlinderSystem(
            system_id=system,
            source_mmcif=get_original_pdb_mmcif(data_dir, system),
            reconstruction_dir=data_dir / "reconstructed_systems" / system,
            canonical_ligand_dir=get_canonical_ligand_dir(data_dir, system),
        )
    except Exception as e:
        LOG.error(
            f"system_save_and_score_representatives: Error in making reference system: {e}"
        )
        return
    saver_input.apply(
        lambda row: system_save_and_score_representative(
            row,
            reference_system=reference_system,
            data_dir=data_dir,
            search_db=search_db,
            output_folder=output_folder,
            overwrite=overwrite,
        ),
        axis=1,
    )


def save_linked_structures(
    links_file: Path,
    data_dir: Path,
    search_db: str,
    output_folder: Path,
    num_threads: int = 8,
    overwrite: bool = False,
) -> None:
    """
    Saves superposed linked structure files

    Parameters
    ----------
    links_file: Path
        The parquet file generated by filter_and_sort_scores
    data_dir: Path
    output_folder: Path
        Path to save superposed files, in the form <output_folder>/<system_id[1:3]>/<system_id>/<linked_system_id>/superposed.cif
    num_threads: int
        Number of processes to use #TODO rename
    overwrite: bool
        Skips existing files if False
    """
    links = pd.read_parquet(links_file)
    with multiprocessing.get_context("spawn").Pool(num_threads) as p:
        p.starmap(
            system_save_and_score_representatives,
            [
                (system, group, data_dir, search_db, output_folder, overwrite)
                for system, group in links.groupby("reference_system_id")
            ],
        )
