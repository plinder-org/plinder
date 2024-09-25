# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import pandas as pd
from ost import io, mol

from plinder.core import PlinderSystem, scores
from plinder.core.utils.log import setup_logger
from plinder.data.utils.annotations.save_utils import save_cif_file
from plinder.eval.docking import utils

LOG = setup_logger(__name__)


def get_resolution(cif_file: Path) -> float | None:
    if not cif_file.exists():
        LOG.info(f"no such file {cif_file}")
        return None
    block = gemmi.cif.read(cif_file.as_posix()).sole_block()
    res = block.find_value("_refine.ls_d_res_high")
    if not res:
        res = block.find_value("_em_3d_reconstruction.resolution")
    if res:
        return float(gemmi.cif.as_number(res))
    return None


def get_plddt(cif_file: Path) -> float | None:
    # may only be present in pred_mmseqs so fail gracefully
    if not cif_file.exists():
        LOG.info(f"no such file {cif_file}")
        return None
    block = gemmi.cif.read(str(cif_file.as_posix())).sole_block()
    metric = block.find_value("_ma_qa_metric_global.metric_value")
    if metric:
        return float(gemmi.cif.as_number(metric))
    return None


def superpose_to_system(
    system_mol: mol.EntityHandle,
    target_cif_file: Path,
    save_folder: Path,
    target_chain: str | None = None,
    name_mapping: dict[str, str] | None = None,
) -> None:
    """
    Superpose a target asymmetric unit and chain to a system.
    Score the ligands as if transplanted from the system to the target

    Parameters
    ----------
    system_mol : mol.EntityHandle
        receptor.cif loaded into an EntityHandle
    target_cif_file : Path
        Path to the target asymmetric unit cif file
    save_folder : Path
        Folder to save the superposed target cif and pdb files
    target_chain : str
        Chain of the target asymmetric unit to superpose
    """
    # Load target asymmetric unit and chain
    target_mol, info = io.LoadMMCIF(
        target_cif_file.as_posix(), info=True, fault_tolerant=True
    )
    if target_chain is not None:
        target_mol = mol.CreateEntityFromView(
            target_mol.Select(f"chain='{target_chain}'"), True
        )
    target_mol = mol.CreateEntityFromView(target_mol.Select("water=False"), True)

    # Superpose target to query system
    superposition = mol.alg.Superpose(target_mol, system_mol, match="local-aln")
    LOG.info(f"target_cif {target_cif_file} rmsd: {superposition.rmsd}")
    target_mol.FixTransform()

    if name_mapping is None:
        assert target_chain is not None
        # Rename target_chain to A for PDB format
        target_pdb = target_mol.Copy()
        if target_chain != "A":
            edi = target_pdb.EditXCS(mol.BUFFERED_EDIT)
            edi.RenameChain(target_pdb.FindChain(target_chain), "A")
            edi.UpdateICS()
    else:
        # Rename target system according to its existing name mapping for PDB format
        target_pdb = target_mol.Copy()
        intermediate_names = {}
        edi = target_pdb.EditXCS(mol.BUFFERED_EDIT)
        for i, chain in enumerate(target_pdb.GetChainList()):
            intermediate_names[f"T{i}"] = chain.name
            edi.RenameChain(chain, f"T{i}")
        edi.UpdateICS()
        for i, chain in enumerate(target_pdb.GetChainList()):
            edi.RenameChain(chain, name_mapping[intermediate_names[chain.name]])

    # Save superposed target cif and pdb
    cif_file = save_folder / "superposed.cif"
    save_cif_file(target_mol, info, cif_file.stem, cif_file)
    pdb_file = save_folder / "superposed.pdb"
    io.SavePDB(target_pdb, pdb_file.as_posix())


@dataclass
class LinkedStructureConfig:
    num_per_system: (
        int
    ) = 5  # Maximum number of apo/pred/cross structures to keep per system
    filter_criteria: dict[str, int] = field(
        default_factory=lambda: {
            "pocket_fident": 95,
            "protein_fident_weighted_sum": 95,
            "protein_fident_qcov_weighted_sum": 80,
            "pocket_lddt": 20,
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
    def get_system_ligand_files(system_id: str) -> list[str]:
        system_folder = get_cif_file(data_dir, "holo", system_id).parent
        return [
            (system_folder / "ligand_files" / f"{c}.sdf").as_posix()
            for c in system_id.split("__")[-1].split("_")
        ]

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
        links["target_id"] = links["target_system"].map(lambda x: x.split("_")[0])
        links.to_parquet(
            output_file.parent / f"{output_file.stem}_intermediate.parquet", index=False
        )

    targets = set(links["target_id"])
    score_dir = search_db
    if search_db == "holo":
        score_dir = "apo"  # always get resolution from ingest
    target_files = [get_cif_file(data_dir, score_dir, x) for x in targets]
    del links

    sort_scores = {}
    if search_db == "pred":
        func = get_plddt
        ascending = False
    else:
        func = get_resolution
        ascending = True
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
    links["ligand_files"] = links["reference_system_id"].apply(get_system_ligand_files)
    links.to_parquet(output_file, index=False)


def get_cif_file(data_dir: Path, search_db: str, system: str) -> Path:
    if search_db == "apo":
        return (
            data_dir
            / "ingest"
            / system[1:3]
            / f"pdb_0000{system[:4]}"
            / f"pdb_0000{system[:4]}_xyz-enrich.cif.gz"
        )
    elif search_db == "pred":
        return (
            data_dir
            / "dbs"
            / "alphafold"
            / f"AF-{system.split('_')[0]}-F1-model_v4.cif"
        )
    elif search_db == "holo":
        return data_dir / "raw_entries" / system[1:3] / system / "receptor.cif"
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
    name_mapping, target_chain = None, None
    if search_db == "holo":
        name_mapping_file = target_cif_file.parent / "chain_mapping.json"
        if not name_mapping_file.exists():
            LOG.error(
                f"get_transplanted_ligand_scores_system: {name_mapping_file} does not exist"
            )
            return False
        with open(name_mapping_file) as f:
            name_mapping = json.load(f)
    else:
        target_chain = link.id.split("_")[-1]
    try:
        superpose_to_system(
            system_mol=reference_system.receptor_entity,
            target_cif_file=target_cif_file,
            save_folder=save_folder,
            name_mapping=name_mapping,
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
            link.ligand_files,
            reference_system,
            score_protein=True,
            score_posebusters=True,
            score_posebusters_full_report=True,
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
        reference_system = PlinderSystem(system_id=system)
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
