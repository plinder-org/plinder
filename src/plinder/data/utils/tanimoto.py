# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist

from plinder.core.utils import schemas
from plinder.core.utils.log import setup_logger
from plinder.data import smallmolecules

if TYPE_CHECKING:
    from plinder.data.utils.annotations.aggregate_annotations import Entry

LOG = setup_logger(__name__)


def get_ecfp_fingerprint(
    smiles: str, radius: int, nbits: int
) -> Optional[np.ndarray[int, Any]]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        return np.array(fp)
    except:
        return None


def load_ligands_from_entry(
    *,
    entry: "Entry",
    ccd_col: str = "unique_ccd_code",
    smiles_col: str = "rdkit_canonical_smiles",
) -> Optional[pd.DataFrame]:
    """
    Load ligands from entries for tanimoto similarity calculations

    Parameters
    ----------
    entry : Entry
        annotation entry
    ccd_col : str
        entry attribute where ccd code is fetched
    smiles_col : str
        entry attribute where smiles is fetched
    """
    ligands = []
    for system_id, system in entry.systems.items():
        if system.system_type != "holo":
            LOG.info(f"skipping {system_id} as it's not a holo system")
            continue
        for ligand in system.ligands:
            smiles = getattr(ligand, smiles_col)
            ccd_code = getattr(ligand, ccd_col)
            if not smiles:
                continue
            ligand_id = f"{system_id}__{ligand.biounit_id}.{ligand.asym_id}"
            ligands.append(
                (
                    entry.pdb_id,
                    system_id,
                    smiles,
                    ccd_code,
                    ligand_id,
                )
            )
    LOG.info(f"{entry.pdb_id} loaded {len(ligands)} ligands")
    if not len(ligands):
        return None
    df = (
        pd.DataFrame(
            ligands,
            columns=[
                "pdb_id",
                "system_id",
                "ligand_rdkit_canonical_smiles",
                "ligand_ccd_code",
                "ligand_id",
            ],
        )
        .dropna(subset="ligand_id")
        .drop_duplicates(subset=["ligand_id"])
        .reset_index(drop=True)
        .sort_values(by="ligand_id")
    )
    LOG.info(f"after deduplication {len(df.index)} entries")
    # aggregate multiple ligands into one:
    df["inchikeys"] = df["ligand_rdkit_canonical_smiles"].map(
        smallmolecules.smil2inchikey
    )
    return df


def compute_ligand_fingerprints(
    *, data_dir: Path, split_char: str = "__", radius: int = 2, nbits: int = 1024
) -> None:
    ligands = (
        pd.read_parquet(data_dir / "ligands").drop_duplicates().reset_index(drop=True)
    )
    for col in ligands.columns:
        nunique = ligands[col].nunique()
        LOG.info(f"compute_ligand_fingerprints: unique {col}={nunique}")

    ligands_unique = (
        ligands[ligands["inchikeys"] != ""][
            ["inchikeys", "ligand_rdkit_canonical_smiles"]
        ]
        .drop_duplicates(subset="inchikeys")
        .sort_values(by="inchikeys")
        .reset_index(drop=True)
        .sort_values(by=["inchikeys"])
    )
    ids, _ = pd.factorize(ligands_unique["inchikeys"])
    ids_str = [f"l{int(i)}" for i in ids]
    ligands_unique["number_id_by_inchikeys"] = ids
    ligands_unique["number_id_by_inchikeys_withl"] = ids_str
    missing_inchi = ligands_unique[ligands_unique["number_id_by_inchikeys"] <= -1]
    LOG.info(f"ligands missing inchikeys: {len(missing_inchi.index)}")
    ligands_unique = ligands_unique[ligands_unique["number_id_by_inchikeys"] > -1]
    output_dir = data_dir / "fingerprints"
    output_dir.mkdir(exist_ok=True, parents=True)
    ligands_unique.to_parquet(output_dir / "ligands_per_inchikey.parquet", index=False)

    LOG.info(f"computing ECFP fingerprints with radius={radius} nbits={nbits}")
    ligands_unique["ECFP4"] = ligands_unique["ligand_rdkit_canonical_smiles"].apply(
        get_ecfp_fingerprint,
        radius=radius,
        nbits=nbits,
    )
    if ligands_unique["ECFP4"].isnull().any():
        raise ValueError(f"found {ligands_unique['ECFP4'].isnull().sum()}")

    ecfp = (
        np.concatenate(ligands_unique.ECFP4.values).reshape(-1, nbits).astype(np.int8)
    )
    LOG.info(f"writing {output_dir}/ligands_per_inchikey_ecfp4.npy")
    np.save(output_dir / "ligands_per_inchikey_ecfp4.npy", ecfp)

    ligs = pd.merge(
        ligands,
        ligands_unique[
            ["inchikeys", "number_id_by_inchikeys", "number_id_by_inchikeys_withl"]
        ],
        on="inchikeys",
    )
    # free up memory
    ligands_unique = None
    ecfp = None

    # remove nan inchikeys entries
    ligs = ligs[ligs["number_id_by_inchikeys"].notna()]
    ligs["multi_number_id_by_inchikeys"] = ligs["number_id_by_inchikeys"]
    ligs["multi_number_id_by_inchikeys_withl"] = ligs["number_id_by_inchikeys_withl"]

    for aggcol in [
        "multi_number_id_by_inchikeys",
        "multi_number_id_by_inchikeys_withl",
    ]:
        ligs[aggcol] = ligs[aggcol].astype(str)
        agg_df = ligs.groupby("system_id")[aggcol].agg(split_char.join).reset_index()
        del ligs[aggcol]
        ligs = pd.merge(ligs, agg_df, on="system_id", how="left")

    LOG.info(f"writing {output_dir}/ligands_per_system.parquet")
    ligs.to_parquet(output_dir / "ligands_per_system.parquet", index=False)


def ligand_scores(
    *,
    ligand_ids: list[int],
    data_dir: Path,
    output_path: Path,
    number_id_col: str = "number_id_by_inchikeys",
    # TODO: always save all ligands!!
    save_top_k_similar_ligands: int = 5000,
    multiply_by: int = 100,
) -> None:
    fp_dir = data_dir / "fingerprints"
    all_ligs_ids = pd.read_parquet(fp_dir / "ligands_per_inchikey.parquet")
    LOG.info(f"ligand_scores: loaded {len(all_ligs_ids.index)} ligands")
    fingerprints = np.load(fp_dir / "ligands_per_inchikey_ecfp4.npy")
    LOG.info(f"ligand_scores: loaded {fingerprints.shape[0]} fingerprints")

    # make sure the shape match and the index is in order
    if fingerprints.shape[0] != len(all_ligs_ids.index):
        raise ValueError("ligands don't match fingerprints!")
    if all_ligs_ids[number_id_col].min() != 0:
        raise ValueError("inconsistency in ligand ids, no index=zero found!")
    if all_ligs_ids[number_id_col].max() != all_ligs_ids.shape[0] - 1:
        raise ValueError("inconsistency in ligand ids, max index != ids shape!")

    query_ecfp = fingerprints[ligand_ids]
    # target_list = fingerprints
    ecfp_distances = cdist(query_ecfp, fingerprints, metric="jaccard")
    # change to similarity score
    ecfp_distances = np.array(1 - ecfp_distances)

    # For each i and j, the metric dist(u=XA[i], v=XB[j]) is computed
    # and stored in the ith, jth entry of ecfp_distances = cdist
    top_k_indices = np.argsort(ecfp_distances, axis=1)[:, -save_top_k_similar_ligands:][
        ..., ::-1
    ]

    table = []
    for row, cols in enumerate(top_k_indices):
        LOG.info(f"row {ligand_ids[row]}")
        query_list = [ligand_ids[row]] * len(cols)
        target_list = cols
        tani_topk = (np.round(ecfp_distances[row, cols], 2) * multiply_by).astype(int)
        single_table = pa.table(
            [
                pa.array(query_list),
                pa.array(target_list),
                pa.array(tani_topk),
            ],
            schema=schemas.TANIMOTO_SCORE_SCHEMA,
        )
        table.append(single_table)
    table = pa.concat_tables(table)
    LOG.info(f"writing {output_path}")
    pq.write_table(table, output_path)
