# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdRascalMCES
from rdkit.Chem.rdchem import Mol
from rdkit.rdBase import BlockLogs

from plinder.core.structure.smallmols_utils import uncharge_mol


def smiles2inchikey(smiles: str, remove_stereo: bool = False) -> str:
    """Return an InChIKey, falling back to standardized canonical SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    with BlockLogs():
        mol = uncharge_mol(mol)
        if remove_stereo:
            Chem.RemoveStereochemistry(mol)
        inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        inchikey = Chem.CanonSmiles(Chem.MolToSmiles(mol), useChiral=not remove_stereo)
    return str(inchikey)


def get_ecfp_fingerprint(
    smiles: str, radius: int, nbits: int
) -> Optional[np.ndarray[int, Any]]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        return np.array(fp)
    except Exception:
        return None


def mol2morgan_fp(
    mol: Mol | str, radius: int = 2, nbits: int = 2048
) -> DataStructs.ExplicitBitVect:
    """Convert an RDKit molecule or SMILES string to a Morgan fingerprint."""
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    if mol is None:
        raise ValueError("cannot fingerprint an invalid molecule")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return generator.GetFingerprint(mol)


def tanimoto_maxsim_and_argmax(
    long_list: list[Any], test_list: list[Any]
) -> tuple[np.ndarray[float], np.ndarray[int]]:
    """Calculate each test fingerprint's maximum similarity to a reference list."""
    similarity_matrix = [
        DataStructs.BulkTanimotoSimilarity(fp, long_list) for fp in test_list
    ]
    return np.max(similarity_matrix, axis=1) * 100, np.argmax(similarity_matrix, axis=1)


def get_mmp_similarity_dict(
    mmp_path: Path, min_constant_size: int = 5
) -> dict[str, dict[str, float]]:
    if str(mmp_path).casefold().endswith(".parquet"):
        mmp_df = pd.read_parquet(
            mmp_path,
            columns=[
                "ligand_smiles_1",
                "ligand_smiles_2",
                "shared_core_num_heavy_atoms",
                "ligand_1_shared_core_fraction",
                "ligand_2_shared_core_fraction",
            ],
        ).rename(
            columns={
                "ligand_smiles_1": "SMILES1",
                "ligand_smiles_2": "SMILES2",
                "shared_core_num_heavy_atoms": "const_size",
                "ligand_1_shared_core_fraction": "sim_1_to_2",
                "ligand_2_shared_core_fraction": "sim_2_to_1",
            }
        )
    else:
        mmp_df = pd.read_csv(
            mmp_path,
            sep="\t",
            compression="gzip",
            names=["SMILES1", "SMILES2", "id1", "id2", "V1>>V2", "CONSTANT"],
        )
        const_size_map = {
            smarts: Chem.MolFromSmarts(smarts).GetNumHeavyAtoms()
            for smarts in mmp_df.CONSTANT.drop_duplicates()
        }
        mmp_df["const_size"] = mmp_df.CONSTANT.map(const_size_map)
        smiles_size_map = {
            smiles: Chem.MolFromSmiles(smiles).GetNumHeavyAtoms()
            for smiles in set(mmp_df.SMILES1.to_list() + mmp_df.SMILES2.to_list())
        }
        mmp_df["SMILES1_size"] = mmp_df.SMILES1.map(smiles_size_map)
        mmp_df["SMILES2_size"] = mmp_df.SMILES2.map(smiles_size_map)
        mmp_df["sim_1_to_2"] = mmp_df["const_size"] / mmp_df["SMILES1_size"]
        mmp_df["sim_2_to_1"] = mmp_df["const_size"] / mmp_df["SMILES2_size"]
    mmp_df = mmp_df[mmp_df["const_size"] >= min_constant_size]
    smiles_inchikey_map = {
        smiles: smiles2inchikey(smiles, remove_stereo=True)
        for smiles in set(mmp_df.SMILES1.to_list() + mmp_df.SMILES2.to_list())
    }
    mmp_df["inchikey1"] = mmp_df.SMILES1.map(smiles_inchikey_map)
    mmp_df["inchikey2"] = mmp_df.SMILES2.map(smiles_inchikey_map)

    similarities: dict[str, dict[str, float]] = {}
    for inchikey1, inchikey2, similarity12, similarity21 in mmp_df[
        ["inchikey1", "inchikey2", "sim_1_to_2", "sim_2_to_1"]
    ].values:
        forward = similarities.setdefault(inchikey1, {})
        reverse = similarities.setdefault(inchikey2, {})
        forward[inchikey2] = max(forward.get(inchikey2, 0.0), similarity12 * 100)
        reverse[inchikey1] = max(reverse.get(inchikey1, 0.0), similarity21 * 100)
    return similarities


def rdRascalMCES_similarity(mol1: Mol, mol2: Mol, sim_threshold: float = 0.4) -> float:
    rascal_opts = rdRascalMCES.RascalOptions()
    rascal_opts.allBestMCESs = False
    rascal_opts.returnEmptyMCES = True
    rascal_opts.completeAromaticRings = False
    rascal_opts.ringMatchesRingOnly = False
    rascal_opts.maxBondMatchPairs = 1000
    rascal_opts.timeout = 2
    rascal_opts.similarityThreshold = sim_threshold
    result = rdRascalMCES.FindMCES(mol1, mol2, rascal_opts)
    return float(result[0].tier2Sim if result else 0)
