# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    rdFingerprintGenerator,
    rdMHFPFingerprint,
    rdRascalMCES,
)
from rdkit.Chem.rdchem import Mol
from rdkit.rdBase import BlockLogs

# MHFP6 = MinHashed FingerPrint at radius 3 (ECFP diameter 6 equivalent).
# Probst, D. & Reymond, J-L. "A probabilistic molecular fingerprint for big
# data settings." J. Cheminform. 10, 8 (2018).
# https://doi.org/10.1186/s13321-018-0321-8
MHFP6_RADIUS = 3
MHFP6_N_PERMUTATIONS = 2048
# A fixed permutation seed makes the MinHash vectors reproducible across runs
# and machines; changing it invalidates every stored MHFP6 fingerprint.
MHFP6_SEED = 42


def smiles2nonstereo(smiles: str) -> str:
    """Return a stereo-stripped RDKit canonical SMILES — the ligand identity key.

    RDKit always yields a canonical SMILES from a parsed molecule, whereas InChIKey
    generation fails for many CCD molecules (organometallics, exotic valences), so
    canonical SMILES is the robust, uniform identifier. Stereochemistry is removed
    (identity is stereo-insensitive — used for split stratification and MMP
    grouping) and no charge normalization is applied: the CCD representation is
    taken as-is, consistent with the exact-SMILES cofactor/artifact matching.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    with BlockLogs():
        Chem.RemoveStereochemistry(mol)
        return str(Chem.MolToSmiles(mol))


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


@cache
def _mhfp6_encoder() -> rdMHFPFingerprint.MHFPEncoder:
    """Return the shared, seeded MHFP6 encoder (permutations are seed-derived)."""
    return rdMHFPFingerprint.MHFPEncoder(MHFP6_N_PERMUTATIONS, MHFP6_SEED)


def mol2mhfp6(mol: Mol | str) -> NDArray[np.uint32]:
    """Return the MHFP6 MinHash fingerprint of a molecule or SMILES string.

    MHFP6 hashes the set of circular SMILES shingles (radii 1..3) and keeps the
    per-permutation minima, so the fraction of matching positions between two
    fingerprints estimates the Jaccard similarity of their shingle sets
    (Probst & Reymond, J. Cheminform. 2018). Unlike a folded Morgan bit-vector,
    the result is a dense vector of ``MHFP6_N_PERMUTATIONS`` uint32 hashes.
    """
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    if mol is None:
        raise ValueError("cannot fingerprint an invalid molecule")
    encoded = _mhfp6_encoder().EncodeMol(mol, radius=MHFP6_RADIUS)
    return np.asarray(encoded, dtype=np.uint32)


def mhfp6_bulk_jaccard(
    query: NDArray[np.uint32], reference_matrix: NDArray[np.uint32]
) -> NDArray[np.float64]:
    """Estimate one MHFP6 vector's Jaccard similarity to every reference row.

    Each reference row is a MinHash vector; the estimated Jaccard similarity is
    the fraction of permutations at which the two vectors share the same minimum
    hash -- the MinHash counterpart of ``BulkTanimotoSimilarity`` for bit-vectors.
    """
    if reference_matrix.ndim != 2 or reference_matrix.shape[1] != query.shape[0]:
        raise ValueError("query and reference MinHash widths do not match")
    matches = np.count_nonzero(reference_matrix == query, axis=1)
    return matches / reference_matrix.shape[1]


def mhfp6_maxsim_and_argmax(
    long_matrix: NDArray[np.uint32], test_matrix: NDArray[np.uint32]
) -> tuple[NDArray[np.float64], NDArray[np.intp]]:
    """Calculate each test MHFP6 vector's max estimated Jaccard to a reference set."""
    similarity_matrix = np.stack(
        [mhfp6_bulk_jaccard(row, long_matrix) for row in test_matrix]
    )
    return (
        np.max(similarity_matrix, axis=1) * 100,
        np.argmax(similarity_matrix, axis=1),
    )


def get_mmp_similarity_dict(
    mmp_path: Path, min_constant_size: int = 5
) -> dict[str, dict[str, float]]:
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
    mmp_df = mmp_df[mmp_df["const_size"] >= min_constant_size]
    smiles_inchikey_map = {
        smiles: smiles2nonstereo(smiles)
        for smiles in set(mmp_df.SMILES1.to_list() + mmp_df.SMILES2.to_list())
    }
    mmp_df["inchikey1"] = mmp_df.SMILES1.map(smiles_inchikey_map)
    mmp_df["inchikey2"] = mmp_df.SMILES2.map(smiles_inchikey_map)
    smiles_size_map = {
        smiles: Chem.MolFromSmiles(smiles).GetNumHeavyAtoms()
        for smiles in set(mmp_df.SMILES1.to_list() + mmp_df.SMILES2.to_list())
    }
    mmp_df["SMILES1_size"] = mmp_df.SMILES1.map(smiles_size_map)
    mmp_df["SMILES2_size"] = mmp_df.SMILES2.map(smiles_size_map)
    mmp_df["sim_1_to_2"] = mmp_df["const_size"] / mmp_df["SMILES1_size"]
    mmp_df["sim_2_to_1"] = mmp_df["const_size"] / mmp_df["SMILES2_size"]

    similarities: dict[str, dict[str, float]] = {}
    for inchikey1, inchikey2, similarity12, similarity21 in mmp_df[
        ["inchikey1", "inchikey2", "sim_1_to_2", "sim_2_to_1"]
    ].values:
        similarities.setdefault(inchikey1, {})[inchikey2] = similarity12 * 100
        similarities.setdefault(inchikey2, {})[inchikey1] = similarity21 * 100
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
