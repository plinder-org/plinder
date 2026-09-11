# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from functools import cache
from pathlib import Path
from typing import Any, NamedTuple

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

from plinder.core.utils.sanitize import mol_from_smiles

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

    Parsing uses the tolerant :func:`mol_from_smiles` so over-valent ligands
    (boron cages, hypervalent metals) still yield a canonical 2D graph key
    instead of falling back to their raw, non-canonical SMILES.
    """
    mol = mol_from_smiles(smiles)
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
        mol = mol_from_smiles(mol)
    if mol is None:
        raise ValueError("cannot fingerprint an invalid molecule")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return generator.GetFingerprint(mol)


# TODO: legacy - no live callers; superseded by mhfp6_maxsim_and_argmax.
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
        mol = mol_from_smiles(mol)
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


# TODO: legacy - no live callers (tests only); MMP pairs now come from
# make_ligand_mmp_pairs and query_ccd_mmp_pairs.
def _sanitized(smiles: str) -> Mol:
    """Parse without RDKit's strict sanitizer, then sanitize the peppr way."""
    from plinder.core.utils.sanitize import sanitize as peppr_sanitize

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    peppr_sanitize(mol)
    return mol


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
            smiles: _sanitized(smiles).GetNumHeavyAtoms()
            for smiles in set(mmp_df.SMILES1.to_list() + mmp_df.SMILES2.to_list())
        }
        mmp_df["SMILES1_size"] = mmp_df.SMILES1.map(smiles_size_map)
        mmp_df["SMILES2_size"] = mmp_df.SMILES2.map(smiles_size_map)
        mmp_df["sim_1_to_2"] = mmp_df["const_size"] / mmp_df["SMILES1_size"]
        mmp_df["sim_2_to_1"] = mmp_df["const_size"] / mmp_df["SMILES2_size"]
    mmp_df = mmp_df[mmp_df["const_size"] >= min_constant_size]
    smiles_inchikey_map = {
        smiles: smiles2nonstereo(smiles)
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


Centres = dict[int, tuple[Chem.StereoDescriptor, list[int]]]


def _stereo_centres(mol: Mol) -> Centres:
    """Specified tetrahedral centres: descriptor and neighbours, padded to four."""
    return {
        si.centeredOn: (si.descriptor, (list(si.controllingAtoms) + [-1])[:4])
        for si in Chem.FindPotentialStereo(mol)
        if si.type == Chem.StereoType.Atom_Tetrahedral
        and si.specified == Chem.StereoSpecified.Specified
    }


def _opposed_centres(
    first: Centres, second: Centres, partner: dict[int, int], atoms: set[int]
) -> int:
    """Count matched tetrahedral centres of opposite handedness under the mapping.

    Tet_CW/Tet_CCW is relative to each molecule's own neighbour order, so the
    mapping between the two orders is a permutation and the descriptors are
    equal only up to its parity. A centre needs three mapped neighbours; the
    leftover slot on each side (fourth neighbour or implicit H, which RDKit
    orders last) is taken as corresponding.
    """
    opposed = 0
    for i in atoms & first.keys():
        j = partner.get(i)
        if j not in second:
            continue
        (descriptor, neighbours), (other_descriptor, other) = first[i], second[j]
        # i's neighbour order translated into mol2 indices, next to j's own order
        mapped = [partner.get(a, -1) for a in neighbours]
        common = set(mapped) & set(other) - {-1}
        if len(common) < 3:
            continue
        mapped, other = (
            [a if a in common else -1 for a in ns] for ns in (mapped, other)
        )
        order = [other.index(a) for a in mapped]
        odd = sum(order[x] > order[y] for x in range(4) for y in range(x + 1, 4)) % 2
        opposed += (descriptor == other_descriptor) == bool(odd)
    return opposed


class RascalParityMatch(NamedTuple):
    """One Rascal MCES reduced to what the PARITY-like score needs."""

    atoms: dict[int, int]  # largest connected fragment, mol1 -> mol2 atom index
    bonds: int  # matched bonds over all fragments
    opposed: int  # fragment centres of opposite handedness under the mapping
    alternatives: tuple[dict[int, int], ...] = ()  # fragments of the other best MCESs

    def score(self, mol1: Mol, mol2: Mol, *, stereo: bool = True) -> float:
        """PARITY-like similarity in [0, 1]; an opposed centre counts half an atom."""
        atoms = len(self.atoms) - (0.5 * self.opposed if stereo else 0.0)
        return parity_similarity(
            atoms,
            self.bonds,
            mol1.GetNumAtoms() + mol2.GetNumAtoms(),
            mol1.GetNumBonds() + mol2.GetNumBonds(),
        )


def parity_similarity(
    atoms: float, bonds: float, atom_total: int, bond_total: int
) -> float:
    """Mean of the atom and bond Tanimotos; totals are both molecules' counts summed."""
    if bonds <= 0:
        return 0.0
    return (atoms / (atom_total - atoms) + bonds / (bond_total - bonds)) / 2


def rascal_parity_match(
    mol1: Mol,
    mol2: Mol,
    *,
    target: float = 0.3,
    timeout: int = 2,
    all_best: bool = False,
) -> RascalParityMatch:
    """Element-exact Rascal MCES with fragments allowed, pruned to a target bond similarity.

    Pairs whose bond Tanimoto cannot reach ``target`` are pruned before the
    clique search and come back empty. Rascal maximises matched bonds, and the
    MCESs of that size differ in their largest fragment, so the bonds come from
    one MCES and the fragment from a second run that keeps only the largest
    fragment over all of them (``singleLargestFrag``). With ``all_best`` every
    fragment is enumerated instead (a further 1.5x) and the one scoring best is
    taken; the rest follow, best first and then in a fixed order.
    """
    options = rdRascalMCES.RascalOptions()
    # tier screens use Johnson similarity; below (2t/(1+t))^2 no pair reaches bond Tanimoto t
    options.similarityThreshold = (2 * target / (1 + target)) ** 2
    # fewest matched bonds giving bond Tanimoto >= target; Rascal's bound is exclusive
    options.minCliqueSize = max(
        int(np.ceil(target * (mol1.GetNumBonds() + mol2.GetNumBonds()) / (1 + target)))
        - 1,
        0,
    )
    options.returnEmptyMCES = True
    options.completeAromaticRings = False
    options.ringMatchesRingOnly = False
    options.maxBondMatchPairs = 5000
    options.timeout = timeout
    options.allBestMCESs = all_best
    matches: dict[tuple[tuple[int, int], ...], dict[int, int]] = {}
    bonds = 0
    for result in rdRascalMCES.FindMCES(mol1, mol2, options):
        matched = result.bondMatches()
        if not matched:
            continue
        # largest connected fragment of the matched bonds, in mol1 atom indices
        atom_map: dict[int, int] = {}
        core = Chem.PathToSubmol(mol1, [i for i, _ in matched], atomMap=atom_map)
        original = {new: old for old, new in atom_map.items()}
        partner = dict(result.atomMatches())
        fragment = {
            original[i]: partner[original[i]]
            for i in max(Chem.GetMolFrags(core), key=len)
        }
        matches.setdefault(tuple(sorted(fragment.items())), partner)
        bonds = len(matched)
    if not matches:
        return RascalParityMatch({}, 0, 0)
    if not all_best:
        options.singleLargestFrag = True
        best = rdRascalMCES.FindMCES(mol1, mol2, options)
        if best and best[0].bondMatches():
            fragment = dict(best[0].atomMatches())
            matches = {tuple(sorted(fragment.items())): fragment}
    centres_1, centres_2 = _stereo_centres(mol1), _stereo_centres(mol2)
    opposed = {
        key: _opposed_centres(centres_1, centres_2, partner, {i for i, _ in key})
        for key, partner in matches.items()
    }
    ordered = sorted(matches, key=lambda key: (0.5 * opposed[key] - len(key), key))
    return RascalParityMatch(
        dict(ordered[0]),
        bonds,
        opposed[ordered[0]],
        tuple(dict(key) for key in ordered[1:]),
    )


def rascal_parity_score(
    mol1: Mol, mol2: Mol, *, target: float = 0.3, stereo: bool = True, timeout: int = 2
) -> float:
    """PARITY-like similarity from one Rascal MCES, in [0, 1].

    Element-exact MCES with fragments allowed, scored as the mean of two
    Tanimoto terms: atoms of the largest connected fragment over the atom union,
    and all matched bonds over the bond union. The largest fragment carries
    residue order (a shuffled peptide splits into one fragment per residue, a
    substitution does not); the bonds keep credit for matches beyond a
    substituted connecting atom. With ``stereo`` a matched tetrahedral centre of
    opposite handedness counts half an atom.
    """
    match = rascal_parity_match(mol1, mol2, target=target, timeout=timeout)
    return match.score(mol1, mol2, stereo=stereo)
