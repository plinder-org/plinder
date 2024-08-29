# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
# This script contains all funcs related to small molecules.
from __future__ import annotations

from typing import Any, Optional

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdDepictor, rdFingerprintGenerator
from rdkit.Chem.rdchem import Mol

rdDepictor.SetPreferCoordGen(True)


def nha(smi: str) -> int:
    try:
        return int(Chem.MolFromSmiles(smi).GetNumHeavyAtoms())
    except:
        return -1


def mw(smi: str) -> Optional[float]:
    try:
        return round(float(Descriptors.ExactMolWt(Chem.MolFromSmiles(smi))), 2)
    except:
        return round(-1, 2)


def n_ring(smi: str) -> int:
    try:
        return int(Descriptors.RingCount(Chem.MolFromSmiles(smi)))
    except:
        return -1


def n_ro_bonds(smi: str) -> int:
    try:
        return int(Descriptors.NumRotatableBonds(Chem.MolFromSmiles(smi)))
    except:
        return -1


def fcsp3(smi: str) -> Optional[float]:
    try:
        return round(float(Descriptors.FractionCSP3(Chem.MolFromSmiles(smi))), 2)
    except:
        return round(-1, 2)


def NumHAcceptors(smi: str) -> int:
    try:
        return int(Chem.Lipinski.NumHAcceptors(Chem.MolFromSmiles(smi)))
    except:
        return -1


def NumHDonors(smi: str) -> int:
    try:
        return int(Chem.Lipinski.NumHDonors(Chem.MolFromSmiles(smi)))
    except:
        return -1


def smil2inchikey(smiles: str) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        # Chem.RemoveStereochemistry(mol)
        inchikey = Chem.MolToInchiKey(mol)
        return str(inchikey)
    except:
        return None


def neutralize_mol(mol: Mol) -> Optional[Mol]:
    try:
        pattern = Chem.MolFromSmarts(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
        )
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        return mol
    except:
        return None


def smil2nochiral_nocharge(smiles: str) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        Chem.RemoveStereochemistry(mol)
        newmol = neutralize_mol(mol)
        if newmol is None:
            return None
        else:
            smilnew = Chem.MolToSmiles(newmol)
        return str(smilnew)
    except:
        return None


def rdkit_valid(smiles: str) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        smi = Chem.MolToSmiles(mol)
        if smi:
            return True
        else:
            return False
    except:
        return False


def containC(smiles: str) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "C":
                return True
        return False
    except:
        return False


def get_ecfp_fingerprint(
    smiles: str, radius: int = 3, nbits: int = 2048
) -> None | np.ndarray[int, Any]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        return np.array(fp)
    except:
        return None


def mol2morgan_fp(
    mol: Mol | str, radius: int = 2, nbits: int = 2048
) -> DataStructs.ExplicitBitVect:
    """Convert an RDKit molecule to a Morgan fingerprint
    :param mol: RDKit molecule or SMILES str
    :param radius: fingerprint radius
    :param nbits: number of fingerprint bits
    :return: RDKit Morgan fingerprint
    """
    if type(mol) == str:
        mol = Chem.MolFromSmiles(mol)
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    fp = mfpgen.GetFingerprint(mol)
    return fp


def tanimoto_maxsim_matrix(
    fp_list1: list[Any], fp_list2: list[Any]
) -> np.ndarray[float]:
    """Calculate maximum similarity for the fingerprint second list to the first fingerprint lists"""
    similarity_matrix = [
        np.max(DataStructs.BulkTanimotoSimilarity(fp, fp_list1)) for fp in fp_list2
    ]
    return np.array(similarity_matrix) * 100
