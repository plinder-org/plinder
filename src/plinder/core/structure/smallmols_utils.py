# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdDepictor, rdMolDescriptors, rdRascalMCES
from rdkit.Chem.MolStandardize import rdMolStandardize

from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)


def uncharge_mol(mol: Mol) -> Mol:
    """Neutralize formal charges where possible."""
    if sum([at.GetFormalCharge() != 0 for at in mol.GetAtoms()]):
        uncharger = rdMolStandardize.Uncharger(canonicalOrder=True, force=False)
        res = uncharger.uncharge(mol)
        res.UpdatePropertyCache(strict=False)
        return res
    return mol


def generate_input_conformer(
    template_mol: Chem.Mol,
    addHs: bool = False,
    minimize_maxIters: int = -1,
    skip_3d_confgen: bool = False,
) -> Chem.Mol:
    _mol = copy.deepcopy(template_mol)
    # need to add Hs to generate sensible conformers
    _mol = Chem.AddHs(_mol)

    if skip_3d_confgen:
        confid = -1
    else:
        # try embedding molecule using ETKDGv2 (default)
        confid = AllChem.EmbedMolecule(
            _mol,
            # ps,
            useRandomCoords=True,
            useBasicKnowledge=True,
            maxAttempts=100,
            randomSeed=42,
        )
        if confid != -1:
            if minimize_maxIters > 0:
                # molecule successfully embedded - minimize
                success = AllChem.MMFFOptimizeMolecule(_mol, maxIters=minimize_maxIters)
                # 0 if the optimization converged,
                # -1 if the forcefield could not be set up,
                # 1 if more iterations are required.
                if success == 1:
                    log.info(
                        f"generate_conformer: MMFFOptimizeMolecule - more iterations are required, doubling the steps (2x {minimize_maxIters})"
                    )
                    # extend optimization to double the steps (extends by the same amount)
                    AllChem.MMFFOptimizeMolecule(_mol, maxIters=minimize_maxIters)
                elif success == -1:
                    log.warning(
                        "generate_conformer: MMFFOptimizeMolecule - the forcefield could not be set up"
                    )
        else:
            # this means EmbedMolecule failed
            log.warning(
                "generate_conformer: default EmbedMolecule - failed, trying using useBasicKnowledge=False"
            )
            # try less optimal approach
            confid = AllChem.EmbedMolecule(
                _mol,
                useRandomCoords=True,
                useBasicKnowledge=False,
                maxAttempts=100,
                randomSeed=42,
            )

    if confid == -1:
        # if 3D confgen fails or skipped
        log.warning(
            "generate_conformer: using 2D (rdDepictor.Compute2DCoords) instead 3D"
        )
        confid = rdDepictor.Compute2DCoords(_mol)

    # verify that mol has conformers
    if _mol.GetNumConformers() == 0:
        raise ValueError("Could not generate conformer")

    if not addHs:
        # remove Hs if they should not be kept
        _mol = Chem.RemoveAllHs(_mol, sanitize=False)

    return _mol


def match_ligands(
    input_smiles: str,
    resolved_sdf: str | Path,
) -> tuple[Chem.Mol, Chem.Mol, tuple[NDArray, NDArray]]:
    template_mol = Chem.MolFromSmiles(input_smiles)
    resolved_mol = Chem.MolFromMolFile(resolved_sdf.__str__())
    atom_order_stacks = get_template_to_mol_matches(template_mol, resolved_mol)
    return template_mol, resolved_mol, atom_order_stacks


def get_template_to_mol_matches(
    template: Chem.Mol, mol: Chem.Mol
) -> tuple[NDArray, NDArray]:
    """
    Function that works a lot like get_matched_template but can better deal with fragmented molecules
    """
    rascal_opts = rdRascalMCES.RascalOptions()
    rascal_opts.similarityThreshold = 0.1
    rascal_opts.allBestMCESs = True
    rascal_opts.returnEmptyMCES = True
    rascal_opts.completeAromaticRings = False
    rascal_opts.ringMatchesRingOnly = False
    rascal_opts.maxBondMatchPairs = 5000
    rascal_opts.timeout = 20

    results = rdRascalMCES.FindMCES(mol, template, rascal_opts)
    atom_matches = np.array(results[0].atomMatches())
    bond_matches = np.array(results[0].bondMatches())

    numHA_template = rdMolDescriptors.CalcNumHeavyAtoms(template)
    numHA_mol = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if len(atom_matches) < min(numHA_template, numHA_mol):
        # if not complete molecule is matched
        match_mol = copy.deepcopy(mol)
        ref_mol = copy.deepcopy(template)

        log.warning(
            "get_template_to_mol_matches: could not match template fully - retry with unmatched bonds set as UNSPECIFIED"
        )
        # set all unmatched bonds to UNSPECIFIED to help with the match
        if len(bond_matches):
            [
                b.SetBondType(Chem.BondType.UNSPECIFIED)
                for b in match_mol.GetBonds()
                if not b.GetIdx() in bond_matches[:, 0]
            ]
            [
                b.SetBondType(Chem.BondType.UNSPECIFIED)
                for b in ref_mol.GetBonds()
                if not b.GetIdx() in bond_matches[:, 1]
            ]
            # run again
            results2 = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)
            if len(results2[0].atomMatches()) > len(results[0].atomMatches()):
                results = results2

        # if still not fully matched - attempt one more!
        if len(results[0].atomMatches()) < min(numHA_template, numHA_mol):
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in match_mol.GetBonds()]
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in ref_mol.GetBonds()]
            # run again
            results3 = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)
            if len(results3[0].atomMatches()) > len(results[0].atomMatches()):
                results = results3

    # convert to atom order array stacks
    template_atom_order_stack1 = np.array(
        [[ix_templ for _, ix_templ in match.atomMatches()] for match in results]
    )
    mol_atom_order_stack2 = np.array(
        [[ix_mol for ix_mol, _ in match.atomMatches()] for match in results]
    )
    return template_atom_order_stack1, mol_atom_order_stack2


def compare_stereo_to_template(
    resolved_mol: Mol,
    template_mol: Mol,
) -> bool:
    """Compare per-atom CIP codes between a resolved mol and a template.

    If the resolved mol has fewer atoms than the template (partial
    resolution), the template is trimmed via MCS and CIP codes are
    re-assigned on the trimmed template before comparison.

    Only stereocenters defined in *both* mols are compared.  Achiral
    compounds (no stereocenters in either mol) return True — no conflict.

    Parameters
    ----------
    resolved_mol : Mol
        RDKit Mol with ``AssignStereochemistryFrom3D`` already called.
        Must have PDB residue info on each atom.
    template_mol : Mol
        CCD template Mol with stereo assigned from ideal 3D.

    Returns
    -------
    bool
        True if stereo matches (or achiral), False if any center differs.
    """
    # Build atom name → CIP map from resolved mol
    resolved_cip: dict[str, str] = {}
    for atom in resolved_mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(
                f"Atom {atom.GetIdx()} in resolved mol has no PDB residue info"
            )
        cip = atom.GetPropsAsDict().get("_CIPCode", "")
        if cip:
            resolved_cip[info.GetName().strip()] = cip

    # Trim template if partially resolved
    if resolved_mol.GetNumAtoms() < template_mol.GetNumAtoms():
        try:
            trimmed = get_matched_template(template_mol, resolved_mol)
            Chem.AssignStereochemistry(trimmed, cleanIt=True, force=True)
        except Exception:
            trimmed = template_mol
    else:
        trimmed = template_mol

    # Compare CIP codes where both sides are defined
    for atom in trimmed.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(
                f"Atom {atom.GetIdx()} in template lost PDB residue info after trimming"
            )
        template_cip = atom.GetPropsAsDict().get("_CIPCode", "")
        if not template_cip:
            continue
        atom_name = info.GetName().strip()
        resolved_cip_val = resolved_cip.get(atom_name, "")
        if not resolved_cip_val:
            continue
        if resolved_cip_val != template_cip:
            return False

    # No mismatches found (including achiral — no stereocenters = no conflict)
    return True


# below functions used for data ingest
def mol_assigned_bond_orders_by_template(template_mol: Mol, mol: Mol) -> Mol:
    try:
        fixed_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    except Exception as e:
        log.warning(
            f"mol_assigned_bond_orders_by_template: {e} - try get_matched_template"
        )
        template_mol = get_matched_template(template_mol, mol)
        fixed_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    return fixed_mol


def _remove_unmatched(
    mol: Chem.Mol, matched_atoms: NDArray, matched_bonds: NDArray
) -> Chem.Mol:
    """Remove atoms and bonds in mol whose indices are not in match.
    Parameters
    ----------
    mol : Chem.Mol
        the mol to be modified
    match_bonds : NDArray
        indices that are matches and should not be removed
    matched_bonds : NDArray
        indices that are matches and should not be removed
    Returns
    -------
    Chem.Mol
        the mol with unmatched atoms removed
    """
    res = Chem.RWMol(mol)
    bonds_to_remove = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in mol.GetBonds()
        if b.GetIdx() not in matched_bonds
    ]
    atoms_to_remove = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in matched_atoms
    ]
    res.BeginBatchEdit()
    for atom_idx in atoms_to_remove:
        neighbors = res.GetAtomWithIdx(atom_idx).GetNeighbors()
        for neighbor in neighbors:
            bonds_to_remove.append((atom_idx, neighbor.GetIdx()))
            # if atom to be removed is neighbour to double bond - set that stereo bond to undefined
            [
                nb.SetStereo(Chem.rdchem.BondStereo.STEREONONE)
                for nb in neighbor.GetBonds()
                if nb.GetBondType() == Chem.rdchem.BondType.DOUBLE
            ]
        res.RemoveAtom(atom_idx)
    # now remove all unmatched bonds!
    for bond_idx1, bond_idx2 in bonds_to_remove:
        res.RemoveBond(bond_idx1, bond_idx2)

    res.CommitBatchEdit()
    res = Chem.Mol(res)
    try:
        Chem.SanitizeMol(res)
    except Exception:
        pass
    [a.SetNumRadicalElectrons(0) for a in res.GetAtoms()]
    return res


def get_matched_template(template: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
    """Trim template to the MCS with mol using Rascal MCES.

    Handles fragmented molecules and unmatched bonds correctly.
    """
    rascal_opts = rdRascalMCES.RascalOptions()
    rascal_opts.similarityThreshold = 0.1
    rascal_opts.allBestMCESs = False
    rascal_opts.returnEmptyMCES = True
    rascal_opts.completeAromaticRings = False
    rascal_opts.ringMatchesRingOnly = False
    rascal_opts.maxBondMatchPairs = 5000
    rascal_opts.timeout = 20

    result = rdRascalMCES.FindMCES(mol, template, rascal_opts)[0]
    atom_matches = np.array(result.atomMatches())
    bond_matches = np.array(result.bondMatches())

    numHA_template = rdMolDescriptors.CalcNumHeavyAtoms(template)
    numHA_mol = rdMolDescriptors.CalcNumHeavyAtoms(mol)

    if len(atom_matches) < min(numHA_template, numHA_mol):
        # if not complete molecule is matched
        match_mol = copy.deepcopy(mol)
        ref_mol = copy.deepcopy(template)

        log.warning(
            "get_matched_template: could not match template fully - retry with unmatched bonds set as UNSPECIFIED"
        )
        # set all unmatched bonds to UNSPECIFIED to help with the match
        if len(bond_matches):
            [
                b.SetBondType(Chem.BondType.UNSPECIFIED)
                for b in match_mol.GetBonds()
                if not b.GetIdx() in bond_matches[:, 0]
            ]
            [
                b.SetBondType(Chem.BondType.UNSPECIFIED)
                for b in ref_mol.GetBonds()
                if not b.GetIdx() in bond_matches[:, 1]
            ]
            # run again
            result = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)[0]

        # if still not fully matched - attempt one more!
        if len(result.atomMatches()) < min(numHA_template, numHA_mol):
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in match_mol.GetBonds()]
            [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in ref_mol.GetBonds()]
            # run again
            result2 = rdRascalMCES.FindMCES(match_mol, ref_mol, rascal_opts)[0]
            if len(result2.atomMatches()) > len(result.atomMatches()):
                result = result2

    # used to remove all atoms and bonds from the ref that are not matched
    atom_map_template = np.array([j for i, j in result.atomMatches()])
    bond_map_template = np.array([j for i, j in result.bondMatches()])
    if len(atom_map_template) == 0:
        raise ValueError("get_matched_template: cannot match mol to template")

    # Removes unmatched atoms and bonds from the template
    matched_template_mol = _remove_unmatched(
        template, atom_map_template, bond_map_template
    )
    return matched_template_mol
