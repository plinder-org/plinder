# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy

import numpy as np
from numpy.typing import NDArray
from openbabel import pybel
from ost import conop, io
from ost import mol as omol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdFMCS import FindMCS

from plinder.core.utils.log import setup_logger
from plinder.data.common.constants import BASE_DIR

LOG = setup_logger(__name__)
COMPOUND_LIB = conop.GetDefaultLib()
PRD_LIB = conop.CompoundLib.Load(
    str(BASE_DIR / "data/utils/annotations/static_files/prdcc.chemlib")
)


# def process_organometalics(mol: Mol) -> Mol:
#     """Converts covalent organometallic bond to dative bond.
#     This is helpful in cases where covalent metal bonding to
#     ligands is causing rdkit to try errors

#     Parameters
#     ----------
#     mol : Chem.rdchem.Mol

#     Returns
#     -------
#     Chem.rdchem.Mol
#         Mol with covalent organometallic bond
#         replaced with dative bond
#     """
#     metals = ["Ti", "Al", "Mo", "Ru", "Co", "Rh", "Ir", "Ni", "Zr", "Hf", "W", "Fe"]
#     for bond in mol.GetBonds():
#         if (bond.GetEndAtom().GetSymbol()) in metals or (
#             bond.GetBeginAtom().GetSymbol() in metals
#         ):
#             logging.info("found metal-ligand bond")
#             logging.info("original type: " + str(bond.GetBondType()))
#             btype = Chem.rdchem.BondType.DATIVE
#             bond.SetBondType(btype)
#             logging.info(
#                 "changed to: " + str(mol.GetBonds()[bond.GetIdx()].GetBondType())
#             )
#             try:
#                 sanitize_mol(mol)
#             except ValueError:
#                 logging.error("Sanitization failed: {e}")
#     return mol


def sanitize_mol(mol: Mol) -> None:
    """Santitize while keeping hydrogen as is.

    Parameters
    ----------
    mol : Chem.rdchem.Mol

    Returns
    -------
    None
        Sanitizes molecule in place
    """
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        Chem.SanitizeMol(
            mol,
            # sanitize all but keep hydrogens as is
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
        )


def params_removeHs(mol: Chem.Mol) -> Chem.Mol:
    params = Chem.rdmolops.RemoveHsParameters()
    params.removeIsotopes = True
    params.removeDegreeZero = True
    params.removeHigherDegrees = True
    params.removeOnlyHNeighbors = True
    params.removeNontetrahedralNeighbors = True
    return Chem.rdmolops.RemoveHs(mol, params, sanitize=False)


def explicit_H_remover(mol: Mol, remove_hydrogens: list[int]) -> Mol:
    """removes all H atoms in the list and all bonds to those hydrogens"""
    res = Chem.RWMol(mol)
    res.BeginBatchEdit()
    for aid in remove_hydrogens:
        neighbors = res.GetAtomWithIdx(aid).GetNeighbors()
        for neighbor in neighbors:
            res.RemoveBond(aid, neighbor.GetIdx())
        res.RemoveAtom(aid)
    res.CommitBatchEdit()
    return res


def fix_valency_issues(mol: Mol) -> Mol:
    """Fix valency issues with rdkit mol and return sanitized.
    Deals with cases like:
    Removed hydrogens if there is an issue with their valence!
    Explicit valence for atom # X N, 4, is greater than permitted
    Explicit valence for atom # X O, 3, is greater than permitted
    # NOT: Explicit valence for atom # X C, 5, is greater than permitted
    # skipped C - as we don't like Texas carbons :)

    Parameters
    ----------
    mol : Chem.rdchem.Mol

    Returns
    -------
    Chem.rdchem.Mol
        Sanitized Mol with valency issues fixed
    """
    max_explicit_valency_per_element = {
        # 6: 4,
        7: 3,
        8: 2,
        # 1: 1,
    }
    mol.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(mol)
    if not ps:
        # if no problems - just sanitize and return
        sanitize_mol(mol)
        return mol

    # quick scan if the issue with hydrogens - see if needed to remove
    delete_hydrogens = set()
    ps = Chem.DetectChemistryProblems(mol)
    for p in ps:
        if p.GetType() == "AtomValenceException":
            at = mol.GetAtomWithIdx(p.GetAtomIdx())
            atm_no = at.GetAtomicNum()
            if atm_no == 1:
                delete_hydrogens.add(p.GetAtomIdx())
            elif atm_no == 6:
                delete_hydrogens |= {
                    nat.GetIdx() for nat in at.GetNeighbors() if nat.GetAtomicNum() == 1
                }
    # remove explicit hydrogents when some are causing issues
    if delete_hydrogens:
        LOG.warn(
            f"fix_valency_issues: found issues with H atoms {delete_hydrogens} - will try removing these atoms explicitly!"
        )
        mol = explicit_H_remover(mol, list(delete_hydrogens))
        # scan again for remaining problems
        ps = Chem.DetectChemistryProblems(mol)

    # deal with remainng issues, if any
    for p in ps:
        if p.GetType() == "AtomValenceException":
            at = mol.GetAtomWithIdx(p.GetAtomIdx())
            atm_no = at.GetAtomicNum()
            formal_charge = at.GetFormalCharge()
            valency = at.GetExplicitValence()
            elem_max_explicit_valency = max_explicit_valency_per_element[atm_no]
            expected_charge = valency - elem_max_explicit_valency
            if expected_charge > formal_charge:
                # Fix Explicit valence issue
                at.SetFormalCharge(expected_charge)
        if p.GetType() == "KekulizeException":
            # hack: only works for nitrogens with missing explicit Hs
            for atidx in p.GetAtomIndices():
                at = mol.GetAtomWithIdx(atidx)
                # set one of the nitrogens with two bonds in a ring system as "[nH]"
                if at.GetAtomicNum() == 7 and at.GetDegree() == 2:
                    at.SetNumExplicitHs(1)
                    break
    sanitize_mol(mol)
    return mol


def assign_bond_from_smiles(smiles: str, mol: Mol) -> Mol:
    """Assign bonds from a list smiles.

    The goal of this bond assigner is to capture all
    ligands, including the ones with subunits.

    Parameters
    ----------
    smiles : str
        smiles string of template for bond assignment
    mol : Chem.rdchem.Mol
        Mol that needs bond assignment

    Returns
    -------
    Chem.rdchem.Mol
        Mol bonds assigned
    """
    try:
        # Iterative assign bonds to subunits; this captures multisubunit ligands
        template = AllChem.MolFromSmiles(smiles)
        return AllChem.AssignBondOrdersFromTemplate(template, mol)

    except ValueError:
        return None


def get_element_count(mol: Mol) -> dict[int, int]:
    atomic_count: dict[int, int] = {}
    for atom in mol.GetAtoms():
        elem = atom.GetAtomicNum()
        if elem == 1:
            continue
        atomic_count.setdefault(elem, 0)
        atomic_count[elem] += 1
    return atomic_count


def remove_unmatched_atoms(mol: Chem.Mol, match: NDArray) -> Chem.Mol:
    """Remove atoms in mol whose indices are not in match.
    Parameters
    ----------
    mol : Chem.Mol
        the mol to be modified
    match : NDArray
        indices that are matches and should not be removed
    Returns
    -------
    Chem.Mol
        the mol with unmatched atoms removed
    """
    res = Chem.RWMol(mol)
    atoms_to_remove = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in match]
    res.BeginBatchEdit()
    for atom_idx in atoms_to_remove:
        neighbors = res.GetAtomWithIdx(atom_idx).GetNeighbors()
        for neighbor in neighbors:
            res.RemoveBond(atom_idx, neighbor.GetIdx())
        res.RemoveAtom(atom_idx)
    res.CommitBatchEdit()
    res = Chem.Mol(res)
    try:
        Chem.SanitizeMol(res)
    except:
        pass
    [a.SetNumRadicalElectrons(0) for a in res.GetAtoms()]
    return res


def get_matched_template(template: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
    """
    Perform MCS matching between a (subject) mol and a template; and return the matched template with
    the bond orders of the template. Used to assign bond orders in
    `safe_mol_from_pdb_assign_bond_orders`. Known limitation: if the template has
    double/triple bonds and the mol doesn't (because it's read from PDB), this leads to
    removing all atoms that don't match, incl. the ones bound via e.g., a double bond. This is
    only a problem if we need to use this fallback option because previous attempts in
    `safe_mol_from_pdb_assign_bond_orders` have failed.
    Returns
    -------
    Chem.Mol
        the matching template with bond orders from template
    """
    # set all bonds to unspecified to help with the match
    match_mol = copy.deepcopy(mol)
    [b.SetBondType(Chem.BondType.UNSPECIFIED) for b in match_mol.GetBonds()]

    mcs = FindMCS(
        [match_mol, template],
        completeRingsOnly=False,
        ringMatchesRingOnly=False,
        timeout=10,
    )
    patt = Chem.MolFromSmarts(mcs.smartsString)
    atom_map_template = np.array(template.GetSubstructMatch(patt))
    # remove all atoms from the ref that are not in the MCS --> use this as template for
    # bond orders
    matched_template_mol = remove_unmatched_atoms(template, atom_map_template)
    return matched_template_mol


def mol_assigned_bond_orders_by_template(template_mol: Mol, mol: Mol) -> Mol:
    try:
        # Assign bonds according to template smiles!
        fixed_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    except Exception as e:
        # raise AssertionError(f"mol_assigned_bond_orders_by_template: {e}")
        # update template in case fully resovled mol but bonding is an issue
        LOG.warn(
            f"mol_assigned_bond_orders_by_template: {e} - try get_matched_template"
        )
        template_mol = get_matched_template(template_mol, mol)
        fixed_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    return fixed_mol


def make_rdkit_compatible_mol(mol: Mol) -> Mol | None:
    """Extract rdk mol from mmcif.

    Parameters
    ----------
    mol : Chem.rdchem.Mol

    Returns
    -------
    Chem.rdchem.Mol | None
        Mol of relevant molecule
    """
    try:
        sanitize_mol(mol)
    except:
        try:
            # fix N, O, C, H valency issues and then sanitize
            mol = fix_valency_issues(mol)
        except Exception:
            mol = None
    return mol


def ligand_ost_ent_to_rdkit_mol(
    ent: omol.EntityHandle,
    ligand_smiles: str | None = None,
    ligand_num_unresolved_heavy_atoms: int = 0,
) -> Mol:
    new_chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    edi = ent.EditXCS(omol.BUFFERED_EDIT)
    for i, chain in enumerate(ent.GetChainList()):
        edi.RenameChain(chain, f"{new_chains[i]}")
    edi.UpdateICS()

    pdbstring = io.EntityToPDBStr(ent).strip()
    # NOTE: rdkit's Chem.MolFromPDBBlock does not read connect records
    # work around via openbabel bond perception
    sdfstring = pybel.readstring("pdb", pdbstring).write("sdf")
    rdkit_mol = Chem.MolFromMolBlock(sdfstring, sanitize=False)  # , removeHs=True,
    # removeHs does not work when sanitize is False
    rdkit_mol = params_removeHs(rdkit_mol)

    if ligand_smiles:
        try:
            rdkit_mol = make_rdkit_compatible_mol(rdkit_mol)
            # if smiles match - no fix is needed!
            if Chem.CanonSmiles(ligand_smiles) == Chem.CanonSmiles(
                Chem.MolToSmiles(rdkit_mol)
            ):
                return rdkit_mol
            else:
                raise AssertionError("SMILES do not match reference - will try fixing")
        except Exception as e:
            LOG.warn(f"ligand_ost_ent_to_rdkit_mol: {e}")
        try:
            # another try via OST SDF
            # open structure output singly bonded SDF that can be adjusted with template
            # first - try to get a fixed molecule
            sdfstring_ost = io.EntityToSDFStr(ent).strip()
            rdkit_mol_tmp = Chem.MolFromMolBlock(
                sdfstring_ost,
                sanitize=False,
                # removeHs=True,
            )
            # Note: removeHs does not work when sanitize is False
            rdkit_mol_tmp = params_removeHs(rdkit_mol_tmp)
            # attempt to pre-emptively fix some valency problems
            try:
                rdkit_mol_tmp = fix_valency_issues(rdkit_mol_tmp)
            except Exception:
                LOG.warn(
                    "fix_valency_issues: failed before mol_assigned_bond_orders_by_template"
                )
            # get template from smiles
            template = Chem.MolFromSmiles(ligand_smiles)
            if ligand_num_unresolved_heavy_atoms > 0:
                # update template
                template = get_matched_template(template, rdkit_mol_tmp)
            try:
                # Assign bonds by template
                rdkit_mol = mol_assigned_bond_orders_by_template(
                    template, rdkit_mol_tmp
                )
            except ValueError as e:
                LOG.error(
                    f"template bonds could not be assigned: {e}; "
                    f"template_smiles: {ligand_smiles}"
                )
                raise ValueError(
                    "cannot assign bonds by SMILES, use OpenBabel inference instead"
                )
        except Exception as e:
            LOG.warn(f"ligand_ost_ent_to_rdkit_mol: {e}")
    try:
        # Fix issues if any
        rdkit_mol = make_rdkit_compatible_mol(rdkit_mol)
        if rdkit_mol is None:
            raise ValueError("make_rdkit_compatible_mol: returned None")
        elif len(Chem.MolToSmiles(rdkit_mol).split(".")) > 1:
            raise ValueError(
                f"rdkit_mol seems fragmented: molecule is not connected: {Chem.MolToSmiles(rdkit_mol)}"
            )
    except Exception as e:
        LOG.error(f"ligand_ost_ent_to_rdkit_mol: could not fix: {e}")
    return rdkit_mol


def set_smiles_from_ligand_ost(ent: omol.EntityHandle) -> str:
    residues = [res.name for res in ent.residues]
    if len(residues) == 1:
        resname = residues[0]
        if resname.startswith("PRD_"):
            # TODO - need to make this line used!
            # currently, PRD_entries are not mapped to residue names and are more than one residue!
            mol = PRD_LIB.FindCompound(resname)
        else:
            mol = COMPOUND_LIB.FindCompound(resname)
        if mol is not None:
            try:
                rdkit_mol = Chem.MolFromSmiles(str(mol.smiles), sanitize=False)
                rdkit_mol = make_rdkit_compatible_mol(rdkit_mol)
                return str(Chem.MolToSmiles(rdkit_mol))
            except Exception:
                LOG.warn(
                    "set_smiles_from_ligand_ost: CCD smiles could not be loaded by rdkit, moving to fix"
                )
    rdkit_mol = ligand_ost_ent_to_rdkit_mol(ent)
    return str(Chem.MolToSmiles(rdkit_mol))
