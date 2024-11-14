# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
import gzip
from pathlib import Path
from typing import Any, Union

import numpy as np
from biotite import TextFile
from biotite.structure import get_residues
from biotite.structure.atoms import AtomArray, AtomArrayStack
from biotite.structure.io.pdbx import get_structure
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdDepictor, rdMolDescriptors, rdRascalMCES
from rdkit.Chem.rdFMCS import FindMCS

from plinder.core.structure.vendored import (
    _convert_resn_to_sequence_and_numbering,
    _get_structure_and_res_info,
    align_sequences,
    apply_mask,
)
from plinder.core.utils import constants as pc
from plinder.core.utils.log import setup_logger

log = setup_logger(__name__)

DOCKQ_BACKBONE_ATOMS = ["C", "CA", "N", "O"]

THREE_LETTER_AA = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "SER",
    "THR",
    "VAL",
    "TYR",
    "TRP",
    "PRO",
    "HIS",
    "PHE",
]

BINDING_SITE_METALS = [
    "MG",
    "K",
    "MN",
    "NA",
    "ZN",
    "MG",
    "CA",
    "CD",
    "FE",
    "CU",
    "4MO",
    "CO",
]

_AtomArrayOrStack = Union[AtomArray, AtomArrayStack]


def biotite_ciffile() -> TextFile:
    from biotite.structure.io.pdbx import CIFFile

    return CIFFile


def atom_array_from_cif_file(
    structure: Path | _AtomArrayOrStack, use_author_fields: bool = True
) -> AtomArray | None:
    # TODO: this conversion may be unnecessary as the annotation is
    # Path | _AtomArrayOrStack
    if isinstance(structure, str):
        structure = Path(structure)

    if isinstance(structure, Path):
        reader = biotite_ciffile()
        try:
            if structure.suffix == ".gz":
                with gzip.open(str(structure), "rt", encoding="utf-8") as f:
                    mod = reader.read(f)
            else:
                mod = reader.read(structure)
            arr = get_structure(
                mod, model=1, use_author_fields=use_author_fields, include_bonds=True
            )  # noqa
            return arr
        except Exception as e:
            log.error(f"Unable to parse {structure}! {str(e)}")
    return structure


def generate_input_conformer(
    template_mol: Chem.Mol, addHs: bool = False, minimize_maxIters: int = -1
) -> Chem.Mol:
    _mol = copy.deepcopy(template_mol)
    # need to add Hs to generate sensible conformers
    _mol = Chem.AddHs(_mol)
    # ps = AllChem.ETKDGv2()
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
            "generate_conformer: default EmbedMolecule - failed, try using useBasicKnowledge=False"
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
            # if that still fails - try generating just 2D conformer
            log.warning(
                "generate_conformer: EmbedMolecule - failed, trying rdDepictor.Compute2DCoords instead"
            )
            confid = rdDepictor.Compute2DCoords(_mol)

    # verify that mol has conformers
    if _mol.GetNumConformers() == 0:
        raise ValueError("Could not generate conformer")

    if not addHs:
        # remove Hs if they should not be kept
        _mol = params_removeHs(_mol)

    return _mol


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
    params.removeDummyNeighbors = True
    params.removeNontetrahedralNeighbors = True
    params.removeDefiningBondStereo = True
    params.removeWithWedgedBond = True
    params.showWarnings = True
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
        log.warning(
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


def remove_unmatched_atoms_and_bonds(
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
    except:
        pass
    [a.SetNumRadicalElectrons(0) for a in res.GetAtoms()]
    return res


def get_matched_template_v2(template: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
    """
    Function that works a lot like get_matched_template but can better deal with fragmented molecules
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
            "get_matched_template_v2: could not match template fully - retry with unmatched bonds set as UNSPECIFIED"
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
        raise ValueError("get_matched_template_v2: cannot match mol to template")

    # Removes unmatched atoms and bonds from the template
    matched_template_mol = remove_unmatched_atoms_and_bonds(
        template, atom_map_template, bond_map_template
    )
    return matched_template_mol


def mol_assigned_bond_orders_by_template(template_mol: Mol, mol: Mol) -> Mol:
    try:
        # Assign bonds according to template smiles!
        fixed_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    except Exception as e:
        # raise AssertionError(f"mol_assigned_bond_orders_by_template: {e}")
        # update template in case fully resovled mol but bonding is an issue
        log.warning(
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


def match_ligands(
    input_smiles: str,
    resolved_sdf: str | Path,
) -> tuple[Chem.Mol, Chem.Mol, tuple[_AtomArrayOrStack, _AtomArrayOrStack]]:
    template_mol = Chem.MolFromSmiles(input_smiles)
    resolved_mol = Chem.MolFromMolFile(resolved_sdf.__str__())
    atom_order_stacks = get_template_to_mol_matches(template_mol, resolved_mol)
    return template_mol, resolved_mol, atom_order_stacks


def get_template_to_mol_matches(
    template: Chem.Mol, mol: Chem.Mol
) -> tuple[_AtomArrayOrStack, _AtomArrayOrStack]:
    """
    Function that works a lot like get_matched_template but can better deal with fragmented molecules
    """
    rascal_opts = rdRascalMCES.RascalOptions()
    rascal_opts.similarityThreshold = 0.1
    rascal_opts.allBestMCESs = True
    rascal_opts.returnEmptyMCES = True
    rascal_opts.completeAromaticRings = False
    rascal_opts.ringMatchesRingOnly = False
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


def get_residue_index_mapping_mask(
    ref_seqs: dict[str, str], subject_arr: _AtomArrayOrStack
) -> dict[str, list[int]]:
    mask_map = {}
    for reference_chain, ref_seq in ref_seqs.items():
        # refernece and subject chain are the same
        subject_chain = reference_chain
        subject_mask = subject_arr.chain_id == subject_chain
        subject_ch_arr = apply_mask(subject_arr, subject_mask)

        subj_info = _get_structure_and_res_info(subject_ch_arr)
        subj_seq, subj_numbering = _convert_resn_to_sequence_and_numbering(subj_info)
        alignments = align_sequences(
            ref_seq, subj_seq, ref_numbering=None, subject_numbering=subj_numbering
        )
        (
            _,
            _,
            ref_numbering_mapped,
            _,
        ) = alignments
        mask = np.zeros(len(ref_seq))
        for i in range(len(mask)):
            if (i + 1) in ref_numbering_mapped:
                mask[i] = 1
        mask_map[reference_chain] = mask
    return mask_map


def get_ligand_atom_index_mapping_mask(
    ref_mol: Mol, matching_indices: tuple[int, ...]
) -> NDArray[np._int]:
    mask = np.zeros(len(ref_mol.GetAtoms()))
    for atm_idx in range(len(mask)):
        if atm_idx in matching_indices:
            mask[atm_idx] = 1
    return mask


def make_one_hot_atom_features(atom_name: list[str]) -> list[int]:
    allowed_atom_names = ["C", "N", "O", "S", "P"]
    striped_atom_name = "".join(filter(lambda x: not x.isdigit(), atom_name))[0]
    return [1 if striped_atom_name == atm else 0 for atm in allowed_atom_names]


def _convert_pdb_atom_name_to_elem_symbol(atom_name: str) -> str:
    return "".join(filter(lambda x: not x.isdigit(), atom_name))[0]


def get_per_residue_mask(
    residue_reference_atom_list: list[str], atom_list: list[str]
) -> list[int]:
    res_mask = [1 if i in atom_list else 0 for i in residue_reference_atom_list]
    return res_mask


def get_per_residue_atoms(
    atom_array: _AtomArrayOrStack, resi: int, resn: str
) -> NDArray[np.str_]:
    return atom_array[
        (atom_array.res_id == resi) & (atom_array.res_name == resn)
    ].atom_name


def make_atom_mask(
    atom_array: _AtomArrayOrStack, seq_res: str, seq_mask: list[int]
) -> list[int]:
    seq_res_three_aa = [pc.ONE_TO_THREE[aa] for aa in seq_res]
    resi, resn = get_residues(atom_array)
    residue_tuple = list(zip(tuple(resi), tuple(resn)))

    atom_mask = []
    resolved_residue_start = 0
    for seq_, mask in zip(seq_res_three_aa, seq_mask):
        if mask == 0:
            atom_mask.append([0 for i in range(len(pc.ORDERED_AA_FULL_ATOM[seq_]))])
        else:
            resi, resn = residue_tuple[resolved_residue_start]
            atom_mask.append(
                get_per_residue_mask(
                    pc.ORDERED_AA_FULL_ATOM[resn],
                    get_per_residue_atoms(atom_array, resi, resn),
                )
            )
    return [atm for res in atom_mask for atm in res]


def _stack_atom_array_features(
    atom_arr: _AtomArrayOrStack,
    atom_arr_feat: str,
    chain_order_list: list[str] | None,
) -> list[NDArray[np.int_ | np.str_ | np.float_]]:
    assert chain_order_list is not None
    return [
        getattr(atom_arr[atom_arr.chain_id == chain], atom_arr_feat)
        for chain in chain_order_list
    ]


def _stack_ligand_feat(
    feat_dict: dict[str, Any], chain_order_list: list[str] | None
) -> list[list[list[int]]]:
    assert chain_order_list is not None
    return [feat_dict[chain] for chain in chain_order_list]


def _one_hot_encode_stack(
    stack: list[NDArray],
    feature_dict: dict[str, int],
    unknown_name_filler: str,
) -> list[NDArray]:
    feat_array = []
    unknown_name_filler_value = feature_dict[unknown_name_filler]
    for per_chain_feat in stack:
        feat_array_by_chain = np.zeros(
            (len(per_chain_feat), len(set(list(feature_dict.values()))))
        )
        for index, value in enumerate(per_chain_feat):
            feat_array_by_chain[
                index, feature_dict.get(value, unknown_name_filler_value)
            ] = 1.0
        feat_array.append(feat_array_by_chain)
    return feat_array


def _sequence_full_atom_type_array(
    input_sequences: dict[str, str]
) -> dict[str, NDArray]:
    """Resolved sequence full atom features."""
    seq_atom_dict = {}
    for chain, sequence in input_sequences.items():
        feat = []
        for res in sequence:
            for atom in pc.ORDERED_AA_FULL_ATOM[pc.ONE_TO_THREE[res]]:
                feat.append(_convert_pdb_atom_name_to_elem_symbol(atom))
        seq_atom_dict[chain] = np.array(feat)
    return seq_atom_dict
