# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Mol,
    rdDepictor,
    rdDistGeom,
    rdMolDescriptors,
    rdRascalMCES,
)

from plinder.core.utils.log import setup_logger

# TODO(peppr): peppr_sanitize is the vendored local copy of peppr.sanitize
# (see plinder.core.utils.sanitize); revert to `from peppr import sanitize`
# once the over-valence fixes land in a released peppr.
from plinder.core.utils.sanitize import sanitize as peppr_sanitize

log = setup_logger(__name__)


def generate_input_conformer(
    template_mol: Chem.Mol,
    addHs: bool = False,
    minimize_maxIters: int = -1,
    skip_3d_confgen: bool = False,
) -> Chem.Mol:
    _mol = copy.deepcopy(template_mol)
    # need to add Hs to generate sensible conformers
    _mol = Chem.AddHs(_mol)
    # peppr_sanitize resolves valences, rings and hybridization — including
    # over-valent centres (boron cages, metals) via dative bonds — so the
    # embedder has the chemistry it needs. Paired with
    # ``embedFragmentsSeparately=False`` below, this is the rdkit#8653 work-around
    # that lets EmbedMolecule skip its internal strict-valence sanitize, which
    # would otherwise reject those molecules with an AtomValenceException.
    try:
        peppr_sanitize(_mol)
    except Exception as exc:  # embed anyway; RDKit may still cope
        log.warning(f"generate_input_conformer: peppr_sanitize failed ({exc})")

    def _embed(use_basic_knowledge: bool) -> int:
        params = rdDistGeom.ETKDGv3()
        params.useRandomCoords = True
        params.useBasicKnowledge = use_basic_knowledge
        params.randomSeed = 42
        params.maxIterations = 100
        # skip RDKit's internal (strict) sanitize during embedding; the
        # chemistry it needs was supplied by peppr_sanitize above (rdkit#8653).
        params.embedFragmentsSeparately = False
        return int(rdDistGeom.EmbedMolecule(_mol, params))

    if skip_3d_confgen:
        confid = -1
    else:
        # try embedding molecule using ETKDGv3
        confid = _embed(use_basic_knowledge=True)
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
            confid = _embed(use_basic_knowledge=False)

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
    template_mol = Chem.MolFromSmiles(input_smiles, sanitize=False)
    try:
        peppr_sanitize(template_mol)
    except Exception as exc:  # embed anyway; RDKit may still cope
        log.warning(f"template_mol: peppr_sanitize failed ({exc})")
    resolved_mol = Chem.MolFromMolFile(resolved_sdf.__str__(), sanitize=False)
    try:
        peppr_sanitize(resolved_mol)
    except Exception as exc:  # embed anyway; RDKit may still cope
        log.warning(f"resolved_mol: peppr_sanitize failed ({exc})")
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
) -> bool | None:
    """Compare resolved 3D stereo against a template by CIP code.

    The resolved mol's *own bonds are not trusted*: aromatic ring bonds can
    arrive order-unspecified (which blocks CIP labelling entirely), and
    even with the atom chiral tags in place there is no way to compute a
    ``_CIPCode`` on such a mol. Instead we transplant the resolved 3D
    coordinates onto the *template* graph — which carries correct bond
    orders and matching PDB atom names — perceive tetrahedral chirality
    there via :func:`Chem.AssignAtomChiralTagsFromStructure` (which keeps
    quaternary centers that legacy perception drops), and CIP-label both the
    transplanted probe and the template with the **same** algorithm
    (:func:`Chem.AssignCIPLabels`) so the R/S codes are directly comparable.

    Only stereocenters that are (a) present in the template with a defined
    CIP code and (b) fully resolved — the center *and* all its immediate
    neighbors have coordinates — are compared. An unresolved neighbor makes
    the perceived 3D chirality meaningless, so such centers are skipped.

    Parameters
    ----------
    resolved_mol : Mol
        RDKit Mol with a 3D conformer and PDB residue info on every atom.
        Bond orders may be incomplete — only coordinates and atom names are
        read from it.
    template_mol : Mol
        Template Mol (user SMILES or CCD) with correct bond orders, PDB atom
        names, and reference stereo (from SMILES parity or ideal 3D).

    Returns
    -------
    bool | None
        True if all comparable centers match (or none are chiral), False if
        any center differs, None if no atom could be resolved onto the
        template.
    """
    from rdkit.Geometry import Point3D

    # Resolved 3D coordinates keyed by PDB atom name.
    conf_r = resolved_mol.GetConformer()
    coords: dict[str, Point3D] = {}
    for atom in resolved_mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(
                f"Atom {atom.GetIdx()} in resolved mol has no PDB residue info"
            )
        coords[info.GetName().strip()] = conf_r.GetAtomPosition(atom.GetIdx())

    # Transplant those coordinates onto the template graph. Unresolved atoms
    # are parked at the origin and excluded from comparison below.
    probe = Chem.Mol(template_mol)
    conf = Chem.Conformer(probe.GetNumAtoms())
    resolved_names: set[str] = set()
    for atom in probe.GetAtoms():
        info = atom.GetPDBResidueInfo()
        name = info.GetName().strip() if info is not None else None
        if name is not None and name in coords:
            conf.SetAtomPosition(atom.GetIdx(), coords[name])
            resolved_names.add(name)
        else:
            conf.SetAtomPosition(atom.GetIdx(), Point3D(0.0, 0.0, 0.0))
    if not resolved_names:
        return None
    conf.Set3D(True)
    probe.RemoveAllConformers()
    probe.AddConformer(conf, assignId=True)

    # Perceive chirality from the transplanted geometry, then CIP-label the
    # probe and the template with the same labeller for comparable R/S codes.
    Chem.AssignAtomChiralTagsFromStructure(probe)
    Chem.AssignCIPLabels(probe)
    ref = Chem.Mol(template_mol)
    Chem.AssignCIPLabels(ref)

    # Resolved CIP by name, only for fully-resolved centers (an unresolved
    # neighbor makes the perceived 3D tag meaningless).
    resolved_cip: dict[str, str] = {}
    for atom in probe.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None or info.GetName().strip() not in resolved_names:
            continue
        if any(
            n.GetPDBResidueInfo() is None
            or n.GetPDBResidueInfo().GetName().strip() not in resolved_names
            for n in atom.GetNeighbors()
        ):
            continue
        cip = atom.GetPropsAsDict().get("_CIPCode", "")
        if cip:
            resolved_cip[info.GetName().strip()] = cip

    # Compare where the template defines a center and the resolved side has one.
    for atom in ref.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        template_cip = atom.GetPropsAsDict().get("_CIPCode", "")
        if not template_cip:
            continue
        resolved_cip_val = resolved_cip.get(info.GetName().strip(), "")
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
