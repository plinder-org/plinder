# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import copy
from pathlib import Path

import biotite.structure as struc
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
        # resanitize after RemoveAllHs
        peppr_sanitize(_mol)

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
    """Check resolved 3D stereo against a template — tetrahedral R/S and E/Z.

    Atoms are matched by PDB name. Template stereo is read from its tags via
    :func:`Chem.FindPotentialStereo`; a template that ships ideal coordinates (CCD)
    first has its atom + bond stereo perceived from them via
    :func:`Chem.AssignStereochemistryFrom3D`, while SMILES templates already carry
    tags — no conformer is generated either way. Each tetrahedral center's descriptor
    (Tet_CW/Tet_CCW) fixes the sign of the signed volume of its ``controllingAtoms``;
    each double bond's descriptor (Bond_Cis/Bond_Trans) fixes the dihedral of its
    reference substituents. Both are checked against the resolved coordinates only —
    no CIP priority, no conformer generation and no empty filled coordinates.

    Returns True if every fully-resolved stereo element matches (or the template has
    none), False on any mismatch, None if nothing could be matched. Elements with
    unresolved atoms or degenerate geometry are skipped.
    """

    def _name(atom: Chem.Atom) -> str | None:
        info = atom.GetPDBResidueInfo()
        return info.GetName().strip() if info is not None else None

    # perceive full stereo (atom + bond) from a template's ideal coordinates (CCD) so
    # double-bond E/Z is available; SMILES templates already carry their stereo tags
    if template_mol.GetNumConformers() > 0:
        template_mol = Chem.Mol(template_mol)
        Chem.AssignStereochemistryFrom3D(template_mol)

    conf = resolved_mol.GetConformer()
    resolved_pos: dict[str, NDArray] = {
        n: np.array(conf.GetAtomPosition(a.GetIdx()))
        for a in resolved_mol.GetAtoms()
        if (n := _name(a)) is not None
    }
    template_names = {n for a in template_mol.GetAtoms() if (n := _name(a)) is not None}
    if not resolved_pos.keys() & template_names:
        return None  # nothing maps onto the template

    for si in Chem.FindPotentialStereo(template_mol):
        if si.specified != Chem.StereoSpecified.Specified:
            continue

        if si.type == Chem.StereoType.Atom_Tetrahedral:
            center = _name(template_mol.GetAtomWithIdx(si.centeredOn))
            names = [
                _name(template_mol.GetAtomWithIdx(i))
                for i in list(si.controllingAtoms)[:3]
            ]
            if center is None or center not in resolved_pos or len(names) < 3:
                continue
            if any(n is None or n not in resolved_pos for n in names):
                continue  # center or reference neighbours not resolved -> undefined
            order = [n for n in names if n is not None]  # resolved; list[str] for mypy
            # Tet_CW -> -1, Tet_CCW -> +1 signed volume of controllingAtoms (see tests)
            ref = -1.0 if si.descriptor == Chem.StereoDescriptor.Tet_CW else 1.0
            u, v, w = (resolved_pos[n] - resolved_pos[center] for n in order)
            res = float(np.sign(np.dot(u, np.cross(v, w))))
            if res != 0 and res != ref:  # resolved handedness opposes the template
                return False

        elif si.type == Chem.StereoType.Bond_Double:
            ctrl = list(si.controllingAtoms)
            if len(ctrl) < 3:
                continue
            bond = template_mol.GetBondWithIdx(si.centeredOn)
            # controllingAtoms = [begin_ref, begin_ref, end_ref, end_ref]; the dihedral
            # of one reference per end (ref_a1 - a1 - a2 - ref_a2) fixes cis/trans
            quad = [
                _name(template_mol.GetAtomWithIdx(i))
                for i in (
                    ctrl[0],
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    ctrl[2],
                )
            ]
            if any(n is None or n not in resolved_pos for n in quad):
                continue  # a reference or bond atom not resolved -> undefined
            p = [resolved_pos[n] for n in quad if n is not None]  # 4 pts; list for mypy
            dih = struc.dihedral(p[0], p[1], p[2], p[3])  # radians; nan if degenerate
            if np.isnan(dih):
                continue
            # Bond_Cis -> ~0, Bond_Trans -> ~180 (biotite dihedral convention; see tests)
            expected_cis = si.descriptor == Chem.StereoDescriptor.Bond_Cis
            if (abs(dih) < np.pi / 2) != expected_cis:  # resolved E/Z opposes template
                return False
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


# Categorical ligand-atom descriptors, one vocabulary per column. A value
# outside its vocabulary maps to a trailing "other" bucket (index
# ``len(vocabulary)``), so the encoding is stable across RDKit versions and
# exotic chemistry.
LIGAND_ATOM_FEATURE_NAMES = (
    "atomic_number",
    "chiral_tag",
    "total_degree",
    "formal_charge",
    "implicit_valence",
    "total_hydrogens",
    "radical_electrons",
    "hybridization",
    "is_aromatic",
    "ring_count",
    "in_ring_3",
    "in_ring_4",
    "in_ring_5",
    "in_ring_6",
    "in_ring_7",
    "in_ring_8",
)
_ATOMIC_NUMBERS = tuple(range(1, 119))
_CHIRAL_TAGS = (
    "CHI_UNSPECIFIED",
    "CHI_TETRAHEDRAL_CW",
    "CHI_TETRAHEDRAL_CCW",
    "CHI_OTHER",
)
_TOTAL_DEGREES = tuple(range(11))
_FORMAL_CHARGES = tuple(range(-5, 6))
_IMPLICIT_VALENCES = tuple(range(7))
_HYDROGEN_COUNTS = tuple(range(9))
_RADICAL_ELECTRONS = tuple(range(5))
_HYBRIDIZATIONS = ("SP", "SP2", "SP3", "SP3D", "SP3D2")
_RING_COUNTS = tuple(range(7))
_RING_SIZES = (3, 4, 5, 6, 7, 8)


def _bucket(vocabulary: tuple[int, ...] | tuple[str, ...], value: int | str) -> int:
    """Return the index of ``value`` in ``vocabulary`` or the "other" bucket."""
    try:
        return vocabulary.index(value)
    except ValueError:
        return len(vocabulary)


def ligand_atom_features(mol: Chem.Mol) -> NDArray[np.int64]:
    """Encode every atom of ``mol`` as categorical vocabulary indices.

    Columns follow :data:`LIGAND_ATOM_FEATURE_NAMES`: atomic number, chiral
    tag, total degree, formal charge, implicit valence, total hydrogen count,
    radical electrons, hybridization, aromaticity, number of rings the atom
    belongs to, and membership in rings of size 3 to 8.  Chiral tags other
    than the two tetrahedral ones collapse onto ``CHI_OTHER``.  Only graph
    properties are used, so a 2D template without a conformer is sufficient.

    Parameters
    ----------
    mol : Chem.Mol
        Ligand with explicit atoms in the order the features should follow.

    Returns
    -------
    NDArray[np.int64]
        Array of shape ``(n_atoms, len(LIGAND_ATOM_FEATURE_NAMES))``.
    """
    ring_info = mol.GetRingInfo()
    features = np.empty(
        (mol.GetNumAtoms(), len(LIGAND_ATOM_FEATURE_NAMES)), dtype=np.int64
    )
    for row, atom in zip(features, mol.GetAtoms()):
        index = atom.GetIdx()
        chiral_tag = str(atom.GetChiralTag())
        if chiral_tag not in _CHIRAL_TAGS:
            chiral_tag = "CHI_OTHER"
        row[:] = (
            _bucket(_ATOMIC_NUMBERS, atom.GetAtomicNum()),
            _CHIRAL_TAGS.index(chiral_tag),
            _bucket(_TOTAL_DEGREES, atom.GetTotalDegree()),
            _bucket(_FORMAL_CHARGES, atom.GetFormalCharge()),
            _bucket(_IMPLICIT_VALENCES, atom.GetImplicitValence()),
            _bucket(_HYDROGEN_COUNTS, atom.GetTotalNumHs()),
            _bucket(_RADICAL_ELECTRONS, atom.GetNumRadicalElectrons()),
            _bucket(_HYBRIDIZATIONS, str(atom.GetHybridization())),
            int(atom.GetIsAromatic()),
            _bucket(_RING_COUNTS, ring_info.NumAtomRings(index)),
            *(int(ring_info.IsAtomInRingOfSize(index, size)) for size in _RING_SIZES),
        )
    return features
