# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from openbabel import pybel
from ost import conop, io
from ost import mol as omol
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from plinder.core.structure.smallmols_utils import (
    fix_valency_issues,
    get_matched_template,
    get_matched_template_v2,
    make_rdkit_compatible_mol,
    mol_assigned_bond_orders_by_template,
    params_removeHs,
    uncharge_mol,
)
from plinder.core.utils.constants import BASE_DIR
from plinder.core.utils.log import setup_logger

LOG = setup_logger(__name__)
COMPOUND_LIB = conop.GetDefaultLib()
PRD_LIB = conop.CompoundLib.Load(
    str(BASE_DIR / "data/utils/annotations/static_files/prdcc.chemlib")
)


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
            LOG.warning(f"ligand_ost_ent_to_rdkit_mol: {e}")
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
                LOG.warning(
                    "fix_valency_issues: failed before mol_assigned_bond_orders_by_template"
                )
            # get template from smiles
            template = Chem.MolFromSmiles(ligand_smiles)
            if ligand_num_unresolved_heavy_atoms > 0:
                # update template
                # TODO: could be replaced by get_matched_template_v2
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
            LOG.warning(f"ligand_ost_ent_to_rdkit_mol: {e}")
    try:
        # Fix issues if any
        rdkit_mol = make_rdkit_compatible_mol(rdkit_mol)
        # run uncharger - for consistent protonation
        rdkit_mol = uncharge_mol(rdkit_mol)
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
                # run uncharger - for consistent protonation
                rdkit_mol = uncharge_mol(rdkit_mol)
                return str(Chem.MolToSmiles(rdkit_mol))
            except Exception:
                LOG.warning(
                    "set_smiles_from_ligand_ost: CCD smiles could not be loaded by rdkit, moving to fix"
                )
    rdkit_mol = ligand_ost_ent_to_rdkit_mol(ent)
    return str(Chem.MolToSmiles(rdkit_mol))


def set_smiles_from_ligand_ost_v2(ent: omol.EntityHandle) -> tuple[str, str]:
    input_smiles = ""
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
                template_mol = Chem.MolFromSmiles(str(mol.smiles), sanitize=False)
                template_mol = make_rdkit_compatible_mol(template_mol)
                input_smiles = str(Chem.MolToSmiles(template_mol))
            except Exception:
                LOG.warning(
                    "set_smiles_from_ligand_ost_v2: CCD smiles could not be loaded by RDKit, moving to fix"
                )
    resolved_mol = ligand_ost_ent_to_rdkit_mol(ent)
    if input_smiles:
        matched_template = get_matched_template_v2(template_mol, resolved_mol)
        matched_smiles = Chem.CanonSmiles(Chem.MolToSmiles(matched_template))
    else:
        matched_smiles = Chem.CanonSmiles(str(Chem.MolToSmiles(resolved_mol)))
        # when no reference - define this as a match and the ground truth
        input_smiles = matched_smiles
    return input_smiles, matched_smiles
