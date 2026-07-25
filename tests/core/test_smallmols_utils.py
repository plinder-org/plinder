# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
import pytest
from rdkit import Chem


@pytest.mark.parametrize(
    ["smiles", "num_problems"],
    [
        ["CC(=O)OCCN(C)(C)C", 0],  # AtomValenceException
        ["c1ccnc1", 0],  # KekulizeException
    ],
)
def test_valence_issue_handling(smiles, num_problems):
    from peppr import sanitize as peppr_sanitize

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    peppr_sanitize(mol)
    problems = Chem.DetectChemistryProblems(mol)
    assert len(problems) == num_problems


@pytest.mark.parametrize(
    ["smiles", "num_charged_atoms"],
    [
        ["CC(=O)OCCN(C)(C)C", 0],
        ["[O-]C(=[O])CC[NH+](C)(C)", 0],
        ["[O-]C(=[OH+])CC[NH+](C)(C)", 0],
        ["OC(=[O])CC[N+](C)(C)(C)", 1],
        ["[O-]C(=[O])CC[N+](C)(C)(C)", 2],
        ["[O-]C(=[O])C.C[N+](C)(C)(C)", 2],
    ],
)
def test_uncharge_mol(smiles, num_charged_atoms):
    from plinder.core.structure.smallmols_utils import uncharge_mol

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = uncharge_mol(mol)
    assert (
        sum([at.GetFormalCharge() != 0 for at in mol.GetAtoms()]) == num_charged_atoms
    )


@pytest.mark.parametrize(
    ["smiles", "inchikey", "remove_stereo"],
    [
        ["CC/C=C/Cl", "DUDKKPVINWLFBI-ONEGZZNKSA-N", False],
        ["CC/C=C\\Cl", "DUDKKPVINWLFBI-ARJAWSKDSA-N", False],
        ["CC/C=C/Cl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", True],
        ["CC/C=C\\Cl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", True],
        ["CCC=CCl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", False],
        ["CCC=CCl", "DUDKKPVINWLFBI-UHFFFAOYSA-N", True],
        ["C[C@@](F)(Cl)CBr", "REKDFINPOZVXJS-VKHMYHEASA-N", False],
        ["C[C@](F)(Cl)CBr", "REKDFINPOZVXJS-GSVOUGTGSA-N", False],
        ["C[C@](F)(Cl)CBr", "REKDFINPOZVXJS-UHFFFAOYSA-N", True],
        ["C[C@@](F)(Cl)CBr", "REKDFINPOZVXJS-UHFFFAOYSA-N", True],
        ["CC(F)(Cl)CBr", "REKDFINPOZVXJS-UHFFFAOYSA-N", True],
    ],
)
def test_inchikey(smiles, inchikey, remove_stereo):
    from plinder.core.structure.smallmols_similarity import smiles2inchikey

    assert inchikey == smiles2inchikey(smiles, remove_stereo=remove_stereo)


def test_load_ligands_from_index_uses_proper_holo_ligand_rows():
    from plinder.data.annotations.get_similarity_scores import (
        load_ligands_from_index,
    )

    annotation = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "2def"],
            "system_id": ["1abc__1__A", "1abc__1__A", "2def__1__B"],
            "system_type": ["holo", "holo", "apo"],
            "ligand_is_proper": [True, False, True],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCO", "CCN"],
            "ligand_unique_ccd_code": ["LIG", "LIG", "OTH"],
            "ligand_id": ["1abc__1__A__1.C", "1abc__1__A__1.C", "2def__1__B__1.D"],
            "ligand_asym_id": ["C", "C", "D"],
        }
    )

    ligands = load_ligands_from_index(annotation=annotation)

    assert ligands["ligand_id"].tolist() == ["1abc__1__A__1.C"]
    assert ligands["pdb_id"].tolist() == ["1abc"]
    assert ligands["ligand_ccd_code"].tolist() == ["LIG"]
    assert "inchikeys" not in ligands.columns


def test_compare_stereo_to_template():
    """Test compare_stereo_to_template: match, mismatch, achiral."""
    import biotite.structure as struc
    import biotite.structure.info as bt_info
    from biotite.interface import rdkit as rdkit_interface
    from peppr import sanitize as peppr_sanitize
    from plinder.core.structure.atoms import is_hydrogen_isotope
    from plinder.core.structure.smallmols_utils import compare_stereo_to_template

    # Build a CCD mol with stereo (NAG — chiral sugar)
    ref = bt_info.residue("NAG")
    ref_heavy = ref[~is_hydrogen_isotope(ref.element)]
    ref_heavy.bonds = struc.connect_via_residue_names(ref_heavy)
    template = rdkit_interface.to_mol(ref_heavy)
    peppr_sanitize(template)
    Chem.AssignAtomChiralTagsFromStructure(template)

    # Resolved mol = same as template (exact match)
    resolved = rdkit_interface.to_mol(ref_heavy)
    peppr_sanitize(resolved)
    Chem.AssignAtomChiralTagsFromStructure(resolved)
    assert compare_stereo_to_template(resolved, template) is True

    # Invert the 3D geometry (improper reflection through the x=0 plane) →
    # enantiomer → mismatch. compare_stereo_to_template judges stereo from
    # coordinates, so flipping a chiral *tag* without moving atoms would be a
    # no-op; the geometry is the source of truth, so we must move atoms.
    from rdkit.Geometry import Point3D

    flipped = Chem.Mol(resolved)
    conf = flipped.GetConformer()
    for i in range(flipped.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(-p.x, p.y, p.z))
    assert compare_stereo_to_template(flipped, template) is False

    # Achiral mol (DMS — no stereocenters)
    ref_dms = bt_info.residue("DMS")
    ref_dms_heavy = ref_dms[~is_hydrogen_isotope(ref_dms.element)]
    ref_dms_heavy.bonds = struc.connect_via_residue_names(ref_dms_heavy)
    dms_mol = rdkit_interface.to_mol(ref_dms_heavy)
    peppr_sanitize(dms_mol)
    dms_template = rdkit_interface.to_mol(ref_dms_heavy)
    peppr_sanitize(dms_template)
    assert (
        compare_stereo_to_template(dms_mol, dms_template) is True
    )  # achiral = no conflict


def test_sequences_match_core():
    """Test sequence matching for binding affinity validation."""
    from plinder.data.annotations.protein_utils import sequences_match_core

    assert sequences_match_core("ABCDEFGH", "ABCDEFGH") is True
    assert sequences_match_core("MHHHHHABCDEFGH", "ABCDEFGH") is True
    assert sequences_match_core("ABCDEFGHLEVLFQ", "ABCDEFGH") is True
    assert sequences_match_core("ABXDEFGH", "ABCDEFGH") is False
    assert sequences_match_core("BCDEFG", "ABCDEFGH") is True
    assert sequences_match_core("", "ABCDEFGH") is False
    assert sequences_match_core("AB", "ABCDEFGHIJKLMNOP") is False


def test_matched_templates():
    from plinder.core.structure.smallmols_utils import (
        get_matched_template,
        mol_assigned_bond_orders_by_template,
    )

    mol1 = Chem.MolFromSmiles("FC(Cl)(Br)C.CNCC1CCCCC1.CCC(OC)O")
    template = Chem.MolFromSmiles("F[C@@](Br)(Cl)CCCNCc1cc(C(=O)N/C=C/C(OC)=O)ccc1")
    matched_template = get_matched_template(template, mol1)
    fixed_mol = mol_assigned_bond_orders_by_template(matched_template, mol1)
    fixed_mol_SMILES = Chem.CanonSmiles(Chem.MolToSmiles(fixed_mol))
    assert fixed_mol_SMILES.count("=") >= 2
    assert fixed_mol_SMILES == "C=CC(=O)OC.CC(F)(Cl)Br.CNCc1ccccc1"
