# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import numpy as np
import pandas as pd
import pytest
from rdkit import Chem


@pytest.mark.parametrize(
    ["smiles", "num_problems"],
    [
        # positive controls: peppr repairs the issue, so no problems remain
        ["CC(=O)OCCN(C)(C)C", 0],  # AtomValenceException, fixed via onium charge
        ["c1ccnc1", 0],  # KekulizeException, fixed by adding [nH]
        # negative controls: peppr TOLERATES main-group over-valence (does not
        # raise) but cannot neutralize it, so one residual problem is still found
        ["[B]12[B]3[B]4[B]1[B]5[B]6[B]2[B]3[B]45C6", 1],  # boron cage: over-valent B
        ["[Be](C)(C)(C)C", 1],  # over-valent beryllium centre, no 2-centre-bond fix
    ],
)
def test_valence_issue_handling(smiles, num_problems):
    # TODO(peppr): use the vendored sanitize the pipeline uses (boron/main-group
    # over-valence fixes not yet in a released peppr). Revert to
    # `from peppr import sanitize as peppr_sanitize` once upstream.
    from plinder.core.utils.sanitize import sanitize as peppr_sanitize

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    peppr_sanitize(mol)
    problems = Chem.DetectChemistryProblems(mol)
    assert len(problems) == num_problems


@pytest.mark.parametrize(
    "smiles",
    [
        "[B]12[B]3[B]4[B]1[B]5[B]6[B]2[B]3[B]45C6",  # boron cage
        "[Be](C)(C)(C)C",  # over-valent beryllium
    ],
)
def test_mol_from_smiles_tolerates_overvalence(smiles):
    """mol_from_smiles keeps over-valent main-group SMILES that strict rejects.

    Strict ``Chem.MolFromSmiles`` returns None for these, which would silently
    drop the ligand from descriptors/fingerprints; the tolerant helper keeps it.
    """
    from plinder.core.utils.sanitize import mol_from_smiles

    assert Chem.MolFromSmiles(smiles) is None  # strict rejects
    assert mol_from_smiles(smiles) is not None  # tolerant keeps


@pytest.mark.parametrize(
    ["smiles", "expected_stereo"],
    [
        ["O=C(O)/C=C/C(=O)O", Chem.BondStereo.STEREOE],  # fumarate (trans, E)
        ["O=C(O)/C=C\\C(=O)O", Chem.BondStereo.STEREOZ],  # maleate (cis, Z)
    ],
)
def test_mol_from_smiles_perceives_double_bond_ez(smiles, expected_stereo):
    """mol_from_smiles perceives double-bond E/Z, matching Chem.MolFromSmiles.

    The sanitize=False parse skips RDKit's assignStereochemistry finalization, so
    the helper re-runs it — otherwise fumarate and maleate collapse to one
    identity (their C=C would be Unspecified).
    """
    from plinder.core.utils.sanitize import mol_from_smiles

    # the stereogenic C=C (there are also carbonyls, which stay STEREONONE)
    ez = [
        b.GetStereo()
        for b in mol_from_smiles(smiles).GetBonds()
        if b.GetBondType() == Chem.BondType.DOUBLE
        and b.GetStereo() != Chem.BondStereo.STEREONONE
    ]
    assert ez == [expected_stereo]


@pytest.mark.parametrize(
    ["s1", "s2"],
    [
        # E/Z isomers
        ["CC/C=C/Cl", "CC/C=C\\Cl"],
        # R/S enantiomers
        ["C[C@](F)(Cl)CBr", "C[C@@](F)(Cl)CBr"],
    ],
)
def test_smiles2nonstereo(s1, s2):
    # uncharge / InChIKey were retired: ligand identity is now a stereo-stripped
    # RDKit canonical SMILES (robust where InChIKey fails on organometallics).
    from plinder.core.structure.smallmols_similarity import smiles2nonstereo

    out = smiles2nonstereo(s1)
    # a valid canonical SMILES with no stereo markers
    assert out and "@" not in out and "/" not in out and "\\" not in out
    # stereoisomers collapse to one identity
    assert smiles2nonstereo(s1) == smiles2nonstereo(s2)


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
            "ligand_smiles": ["CCO", "CCO", "CCN"],
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
    from biotite.structure import filter_heavy
    from plinder.core.structure.smallmols_utils import compare_stereo_to_template

    # TODO(peppr): use the vendored sanitize the pipeline uses (boron/main-group
    # over-valence fixes not yet in a released peppr). Revert to
    # `from peppr import sanitize as peppr_sanitize` once upstream.
    from plinder.core.utils.sanitize import sanitize as peppr_sanitize

    # Build a CCD mol with stereo (NAG — chiral sugar)
    ref = bt_info.residue("NAG")
    ref_heavy = ref[filter_heavy(ref)]
    ref_heavy.bonds = struc.connect_via_residue_names(ref_heavy)
    template = rdkit_interface.to_mol(ref_heavy)
    peppr_sanitize(template)

    # Resolved mol = same as template (exact match). No manual chiral-tag
    # assignment: compare_stereo_to_template perceives the template's stereo
    # itself (from its coords) and reads only the resolved mol's coordinates.
    resolved = rdkit_interface.to_mol(ref_heavy)
    peppr_sanitize(resolved)
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
    ref_dms_heavy = ref_dms[filter_heavy(ref_dms)]
    ref_dms_heavy.bonds = struc.connect_via_residue_names(ref_dms_heavy)
    dms_mol = rdkit_interface.to_mol(ref_dms_heavy)
    peppr_sanitize(dms_mol)
    dms_template = rdkit_interface.to_mol(ref_dms_heavy)
    peppr_sanitize(dms_template)
    assert (
        compare_stereo_to_template(dms_mol, dms_template) is True
    )  # achiral = no conflict


def test_compare_stereo_to_template_ez():
    """compare_stereo_to_template detects cis/trans (E/Z) double-bond mismatches."""
    from plinder.core.structure.smallmols_utils import compare_stereo_to_template
    from rdkit.Chem import AllChem

    names = ["C1", "C2", "C3", "CL"]

    def _mol(smiles, embed):
        # stamp PDB names positionally so template and resolved map to each other
        m = Chem.RemoveHs(Chem.MolFromSmiles(smiles), sanitize=False)
        for atom, nm in zip(m.GetAtoms(), names):
            info = Chem.AtomPDBResidueInfo()
            info.SetName(nm)
            info.SetResidueName("LIG")
            info.SetResidueNumber(1)
            atom.SetMonomerInfo(info)
        if embed:
            mh = Chem.AddHs(m)
            AllChem.EmbedMolecule(mh, randomSeed=1)
            m = Chem.RemoveHs(mh)
        return m

    trans_template = _mol("C/C=C/Cl", embed=False)  # SMILES tags, no conformer
    match = compare_stereo_to_template(_mol("C/C=C/Cl", embed=True), trans_template)
    flip = compare_stereo_to_template(_mol("C/C=C\\Cl", embed=True), trans_template)
    assert match is True  # same E/Z geometry
    assert flip is False  # inverted (cis vs trans) -> mismatch


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


def test_mhfp6_fingerprint_is_deterministic_and_sized():
    from plinder.core.structure.smallmols_similarity import (
        MHFP6_N_PERMUTATIONS,
        mol2mhfp6,
    )

    fp_a = mol2mhfp6("c1ccccc1O")
    fp_b = mol2mhfp6("c1ccccc1O")
    assert fp_a.shape == (MHFP6_N_PERMUTATIONS,)
    assert str(fp_a.dtype) == "uint32"
    # a fixed permutation seed makes the MinHash vector reproducible
    assert (fp_a == fp_b).all()
    with pytest.raises(ValueError):
        mol2mhfp6("not a molecule")


def test_mhfp6_bulk_jaccard_counts_matching_minima():
    import numpy as np
    from plinder.core.structure.smallmols_similarity import (
        MHFP6_N_PERMUTATIONS,
        mhfp6_bulk_jaccard,
        mol2mhfp6,
    )

    # self-similarity is exact; a distinct molecule scores strictly below 1
    query = mol2mhfp6("c1ccccc1O")
    other = mol2mhfp6("c1ccccc1N")
    similarities = mhfp6_bulk_jaccard(query, np.stack([query, other]))
    assert similarities[0] == 1.0
    assert 0.0 <= similarities[1] < 1.0

    # the estimate is by definition the fraction of permutations whose minima
    # agree: break a controlled 1/4 of the positions and expect exactly 0.75.
    # (RDKit's MHFPEncoder.Distance is not used here — the C++ argument type it
    # accepts is registered differently across RDKit versions.)
    base = np.arange(MHFP6_N_PERMUTATIONS, dtype=np.uint32)
    changed = base.copy()
    changed[: MHFP6_N_PERMUTATIONS // 4] = np.uint32(2**32 - 1)
    assert mhfp6_bulk_jaccard(base, np.stack([changed]))[0] == 0.75


def test_is_chemical_cluster_metric_covers_ecfp4_and_mhfp6():
    from plinder.core.scores.metrics import (
        CHEMICAL_CLUSTER_METRICS,
        is_chemical_cluster_metric,
    )

    assert set(CHEMICAL_CLUSTER_METRICS) == {
        "tanimoto_similarity_ecfp4_1024",
        "jaccard_similarity_mhfp6_2048",
    }
    assert is_chemical_cluster_metric("jaccard_similarity_mhfp6_2048")
    assert is_chemical_cluster_metric("tanimoto_similarity_ecfp4_1024")
    assert not is_chemical_cluster_metric("pocket_qcov")


def test_ecfp4_vs_mhfp6_on_sequence_isomer_edge_cases():
    """ECFP4 and MHFP6 both handle monomer-order isomers the same way.

    Peptides and nucleic acids are built with RDKit's sequence->mol conversion.
    Two edge cases are probed: a genuine single-monomer change (should stay
    similar) versus a sequence permutation of the *same* monomers (a distinct
    molecule that a local-substructure fingerprint struggles to tell apart).
    """
    import numpy as np
    from plinder.core.structure.smallmols_similarity import (
        mhfp6_bulk_jaccard,
        mol2mhfp6,
        mol2morgan_fp,
    )
    from rdkit import DataStructs

    def ecfp4(m1, m2):
        return DataStructs.TanimotoSimilarity(
            mol2morgan_fp(m1, radius=2, nbits=1024),
            mol2morgan_fp(m2, radius=2, nbits=1024),
        )

    def mhfp6(m1, m2):
        return float(mhfp6_bulk_jaccard(mol2mhfp6(m1), np.stack([mol2mhfp6(m2)]))[0])

    # ---- peptides (flavor=0: L-amino-acid chain) ----
    gya = Chem.MolFromSequence("GYA", flavor=0)  # Gly-Tyr-Ala
    gfa = Chem.MolFromSequence("GFA", flavor=0)  # Gly-Phe-Ala (Tyr->Phe: one -OH)
    agf = Chem.MolFromSequence("AGF", flavor=0)  # same residues, permuted order
    assert all(m is not None for m in (gya, gfa, agf))
    # GFA and AGF are genuinely different molecules, not the same input twice
    assert Chem.MolToSmiles(gfa) != Chem.MolToSmiles(agf)

    # the single-residue change (GYA/GFA) is more similar than the permutation
    # (GFA/AGF) under BOTH fingerprints, by a clear margin
    ecfp_margin = ecfp4(gya, gfa) - ecfp4(gfa, agf)
    mhfp_margin = mhfp6(gya, gfa) - mhfp6(gfa, agf)
    assert ecfp4(gya, gfa) > ecfp4(gfa, agf)
    assert mhfp6(gya, gfa) > mhfp6(gfa, agf)
    assert ecfp_margin > 0.1 and mhfp_margin > 0.1
    # MHFP6 does not separate this edge case better: the margins are comparable
    assert abs(mhfp_margin - ecfp_margin) < 0.05

    # ---- nucleic acids (flavor=6: DNA) ----
    gac = Chem.MolFromSequence("GAC", flavor=6)
    ggc = Chem.MolFromSequence("GGC", flavor=6)  # A->G: single-base change
    acg = Chem.MolFromSequence("ACG", flavor=6)  # same bases, permuted order
    assert all(m is not None for m in (gac, ggc, acg))
    # distinct molecules with identical base composition
    assert Chem.MolToSmiles(gac) != Chem.MolToSmiles(acg)

    # bases sit ~10+ bonds apart along the backbone, far beyond ECFP4 radius 2
    # and MHFP6 radius 3, so neither fingerprint can see base ORDER: a pure
    # permutation collides at ~1.0 for BOTH, and even scores at least as high as
    # a real single-base substitution. MHFP6 does not fix this blind spot.
    assert ecfp4(gac, acg) > 0.99
    assert mhfp6(gac, acg) > 0.99
    assert ecfp4(gac, acg) >= ecfp4(gac, ggc)
    assert mhfp6(gac, acg) >= mhfp6(gac, ggc)


def test_ligand_atom_features_are_hand_checkable_vocabulary_indices():
    from plinder.core.structure.smallmols_utils import (
        LIGAND_ATOM_FEATURE_NAMES,
        ligand_atom_features,
    )

    benzene = ligand_atom_features(Chem.MolFromSmiles("c1ccccc1"))

    assert benzene.shape == (6, len(LIGAND_ATOM_FEATURE_NAMES))
    assert benzene.dtype == np.int64
    # Aromatic carbon: Z=6 -> 5, unspecified chirality, degree 3, charge 0 -> 5,
    # one implicit H, one H, no radicals, SP2 -> 1, aromatic, in one ring, and
    # only the six-membered ring flag set.
    assert benzene.tolist() == [[5, 0, 3, 5, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0]] * 6

    alanine = ligand_atom_features(Chem.MolFromSmiles("C[C@H](N)C(=O)O"))
    assert alanine.tolist() == [
        [5, 0, 4, 5, 3, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],  # methyl carbon
        [5, 2, 4, 5, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],  # CCW alpha carbon
        [6, 0, 3, 5, 2, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],  # amine nitrogen
        [5, 0, 3, 5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # carbonyl carbon
        [7, 0, 1, 5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # carbonyl oxygen
        [7, 0, 2, 5, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # hydroxyl oxygen
    ]

    cyclopropane = ligand_atom_features(Chem.MolFromSmiles("C1CC1"))
    assert (
        cyclopropane.tolist() == [[5, 0, 4, 5, 2, 2, 0, 2, 0, 1, 1, 0, 0, 0, 0, 0]] * 3
    )

    # Hexafluorophosphate: P is Z=15 -> 14, degree 6, charge -1 -> 4, SP3D2 -> 4.
    phosphorus = ligand_atom_features(Chem.MolFromSmiles("F[P-](F)(F)(F)(F)F"))[1]
    assert phosphorus.tolist() == [14, 0, 6, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0]


def test_ligand_atom_features_bucket_out_of_vocabulary_values():
    from plinder.core.structure.smallmols_utils import ligand_atom_features

    # A dummy atom has Z=0 and no hybridization: both land in the "other"
    # bucket after the 118 elements and the five hybridizations.
    dummy = ligand_atom_features(Chem.MolFromSmiles("*C"))[0]
    assert dummy.tolist() == [118, 0, 1, 5, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0]

    ammonium = ligand_atom_features(Chem.MolFromSmiles("[NH4+]"))
    assert ammonium.tolist() == [[6, 0, 4, 6, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0]]

    # Chiral tags outside the tetrahedral pair collapse onto CHI_OTHER.
    exotic = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
    exotic.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL)
    assert ligand_atom_features(exotic)[1, 1] == 3

    assert ligand_atom_features(Chem.Mol()).shape == (0, 16)
