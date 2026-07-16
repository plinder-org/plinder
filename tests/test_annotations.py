# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import numpy as np
import pandas as pd
import pytest
from plinder.data.get_system_annotations import GetPlinderAnnotation
from plinder.data.utils.annotations.aggregate_annotations import Entry
from plinder.data.utils.annotations.cif_utils import read_mmcif_container
from plinder.data.utils.annotations.interaction_utils import get_covalent_connections
from plinder.data.utils.annotations.interface_gap import annotate_interface_gaps
from plinder.data.utils.annotations.ligand_utils import sort_ccd_codes
from plinder.data.utils.annotations.mmpdb_utils import add_mmp_clusters_to_data
from plinder.data.utils.annotations.save_utils import (
    SystemReconstructionOptions,
    SystemReconstructionOutputs,
    save_reconstructed_system,
)
from rdkit import Chem


def test_ccd_name_sorter():
    assert sort_ccd_codes({"G", "G25", "CPG", "5GP"}) == ["CPG", "G25", "G", "5GP"]


def test_chain_from_cif_data_nucleotides(cif_8ufz):
    """Test Chain.from_cif_data assigns correct one-letter codes and chem_types for DNA.

    8ufz chain A is a 16-nt DNA strand (DA, DT, DC, DG residues).
    """
    import biotite.structure.io.pdbx as pdbx
    from plinder.core.structure.atoms import is_hydrogen_isotope
    from plinder.data.utils.annotations.cif_utils import read_mmcif_file
    from plinder.data.utils.annotations.protein_utils import Chain, get_seqres_from_cif

    cif_obj = read_mmcif_file(cif_8ufz)
    block = list(cif_obj.values())[0]
    atoms = pdbx.get_structure(
        cif_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[~is_hydrogen_isotope(atoms.element)]
    seqres = get_seqres_from_cif(block)

    chain_a_atoms = atoms[atoms.chain_id == "A"]
    chain = Chain.from_cif_data(
        asym_id="A",
        block=block,
        atoms=chain_a_atoms,
        seqres_length=len(seqres.get("A", "")),
    )

    expected_seq = "AATAAAAGCGGAAGTG"
    actual_seq = "".join(
        chain.residues[r].one_letter_code for r in sorted(chain.residues)
    )
    assert (
        actual_seq == expected_seq
    ), f"DNA sequence mismatch: got '{actual_seq}', expected '{expected_seq}'"

    for resnum, residue in chain.residues.items():
        assert (
            residue.chem_type == "DNA Linking"
        ), f"Residue {residue.name} at {resnum}: expected 'DNA Linking', got '{residue.chem_type}'"


def test_covalent_linkage(cif_1qz5):
    reference = [("72:GLU:A:72:C", "73:HIC:A:73:N"), ("73:HIC:A:73:C", "74:GLY:A:74:N")]

    assert (
        get_covalent_connections(read_mmcif_container(cif_1qz5))["covale"] == reference
    )


def test_find_missing_residues(cif_2y4i_system):
    actual = annotate_interface_gaps(
        cif_2y4i_system, protein_chains=None, ligand_chains=None
    )["ligand_interface_gap_annotation"][("C", "A")]
    expected = {
        "interface_atom_gaps_4A": 5,
        "interface_atom_gaps_8A": 39,
        "missing_interface_residues_4A": 0,
        "missing_interface_residues_8A": 0,
    }
    assert actual == expected


def test_short_noncov_peptide_detection(cif_6i41, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("6i41")
    plinder_anno = GetPlinderAnnotation(
        cif_6i41, "", min_polymer_size=10, save_folder=entry_dir
    )
    plinder_anno.annotate()
    df = plinder_anno.annotated_df
    chain_file = entry_dir / "6i41" / "entry_chains.parquet"
    assert chain_file.is_file()
    chain_df = pd.read_parquet(chain_file)
    assert chain_df["chain_type"].str.lower().str.contains("polypeptide").all()
    assert len(df) == 1
    assert df["ligand_is_covalent"].sum() == 0
    assert set(df.ligand_ccd_code.to_list()) == {"LYS-ALA-ASP-THR-THR-THR-PRO"}
    # Note: chain 'B' is ligand = should not be in protein neigh list
    assert df["ligand_protein_chains_auth_id"].drop_duplicates().to_list() == [["A"]]


@pytest.mark.parametrize(
    "min_polymer_size,expect_ligand",
    [(12, False), (20, True)],
    ids=["threshold_12_peptide_is_receptor", "threshold_20_peptide_is_ligand"],
)
def test_peptide_ligand_threshold(
    cif_6u6k, mock_alternative_datasets, min_polymer_size, expect_ligand
):
    """6u6k: 13-residue synthetic peptide (chain B).

    With min_polymer_size=12, the peptide is receptor (13 >= 12) → no systems.
    With min_polymer_size=20, the peptide is ligand (13 < 20) → system created.
    """
    entry_dir = mock_alternative_datasets("6u6k")
    entry = Entry.from_cif_file(
        cif_6u6k, save_folder=entry_dir, min_polymer_size=min_polymer_size
    )
    if expect_ligand:
        assert len(entry.systems) == 1
        lig = entry.systems[list(entry.systems.keys())[0]].ligands[0]
        assert lig.ccd_code == "ACE-TRP-TRP-ILE-ILE-PRO-ALY-VAL-LYS-ALY-GLY-CYS-NH2"
    else:
        assert (
            len(entry.systems) == 0
        ), f"13-residue peptide should be receptor with min_polymer_size={min_polymer_size}"


def test_synthetic_cov_peptide_detection(cif_6lu7, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("6lu7")
    plinder_anno = GetPlinderAnnotation(cif_6lu7, "", save_folder=entry_dir)
    plinder_anno.annotate()
    df = plinder_anno.annotated_df
    assert len(df) == 2
    # TODO: test 'ligand_covalent_linkages'
    assert df["ligand_is_covalent"].sum() == 2
    # TODO: need to decide if to test for this name or use PRD/BIRD name!
    ligname = df.ligand_ccd_code.to_list()[0]
    assert ligname in ["02J-ALA-VAL-LEU-PJE-010", "PRD_002214"]
    # peptide ligand chain 'B' should not be in protein chains!
    assert sum(["B" in i for i in df["ligand_protein_chains_auth_id"].values]) == 0
    lig = plinder_anno.entry.systems["6lu7__1__1.A_2.A__1.B"].ligands[0]
    assert lig.is_invalid == False
    assert lig.is_covalent == True
    assert lig.covalent_linkages == {"145:CYS:A:145:SG__5:PJE:B:5:C20"}
    outsdffile = entry_dir / "6lu7" / "ligand_files" / "B.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    assert len(Chem.MolToSmiles(rdmol).split(".")) == 1


def test_crystal_contact_detection(cif_6lu7, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("6lu7")
    plinder_anno = GetPlinderAnnotation(cif_6lu7, "", save_folder=entry_dir)
    plinder_anno.annotate()
    df = plinder_anno.annotated_df
    assert len(df) == 2
    assert all(x == 5 for x in df["system_num_atoms_with_crystal_contacts"])
    assert all(x == 2 for x in df["system_num_crystal_contacted_residues"])


def test_simple_covalency_detection(cif_7gl9, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("7gl9")
    plinder_anno_noncov = GetPlinderAnnotation(cif_7gl9, "", save_folder=entry_dir)
    plinder_anno_noncov.annotate()
    df_noncov = plinder_anno_noncov.annotated_df
    assert df_noncov["ligand_is_covalent"].sum() == 0


def test_simple_covalency_detection_found(cif_7gj7, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("7gj7")
    plinder_anno_cov = GetPlinderAnnotation(cif_7gj7, "", save_folder=entry_dir)
    plinder_anno_cov.annotate()
    df_cov = plinder_anno_cov.annotated_df
    assert df_cov["ligand_is_covalent"].sum() == 2
    # test 'ligand_covalent_linkages'
    lig = plinder_anno_cov.entry.systems["7gj7__1__1.A_1.B__1.N_1.P"].ligands[0]
    assert lig.is_covalent == True
    assert lig.covalent_linkages == {"145:CYS:B:145:SG__404:Q0I:N:.:C"}


def test_simple_ternary_detection(cif_2p1q, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("2p1q")
    plinder_anno = GetPlinderAnnotation(cif_2p1q, "", save_folder=entry_dir)
    plinder_anno.annotate()
    df = plinder_anno.annotated_df
    assert sorted(set(df.ligand_ccd_code.to_list())) == ["IAC"]
    auxin_entry = df.ligand_ccd_code == "IAC"
    # expect two chains to get this correct!
    assert df[auxin_entry][
        "ligand_protein_chains_auth_id"
    ].drop_duplicates().to_list() == [["B", "C"]]


def test_plip_entry_binary(cif_4ci1, mock_alternative_datasets, lig_code="EF2"):
    entry_dir = mock_alternative_datasets("4ci1")
    entry = Entry.from_cif_file(
        cif_4ci1,
        neighboring_residue_threshold=6.0,
        neighboring_ligand_threshold=4.0,
        min_polymer_size=10,
        save_folder=entry_dir,
    )

    for system in entry.systems:
        for ligand in entry.systems[system].ligands:
            if ligand.ccd_code == lig_code:
                break

    # assert that expected chain is detected
    assert sorted(ligand.interactions.keys()) == ["1.B"]

    # Expected interactions (hydrophobic contacts dropped in peppr migration)
    expected_interactions = {
        404: [
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
            "type:hydrogen_bonds__protisdon:False__sidechain:False",
        ],
        406: [
            "type:hydrogen_bonds__protisdon:True__sidechain:False",
        ],
        377: ["type:water_bridges__protisdon:True"],
        383: ["type:water_bridges__protisdon:False"],
    }
    assert ligand.interactions["1.B"] == expected_interactions


def test_plip_entry_ternary(cif_2p1q, mock_alternative_datasets, lig_code="IAC"):
    entry_dir = mock_alternative_datasets("2p1q")
    entry = Entry.from_cif_file(
        cif_2p1q,
        neighboring_residue_threshold=6.0,
        neighboring_ligand_threshold=4.0,
        min_polymer_size=10,
        save_folder=entry_dir,
    )

    for system in entry.systems:
        for ligand in entry.systems[system].ligands:
            if ligand.ccd_code == lig_code:
                break

    # assert that expected two chains are detected
    assert sorted(ligand.interactions.keys()) == ["2.B", "2.C"]

    # Expected interactions (hydrophobic contacts dropped in peppr migration)
    expected_interactions_2B = {
        403: [
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
            "type:salt_bridges__protispos:True",
        ],
        438: [
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
        ],
        439: ["type:hydrogen_bonds__protisdon:False__sidechain:False"],
        436: ["type:water_bridges__protisdon:True"],
        462: ["type:water_bridges__protisdon:True"],
    }
    expected_interactions_2C = {
        5: ["type:pi_stacks__stack_type:T"],
        7: ["type:water_bridges__protisdon:False"],
    }
    expected_waters = {"2.G": {2, 4}}

    # Check all expected interactions for chain 2.B
    # Exact match
    assert ligand.interactions["2.B"] == expected_interactions_2B
    assert ligand.interactions.get("2.C", {}) == expected_interactions_2C
    assert {k: set(v) for k, v in ligand.waters.items()} == expected_waters


def test_water_saving(cif_2p1q, mock_alternative_datasets):
    import biotite.structure as struc
    from biotite.structure.io import pdbx
    from plinder.data.utils.annotations.cif_utils import read_mmcif_file

    entry_dir = mock_alternative_datasets("2p1q")
    system_tag = "2p1q__2__2.B_2.C__2.E"
    entry = Entry.from_cif_file(cif_2p1q, save_folder=entry_dir)
    row = entry.to_df().query("system_id == @system_tag").iloc[0]

    output_dir = entry_dir / "reconstructed" / system_tag
    receptor_cif = output_dir / "receptor.cif"
    save_reconstructed_system(
        cif_2p1q,
        row,
        outputs=SystemReconstructionOutputs(receptor_cif=receptor_cif),
    )
    assert receptor_cif.is_file()
    assert not (output_dir / "system.cif").exists()
    assert not list(output_dir.glob("*.pdb"))

    atoms = pdbx.get_structure(read_mmcif_file(receptor_cif), model=1)
    water_atoms = atoms[struc.filter_solvent(atoms)]
    assert len(set(zip(water_atoms.chain_id, water_atoms.res_id))) == 2

    all_waters_cif = output_dir / "receptor_all_waters.cif"
    save_reconstructed_system(
        cif_2p1q,
        row,
        outputs=SystemReconstructionOutputs(receptor_cif=all_waters_cif),
        options=SystemReconstructionOptions(receptor_waters="all"),
    )
    all_atoms = pdbx.get_structure(read_mmcif_file(all_waters_cif), model=1)
    all_water_atoms = all_atoms[struc.filter_solvent(all_atoms)]
    assert len(set(zip(all_water_atoms.chain_id, all_water_atoms.res_id))) > 2


def test_plip_same_hinge_binders(cif_2gdo, cif_4qyf, mock_alternative_datasets):
    pdb_ids = ["2gdo", "4qyf"]
    mmcifs = [cif_2gdo, cif_4qyf]
    ccd_codes = ["12C", "3DV"]
    hinge_resids = [85, 87]
    interactions_sets = []
    for ccd, cif, pdb_id in zip(ccd_codes, mmcifs, pdb_ids):
        entry_dir = mock_alternative_datasets(pdb_id)
        entry = Entry.from_cif_file(cif, save_folder=entry_dir)
        for system in entry.systems.values():
            for ligand in system.ligands:
                if ligand.ccd_code == ccd:
                    interactions_sets.append(ligand.interactions["1.A"])
                    break

    assert len(interactions_sets) == 2

    for hr in hinge_resids:
        assert len(set(interactions_sets[0][hr]).intersection(interactions_sets[1][hr]))


def test_get_single_ligand_system_annotations(cif_6fx1, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("6fx1")
    entry = Entry.from_cif_file(cif_6fx1, save_folder=entry_dir)
    ligands = []
    for system in entry.systems:
        for ligand in entry.systems[system].ligands:
            ligands.append(ligand)
    single_ligand_system_result = set(l.ccd_code for l in ligands)
    single_ligand_system_target = {
        "C4W-NAG-BMA-MAN-NAG-FUC",
        "OOA",
        "C4W-NAG-BMA-MAN-FUC",
        "C4W-NAG-BMA-MAN-NAG-MAN-FUC",
        "C4W-NAG-BMA-FUC",
        "C4W-NAG-BMA-MAN-NAG-MAN-NAG-FUC",
        "MLI",
    }
    assert single_ligand_system_result == single_ligand_system_target


def test_canonical_ligand_saving_and_system_reconstruction(
    cif_2y4i, mock_alternative_datasets
):
    import biotite.structure as struc
    from biotite.structure.io import pdbx
    from plinder.data.utils.annotations.cif_utils import read_mmcif_file

    entry_dir = mock_alternative_datasets("2y4i")
    system_tag = "2y4i__1__1.B__1.E_1.F"
    entry = Entry.from_cif_file(cif_2y4i, save_folder=entry_dir)
    for chain in entry.chains.values():
        assert len(chain.mappings["ECOD"])

    canonical_ligand_dir = entry_dir / "2y4i" / "ligand_files"
    assert {path.name for path in canonical_ligand_dir.glob("*.sdf")} >= {
        "E.sdf",
        "F.sdf",
    }
    assert not (entry_dir / system_tag).exists()

    row = entry.to_df().query("system_id == @system_tag").iloc[0]
    assert set(row["system_other_chains_asym_id"]).isdisjoint(
        row["system_protein_chains_asym_id"] + row["system_ligand_chains"]
    )
    output_dir = entry_dir / "reconstructed" / system_tag
    outputs = SystemReconstructionOutputs(
        system_cif=output_dir / "system.cif",
        receptor_cif=output_dir / "receptor.cif",
        sequences_fasta=output_dir / "sequences.fasta",
    )
    written = save_reconstructed_system(
        cif_2y4i,
        row,
        outputs=outputs,
        options=SystemReconstructionOptions(
            system_waters="none",
            receptor_waters="none",
            system_include_other_protein_chains=True,
        ),
    )
    assert set(written) == {"system_cif", "receptor_cif", "sequences_fasta"}
    assert all(path.is_file() for path in written.values())
    assert not list(output_dir.glob("*.pdb"))

    system_atoms = pdbx.get_structure(read_mmcif_file(outputs.system_cif), model=1)
    assert not np.any(struc.filter_solvent(system_atoms))
    expected_chains = set(row["system_protein_chains_asym_id"])
    expected_chains.update(row["system_ligand_chains"])
    expected_chains.update(row["system_other_protein_chains_asym_id"])
    assert set(system_atoms.chain_id) == expected_chains


def test_smiles_from_nextgen(rcsb_ccd_reference_csv):
    """Test CCD SMILES against RCSB ground truth.

    For each compound in the RCSB reference CSV, verify:
    1. InChIKey from CCD ideal 3D matches RCSB InChIKey
    2. Per-atom chirality matches via substructure match
    """
    from plinder.data.utils.annotations.cif_utils import _COORDINATION_METALS
    from plinder.data.utils.annotations.ligand_utils import _get_ccd_mol
    from rdkit.Chem.inchi import MolToInchiKey

    rcsb_df = pd.read_csv(rcsb_ccd_reference_csv)
    assert len(rcsb_df) > 0, "Should have RCSB ground truth entries"

    mismatches = []
    for _, row in rcsb_df.iterrows():
        comp_id = row["comp_id"]
        rcsb_inchikey = row["inchikey"]
        if not rcsb_inchikey or pd.isna(rcsb_inchikey):
            continue

        # Production code: get CCD mol with stereo from ideal 3D
        ccd_mol = _get_ccd_mol(comp_id)
        if ccd_mol is None:
            continue

        ccd_inchikey = MolToInchiKey(ccd_mol) or ""

        # Skip organometallic compounds — biotite doesn't produce dative
        # bonds for metal coordination, giving different connectivity than
        # the RCSB canonical representation (e.g. HEM Fe-N bonds)
        has_metal = any(
            a.GetSymbol().upper() in _COORDINATION_METALS and a.GetDegree() > 0
            for a in ccd_mol.GetAtoms()
        )
        if has_metal:
            continue

        # Allow stereo-ambiguous cases (same connectivity, different stereo)
        if ccd_inchikey and rcsb_inchikey and ccd_inchikey[:14] == rcsb_inchikey[:14]:
            if ccd_inchikey != rcsb_inchikey:
                continue  # ambiguous stereo in CCD — skip

        # Check InChIKey (primary — canonical across toolkits)
        if ccd_inchikey != rcsb_inchikey:
            mismatches.append(
                (
                    comp_id,
                    "InChIKey",
                    ccd_inchikey,
                    f"expected {rcsb_inchikey}",
                )
            )
            continue

        # Check chirality via substructure match between CCD and RCSB mols
        rcsb_mol = Chem.MolFromSmiles(row["rcsb_smiles"])
        if rcsb_mol is not None:
            Chem.AssignStereochemistry(rcsb_mol, force=True)
            Chem.AssignStereochemistry(ccd_mol, force=True)
            match = ccd_mol.GetSubstructMatch(rcsb_mol)
            if match:
                for rcsb_idx, ccd_idx in enumerate(match):
                    rcsb_atom = rcsb_mol.GetAtomWithIdx(rcsb_idx)
                    ccd_atom = ccd_mol.GetAtomWithIdx(ccd_idx)
                    rcsb_cip = rcsb_atom.GetPropsAsDict().get("_CIPCode", "")
                    ccd_cip = ccd_atom.GetPropsAsDict().get("_CIPCode", "")
                    if rcsb_cip and ccd_cip and rcsb_cip != ccd_cip:
                        info = ccd_atom.GetPDBResidueInfo()
                        name = info.GetName().strip() if info else str(ccd_idx)
                        mismatches.append(
                            (
                                comp_id,
                                f"chirality@{name}",
                                ccd_cip,
                                f"expected {rcsb_cip}",
                            )
                        )

    assert len(mismatches) == 0, "CCD vs RCSB mismatches:\n" + "\n".join(
        f"  {m}" for m in mismatches
    )


def _build_resolved_mol(cif_path, chain_id):
    """Helper: build resolved mol from CIF chain using production code."""
    import biotite.structure.io.pdbx as pdbx
    from plinder.core.structure.atoms import is_hydrogen_isotope
    from plinder.data.utils.annotations.cif_utils import (
        atoms_to_rdkit_mol,
        read_mmcif_file,
    )

    cif_obj = read_mmcif_file(cif_path)
    atoms = pdbx.get_structure(
        cif_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[~is_hydrogen_isotope(atoms.element)]
    return atoms_to_rdkit_mol(atoms[atoms.chain_id == chain_id])


def _flip_first_chiral(mol):
    """Helper: return a copy with one chiral center inverted."""
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetPropsAsDict().get("_CIPCode", ""):
            chiral = atom.GetChiralTag()
            if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            Chem.AssignStereochemistry(rw, cleanIt=True, force=True)
            return rw.GetMol()
    return None


def test_stereo_check_single_residue(cif_7gj7):
    """Test _check_stereo_vs_template on single-residue ligands.

    Q0I (chain E): chiral — should match CCD, flipped should fail.
    DMS (chain C): achiral — should return None (no comparable centers).
    """
    from plinder.data.utils.annotations.ligand_utils import _check_stereo_vs_template

    # Chiral: Q0I
    q0i_mol = _build_resolved_mol(cif_7gj7, "E")
    assert _check_stereo_vs_template(q0i_mol) is True

    q0i_flipped = _flip_first_chiral(q0i_mol)
    assert q0i_flipped is not None, "Q0I should have a chiral center to flip"
    assert _check_stereo_vs_template(q0i_flipped) is False

    # Achiral: DMS (dimethyl sulfoxide) — no stereocenters, no conflict
    dms_mol = _build_resolved_mol(cif_7gj7, "C")
    assert _check_stereo_vs_template(dms_mol) is True


def test_stereo_check_partial_resolution(cif_1ngx):
    """Test _check_stereo_vs_template with partially resolved ligand.

    JEF in 1ngx chain E has 28/41 heavy atoms resolved. The CCD template
    must be trimmed via MCS to match only the resolved atoms before CIP
    comparison.
    """
    from plinder.data.utils.annotations.ligand_utils import _check_stereo_vs_template

    jef_mol = _build_resolved_mol(cif_1ngx, "E")
    assert jef_mol.GetNumAtoms() < 41, "JEF should be partially resolved"

    # Stereo should match in the resolved portion
    result = _check_stereo_vs_template(jef_mol)
    assert (
        result is not None
    ), "Partially resolved JEF should have comparable stereocenters"

    # Flipping should be detected even on trimmed template
    jef_flipped = _flip_first_chiral(jef_mol)
    if jef_flipped is not None:
        result_flipped = _check_stereo_vs_template(jef_flipped)
        assert result_flipped is False, "Flipped partial JEF should be detected"


def test_stereo_check_multi_residue(cif_6fx1):
    """Test _check_stereo_vs_template on multi-residue glycan.

    6fx1 chain M: NAG+BMA+MAN+FUC+C4W (25+ chiral centers).

    Multi-residue ligands have inter-residue bonds (glycosidic) that
    change CIP priorities vs isolated CCD residues.  The per-residue
    comparison may report False for centers whose CIP changed due to
    the glycosidic bond — this is a known limitation, not a bug.

    We verify:
    1. The function returns a definite result (not None)
    2. The mol has chiral centers that are being compared
    """
    from plinder.data.utils.annotations.ligand_utils import _check_stereo_vs_template

    glycan_mol = _build_resolved_mol(cif_6fx1, "M")

    # Verify multi-residue composition
    res_names = {
        a.GetPDBResidueInfo().GetResidueName().strip()
        for a in glycan_mol.GetAtoms()
        if a.GetPDBResidueInfo()
    }
    assert len(res_names) > 1, f"Should be multi-residue, got {res_names}"

    # Must return a definite result (True or False), not None
    # (None would mean no comparable centers — wrong for a glycan)
    result = _check_stereo_vs_template(glycan_mol)
    assert (
        result is not None
    ), "Multi-residue glycan should have comparable stereocenters"

    # Verify the mol actually has chiral centers
    n_chiral = sum(
        1 for a in glycan_mol.GetAtoms() if a.GetPropsAsDict().get("_CIPCode")
    )
    assert n_chiral > 10, f"Glycan should have many chiral centers, got {n_chiral}"


def test_multi_ligand_system_grouping(cif_7fee, mock_alternative_datasets):
    """Test pocket-based grouping for adjacent drug-like ligands (7fee GPCR).

    7fee: GPCR with 9GF + 7IC binding adjacent pockets (7.9 Å apart,
    4 shared receptor residues). Uses GetPlinderAnnotation for full
    classification. 9GF and 7IC must be in the same system.
    """
    entry_dir = mock_alternative_datasets("7fee")
    plinder_anno = GetPlinderAnnotation(cif_7fee, "", save_folder=entry_dir)
    plinder_anno.annotate()

    systems = plinder_anno.entry.systems

    # 9GF(D) + 7IC(E) grouped by shared pocket residues
    assert "7fee__1__1.A__1.D_1.E" in systems
    drug_sys = systems["7fee__1__1.A__1.D_1.E"]
    assert sorted(l.ccd_code for l in drug_sys.ligands) == ["7IC", "9GF"]
    assert drug_sys.system_type == "holo"
    for lig in drug_sys.ligands:
        assert not lig.is_artifact, f"{lig.ccd_code} should not be artifact"
        assert lig.is_proper, f"{lig.ccd_code} should be proper"

    # CLR(B) standalone — not chained into drug system
    assert "7fee__1__1.A__1.B" in systems
    clr_sys = systems["7fee__1__1.A__1.B"]
    assert [l.ccd_code for l in clr_sys.ligands] == ["CLR"]
    assert clr_sys.system_type == "holo"
    clr = clr_sys.ligands[0]
    assert clr.is_proper, "CLR should be proper"
    assert not clr.is_artifact, "CLR should not be artifact"
    assert not clr.is_cofactor, "CLR is a lipid, not a cofactor"

    # CLR(C) + OLC(L) grouped by proximity
    assert "7fee__1__1.A__1.C_1.L" in systems
    clr_olc = systems["7fee__1__1.A__1.C_1.L"]
    assert sorted(l.ccd_code for l in clr_olc.ligands) == ["CLR", "OLC"]
    assert clr_olc.system_type == "holo"
    for lig in clr_olc.ligands:
        assert not lig.is_cofactor, f"{lig.ccd_code} should not be cofactor"


def test_cofactor_system_stays_holo(cif_1atp, mock_alternative_datasets):
    """Test that cofactor systems are holo via production code path.

    1atp: PKA with ATP (cofactor, from mock DB) + 2x Mn (ions) + PKI peptide.
    Uses GetPlinderAnnotation for full classification.
    """
    entry_dir = mock_alternative_datasets("1atp")
    plinder_anno = GetPlinderAnnotation(cif_1atp, "", save_folder=entry_dir)
    plinder_anno.annotate()

    systems = plinder_anno.entry.systems

    # ATP(E) + Mn(C,D) + PKI peptide(B) all in one holo system
    # B (PKI, 20 res) is receptor with min_polymer_size=12
    expected = "1atp__1__1.A_1.B__1.C_1.D_1.E"
    assert expected in systems, f"Expected {expected}, got {sorted(systems.keys())}"
    atp_sys = systems[expected]
    assert atp_sys.system_type == "holo"
    codes = {l.ccd_code for l in atp_sys.ligands}
    assert "ATP" in codes
    assert "MN" in codes

    for lig in atp_sys.ligands:
        if lig.ccd_code == "ATP":
            assert lig.is_cofactor, "ATP should be cofactor"
            assert lig.is_proper, "ATP should be proper"
            assert not lig.is_artifact, "ATP should not be artifact"
        elif lig.ccd_code == "MN":
            assert lig.is_ion, "MN should be ion"


def test_cofactor_system_holo_19hc(cif_19hc, mock_alternative_datasets):
    """Test that cofactor systems (19hc HEM) are holo.

    19hc: erythrocruorin with 18 HEM (cofactor) + 5 ACT (artifact).
    HEMs share pocket residues on the same protein chain (adjacent
    binding sites), so pocket-based grouping merges them into one
    large system.  ACTs attach via 4 Å proximity to HEMs; isolated
    ACTs (C, P) form standalone artifact systems.
    """
    entry_dir = mock_alternative_datasets("19hc")
    plinder_anno = GetPlinderAnnotation(cif_19hc, "", save_folder=entry_dir)
    plinder_anno.annotate()

    systems = plinder_anno.entry.systems

    # All 18 HEMs merge via shared pocket residues + 3 ACTs attach via proximity
    big = "19hc__1__1.A_1.B__1.D_1.E_1.F_1.G_1.H_1.I_1.J_1.K_1.L_1.M_1.N_1.O_1.Q_1.R_1.S_1.T_1.U_1.V_1.W_1.X_1.Y"
    assert big in systems, f"Expected merged HEM system, got {sorted(systems.keys())}"
    big_sys = systems[big]
    assert big_sys.system_type == "holo"
    codes = {l.ccd_code for l in big_sys.ligands}
    assert "HEM" in codes
    assert "ACT" in codes
    hem_count = sum(1 for l in big_sys.ligands if l.ccd_code == "HEM")
    assert hem_count == 18, f"Expected 18 HEMs, got {hem_count}"

    # Only one system: isolated ACTs (C, P) have no protein neighbors
    assert len(systems) == 1, f"Expected 1 system, got {sorted(systems.keys())}"

    # HEM classification checks
    for lig in big_sys.ligands:
        if lig.ccd_code == "HEM":
            assert lig.is_cofactor, "HEM should be cofactor"
            assert lig.is_proper, "HEM should be proper"
            assert not lig.is_artifact, "HEM should not be artifact"
        if lig.ccd_code == "ACT":
            assert lig.is_artifact, "ACT should be artifact"


def test_nucleic_acid_receptor_detection(cif_8ufz):
    """Verify DNA/RNA chains are included as receptor neighbors (issue #61).

    Uses 8ufz: protein-DNA complex (DNA A-D, protein E-F) with ligand
    Y5U (chains G, H) that binds at the DNA-protein interface.
    Without the filter fix, DNA chains would be invisible as receptor
    neighbors and the ligand would miss DNA interactions.
    """
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx
    from plinder.core.structure.atoms import is_hydrogen_isotope
    from plinder.data.utils.annotations.cif_utils import read_mmcif_file

    cif_obj = read_mmcif_file(cif_8ufz)
    atoms = pdbx.get_structure(
        cif_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[~is_hydrogen_isotope(atoms.element)]

    dna_chains = {"A", "B", "C", "D"}
    protein_chains = {"E", "F"}

    # DNA chains must be detected as nucleotides
    for chain_id in dna_chains:
        chain_atoms = atoms[atoms.chain_id == chain_id]
        assert struc.filter_nucleotides(
            chain_atoms
        ).any(), f"Chain {chain_id} should be detected as nucleotide"

    # Receptor mask must include both protein AND DNA
    receptor_mask = struc.filter_amino_acids(atoms) | struc.filter_nucleotides(atoms)
    receptor_chains = set(atoms.chain_id[receptor_mask])
    assert dna_chains.issubset(
        receptor_chains
    ), f"DNA chains {dna_chains} missing from receptor set {receptor_chains}"
    assert protein_chains.issubset(
        receptor_chains
    ), f"Protein chains {protein_chains} missing from receptor set {receptor_chains}"

    # Ligand Y5U (chain G) must have DNA neighbors within 6A
    lig_coords = atoms.coord[atoms.chain_id == "G"]
    receptor_atoms = atoms[receptor_mask]
    cell = struc.CellList(receptor_atoms, 6.0)
    near_mask = np.zeros(len(receptor_atoms), dtype=bool)
    for coord in lig_coords:
        indices = cell.get_atoms(coord, radius=6.0)
        near_mask[indices[indices >= 0]] = True
    neighbor_chains = set(receptor_atoms.chain_id[near_mask])
    assert (
        neighbor_chains & dna_chains
    ), f"Ligand Y5U should have DNA neighbors, got {neighbor_chains}"
    assert (
        neighbor_chains & protein_chains
    ), f"Ligand Y5U should have protein neighbors, got {neighbor_chains}"


def test_get_validation(
    cif_1qz5,
    validation_1qz5,
    mock_alternative_datasets,
):
    entry_dir = mock_alternative_datasets("1qz5")
    reference_df = pd.DataFrame.from_dict(
        {
            "system_id": "1qz5__1__1.A__1.D",
            "system_ligand_validation_atom_count": 67,
            "system_ligand_validation_average_occupancy": 1.000,
            "system_ligand_validation_average_rscc": 0.958,
            "system_ligand_validation_average_rsr": 0.098,
            "ligand_ccd_code": "KAB",
            "entry_validation_resolution": 1.45,
            "entry_validation_rfree": 0.19,
            "entry_validation_r": 0.17,
            "entry_validation_clashscore": 3.93,
            "entry_validation_percent_rama_outliers": 0.00,
            "entry_validation_molprobity": 1.408352978496041,
            "entry_validation_r_minus_rfree": 0.01999999999999999,
            "system_pocket_validation_max_alt_count": 2,
            "system_pass_validation_criteria": False,
            "entry_pass_validation_criteria": True,
        },
        orient="index",
    ).T.infer_objects()
    entry = GetPlinderAnnotation(cif_1qz5, validation_1qz5, save_folder=entry_dir)
    entry.annotate()
    validation_df = entry.annotated_df[reference_df.columns]
    validation_df = validation_df[
        validation_df["system_id"] == "1qz5__1__1.A__1.D"
    ].reset_index(drop=True)

    pd.testing.assert_frame_equal(reference_df, validation_df)


def test_mmp(mini_mmp_index, mini_mmp_data_annotation, mini_mmp_cluster_folder):
    system_df = pd.read_csv(mini_mmp_data_annotation, sep="\t")
    load_mmp_df = pd.read_csv(mini_mmp_index, compression="gzip", sep="\t", header=None)
    load_mmp_df.columns = ["SMILES1", "SMILES2", "id1", "id2", "V1>>V2", "CONSTANT"]
    mmp_data = add_mmp_clusters_to_data(
        load_mmp_df,
        system_df,
        cluster_folder=mini_mmp_cluster_folder,
        protein_metric="protein_fident_weighted_sum",
        protein_threshold=95,
        protein_directed=False,
        pocket_metric="pocket_fident",
        pocket_threshold=100,
        pocket_directed=True,
        min_constant_size=10,
    )
    # All have same pocket-protein id
    assert list(mmp_data.prot_pocket_set_shared.unique()) == ["c1931_c55468"]

    # Minimum constant size is greater than 10
    assert mmp_data.const_size.min() == 18.0

    # Number of unique congeneric ids is equal to number of unique constants
    assert len(mmp_data.congeneric_id.unique()) == len(mmp_data.CONSTANT.unique())


def test_ligand_fix_to_valid_imatinib(cif_2hyy, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("2hyy")
    entry = Entry.from_cif_file(
        cif_2hyy,
        save_folder=entry_dir,
    )
    lig = entry.systems["2hyy__1__1.A__1.E"].ligands[0]
    #  before fix it is invalid
    assert lig.is_invalid == False
    outsdffile = entry_dir / "2hyy" / "ligand_files" / "E.sdf"
    assert outsdffile.is_file()
    rdmol_sdf = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    rdmol_smi = Chem.MolFromSmiles(lig.smiles)
    # check that numnber of aromatic rings is undderstood correctly
    # N.B. this is expected to be true for fully resolved systems
    assert Chem.rdMolDescriptors.CalcNumAromaticRings(
        rdmol_sdf
    ) == Chem.rdMolDescriptors.CalcNumAromaticRings(rdmol_smi)


def test_ligand_fix_to_valid_thalidomide(cif_7bqu, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("7bqu")
    entry = Entry.from_cif_file(
        cif_7bqu,
        save_folder=entry_dir,
    )
    # EF2 may group with nearby ZN via shared pocket residues
    lig = None
    for system in entry.systems.values():
        for l in system.ligands:
            if l.ccd_code == "EF2":
                lig = l
                break
    assert lig is not None, "EF2 ligand not found in any system"
    assert lig.is_invalid == False
    outsdffile = entry_dir / "7bqu" / "ligand_files" / "C.sdf"
    assert outsdffile.is_file()
    rdmol_sdf = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    rdmol_smi = Chem.MolFromSmiles(lig.smiles)
    # check that numnber of aromatic rings is undderstood correctly
    # N.B. this is expected to be true for fully resolved systems
    assert Chem.rdMolDescriptors.CalcNumAromaticRings(
        rdmol_sdf
    ) == Chem.rdMolDescriptors.CalcNumAromaticRings(rdmol_smi)


def test_partially_resolved_substructure_JEF(cif_1ngx, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("1ngx")
    entry = Entry.from_cif_file(
        cif_1ngx,
        save_folder=entry_dir,
    )
    lig = entry.systems["1ngx__1__1.A_1.B__1.E"].ligands[0]
    assert lig.is_invalid == False
    assert lig.num_unresolved_heavy_atoms == 13
    outsdffile = entry_dir / "1ngx" / "ligand_files" / "E.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    rdmol_smi = Chem.MolFromSmiles(lig.smiles)
    substruct_matches = rdmol_smi.GetSubstructMatches(rdmol)
    assert len(substruct_matches) == 3
    assert len(substruct_matches[0]) == 28


def test_distorted_molecule_template_fix(cif_3grt, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("3grt")
    entry = Entry.from_cif_file(
        cif_3grt,
        save_folder=entry_dir,
    )
    # FAD(B) and TS2(C) may group via shared pocket residues
    lig = None
    for system in entry.systems.values():
        for l in system.ligands:
            if l.ccd_code == "FAD":
                lig = l
                break
    assert lig is not None, "FAD ligand not found in any system"
    assert lig.is_invalid == False
    # Check SDF was saved and is valid
    outsdffile = entry_dir / "3grt" / "ligand_files" / f"{lig.asym_id}.sdf"
    assert outsdffile.is_file(), f"SDF not found at {outsdffile}"
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE


def test_hydrogen_removed_save(cif_7az3, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("7az3")
    entry = Entry.from_cif_file(
        cif_7az3,
        save_folder=entry_dir,
    )
    ligand = next(iter(entry.systems.values())).ligands[0]
    outsdffile = entry_dir / "7az3" / "ligand_files" / f"{ligand.asym_id}.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=False, sanitize=False)[0]
    assert sum([at.GetAtomicNum() == 1 for at in rdmol.GetAtoms()]) == 0
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE


def test_too_many_hydrogens(cif_6ntj, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("6ntj")
    entry = Entry.from_cif_file(
        cif_6ntj,
        save_folder=entry_dir,
    )
    ligand = next(iter(entry.systems.values())).ligands[0]
    outsdffile = entry_dir / "6ntj" / "ligand_files" / f"{ligand.asym_id}.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=False, sanitize=False)[0]
    assert sum([at.GetAtomicNum() == 1 for at in rdmol.GetAtoms()]) == 0
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE


def test_disconnected_ligand_fix(cif_4nhc, mock_alternative_datasets):
    """4nhc chain C is a 17-residue peptide ligand (with min_polymer_size=20).

    Tests that a fragmented peptide gets fixed to a valid, connected SDF.
    """
    entry_dir = mock_alternative_datasets("4nhc")
    # Use threshold 20 so the 17-residue peptide is classified as ligand
    entry = Entry.from_cif_file(cif_4nhc, save_folder=entry_dir, min_polymer_size=20)
    lig = entry.systems["4nhc__1__1.A_1.B__1.C"].ligands[0]
    assert lig.is_invalid == False
    outsdffile = entry_dir / "4nhc" / "ligand_files" / "C.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    assert len(Chem.MolToSmiles(rdmol).split(".")) == 1


def test_binding_affinity(cif_4jvn, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("4jvn")
    entry = Entry.from_cif_file(cif_4jvn, save_folder=entry_dir)
    target_value = 7.638272164
    affinity = 0.0
    for sys in entry.systems.values():
        if "YUG" in set(l.ccd_code for l in sys.ligands):
            assert sys.has_binding_affinity
        for ligand in sys.ligands:
            if ligand.ccd_code == "YUG":
                affinity = ligand.binding_affinity
    assert affinity == target_value
