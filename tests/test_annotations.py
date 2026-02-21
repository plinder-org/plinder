# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
from rdkit import Chem

from plinder.data.get_system_annotations import GetPlinderAnnotation
from plinder.data.utils.annotations.aggregate_annotations import Entry
from plinder.data.utils.annotations.interaction_utils import get_covalent_connections
from plinder.data.utils.annotations.interface_gap import annotate_interface_gaps
from plinder.data.utils.annotations.ligand_utils import (
    get_smiles_from_cif,
    sort_ccd_codes,
)
from plinder.data.utils.annotations.mmpdb_utils import add_mmp_clusters_to_data
from plinder.data.utils.annotations.protein_utils import read_mmcif_container


def test_ccd_name_sorter():
    assert sort_ccd_codes({"G", "G25", "CPG", "5GP"}) == ["CPG", "G25", "G", "5GP"]

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
    assert len(df) == 1
    assert df["ligand_is_covalent"].sum() == 0
    assert set(df.ligand_ccd_code.to_list()) == {"LYS-ALA-ASP-THR-THR-THR-PRO"}
    # Note: chain 'B' is ligand = should not be in protein neigh list
    assert df["ligand_protein_chains_auth_id"].drop_duplicates().to_list() == [["A"]]

def test_synthetic_noncov_peptide_detection(cif_6u6k, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("6u6k")
    plinder_anno = GetPlinderAnnotation(cif_6u6k, "", save_folder=entry_dir)
    plinder_anno.annotate()
    df = plinder_anno.annotated_df
    assert len(df) == 1
    assert df["ligand_is_covalent"].sum() == 0
    assert set(df.ligand_ccd_code.to_list()) == {
        "ACE-TRP-TRP-ILE-ILE-PRO-ALY-VAL-LYS-ALY-GLY-CYS-NH2"
    }

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
    outsdffile = entry_dir / "6lu7__1__1.A_2.A__1.B/ligand_files/1.B.sdf"
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

    # 10 PLIPs detected:
    # consistent with SWISSMODEL as of 2024-04-18
    # https://swissmodel.expasy.org/templates/4ci1

    # expected_interactions = {
    #     404: ['type:hydrogen_bonds__donortype:Nam__acceptortype:O2__protisdon:False__sidechain:False',
    #           'type:hydrogen_bonds__donortype:Nar__acceptortype:O2__protisdon:True__sidechain:True'],
    #     406: ['type:hydrogen_bonds__donortype:Nam__acceptortype:O2__protisdon:True__sidechain:False', 'type:hydrophobic_contacts'],
    #     377: ['type:hydrogen_bonds__donortype:Nam__acceptortype:O2__protisdon:True__sidechain:True', 'type:hydrophobic_contacts'],
    #     412: ['type:hydrophobic_contacts', 'type:hydrophobic_contacts'],
    #     426: ['type:hydrophobic_contacts'],
    #     428: ['type:hydrophobic_contacts']
    # }
    expected_interactions = {
        404: [
            "type:hydrogen_bonds__protisdon:False__sidechain:False",
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
        ],
        406: [
            "type:hydrogen_bonds__protisdon:True__sidechain:False",
            "type:hydrophobic_contacts",
        ],
        377: [
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
            "type:hydrophobic_contacts",
        ],
        412: ["type:hydrophobic_contacts", "type:hydrophobic_contacts"],
        426: ["type:hydrophobic_contacts"],
        428: ["type:hydrophobic_contacts"],
    }
    # get if the count is right
    assert len(ligand.interactions["1.B"]) == len(expected_interactions)
    # exact report matching
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

    # 12 PLIPs detected:
    # consistent with SWISSMODEL as of 2024-04-18
    # https://swissmodel.expasy.org/templates/2p1q.1
    # expected_interactions_2B =  {
    #     438: ['type:hydrogen_bonds__donortype:O.co2__acceptortype:O3__protisdon:False__sidechain:True',
    #           'type:hydrogen_bonds__donortype:O3__acceptortype:O.co2__protisdon:True__sidechain:True'],
    #     79: ['type:hydrophobic_contacts', 'type:hydrophobic_contacts'],
    #     464: ['type:hydrophobic_contacts'],
    #     403: ['type:water_bridges__donortype:Ng+__acceptortype:O.co2__protisdon:True',
    #           'type:water_bridges__donortype:Ng+__acceptortype:O.co2__protisdon:True',
    #           'type:salt_bridges__lig_group:carboxylate__protispos:True'],
    #     78: ['type:salt_bridges__lig_group:carboxylate__protispos:True']
    # }
    expected_interactions_2B = {
        438: [
            "type:hydrogen_bonds__protisdon:True__sidechain:True",
        ],
        439: ["type:hydrogen_bonds__protisdon:False__sidechain:False"],
        79: ["type:hydrophobic_contacts", "type:hydrophobic_contacts"],
        464: ["type:hydrophobic_contacts"],
        403: [
            "type:water_bridges__protisdon:True",
            "type:water_bridges__protisdon:True",
            "type:salt_bridges__protispos:True",
        ],
        78: ["type:salt_bridges__protispos:True"],
    }
    expected_interactions_2C = {
        7: ["type:hydrophobic_contacts", "type:water_bridges__protisdon:False"],
        5: ["type:pi_stacks__stack_type:T"],
    }

    expected_waters = {"2.G": {66, 4, 2}}

    # get if the count is right
    assert len(ligand.interactions["2.B"]) == len(expected_interactions_2B)
    assert len(ligand.interactions["2.C"]) == len(expected_interactions_2C)

    # exact report matching
    assert ligand.interactions["2.B"] == expected_interactions_2B
    assert ligand.interactions["2.C"] == expected_interactions_2C

    # waters
    assert {k: set(v) for k, v in ligand.waters.items()} == expected_waters

def test_water_saving(cif_2p1q, mock_alternative_datasets):
    from ost import io

    entry_dir = mock_alternative_datasets("2p1q")
    system_tag = "2p1q__2__2.B_2.C__2.E"
    Entry.from_cif_file(cif_2p1q, save_folder=entry_dir)
    for filename in [
        "sequences.fasta",
        "receptor.pdb",
        "system.cif",
        "receptor.cif",
        "chain_mapping.json",
        "water_mapping.json",
    ]:
        assert (entry_dir / system_tag / filename).exists()
    assert (entry_dir / system_tag / "ligand_files" / "2.E.sdf").exists()
    ent = io.LoadPDB(str(entry_dir / system_tag / "receptor.pdb"))
    assert len(ent.FindChain("_").residues) == 3

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

def test_system_saving(cif_2y4i, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("2y4i")
    system_tag = "2y4i__1__1.B__1.E_1.F"
    entry = Entry.from_cif_file(cif_2y4i, save_folder=entry_dir)
    for chain in entry.chains.values():
        assert len(chain.mappings["ECOD"])
    for filename in [
        "sequences.fasta",
        "receptor.pdb",
        "system.cif",
        "receptor.cif",
        "chain_mapping.json",
    ]:
        assert (entry_dir / system_tag / filename).exists()
    for chain in ["1.E", "1.F"]:
        assert (entry_dir / system_tag / "ligand_files" / f"{chain}.sdf").exists()

def test_smiles_from_nextgen(test_dir, smiles_sample_csv):
    from ost import io

    results = []
    pdbids = ["1ppc", "6fx1", "6m92", "2dty", "7gj7", "2e84", "6u6k"]
    for pdbid in pdbids:
        cif_file = test_dir / f"xx/pdb_0000{pdbid}/pdb_0000{pdbid}_xyz-enrich.cif.gz"
        data = read_mmcif_container(cif_file)
        ent = io.LoadMMCIF(str(cif_file))
        pdbid = cif_file.stem.split("_")[1].split("0000")[-1]
        result = get_smiles_from_cif(data, ent)
        result = [(pdbid, k, v) for k, v in result.items()]
        results.extend(result)
    result_df = pd.DataFrame(results, columns=["pdbid", "chain", "smiles"])
    result_df = result_df.sort_values(by=["pdbid", "chain"]).reset_index(drop=True)
    target_df = pd.read_csv(smiles_sample_csv)
    target_df = target_df.sort_values(by=["pdbid", "chain"]).reset_index(drop=True)
    # Canonicalize SMILES to absorb differences across OST versions
    for df in [result_df, target_df]:
        df["smiles"] = df["smiles"].apply(
            lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
            if Chem.MolFromSmiles(s) is not None
            else s
        )
    pd.testing.assert_frame_equal(result_df, target_df)

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
    outsdffile = entry_dir / "2hyy__1__1.A__1.E/ligand_files/1.E.sdf"
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
    lig = entry.systems["7bqu__1__1.A_1.B__1.C"].ligands[0]
    assert lig.is_invalid == False
    outsdffile = entry_dir / "7bqu__1__1.A_1.B__1.C/ligand_files/1.C.sdf"
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
    outsdffile = entry_dir / "1ngx__1__1.A_1.B__1.E/ligand_files/1.E.sdf"
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
    lig = entry.systems["3grt__1__1.A_2.A__1.B"].ligands[0]
    assert lig.is_invalid == False
    outsdffile = entry_dir / "3grt__1__1.A_2.A__1.B/ligand_files/1.B.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE

def test_hydrogen_removed_save(cif_7az3, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("7az3")
    entry = Entry.from_cif_file(
        cif_7az3,
        save_folder=entry_dir,
    )
    pli_entry = list(entry.systems.keys())[0]
    outsdffile = entry_dir / pli_entry / f"ligand_files/{pli_entry[-3:]}.sdf"
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
    pli_entry = list(entry.systems.keys())[0]
    outsdffile = entry_dir / pli_entry / f"ligand_files/{pli_entry[-3:]}.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=False, sanitize=False)[0]
    assert sum([at.GetAtomicNum() == 1 for at in rdmol.GetAtoms()]) == 0
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE

def test_disconnected_ligand_fix(cif_4nhc, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("4nhc")
    entry = Entry.from_cif_file(cif_4nhc, save_folder=entry_dir, skip_posebusters=True)
    lig = entry.systems["4nhc__1__1.A_1.B__1.C"].ligands[0]
    assert lig.is_invalid == False
    outsdffile = entry_dir / "4nhc__1__1.A_1.B__1.C/ligand_files/1.C.sdf"
    assert outsdffile.is_file()
    rdmol = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    assert Chem.SanitizeMol(rdmol) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    assert len(Chem.MolToSmiles(rdmol).split(".")) == 1

def test_binding_affinity(cif_4jvn, mock_alternative_datasets):
    entry_dir = mock_alternative_datasets("4jvn")
    entry = Entry.from_cif_file(cif_4jvn, save_folder=entry_dir, skip_posebusters=True)
    target_value = 7.638272164
    affinity = 0.0
    for sys in entry.systems.values():
        if "YUG" in set(l.ccd_code for l in sys.ligands):
            assert sys.has_binding_affinity
        for ligand in sys.ligands:
            if ligand.ccd_code == "YUG":
                affinity = ligand.binding_affinity
    assert affinity == target_value
