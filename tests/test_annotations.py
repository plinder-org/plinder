# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from types import SimpleNamespace

import biotite.structure as struc
import numpy as np
import pandas as pd
import pytest
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import read_mmcif_container
from plinder.data.annotations.get_ligand_validation import EntryValidation
from plinder.data.annotations.interaction_utils import (
    extract_ligand_links_to_neighbouring_chains,
    get_covalent_connections,
)
from plinder.data.annotations.ligand_utils import (
    BiounitSpatialIndex,
    classify_ligand_polymer_classes,
    get_water_chain_ids,
    is_excluded_mol,
    is_known_artifact_ligand,
)
from plinder.data.annotations.mmpdb_utils import add_mmp_clusters_to_data
from plinder.data.annotations.protein_utils import get_receptor_type
from plinder.data.annotations.save_utils import (
    SystemReconstructionOptions,
    SystemReconstructionOutputs,
    save_ligands,
    save_reconstructed_system,
)
from plinder.data.get_system_annotations import GetPlinderAnnotation
from rdkit import Chem


@pytest.mark.parametrize(
    ("smiles", "expected_true"),
    [
        ("OC1OC(O)C(O)C(O)C1O", "is_monosaccharide"),
        (
            "OC1OC(O)C(O)C(O)C1OC2OC(O)C(O)C(O)C2O",
            "is_oligosaccharide",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O",
            "is_mononucleotide",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O[C@@H]2[C@@H](O)"
            "[C@@H](n3cnc4c(N)ncnc43)O[C@@H]2CO)[C@@H](O)[C@H]1O",
            "is_oligonucleotide",
        ),
        ("NCC(=O)O", "is_monopeptide"),
        ("NCC(=O)NCC(=O)O", "is_oligopeptide"),
    ],
)
def test_classify_ligand_polymer_classes_from_structure(
    smiles: str, expected_true: str
) -> None:
    classes = classify_ligand_polymer_classes(smiles)

    assert classes[expected_true]
    family = expected_true.removeprefix("is_mono").removeprefix("is_oligo")
    opposite = (
        f"is_oligo{family}"
        if expected_true.startswith("is_mono")
        else f"is_mono{family}"
    )
    assert not classes[opposite]


def test_ligand_polymer_units_are_not_summed_across_disconnected_fragments() -> None:
    classes = classify_ligand_polymer_classes("NCC(=O)O.NCC(=O)O")

    assert classes["is_monopeptide"]
    assert not classes["is_oligopeptide"]


def test_get_water_chain_ids_requires_all_chain_atoms_to_be_solvent():
    atoms = struc.AtomArray(6)
    atoms.chain_id = np.array(["A", "A", "W", "W", "M", "M"])
    atoms.res_name = np.array(["ALA", "ALA", "HOH", "HOH", "HOH", "ALA"])

    assert get_water_chain_ids(atoms) == {"W"}


def test_save_ligands_falls_back_to_biotite_without_rdkit_sanitization(
    tmp_path, monkeypatch
):
    atoms = struc.AtomArray(4)
    atoms.chain_id = np.array(["A", "A", "B", "B"])
    atoms.res_id = np.array([1, 1, 2, 2])
    atoms.res_name = np.array(["LIG", "LIG", "OTH", "OTH"])
    atoms.atom_name = np.array(["C1", "N1", "O1", "C1"])
    atoms.element = np.array(["C", "N", "O", "C"])
    atoms.coord = np.array(
        [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0], [5.0, 0.0, 0.0], [6.2, 0.0, 0.0]]
    )
    atoms.bonds = struc.BondList(len(atoms))
    atoms.bonds.add_bond(0, 1, struc.BondType.AROMATIC_SINGLE)
    atoms.bonds.add_bond(2, 3, struc.BondType.DOUBLE)
    monkeypatch.setattr(
        "plinder.data.annotations.cif_utils.atoms_to_rdkit_mol",
        lambda _atoms: (_ for _ in ()).throw(ValueError("cannot sanitize")),
    )

    save_ligands(atoms, ["A"], tmp_path)

    saved = struc.io.mol.get_structure(struc.io.mol.SDFile.read(tmp_path / "A.sdf"))
    assert saved.array_length() == 2
    assert np.allclose(saved.coord, atoms.coord[:2])
    assert not (tmp_path / "B.sdf").exists()


def test_entry_validation_accepts_missing_optional_wwpdb_metrics():
    entry = {
        "PDB-resolution": "2.5",
        "PDB-Rfree": "0.25",
        "PDB-R": "0.20",
        "clashscore": None,
    }
    validation = SimpleNamespace(
        getValidationXML=lambda: SimpleNamespace(getEntry=lambda: entry),
        getReflectionsResolution=lambda: 2.4,
        getMeanIOverSigIObs=lambda: ".",
        countAtoms=lambda: 100,
        calcMolProbityOverallScore=lambda: 1.5,
        calcMeanStructureBFactor=lambda: 20.0,
        calcMedianStructureBFactor=lambda: 18.0,
        getResolution=lambda: 2.5,
    )

    result = EntryValidation.from_entry(validation)

    assert result.clashscore is None
    assert result.meanI_over_sigI_obs is None


def test_crystal_contact_fraction_is_undefined_without_heavy_atoms():
    from plinder.data.annotations.aggregate_annotations import System
    from plinder.data.annotations.ligand_utils import Ligand

    ligand = Ligand(num_heavy_atoms=0)
    system = System(
        pdb_id="1abc",
        biounit_id="1",
        receptor_type="protein",
        ligands=[ligand],
    )

    assert ligand.fraction_atoms_with_crystal_contacts is None
    assert system.fraction_atoms_with_crystal_contacts is None


def test_known_artifact_preflight_is_conservative(monkeypatch):
    monkeypatch.setattr(
        "plinder.data.annotations.ligand_utils._get_ccd_smiles",
        lambda code: {"OHX": "[OH-]", "LIG": "CCNCC"}.get(code),
    )

    assert is_known_artifact_ligand(["GOL"], {"GOL"})
    # UNX is a dummy placeholder (DUM was retired -> UNX; obsolete, never ingested)
    assert is_known_artifact_ligand(["UNX"], set())
    assert is_known_artifact_ligand(["OHX"], set())
    assert not is_known_artifact_ligand(["LIG"], set())
    assert not is_known_artifact_ligand(["OHX", "OHX"], set())


def test_biounit_spatial_index_expands_hits_to_complete_residues():
    atoms = struc.AtomArray(6)
    atoms.chain_id = np.array(["1.A", "1.A", "1.B", "1.B", "1.A", "1.A"])
    atoms.res_id = np.array([1, 1, 1, 1, 2, 2])
    atoms.coord = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [11.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [21.0, 0.0, 0.0],
        ]
    )
    atoms.bonds = struc.BondList(len(atoms))
    atoms.bonds.add_bond(0, 1, struc.BondType.SINGLE)
    atoms.bonds.add_bond(1, 2, struc.BondType.DOUBLE)
    atoms.bonds.add_bond(4, 5, struc.BondType.TRIPLE)
    spatial_index = BiounitSpatialIndex.from_atoms(atoms, max_radius=3.0)

    assert np.array_equal(
        spatial_index.atom_indices_for_chain("1.A"),
        np.array([0, 1, 4, 5]),
    )
    assert np.array_equal(
        spatial_index.complete_residue_indices_near(atoms.coord[[0]], radius=0.5),
        np.array([0, 1]),
    )
    selected_indices = np.array([0, 1, 4, 5])
    expected = atoms[selected_indices]
    actual = spatial_index.take_atoms(atoms, selected_indices, include_bonds=True)
    assert np.array_equal(actual.coord, expected.coord)
    for category in atoms.get_annotation_categories():
        assert np.array_equal(
            actual.get_annotation(category), expected.get_annotation(category)
        )
    assert actual.bonds is not None
    assert expected.bonds is not None
    assert np.array_equal(actual.bonds.as_array(), expected.bonds.as_array())
    assert (
        spatial_index.take_atoms(atoms, selected_indices, include_bonds=False).bonds
        is None
    )


def test_deferred_ions_are_retained_only_when_they_can_join_a_primary_system():
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain

    def chain(asym_id: str, chain_type: str) -> Chain:
        return Chain(
            asym_id=asym_id,
            auth_id=asym_id,
            entity_id=asym_id,
            chain_type_str=chain_type,
            residues={},
            length=100,
            num_unresolved_residues=100,
        )

    entry = Entry(
        pdb_id="1abc",
        chains={
            "A": chain("A", "polypeptide(L)"),
            "B": chain("B", "non-polymer"),
            "C": chain("C", "non-polymer"),
            "D": chain("D", "non-polymer"),
        },
        ligand_like_chains={
            "B": "non-polymer",
            "C": "non-polymer",
            "D": "non-polymer",
        },
    )
    primary = Ligand(
        pdb_id="1abc",
        biounit_id="1",
        asym_id="B",
        instance=1,
        ccd_code="LIG",
        plip_type="SMALLMOLECULE",
        bird_id="",
        centroid=[0.0, 0.0, 1.0],
        smiles="CCNCC",
        residue_numbers=[1],
        neighboring_residues={"1.A": [1, 2, 3]},
    )
    atoms = struc.AtomArray(6)
    atoms.chain_id = np.array(["1.A", "1.A", "1.A", "1.B", "1.C", "1.D"])
    atoms.res_id = np.array([1, 2, 3, 1, 1, 1])
    atoms.coord = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 2.0],
            [50.0, 0.0, 0.0],
        ]
    )
    spatial_index = BiounitSpatialIndex.from_atoms(atoms, max_radius=10.0)

    assert entry._connected_deferred_ligand_chains(
        atoms,
        spatial_index,
        [primary],
        {"1.C", "1.D"},
        min_shared_pocket_members=3,
        neighboring_residue_threshold=6.0,
        interaction_search_threshold=10.0,
    ) == {"1.C"}


@pytest.mark.parametrize(
    "chain_types, expected",
    [
        (["polypeptide(L)"], "protein"),
        (["polydeoxyribonucleotide"], "dna"),
        (["polyribonucleotide"], "rna"),
        (["polypeptide(L)", "polyribonucleotide"], "protein+rna"),
        (
            ["polydeoxyribonucleotide/polyribonucleotide hybrid"],
            "dna+rna",
        ),
    ],
)
def test_receptor_type_classification(chain_types, expected):
    assert get_receptor_type(chain_types) == expected


def test_chain_from_cif_data_nucleotides(cif_8ufz):
    """Test Chain.from_cif_data assigns correct one-letter codes and chem_types for DNA.

    8ufz chain A is a 16-nt DNA strand (DA, DT, DC, DG residues).
    """
    import biotite.structure.io.pdbx as pdbx
    from biotite.structure import filter_heavy
    from plinder.data.annotations.cif_utils import read_mmcif_file
    from plinder.data.annotations.protein_utils import Chain, get_seqres_from_cif

    cif_obj = read_mmcif_file(cif_8ufz)
    block = list(cif_obj.values())[0]
    atoms = pdbx.get_structure(
        cif_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[filter_heavy(atoms)]
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


# link format: "auth_seq:comp:label_asym:label_seq:atom" (asym id at index 2)
@pytest.mark.parametrize(
    "covalent, ligand_asym_id, neighboring, expected",
    [
        # 2-letter asym (common once an entry has >26 chains) must be found,
        # emitted as receptor__ligand — silently dropped by character iteration.
        (
            {"covale": [("1:LIG:AA:.:C1", "50:CYS:B:50:SG")]},
            "AA",
            {"B"},
            {"50:CYS:B:50:SG__1:LIG:AA:.:C1"},
        ),
        # single-character asym still works
        (
            {"covale": [("1:LIG:C:.:C1", "50:CYS:A:50:SG")]},
            "C",
            {"A"},
            {"50:CYS:A:50:SG__1:LIG:C:.:C1"},
        ),
        # covale not involving the ligand -> ignored
        (
            {"covale": [("1:LIG:AA:.:C1", "50:CYS:B:50:SG")]},
            "ZZ",
            {"B"},
            set(),
        ),
        # ligand-ligand covale (no receptor neighbour) -> ignored
        (
            {"covale": [("1:L1:AA:.:C", "2:L2:AB:.:N")]},
            "AA",
            {"B"},
            set(),
        ),
        # ligand "AB" is not in this A<->C bond; a shared character ("A") must
        # not fabricate a spurious linkage.
        (
            {"covale": [("1:XXX:A:1:C", "2:YYY:C:2:N")]},
            "AB",
            {"C"},
            set(),
        ),
    ],
    ids=[
        "multichar-found",
        "singlechar-found",
        "not-ligand",
        "ligand-ligand",
        "shared-char-false-positive",
    ],
)
def test_extract_ligand_links_matches_whole_asym_ids(
    covalent, ligand_asym_id, neighboring, expected
):
    """Covalent-linkage matching compares whole asym ids, so a shared character
    neither drops a real multi-character linkage nor fabricates a spurious one."""
    assert (
        extract_ligand_links_to_neighbouring_chains(
            covalent, ligand_asym_id, neighboring, link_type="covale"
        )
        == expected
    )


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
    biounit_chain_df = pd.read_parquet(
        entry_dir / "6i41" / "entry_biounit_chains.parquet"
    )
    assert set(biounit_chain_df["chain_role"]) == {
        "ligand",
        "receptor",
        "water",
    }
    assert not biounit_chain_df.duplicated(
        ["entry_pdb_id", "biounit_id", "chain_instance"]
    ).any()
    source_df = pd.read_parquet(entry_dir / "6i41" / "entry_source.parquet")
    assert source_df["entry_pdb_id"].tolist() == ["6i41"]
    assert (
        source_df[["source_mmcif_major_revision", "source_mmcif_minor_revision"]]
        .notna()
        .all(axis=None)
    )
    assert len(df) == 1
    assert df["ligand_is_covalent"].sum() == 0
    assert set(df.ligand_ccd_code.to_list()) == {"LYS-ALA-ASP-THR-THR-THR-PRO"}
    # Note: chain 'B' is ligand = should not be in protein neigh list
    assert df["ligand_protein_chains_auth_id"].drop_duplicates().to_list() == [["A"]]


def test_entry_ignores_external_mappings_for_absent_chains(
    cif_6i41, mock_alternative_datasets, monkeypatch
):
    entry_dir = mock_alternative_datasets("6i41")
    monkeypatch.setattr(
        "plinder.data.annotations.aggregate_annotations.get_chain_external_mappings",
        lambda _block: {"P": []},
    )

    entry = Entry.from_cif_file(cif_6i41, save_folder=entry_dir)

    assert "P" not in entry.chains


@pytest.mark.parametrize(
    "min_polymer_size,expect_ligand",
    [(12, False), (20, True)],
    ids=["threshold_12_peptide_is_receptor", "threshold_20_peptide_is_ligand"],
)
def test_peptide_ligand_threshold(
    cif_6u6k,
    mock_alternative_datasets,
    min_polymer_size,
    expect_ligand,
    monkeypatch,
):
    """6u6k: 13-residue synthetic peptide (chain B).

    With min_polymer_size=12, the peptide is receptor (13 >= 12) → no systems.
    With min_polymer_size=20, the peptide is ligand (13 < 20) → system created.
    """
    entry_dir = mock_alternative_datasets("6u6k")
    stale_ligand_dir = entry_dir / "6u6k" / "ligand_files"
    stale_ligand_dir.mkdir(parents=True)
    (stale_ligand_dir / "stale.sdf").touch()
    if not expect_ligand:
        monkeypatch.setattr(
            "plinder.data.annotations.aggregate_annotations.pdbx.get_structure",
            lambda *_args, **_kwargs: pytest.fail(
                "entries without ligand-like chains must skip bonded loading"
            ),
        )
        monkeypatch.setattr(
            "plinder.data.annotations.aggregate_annotations.pdbx.list_assemblies",
            lambda *_args, **_kwargs: pytest.fail(
                "entries without ligand-like chains must skip assembly generation"
            ),
        )
    entry = Entry.from_cif_file(
        cif_6u6k, save_folder=entry_dir, min_polymer_size=min_polymer_size
    )
    if expect_ligand:
        assert len(entry.systems) == 1
        assert not (stale_ligand_dir / "stale.sdf").exists()
        lig = entry.systems[list(entry.systems.keys())[0]].ligands[0]
        assert lig.ccd_code == "ACE-TRP-TRP-ILE-ILE-PRO-ALY-VAL-LYS-ALY-GLY-CYS-NH2"
    else:
        assert (
            len(entry.systems) == 0
        ), f"13-residue peptide should be receptor with min_polymer_size={min_polymer_size}"
        assert not stale_ligand_dir.exists()


def test_annotation_without_systems_skips_validation_and_normalized_tables(
    monkeypatch, tmp_path
):
    class EmptyEntry:
        pdb_id = "1abc"
        systems: dict[str, object] = {}

        def set_validation(self, *_args, **_kwargs):
            pytest.fail("empty entries must skip validation")

    monkeypatch.setattr(
        Entry,
        "from_cif_file",
        lambda *_args, **_kwargs: EmptyEntry(),
    )
    annotation = GetPlinderAnnotation(
        Path("1abc.cif"),
        Path("1abc_validation.xml.gz"),
        save_folder=tmp_path,
    )

    assert annotation.annotate() is None
    assert not (tmp_path / "1abc").exists()


def test_entry_drops_systems_without_a_proper_ligand() -> None:
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain

    receptor = Chain(
        asym_id="A",
        auth_id="A",
        entity_id="1",
        chain_type_str="polypeptide(L)",
        residues={},
        length=100,
        num_unresolved_residues=100,
    )

    def ligand(asym_id: str, *, is_ion: bool = False, is_artifact: bool = False):
        return Ligand(
            pdb_id="1abc",
            biounit_id="1",
            asym_id=asym_id,
            instance=1,
            ccd_code="NA" if is_ion else "GOL" if is_artifact else "LIG",
            plip_type="SMALLMOLECULE",
            bird_id="",
            centroid=[0.0, 0.0, 0.0],
            smiles="[Na+]" if is_ion else "CCO",
            residue_numbers=[1],
            neighboring_residues={"1.A": [1, 2, 3]},
            is_ion=is_ion,
            is_artifact=is_artifact,
        )

    for nonproper in (
        ligand("B", is_ion=True),
        ligand("B", is_artifact=True),
    ):
        entry = Entry(pdb_id="1abc", chains={"A": receptor})
        entry.set_systems({nonproper.id: nonproper})
        assert entry.systems == {}

    proper = ligand("B")
    ion = ligand("C", is_ion=True)
    entry = Entry(pdb_id="1abc", chains={"A": receptor})
    entry.set_systems({item.id: item for item in (proper, ion)})

    assert len(entry.systems) == 1
    retained = next(iter(entry.systems.values()))
    assert retained.system_type == "holo"
    assert {item.ccd_code for item in retained.ligands} == {"LIG", "NA"}


def test_entry_never_groups_ligands_across_biological_assemblies() -> None:
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain

    receptor = Chain(
        asym_id="A",
        auth_id="A",
        entity_id="1",
        chain_type_str="polypeptide(L)",
        residues={},
        length=100,
        num_unresolved_residues=100,
    )

    def ligand(biounit_id: str, asym_id: str) -> Ligand:
        return Ligand(
            pdb_id="1abc",
            biounit_id=biounit_id,
            asym_id=asym_id,
            instance=1,
            ccd_code="LIG",
            plip_type="SMALLMOLECULE",
            bird_id="",
            centroid=[0.0, 0.0, 0.0],
            smiles="CCO",
            residue_numbers=[1],
            neighboring_residues={"1.A": [1, 2, 3]},
        )

    ligands = [
        ligand(biounit_id, asym_id)
        for biounit_id in ("1", "2")
        for asym_id in ("B", "C")
    ]
    entry = Entry(pdb_id="1abc", chains={"A": receptor})
    entry.set_systems({item.id: item for item in ligands})

    assert len(entry.systems) == 2
    assert {system.biounit_id for system in entry.systems.values()} == {"1", "2"}
    for system in entry.systems.values():
        assert {item.biounit_id for item in system.ligands} == {system.biounit_id}
        assert len(system.ligand_chains) == len(set(system.ligand_chains)) == 2


def test_entry_validation_skips_chains_outside_retained_systems(
    monkeypatch, tmp_path
) -> None:
    from types import SimpleNamespace

    from plinder.data.annotations.aggregate_annotations import System
    from plinder.data.annotations.ligand_utils import Ligand
    from plinder.data.annotations.protein_utils import Chain

    def chain(asym_id: str, chain_type: str) -> Chain:
        return Chain(
            asym_id=asym_id,
            auth_id=asym_id,
            entity_id=asym_id,
            chain_type_str=chain_type,
            residues={},
            length=100,
            num_unresolved_residues=100,
        )

    chains = {
        "A": chain("A", "polypeptide(L)"),
        "B": chain("B", "non-polymer"),
        "R": chain("R", "polypeptide(L)"),
        "U": chain("U", "non-polymer"),
    }
    ligand = Ligand(
        pdb_id="1abc",
        biounit_id="1",
        asym_id="B",
        instance=1,
        ccd_code="LIG",
        plip_type="SMALLMOLECULE",
        bird_id="",
        centroid=[0.0, 0.0, 0.0],
        smiles="CCNCC",
        residue_numbers=[1],
        neighboring_residues={"1.A": [1, 2, 3]},
    )
    system = System(
        pdb_id="1abc",
        biounit_id="1",
        ligands=[ligand],
        receptor_type="protein",
    )
    entry = Entry(
        pdb_id="1abc",
        determination_method="X-RAY DIFFRACTION",
        chains=chains,
        ligand_like_chains={"B": "non-polymer", "U": "non-polymer"},
        systems={system.id: system},
    )
    validation_path = tmp_path / "validation.xml.gz"
    validation_path.touch()
    validated: list[str] = []

    monkeypatch.setattr(
        "plinder.data.annotations.aggregate_annotations.ValidationFactory",
        lambda *_args, **_kwargs: SimpleNamespace(getValidation=lambda: object()),
    )
    monkeypatch.setattr(
        "plinder.data.annotations.aggregate_annotations.EntryValidation.from_entry",
        lambda _doc: SimpleNamespace(r=0.2),
    )
    monkeypatch.setattr(
        Chain,
        "set_validation",
        lambda self, _doc, _thresholds: validated.append(self.asym_id),
    )
    monkeypatch.setattr(System, "set_validation", lambda *_args, **_kwargs: None)

    entry.set_validation(validation_path, Path("source.cif"))

    assert set(validated) == {"A", "B"}


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
    # The SDF holds exactly the ligand's heavy atoms — not truncated, and NOT
    # leaking the covalently-bonded receptor Cys145 atom (the cross-chain
    # A->B covale bond is dropped when selecting the ligand chain).
    assert rdmol.GetNumAtoms() == lig.num_heavy_atoms == 49


# Non-standard/modified residues that make up the enlicitide (MK-0616)
# macrocycle in 10sb. ALA/PRO/THR are standard and intentionally excluded.
_ENLICITIDE_MODIFIED_RESIDUES = {"SIN", "2RA", "A1CHA", "FTR", "0A1", "3WX"}


def test_10sb_modified_residues_preserved(cif_10sb, mock_alternative_datasets):
    """Enlicitide's non-standard residues must survive as their own CCD codes.

    The macrocyclic peptide (chain C of 10sb) is built from modified /
    non-standard residues (SIN, 2RA, A1CHA, FTR, 0A1, 3WX). The pipeline
    must NOT standardize them to parent amino acids — that information is
    part of the ligand's chemical identity.
    """
    from plinder.data.annotations.aggregate_annotations import Entry

    entry_dir = mock_alternative_datasets("10sb")
    entry = Entry.from_cif_file(cif_10sb, save_folder=entry_dir)

    codes: set[str] = set()
    for system in entry.systems.values():
        for lig in system.ligands:
            codes.update(lig.ccd_code.split("-"))

    missing = _ENLICITIDE_MODIFIED_RESIDUES - codes
    assert not missing, f"modified residues standardized away / missing: {missing}"


def test_10sb_auth_residue_numbers_are_author_not_label(
    cif_10sb, mock_alternative_datasets
):
    """Integration: 10sb receptor chain A is label-numbered from 31 but author-
    numbered from 61 (a constant +30 offset). ``Residue.auth_number`` must carry
    the AUTHOR number (what RCSB / PyMOL show); before the fix it duplicated the
    label number, since the atoms were loaded with ``use_author_fields=False``.
    """
    from plinder.data.annotations.aggregate_annotations import Entry

    entry_dir = mock_alternative_datasets("10sb")
    entry = Entry.from_cif_file(cif_10sb, save_folder=entry_dir)

    chain = entry.chains["A"]
    first_num, first_res = min(chain.residues.items())
    assert (first_num, first_res.auth_number) == (31, "61")
    # every residue carries the author number, distinct from its label number
    assert not any(res.auth_number == str(num) for num, res in chain.residues.items())


def test_10sb_covalent_macrocycle_is_single_ligand(cif_10sb, mock_alternative_datasets):
    """Enlicitide is one covalent macrocycle spanning chains C + E + F.

    Chain C (cyclic peptide 21) is covalently bonded to E (A1C8P, 3 bonds)
    and F (GOA, 2 bonds). Ligand chains linked covalently are the same
    molecule, so the system must expose exactly ONE ligand covering all
    three chains — not three separate ligands.
    """
    from plinder.data.annotations.aggregate_annotations import Entry

    entry_dir = mock_alternative_datasets("10sb")
    entry = Entry.from_cif_file(cif_10sb, save_folder=entry_dir)

    # The enlicitide macrocycle is the system whose ligand(s) include FTR.
    macrocycle_systems = [
        s
        for s in entry.systems.values()
        if any("FTR" in lig.ccd_code for lig in s.ligands)
    ]
    assert len(macrocycle_systems) == 1
    ligands = macrocycle_systems[0].ligands

    assert (
        len(ligands) == 1
    ), f"expected 1 merged ligand, got {[lig.ccd_code for lig in ligands]}"
    lig = ligands[0]

    components = set(lig.ccd_code.split("-"))
    # merged ligand pulls in the covalently-bonded non-polymer partners...
    assert {"A1C8P", "GOA"} <= components
    # ...while still keeping the modified peptide residues.
    assert _ENLICITIDE_MODIFIED_RESIDUES <= components
    # All covale bonds are internal to the macrocycle (C<->E, C<->F); none
    # reach PCSK9, so the merged ligand is NOT receptor-covalent. Its
    # internal linkages are part of one molecule, not ligand->receptor links.
    assert lig.is_covalent is False
    # one connected molecule, not fragments
    assert len(lig.smiles.split(".")) == 1

    # The merged ligand must carry ALL of its atoms through the downstream
    # atom-selection and SDF-export paths, not just the primary chain.
    assert lig.member_asym_ids == ["C", "E", "F"]
    # selection references every member instance-chain
    for instance_chain in ("1.C", "1.E", "1.F"):
        assert f"cname='{instance_chain}'" in lig.selection
    # The member chains are SERIALIZED and drive the source-reconstruction
    # grouping, so the rebuilt SDF spans the whole molecule (not just the primary
    # chain — the truncation bug). Feed the REAL serialized rows through the
    # reconstruction helper:
    from plinder.core.index.system import _ligand_sdf_groups

    groups = _ligand_sdf_groups(entry.to_df())
    assert groups["1.C"] == ["1.C", "1.E", "1.F"]

    # the written SDF is one connected molecule spanning all member chains
    from rdkit import Chem

    sdf = entry_dir / "10sb" / "ligand_files" / f"{lig.asym_id}.sdf"
    assert sdf.is_file()
    mol = Chem.SDMolSupplier(str(sdf), removeHs=True)[0]
    assert mol is not None
    assert len(Chem.MolToSmiles(mol).split(".")) == 1
    assert mol.GetNumAtoms() == lig.num_heavy_atoms


@pytest.mark.parametrize(
    "fixture, holo_receptors, apo_proteins",
    [
        # 4ci1: chain B binds ligand D (one holo system); chain A is a
        # second, unliganded copy of the same protein in the ASU.
        ("cif_4ci1", {"B"}, {"A"}),
        # 2p1q: chains B + C form the liganded receptor (ligands D, E);
        # chain A is an apo copy that joins no system.
        ("cif_2p1q", {"B", "C"}, {"A"}),
    ],
    ids=["4ci1-apo-copy", "2p1q-apo-copies"],
)
def test_label_chains_marks_holo_and_apo_protein_chains(
    request, fixture, holo_receptors, apo_proteins
):
    """Entry.label_chains flags each protein chain holo/apo from the systems.

    A protein chain that is the receptor of a holo system is ``holo=True``;
    a protein chain that binds no ligand (e.g. a second, unliganded copy in
    the asymmetric unit) is ``holo=False``. Ground truth read from the CIFs:

    * 4ci1 - chain B binds ligand D; chain A is an apo copy.
    * 2p1q - chains B and C form the liganded receptor; chain A is apo.

    The second half corrupts every asserted flag and re-runs ``label_chains``
    directly, proving the labels are (re)computed by that method from
    ``entry.systems`` alone, not merely a leftover of the ingest pipeline.
    """
    cif = request.getfixturevalue(fixture)
    # No save_folder: data_dir stays None so no artifact reclassification —
    # matches the systems the labels below were read from.
    entry = Entry.from_cif_file(cif)

    # Labels produced by the pipeline (label_chains runs inside _finalize).
    for chain in holo_receptors:
        assert entry.chains[chain].holo is True, f"{chain} should be holo"
    for chain in apo_proteins:
        assert entry.chains[chain].holo is False, f"{chain} should be apo"

    # Isolation: set every asserted flag to the WRONG value, re-run
    # label_chains, and confirm it restores the correct labels. This
    # exercises both branches — holo receptors flipped False must go True,
    # apo chains flipped True must go False.
    for chain in holo_receptors:
        entry.chains[chain].holo = False
    for chain in apo_proteins:
        entry.chains[chain].holo = True
    entry.label_chains()
    for chain in holo_receptors:
        assert entry.chains[chain].holo is True
    for chain in apo_proteins:
        assert entry.chains[chain].holo is False


def test_is_excluded_mol_non_molecule_is_excluded():
    """An empty or unparseable SMILES is not a molecule — exclude it.

    A ligand we fail to characterize (CCD miss + RDKit failure) has no real
    SMILES; excluding it is correct, and an invalid SMILES (RDKit returns
    ``None``) must be excluded rather than crash the classifier.
    """
    assert is_excluded_mol("") is True  # empty -> not a molecule
    assert is_excluded_mol("not a smiles") is True  # unparseable (None) -> excluded
    # A genuine tiny fragment is also excluded (< 5 heavy atoms / < 2 carbons).
    assert is_excluded_mol("O") is True
    # A drug-like molecule is not excluded.
    assert is_excluded_mol("CC(=O)Oc1ccccc1C(=O)O") is False  # aspirin


def test_fill_missing_ccd_bonds():
    """A bond-less residue gets its intra-residue bonds back from the CCD (bt_info).

    Strip a standard component's bonds, then confirm _fill_missing_ccd_bonds
    restores them by matching atom names against the CCD dictionary.
    """
    import biotite.structure as struc
    import plinder.data.annotations.cif_utils as cu

    cu._get_ccd_atomarray.cache_clear()
    atoms = cu._get_ccd_atomarray("ATP")
    n_expected = atoms.bonds.as_array().shape[0]
    assert n_expected > 0

    # Simulate a residue that arrived with no internal bonds.
    stripped = atoms.copy()
    stripped.bonds = struc.BondList(stripped.array_length())
    assert stripped.bonds.as_array().shape[0] == 0

    filled = cu._fill_missing_ccd_bonds(stripped)
    assert filled.bonds.as_array().shape[0] == n_expected

    # Idempotent: a residue that already has bonds is untouched.
    again = cu._fill_missing_ccd_bonds(filled)
    assert again.bonds.as_array().shape[0] == n_expected


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
    from plinder.data.annotations.cif_utils import read_mmcif_file

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
    cif_2y4i, mock_alternative_datasets, monkeypatch
):
    import builtins
    import sys

    import biotite.structure as struc
    from biotite.sequence.io.fasta import FastaFile
    from biotite.structure.io import pdbx
    from plinder.data.annotations.cif_utils import read_mmcif_file

    entry_dir = mock_alternative_datasets("2y4i")
    system_tag = "2y4i__1__1.B__1.E_1.F"
    entry = Entry.from_cif_file(cif_2y4i, save_folder=entry_dir)
    assert entry.chains["A"].length == 319
    assert entry.chains["B"].length == 395
    assert entry.chains["A"].num_unresolved_residues >= 0
    assert entry.chains["B"].num_unresolved_residues >= 0

    canonical_ligand_dir = entry_dir / "2y4i" / "ligand_files"
    assert {path.name for path in canonical_ligand_dir.glob("*.sdf")} >= {
        "E.sdf",
        "F.sdf",
    }
    assert not (entry_dir / system_tag).exists()

    row = entry.to_df().query("system_id == @system_tag").iloc[0]
    repeated_membership_columns = {
        column
        for column in row.index
        if column.startswith(("system_biounit_chains", "system_other_chains"))
        or column.startswith(("system_biounit_non_water", "system_biounit_water"))
        or column.startswith(
            ("system_other_protein_chains", "system_other_ligand_chains")
        )
    }
    assert not repeated_membership_columns
    biounit_chains = entry.biounit_chains_to_df()
    assembly_chains = biounit_chains.query("biounit_id == '1'")
    system_members = set(
        row["system_protein_chains_asym_id"] + row["system_ligand_chains"]
    )
    other_receptor_chains = set(
        assembly_chains.loc[
            assembly_chains["chain_role"].eq("receptor")
            & ~assembly_chains["chain_instance"].isin(system_members),
            "chain_instance",
        ]
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
        biounit_chains=biounit_chains,
        options=SystemReconstructionOptions(
            system_waters="none",
            receptor_waters="none",
            system_include_other_protein_chains=True,
        ),
    )
    assert set(written) == {"system_cif", "receptor_cif", "sequences_fasta"}
    assert all(path.is_file() for path in written.values())
    assert not list(output_dir.glob("*.pdb"))
    sequences = dict(FastaFile.read_iter(outputs.sequences_fasta))
    assert set(sequences) == {"1.B"}
    assert len(sequences["1.B"]) == 395

    system_atoms = pdbx.get_structure(read_mmcif_file(outputs.system_cif), model=1)
    assert not np.any(struc.filter_solvent(system_atoms))
    expected_chains = set(row["system_protein_chains_asym_id"])
    expected_chains.update(row["system_ligand_chains"])
    expected_chains.update(other_receptor_chains)
    assert set(system_atoms.chain_id) == expected_chains

    # FASTA reconstruction is part of the base package and must not import
    # pipeline validation or OpenStructure dependencies.
    sys.modules.pop("plinder.data.annotations.protein_utils", None)
    original_import = builtins.__import__

    def reject_optional_imports(name, *args, **kwargs):
        if name == "ost" or name.startswith(("ost.", "PDBValidation")):
            raise AssertionError(f"unexpected optional import: {name}")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_optional_imports)
    dependency_free_fasta = output_dir / "dependency_free.fasta"
    save_reconstructed_system(
        cif_2y4i,
        row,
        outputs=SystemReconstructionOutputs(
            sequences_fasta=dependency_free_fasta,
        ),
        options=SystemReconstructionOptions(
            system_waters="none",
            receptor_waters="none",
        ),
    )
    assert dict(FastaFile.read_iter(dependency_free_fasta)) == sequences

    annotation = entry.to_df()
    query_calls = []

    def query_one_system(*, columns, splits, filters):
        assert columns == ["*"]
        assert splits == ["*"]
        assert len(filters) == 1
        column, operator, value = filters[0]
        assert operator == "=="
        query_calls.append((column, value))
        return annotation[annotation[column] == value].copy()

    monkeypatch.setattr("plinder.core.index.system.query_index", query_one_system)
    from plinder.core import PlinderSystem

    reconstructed_dir = output_dir / "from_plinder_system"
    plinder_system = PlinderSystem(
        system_id=system_tag,
        source_mmcif=cif_2y4i,
        reconstruction_dir=reconstructed_dir,
        canonical_ligand_dir=canonical_ligand_dir,
        biounit_chains=biounit_chains,
        reconstruction_options=SystemReconstructionOptions(
            system_waters="none",
            receptor_waters="none",
        ),
    )
    assert len(plinder_system.entry) == len(
        annotation[annotation["entry_pdb_id"] == "2y4i"]
    )
    assert len(plinder_system.system) == 2
    assert query_calls == [
        ("entry_pdb_id", "2y4i"),
        ("system_id", system_tag),
    ]

    canonical_sdfs = plinder_system.canonical_ligand_sdfs
    assert set(canonical_sdfs) == {"1.E", "1.F"}
    assert {Path(path).name for path in canonical_sdfs.values()} == {
        "E.sdf",
        "F.sdf",
    }
    receptor_path = Path(plinder_system.receptor_cif)
    assert receptor_path.is_file()
    assert not (reconstructed_dir / "system.cif").exists()
    assert not (reconstructed_dir / "sequences.fasta").exists()

    system_aligned_sdfs = plinder_system.ligand_sdfs
    assert set(system_aligned_sdfs) == {"1.E", "1.F"}
    assert all(Path(path).is_file() for path in system_aligned_sdfs.values())
    assert all(
        reconstructed_dir in Path(path).parents for path in system_aligned_sdfs.values()
    )


def test_smiles_from_nextgen(rcsb_ccd_reference_csv):
    """Validate CCD-derived structures against RCSB ground truth by InChIKey.

    The production ``_get_ccd_mol`` (stereo from ideal 3D) is compared to the RCSB
    reference InChIKey, which is canonical and stereo-inclusive: it covers
    connectivity and stereo, both tetrahedral R/S and double-bond E/Z. Metals and
    stereo-underspecified CCD entries (same skeleton, differing stereo layer) are
    tolerated.
    """
    from plinder.data.annotations.interaction_utils import _COORDINATION_METALS
    from plinder.data.annotations.ligand_utils import _get_ccd_mol
    from rdkit.Chem.inchi import MolToInchiKey

    # keep_default_na=False so the sodium comp_id "NA" reads as a string, not NaN
    rcsb_df = pd.read_csv(rcsb_ccd_reference_csv, keep_default_na=False)
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

        # InChIKey is canonical and stereo-inclusive (connectivity + tetrahedral
        # R/S + double-bond E/Z), so this single check validates the full structure.
        if ccd_inchikey != rcsb_inchikey:
            mismatches.append(
                (
                    comp_id,
                    "InChIKey",
                    ccd_inchikey,
                    f"expected {rcsb_inchikey}",
                )
            )

    assert len(mismatches) == 0, "CCD vs RCSB mismatches:\n" + "\n".join(
        f"  {m}" for m in mismatches
    )


def _build_resolved_mol(cif_path, chain_id):
    """Helper: build resolved mol from CIF chain using production code."""
    import biotite.structure.io.pdbx as pdbx
    from biotite.structure import filter_heavy
    from plinder.data.annotations.cif_utils import (
        atoms_to_rdkit_mol,
        read_mmcif_file,
    )

    cif_obj = read_mmcif_file(cif_path)
    atoms = pdbx.get_structure(
        cif_obj, model=1, use_author_fields=False, include_bonds=True
    )
    atoms = atoms[filter_heavy(atoms)]
    return atoms_to_rdkit_mol(atoms[atoms.chain_id == chain_id])


def _flip_first_chiral(mol):
    """Helper: return a copy whose 3D geometry is mirrored (the enantiomer).

    ``compare_stereo_to_template`` judges stereo from coordinates, so we
    invert the geometry (improper reflection through x=0) rather than a
    chiral tag — flipping a tag without moving atoms would be a no-op.
    Mirroring inverts every stereocenter at once. Returns None if the mol
    has no tetrahedral stereocenter (achiral — nothing to detect).
    """
    from rdkit.Geometry import Point3D

    tetrahedral = (
        Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
    )
    if not any(a.GetChiralTag() in tetrahedral for a in mol.GetAtoms()):
        return None
    flipped = Chem.Mol(mol)
    conf = flipped.GetConformer()
    for i in range(flipped.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        # mirror on x=0 plane
        conf.SetAtomPosition(i, Point3D(-p.x, p.y, p.z))
    return flipped


def test_stereo_check_single_residue(cif_7gj7):
    """Test _check_stereo_vs_template on single-residue ligands.

    Q0I (chain E): chiral — should match CCD, flipped should fail.
    DMS (chain C): achiral — should return None (no comparable centers).
    """
    from plinder.data.annotations.ligand_utils import _check_stereo_vs_template

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
    """Test _check_stereo_vs_template with a partially resolved ligand.

    JEF in 1ngx chain E has 28/41 heavy atoms resolved. compare_stereo_to_template
    checks only the stereo elements whose atoms are all resolved (comparing the
    template's descriptors against the resolved coordinates), so the resolved
    portion still yields a definite match/mismatch.
    """
    from plinder.data.annotations.ligand_utils import _check_stereo_vs_template

    jef_mol = _build_resolved_mol(cif_1ngx, "E")
    assert jef_mol.GetNumAtoms() < 41, "JEF should be partially resolved"

    # The resolved portion must MATCH the template (not merely be comparable).
    assert _check_stereo_vs_template(jef_mol) is True

    # And a mirror of the resolved portion must be detected as a mismatch —
    # JEF is chiral, so the flip is always available (assert, don't skip).
    jef_flipped = _flip_first_chiral(jef_mol)
    assert jef_flipped is not None, "JEF should have a chiral center to flip"
    assert _check_stereo_vs_template(jef_flipped) is False


def test_stereo_check_multi_residue(cif_6fx1):
    """Test _check_stereo_vs_template on multi-residue glycan.

    6fx1 chain M: NAG+BMA+MAN+FUC+C4W (25+ chiral centers).

    The per-residue check compares chiral handedness (signed volume), which
    is independent of CIP priority, so an inter-residue (glycosidic) bond does
    not cause the spurious mismatches that changed CIP priorities once did.

    We verify:
    1. The correctly-resolved glycan MATCHES its template (is True)
    2. A mirror (enantiomer) of the glycan is detected as a mismatch (is False)
    3. The mol has chiral centers that are being compared
    """
    from plinder.data.annotations.ligand_utils import _check_stereo_vs_template

    glycan_mol = _build_resolved_mol(cif_6fx1, "M")

    # Verify multi-residue composition
    res_names = {
        a.GetPDBResidueInfo().GetResidueName().strip()
        for a in glycan_mol.GetAtoms()
        if a.GetPDBResidueInfo()
    }
    assert len(res_names) > 1, f"Should be multi-residue, got {res_names}"

    # The correctly-resolved glycan must MATCH the per-residue templates.
    assert _check_stereo_vs_template(glycan_mol) is True

    # A mirror of the whole glycan inverts every stereocenter at once, so the
    # multi-residue path must flag it — this is what proves the check can
    # actually catch a stereo error across residues, not just "run".
    glycan_flipped = _flip_first_chiral(glycan_mol)
    assert glycan_flipped is not None, "Glycan should have chiral centers to flip"
    assert _check_stereo_vs_template(glycan_flipped) is False

    # Verify the mol actually has chiral centers. atoms_to_rdkit_mol assigns
    # chiral *tags* from 3D (not _CIPCode, which needs a CIP-labelling pass),
    # so count tetrahedral tags.
    tetrahedral = (
        Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
    )
    n_chiral = sum(1 for a in glycan_mol.GetAtoms() if a.GetChiralTag() in tetrahedral)
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


def test_mixed_receptor_type_is_written_to_annotation(cif_8ufz):
    entry = Entry.from_cif_file(cif_8ufz)

    assert entry.systems
    assert {system.receptor_type for system in entry.systems.values()} == {
        "protein+dna"
    }
    assert set(entry.to_df()["system_receptor_type"]) == {"protein+dna"}
    chain_types = entry.chains_to_df().set_index("chain_asym_id")["chain_receptor_type"]
    assert {chain_types[chain] for chain in ["A", "B", "C", "D"]} == {"dna"}
    assert {chain_types[chain] for chain in ["E", "F"]} == {"protein"}


@pytest.mark.parametrize(
    "cif_fixture, pdb_id, ccd_code, sdf_name, expected_aromatic_rings, golden_smiles",
    [
        # Imatinib (STI) in 2hyy, fully resolved: 2 benzene + pyridine + pyrimidine.
        (
            "cif_2hyy",
            "2hyy",
            "STI",
            "E.sdf",
            4,
            "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C",
        ),
        # EF2 (phthalimide + glutarimide) in 7bqu, fully resolved: one aromatic ring.
        (
            "cif_7bqu",
            "7bqu",
            "EF2",
            "C.sdf",
            1,
            "c1ccc2c(c1)C(=O)N(C2=O)[C@H]3CCC(=O)NC3=O",
        ),
    ],
)
def test_ligand_fix_to_valid(
    cif_fixture,
    pdb_id,
    ccd_code,
    sdf_name,
    expected_aromatic_rings,
    golden_smiles,
    mock_alternative_datasets,
    request,
):
    """A fully-resolved ligand is rebuilt to the correct covalent structure.

    Beyond the SDF-vs-SMILES self-agreement (two pipeline artifacts), both are
    checked against INDEPENDENT ground truth: the known aromatic-ring count and
    the reference structure's stereo-free InChIKey skeleton (first block —
    connectivity only, robust to protonation/canonicalization). This catches a
    wrong-but-self-consistent bond perception that an SDF-vs-SMILES check alone
    would pass.
    """
    from rdkit.Chem.inchi import MolToInchiKey
    from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings

    cif_path = request.getfixturevalue(cif_fixture)
    entry_dir = mock_alternative_datasets(pdb_id)
    entry = Entry.from_cif_file(cif_path, save_folder=entry_dir)

    lig = next(
        (
            lig
            for system in entry.systems.values()
            for lig in system.ligands
            if lig.ccd_code == ccd_code
        ),
        None,
    )
    assert lig is not None, f"{ccd_code} ligand not found in any system"
    assert lig.is_invalid is False

    outsdffile = entry_dir / pdb_id / "ligand_files" / sdf_name
    assert outsdffile.is_file()
    rdmol_sdf = Chem.SDMolSupplier(str(outsdffile), removeHs=True)[0]
    rdmol_smi = Chem.MolFromSmiles(lig.smiles)

    # SDF (3D-resolved) and SMILES agree AND match the independently known count.
    assert (
        CalcNumAromaticRings(rdmol_sdf)
        == CalcNumAromaticRings(rdmol_smi)
        == expected_aromatic_rings
    )

    # Independent identity: both artifacts match the reference structure's
    # stereo-free InChIKey skeleton.
    golden_skeleton = MolToInchiKey(Chem.MolFromSmiles(golden_smiles)).split("-")[0]
    assert MolToInchiKey(rdmol_sdf).split("-")[0] == golden_skeleton
    assert MolToInchiKey(rdmol_smi).split("-")[0] == golden_skeleton


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
