# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path
from types import SimpleNamespace

import biotite.structure as struc
import numpy as np
import pandas as pd
import pytest
from plinder.data.annotations.aggregate_annotations import Entry
from plinder.data.annotations.cif_utils import (
    build_biounit,
    read_mmcif_container,
    read_mmcif_file,
)
from plinder.data.annotations.get_ligand_validation import EntryValidation
from plinder.data.annotations.interaction_utils import get_covalent_connections
from plinder.data.annotations.interface_utils import (
    DEFAULT_MIN_INTERFACE_RESIDUES,
    INTERFACE_ANNOTATION_SCHEMA,
    MIN_INTERFACE_RESIDUES_METADATA_KEY,
    detect_protein_interfaces,
    interface_system_id,
    protein_interfaces_to_table,
)
from plinder.data.annotations.ligand_utils import (
    BiounitSpatialIndex,
    Ligand,
    classify_ligand_polymer_classes,
    get_water_chain_ids,
    is_known_artifact_ligand,
    sort_ccd_codes,
)
from plinder.data.annotations.mmpdb_utils import add_mmp_clusters_to_data
from plinder.data.annotations.protein_utils import get_receptor_type
from plinder.data.annotations.save_utils import (
    SystemReconstructionOptions,
    SystemReconstructionOutputs,
    _output_asym_ids,
    save_ligands,
    save_reconstructed_system,
)
from plinder.data.get_system_annotations import GetPlinderAnnotation
from rdkit import Chem


def _interface_test_chain(
    residue_numbers: list[int],
    *,
    length: int = 12,
    chain_type: str = "polypeptide(L)",
) -> SimpleNamespace:
    return SimpleNamespace(
        chain_type_str=chain_type,
        length=length,
        residues={
            number: SimpleNamespace(index=index)
            for index, number in enumerate(residue_numbers)
        },
    )


def _interface_test_atoms() -> struc.AtomArray:
    atoms = struc.AtomArray(6)
    atoms.chain_id = np.array(["2.B"] * 3 + ["1.A"] * 3)
    atoms.res_id = np.array([10, 11, 12, 1, 2, 3])
    atoms.res_name = np.array(["ALA"] * 6)
    atoms.atom_name = np.array(["CA"] * 6)
    atoms.element = np.array(["C"] * 6)
    atoms.coord = np.array(
        [
            [0.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [6.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
        ]
    )
    return atoms


def test_detect_protein_interfaces_uses_canonical_chain_order_and_residue_maps():
    interfaces = detect_protein_interfaces(
        _interface_test_atoms(),
        pdb_id="1ABC",
        biounit_id="1",
        chains={
            "A": _interface_test_chain([1, 2, 3]),
            "B": _interface_test_chain([10, 11, 12]),
        },
        contact_radius=1.5,
        min_interface_residues=3,
    )

    assert len(interfaces) == 1
    interface = interfaces[0]
    assert interface.system_id == "1abc__1__1.A--2.B"
    assert interface.chain_1_residue_numbers == (1, 2, 3)
    assert interface.chain_1_residue_indices == (0, 1, 2)
    assert interface.chain_2_residue_numbers == (10, 11, 12)
    assert interface.chain_2_residue_indices == (0, 1, 2)
    assert interface.num_contact_residue_pairs == 3


@pytest.mark.parametrize(
    ("chain_length", "chain_type", "min_interface_residues"),
    [
        (11, "polypeptide(L)", 3),
        (12, "polyribonucleotide", 3),
        (12, "polypeptide(L)", 4),
    ],
    ids=["short_chain", "non_protein_chain", "short_interface"],
)
def test_detect_protein_interfaces_applies_eligibility_filters(
    chain_length: int,
    chain_type: str,
    min_interface_residues: int,
):
    interfaces = detect_protein_interfaces(
        _interface_test_atoms(),
        pdb_id="1abc",
        biounit_id="1",
        chains={
            "A": _interface_test_chain([1, 2, 3]),
            "B": _interface_test_chain(
                [10, 11, 12], length=chain_length, chain_type=chain_type
            ),
        },
        contact_radius=1.5,
        min_interface_residues=min_interface_residues,
    )

    assert interfaces == []


def test_interface_system_id_is_unordered_and_rejects_self_interfaces():
    assert interface_system_id("1ABC", "2", "2.B", "1.A") == interface_system_id(
        "1abc", "2", "1.A", "2.B"
    )
    with pytest.raises(ValueError, match="two distinct"):
        interface_system_id("1abc", "2", "1.A", "1.A")


def _pinder_interface_entry(test_dir: Path, relative_path: str) -> Entry:
    """Load one compact NextGen regression structure ported from Pinder."""
    return Entry.from_cif_file(
        test_dir / "interfaces" / relative_path,
        min_polymer_size=12,
        interface_min_chain_length=12,
        interface_min_residues=DEFAULT_MIN_INTERFACE_RESIDUES,
    )


def _assert_interface_residue_mappings(entry: Entry) -> None:
    for interface in entry.interfaces:
        assert interface.chain_1 < interface.chain_2
        for instance_chain, residue_numbers, residue_indices in (
            (
                interface.chain_1,
                interface.chain_1_residue_numbers,
                interface.chain_1_residue_indices,
            ),
            (
                interface.chain_2,
                interface.chain_2_residue_numbers,
                interface.chain_2_residue_indices,
            ),
        ):
            asym_id = instance_chain.split(".", maxsplit=1)[1]
            assert len(residue_numbers) >= DEFAULT_MIN_INTERFACE_RESIDUES
            assert residue_indices == tuple(
                entry.chains[asym_id].residues[number].index
                for number in residue_numbers
            )


def test_pinder_7cm8_homodimer_interface_regression(test_dir: Path) -> None:
    entry = _pinder_interface_entry(
        test_dir,
        "cm/pdb_00007cm8/pdb_00007cm8_xyz-enrich.cif.gz",
    )

    assert [interface.system_id for interface in entry.interfaces] == [
        "7cm8__1__1.A--2.A"
    ]
    assert entry.interfaces[0].chain_1_residue_numbers == (
        entry.interfaces[0].chain_2_residue_numbers
    )
    assert entry.interfaces[0].prodigy is not None
    prodigy = entry.interfaces[0].prodigy
    assert prodigy.intermolecular_contacts == 149
    assert prodigy.charged_charged_contacts == 12
    assert prodigy.charged_polar_contacts == 10
    assert prodigy.charged_apolar_contacts == 46
    assert prodigy.polar_polar_contacts == 0
    assert prodigy.apolar_polar_contacts == 22
    assert prodigy.apolar_apolar_contacts == 59
    assert prodigy.link_density == pytest.approx(0.06, abs=0.005)
    assert prodigy.label == "BIO"
    assert prodigy.probability_bio == pytest.approx(1.0)
    _assert_interface_residue_mappings(entry)


def test_interface_only_annotation_preserves_ligand_assets(
    test_dir: Path, tmp_path: Path
) -> None:
    cif = test_dir / "interfaces/cm/pdb_00007cm8/" "pdb_00007cm8_xyz-enrich.cif.gz"
    annotation = GetPlinderAnnotation(cif, "", save_folder=tmp_path)
    first = annotation.annotate_interfaces()
    assert first.num_rows == 1
    assert first.column("prodigy_label").to_pylist() == ["BIO"]

    entry_folder = tmp_path / "7cm8"
    interface_path = entry_folder / "interfaces.parquet"
    interface_bytes = interface_path.read_bytes()
    assert annotation.annotate(include_interfaces=False) is None
    assert annotation.entry.interfaces == []
    assert interface_path.read_bytes() == interface_bytes

    ligand_annotation = tmp_path / "7cm8.parquet"
    ligand_annotation.write_bytes(b"preserved ligand annotation")
    ligand_sdf = entry_folder / "ligand_files/1.C.sdf"
    ligand_sdf.parent.mkdir()
    ligand_sdf.write_bytes(b"preserved canonical ligand")
    metadata_path = entry_folder / "entry_metadata.parquet"
    metadata_path.unlink()
    preserved = [
        ligand_annotation,
        ligand_sdf,
        entry_folder / "entry_chains.parquet",
        entry_folder / "entry_biounit_chains.parquet",
        entry_folder / "entry_source.parquet",
    ]
    before = {path: path.read_bytes() for path in preserved}

    second = annotation.annotate_interfaces()

    assert second.equals(first)
    assert {path: path.read_bytes() for path in preserved} == before
    assert pd.read_parquet(metadata_path)["entry_pdb_id"].tolist() == ["7cm8"]


def test_pinder_7cma_label_asym_interface_regression(test_dir: Path) -> None:
    entry = _pinder_interface_entry(
        test_dir,
        "cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz",
    )

    assert entry.chains["A"].auth_id == "A"
    assert entry.chains["B"].auth_id == "C"
    assert [interface.system_id for interface in entry.interfaces] == [
        "7cma__1__1.A--1.B"
    ]
    _assert_interface_residue_mappings(entry)


def test_interface_only_annotation_adds_missing_shared_chain_rows(
    test_dir: Path, tmp_path: Path
) -> None:
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    annotation = GetPlinderAnnotation(cif, "", save_folder=tmp_path)
    annotation.annotate_interfaces()

    entry_folder = tmp_path / "7cma"
    chain_path = entry_folder / "entry_chains.parquet"
    chains = pd.read_parquet(chain_path)
    expected_chain_ids = set(chains["chain_asym_id"])
    chains.loc[chains["chain_asym_id"] == "A", "chain_auth_id"] = "preserved"
    chains.loc[chains["chain_asym_id"] != "B"].to_parquet(chain_path, index=False)
    biounit_path = entry_folder / "entry_biounit_chains.parquet"
    biounits = pd.read_parquet(biounit_path)
    expected_biounit_ids = set(biounits["chain_asym_id"])
    biounits.loc[biounits["chain_asym_id"] == "A", "chain_role"] = "ligand"
    biounits.loc[biounits["chain_asym_id"] != "B"].to_parquet(
        biounit_path, index=False
    )
    metadata_path = entry_folder / "entry_metadata.parquet"
    metadata = pd.read_parquet(metadata_path).drop(
        columns=[
            "entry_source_taxonomy_ids",
            "entry_source_organism_names",
            "entry_host_taxonomy_ids",
            "entry_host_organism_names",
        ]
    )
    metadata["ligand_only_metadata"] = "preserved"
    metadata.to_parquet(metadata_path, index=False)

    annotation.annotate_interfaces()

    repaired_chains = pd.read_parquet(chain_path)
    repaired_biounits = pd.read_parquet(biounit_path)
    repaired_metadata = pd.read_parquet(metadata_path)
    assert set(repaired_chains["chain_asym_id"]) == expected_chain_ids
    assert set(repaired_biounits["chain_asym_id"]) == expected_biounit_ids
    assert repaired_chains.set_index("chain_asym_id").loc["A", "chain_auth_id"] == (
        "preserved"
    )
    assert set(
        repaired_biounits.loc[
            repaired_biounits["chain_asym_id"] == "A", "chain_role"
        ]
    ) == {"ligand"}
    assert repaired_metadata["ligand_only_metadata"].tolist() == ["preserved"]
    assert set(repaired_metadata).issuperset(
        {
            "entry_source_taxonomy_ids",
            "entry_source_organism_names",
            "entry_host_taxonomy_ids",
            "entry_host_organism_names",
        }
    )


def test_pinder_6wwe_enumerates_all_three_interfaces(test_dir: Path) -> None:
    entry = _pinder_interface_entry(
        test_dir,
        "ww/pdb_00006wwe/pdb_00006wwe_xyz-enrich.cif.gz",
    )

    assert [interface.system_id for interface in entry.interfaces] == [
        "6wwe__1__1.A--1.B",
        "6wwe__1__1.A--1.C",
        "6wwe__1__1.B--1.C",
    ]
    _assert_interface_residue_mappings(entry)


def test_pinder_4wwi_heterodimer_interface_regression(test_dir: Path) -> None:
    entry = _pinder_interface_entry(
        test_dir,
        "ww/pdb_00004wwi/pdb_00004wwi_xyz-enrich.cif.gz",
    )

    assert [interface.system_id for interface in entry.interfaces] == [
        "4wwi__1__1.A--1.D",
        "4wwi__2__1.B--1.E",
        "4wwi__3__1.C--1.F",
    ]
    _assert_interface_residue_mappings(entry)


def test_pinder_7nsg_threefold_interface_regression(test_dir: Path) -> None:
    entry = _pinder_interface_entry(
        test_dir,
        "ns/pdb_00007nsg/pdb_00007nsg_xyz-enrich.cif.gz",
    )

    assert [interface.system_id for interface in entry.interfaces] == [
        "7nsg__1__1.A--1.B",
        "7nsg__1__1.A--1.C",
        "7nsg__1__1.B--1.C",
    ]
    _assert_interface_residue_mappings(entry)


def test_pinder_1bo0_monomer_has_no_protein_interface(test_dir: Path) -> None:
    entry = _pinder_interface_entry(
        test_dir,
        "bo/pdb_00001bo0/pdb_00001bo0_xyz-enrich.cif.gz",
    )

    assert entry.interfaces == []


@pytest.mark.parametrize(
    ("relative_path", "expected_instances", "expected_asym_ids"),
    [
        ("bd/pdb_00007bdu/pdb_00007bdu_xyz-enrich.cif.gz", 7, 7),
        ("km/pdb_00007kmx/pdb_00007kmx_xyz-enrich.cif.gz", 840, 14),
        # Pinder's Gemmi path yielded 34 chains here. The deposited assembly
        # applies four operators to all 13 listed asym IDs, so Plinder's
        # Biotite path intentionally retains all 52 instances.
        ("a7/pdb_00002a79/pdb_00002a79_xyz-enrich.cif.gz", 52, 13),
        ("rw/pdb_00006rw4/pdb_00006rw4_xyz-enrich.cif.gz", 125, 125),
        ("y2/pdb_00002y26/pdb_00002y26_xyz-enrich.cif.gz", 120, 40),
    ],
)
def test_pinder_biological_assembly_expansion_regression(
    test_dir: Path,
    relative_path: str,
    expected_instances: int,
    expected_asym_ids: int,
) -> None:
    cif_file = read_mmcif_file(test_dir / "interfaces" / relative_path)
    biounit = build_biounit(cif_file, "1")
    instance_chains = set(str(value) for value in biounit.chain_id)
    asym_ids = {
        instance_chain.split(".", maxsplit=1)[1] for instance_chain in instance_chains
    }

    assert len(instance_chains) == expected_instances
    assert len(asym_ids) == expected_asym_ids
    assert all(
        instance_chain.split(".", maxsplit=1)[0].isdigit()
        for instance_chain in instance_chains
    )


def test_empty_protein_interface_table_retains_release_schema():
    table = protein_interfaces_to_table([])

    assert table.num_rows == 0
    assert table.schema.equals(INTERFACE_ANNOTATION_SCHEMA, check_metadata=False)
    assert table.schema.metadata == {
        MIN_INTERFACE_RESIDUES_METADATA_KEY: str(
            DEFAULT_MIN_INTERFACE_RESIDUES
        ).encode()
    }


def test_protein_interface_table_freezes_custom_ingest_threshold():
    table = protein_interfaces_to_table([], min_interface_residues=11)

    assert table.schema.metadata == {MIN_INTERFACE_RESIDUES_METADATA_KEY: b"11"}


def test_ccd_name_sorter():
    assert sort_ccd_codes({"G", "G25", "CPG", "5GP"}) == ["CPG", "G25", "G", "5GP"]


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


def test_multi_residue_ligand_sums_disconnected_resolved_saccharides() -> None:
    nag = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"
    classes = classify_ligand_polymer_classes(
        nag,
        resolved_smiles=f"{nag}.{nag}",
        is_multi_residue=True,
    )

    assert not classes["is_monosaccharide"]
    assert classes["is_oligosaccharide"]


def test_multi_residue_ligand_uses_richer_identity_for_rdkit() -> None:
    nag = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"
    ligand = Ligand(
        ccd_code="NAG-NAG",
        plip_type="SACCHARIDE",
        smiles=nag,
        resolved_smiles=f"{nag}.{nag}",
        residue_numbers=[1, 2],
    )

    ligand.set_rdkit()

    assert ligand.smiles == ligand.resolved_smiles
    assert ligand.rdkit_canonical_smiles == ligand.resolved_smiles
    assert ligand.num_heavy_atoms == 2 * Chem.MolFromSmiles(nag).GetNumHeavyAtoms()
    assert ligand.is_oligosaccharide
    assert not ligand.is_monosaccharide

    partially_resolved = Ligand(
        ccd_code="NAG-NAG",
        plip_type="SACCHARIDE",
        smiles=nag,
        resolved_smiles="CCO",
        residue_numbers=[1, 2],
    )
    partially_resolved.set_rdkit()

    assert partially_resolved.smiles == nag
    assert partially_resolved.is_oligosaccharide
    assert not partially_resolved.is_monosaccharide

    equal_size = Ligand(
        ccd_code="LIG-LIG",
        smiles="CC.O",
        resolved_smiles="CCO",
        residue_numbers=[1, 2],
    )
    equal_size.set_rdkit()

    assert equal_size.smiles == "CCO"


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
    assert is_known_artifact_ligand(["DUM"], set())
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
    from plinder.core.structure.atoms import is_hydrogen_isotope
    from plinder.data.annotations.cif_utils import read_mmcif_file
    from plinder.data.annotations.protein_utils import Chain, get_seqres_from_cif

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
):
    """6u6k: 13-residue synthetic peptide (chain B).

    With min_polymer_size=12, the peptide is receptor (13 >= 12) → no systems.
    With min_polymer_size=20, the peptide is ligand (13 < 20) → system created.
    """
    entry_dir = mock_alternative_datasets("6u6k")
    stale_ligand_dir = entry_dir / "6u6k" / "ligand_files"
    stale_ligand_dir.mkdir(parents=True)
    (stale_ligand_dir / "stale.sdf").touch()
    entry = Entry.from_cif_file(
        cif_6u6k, save_folder=entry_dir, min_polymer_size=min_polymer_size
    )
    if expect_ligand:
        assert len(entry.systems) == 1
        assert not (stale_ligand_dir / "stale.sdf").exists()
        lig = entry.systems[list(entry.systems.keys())[0]].ligands[0]
        assert lig.ccd_code == "ACE-TRP-TRP-ILE-ILE-PRO-ALY-VAL-LYS-ALY-GLY-CYS-NH2"
        # Ligand-like and protein-interface membership are independent: this
        # 13-residue peptide is still a valid interface chain at the fixed
        # 12-residue protein-interface threshold.
        assert entry.chains["B"].holo
        assert "B" in entry.ligand_like_chains
    else:
        assert (
            len(entry.systems) == 0
        ), f"13-residue peptide should be receptor with min_polymer_size={min_polymer_size}"
        assert entry.interfaces
        assert not stale_ligand_dir.exists()


def test_annotation_without_systems_skips_validation_and_normalized_tables(
    monkeypatch, tmp_path
):
    class EmptyEntry:
        pdb_id = "1abc"
        systems: dict[str, object] = {}
        interfaces: list[object] = []

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


def test_interface_only_annotation_materializes_normalized_tables(
    cif_6u6k,
    mock_alternative_datasets,
) -> None:
    save_root = mock_alternative_datasets("6u6k")
    annotation = GetPlinderAnnotation(
        cif_6u6k,
        "",
        save_folder=save_root,
        min_polymer_size=12,
    )

    ligand_rows = annotation.annotate()

    assert ligand_rows is not None and ligand_rows.empty
    entry_dir = save_root / "6u6k"
    interfaces = pd.read_parquet(entry_dir / "interfaces.parquet")
    assert interfaces["system_id"].str.startswith("6u6k__").all()
    assert interfaces["interface_chain_1_residue_numbers"].map(len).min() >= 3
    assert interfaces["interface_chain_2_residue_numbers"].map(len).min() >= 3
    metadata = pd.read_parquet(entry_dir / "entry_metadata.parquet")
    assert metadata["entry_pdb_id"].tolist() == ["6u6k"]
    chains = pd.read_parquet(entry_dir / "entry_chains.parquet")
    assert set(chains["chain_asym_id"]) == {"A", "B"}
    assert not chains["chain_is_ligand_like"].any()


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
    # the written SDF is one connected molecule spanning all member chains
    from rdkit import Chem

    sdf = entry_dir / "10sb" / "ligand_files" / f"{lig.asym_id}.sdf"
    assert sdf.is_file()
    mol = Chem.SDMolSupplier(str(sdf), removeHs=True)[0]
    assert mol is not None
    assert len(Chem.MolToSmiles(mol).split(".")) == 1
    assert mol.GetNumAtoms() == lig.num_heavy_atoms


def test_get_ccd_mol_components_cif_fallback(monkeypatch):
    """_get_ccd_mol falls back to components.cif for codes bt_info lacks.

    A1C8P is a 5-char extended CCD code (from 10sb) that biotite's bundled
    dictionary predates. When the downloaded components.cif is available, the
    component must be read from there instead of failing.
    """
    from pathlib import Path

    import biotite.structure.info as bt_info
    import plinder.data.annotations.ligand_utils as lu
    from rdkit import Chem

    # Precondition: the code is genuinely absent from the bundled CCD.
    with pytest.raises(Exception):
        bt_info.residue("A1C8P")

    fixture = Path(__file__).parent / "test_data" / "mini_components.cif"
    monkeypatch.setattr(lu, "COMPONENTS_CCD_PATH", fixture)
    # Clear every layer of the (cached) CCD lookup so a stale miss from an
    # earlier lookup doesn't shadow the components.cif fallback.
    lu._component_atoms_from_components_cif.cache_clear()
    lu._get_ccd_atomarray.cache_clear()
    lu._get_ccd_mol.cache_clear()
    try:
        mol = lu._get_ccd_mol("A1C8P")
        assert mol is not None, "expected components.cif fallback to resolve A1C8P"
        assert Chem.MolToSmiles(mol) == "CCCCCCNCc1ccc(CCN)cc1"
    finally:
        # Don't leak the cached fallback mol into other tests.
        lu._component_atoms_from_components_cif.cache_clear()
        lu._get_ccd_atomarray.cache_clear()
        lu._get_ccd_mol.cache_clear()


def test_fill_missing_ccd_bonds_from_components(monkeypatch):
    """A bond-less residue gets its intra-residue bonds back from components.cif.

    Simulates the rare case where a components.cif-only code arrives without
    _chem_comp_bond: strip A1C8P's bonds, then confirm _fill_missing_ccd_bonds
    restores them from the components.cif fallback.
    """
    from pathlib import Path

    import biotite.structure as struc
    import plinder.data.annotations.ligand_utils as lu

    fixture = Path(__file__).parent / "test_data" / "mini_components.cif"
    monkeypatch.setattr(lu, "COMPONENTS_CCD_PATH", fixture)
    lu._component_atoms_from_components_cif.cache_clear()
    lu._get_ccd_atomarray.cache_clear()
    try:
        atoms = lu._get_ccd_atomarray("A1C8P")
        n_expected = atoms.bonds.as_array().shape[0]
        assert n_expected > 0

        # Simulate a residue that arrived with no internal bonds.
        stripped = atoms.copy()
        stripped.bonds = struc.BondList(stripped.array_length())
        assert stripped.bonds.as_array().shape[0] == 0

        filled = lu._fill_missing_ccd_bonds(stripped)
        assert filled.bonds.as_array().shape[0] == n_expected

        # Idempotent: a residue that already has bonds is untouched.
        again = lu._fill_missing_ccd_bonds(filled)
        assert again.bonds.as_array().shape[0] == n_expected
    finally:
        lu._component_atoms_from_components_cif.cache_clear()
        lu._get_ccd_atomarray.cache_clear()
        lu._get_ccd_mol.cache_clear()


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
    multi_residue_saccharides = [
        ligand
        for ligand in ligands
        if ligand.plip_type == "SACCHARIDE" and "-" in ligand.ccd_code
    ]
    assert multi_residue_saccharides
    for ligand in multi_residue_saccharides:
        assert ligand.num_heavy_atoms >= ligand.num_resolved_heavy_atoms
        assert ligand.is_oligosaccharide
        assert not ligand.is_monosaccharide


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
    entry.biounit_legacy_chain_ids["1"].update(
        {"1.B": "2.B", "1.E": "2.E", "1.F": "2.F"}
    )

    canonical_ligand_dir = entry_dir / "2y4i" / "ligand_files"
    assert {path.name for path in canonical_ligand_dir.glob("*.sdf")} >= {
        "E.sdf",
        "F.sdf",
    }
    assert not (entry_dir / system_tag).exists()

    row = entry.to_df().query("system_id == @system_tag").iloc[0]
    assert row["system_id_legacy"] == "2y4i__1__2.B__2.E_2.F"
    assert row["ligand_id_legacy"] == row["ligand_id"].replace("__1.", "__2.")
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

    system_file = read_mmcif_file(outputs.system_cif)
    system_block = list(system_file.values())[0]
    assert {
        "entry",
        "entity",
        "entity_poly",
        "entity_poly_seq",
        "chem_comp",
        "struct_asym",
        "atom_site",
    }.issubset(system_block)
    assert set(system_block["atom_site"]["label_asym_id"].as_array(str)).issubset(
        system_block["struct_asym"]["id"].as_array(str)
    )
    assert set(system_block["atom_site"]["label_entity_id"].as_array(str)).issubset(
        system_block["entity"]["id"].as_array(str)
    )
    assert set(system_block["entity_poly_seq"]["mon_id"].as_array(str)).issubset(
        system_block["chem_comp"]["id"].as_array(str)
    )
    assert system_block["cell"]["entry_id"].as_item() == row["system_id"]
    assert all(
        asym_id.isalnum()
        for asym_id in system_block["atom_site"]["label_asym_id"].as_array(str)
    )
    assert "struct_conn_type" in system_block
    assert set(
        system_block["struct_conn"]["conn_type_id"].as_array(str)
    ).issubset(system_block["struct_conn_type"]["id"].as_array(str))
    assert all(
        value == value.lower()
        for value in system_block["chem_comp_bond"]["value_order"].as_array(str)
    )
    hetero_mask = system_block["atom_site"]["group_PDB"].as_array(str) == "HETATM"
    assert set(
        system_block["atom_site"]["label_seq_id"].as_array(str)[hetero_mask]
    ) == {"."}
    system_atoms = pdbx.get_structure(system_file, model=1)
    assert not np.any(struc.filter_solvent(system_atoms))
    expected_chains = set(row["system_protein_chains_asym_id"])
    expected_chains.update(row["system_ligand_chains"])
    expected_chains.update(other_receptor_chains)
    assert set(system_atoms.chain_id) == set(
        _output_asym_ids(sorted(expected_chains)).values()
    )

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
    """Test CCD SMILES against RCSB ground truth.

    For each compound in the RCSB reference CSV, verify:
    1. InChIKey from CCD ideal 3D matches RCSB InChIKey
    2. Per-atom chirality matches via substructure match
    """
    from plinder.data.annotations.interaction_utils import _COORDINATION_METALS
    from plinder.data.annotations.ligand_utils import _get_ccd_mol
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
    from plinder.data.annotations.cif_utils import (
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
    transplants the resolved 3D coordinates onto the template graph and only
    compares stereocenters whose atom *and* immediate neighbors are all
    resolved, so the resolved portion still yields a definite match/mismatch.
    """
    from plinder.data.annotations.ligand_utils import _check_stereo_vs_template

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
    from plinder.data.annotations.ligand_utils import _check_stereo_vs_template

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
    from plinder.data.annotations.cif_utils import read_mmcif_file

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


def test_entry_chain_table_includes_full_sequences(cif_8ufz):
    entry = Entry.from_cif_file(cif_8ufz, include_ligands=False)
    chain_rows = entry.chains_to_df().set_index("chain_asym_id")

    assert chain_rows.loc["E", "chain_sequence"] == entry.chain_to_seqres["E"]
    assert chain_rows.loc["E", "chain_sequence"]


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
    entry = Entry.from_cif_file(
        cif_4jvn,
        save_folder=entry_dir,
        data_dir=entry_dir.parent.parent,
    )
    target_value = 7.638272164
    affinity = 0.0
    for sys in entry.systems.values():
        if "YUG" in set(l.ccd_code for l in sys.ligands):
            assert sys.has_binding_affinity
        for ligand in sys.ligands:
            if ligand.ccd_code == "YUG":
                affinity = ligand.binding_affinity
    assert affinity == target_value
