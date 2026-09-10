"""Prepare complete predictions without discarding bad or non-proper poses."""

import gzip
import shutil
from types import SimpleNamespace

import biotite.structure as struc
import numpy as np
import pandas as pd
import pytest
from biotite.structure.io import pdbx
from plinder.data.annotations.cif_utils import (
    MissingBondOrderError,
    _get_ccd_atomarray,
    get_structure_with_altloc,
    read_mmcif_file,
)
from plinder.data.annotations.save_utils import save_cif_file
from plinder.eval.commands import run_openstructure
from plinder.eval.inputs import prepare_prediction, reference_ligands
from rdkit import Chem


@pytest.fixture
def complete_model(test_dir):
    return test_dir / "reconstructed_systems/1avd__1__1.A__1.C/system.cif"


@pytest.mark.parametrize("displacement", [0, 100])
def test_prepare_retains_ligand_coordinates(complete_model, tmp_path, displacement):
    cif = read_mmcif_file(complete_model)
    atoms = get_structure_with_altloc(cif, include_bonds=True)
    ligand_atoms = atoms.chain_id == "1.C"
    atoms.coord[ligand_atoms] += displacement
    model = tmp_path / "model.cif"
    save_cif_file(atoms, "prediction", model, source_block=cif.block)
    prepared = prepare_prediction(model, tmp_path / "prepared")
    assert set(prepared.ligands) == {"1.C"}
    ligand = next(Chem.SDMolSupplier(str(prepared.ligands["1.C"])))
    expected = atoms.coord[ligand_atoms & struc.filter_heavy(atoms)]
    np.testing.assert_allclose(
        ligand.GetConformer().GetPositions(), expected, atol=1e-3
    )
    assert ligand.GetNumBonds() > 0
    receptor = get_structure_with_altloc(read_mmcif_file(prepared.receptor))
    assert set(receptor.chain_id) == {"1.A"}
    assert prepared.receptor_molecule.GetNumAtoms() == len(receptor)
    assert all(
        atom.GetPDBResidueInfo().GetChainId() == "1.A"
        for atom in prepared.receptor_molecule.GetAtoms()
    )


def test_preparation_keeps_ions_and_artifacts(complete_model, tmp_path):
    atoms = get_structure_with_altloc(
        read_mmcif_file(complete_model), include_bonds=True
    )
    # One sodium ion and one ethylene glycol, both far from the receptor.
    extra = struc.AtomArray(5)
    extra.chain_id = ["I", "E", "E", "E", "E"]
    extra.res_name = ["NA", "EDO", "EDO", "EDO", "EDO"]
    extra.atom_name = ["NA", "C1", "C2", "O1", "O2"]
    extra.element = ["Na", "C", "C", "O", "O"]
    extra.res_id = np.ones(5, dtype=int)
    extra.hetero[:] = True
    extra.coord = np.array(
        [[100, 0, 0], [200, 0, 0], [201.5, 0, 0], [198.6, 0, 0], [202.9, 0, 0]]
    )
    extra.bonds = struc.BondList(5, np.array([[1, 2, 1], [1, 3, 1], [2, 4, 1]]))
    model = tmp_path / "model.cif"
    save_cif_file(atoms + extra, "prediction", model)
    prepared = prepare_prediction(model, tmp_path / "prepared")
    assert set(prepared.ligands) == {"1.C", "I", "E"}


def test_prepare_compressed_uppercase_cif(complete_model, tmp_path):
    model = tmp_path / "prediction.CIF.gz"
    with gzip.open(model, "wb") as stream:
        stream.write(complete_model.read_bytes())
    assert prepare_prediction(model, tmp_path / "prepared").ligands


def test_prepare_retains_nucleic_acid_receptors(tmp_path):
    atoms = struc.AtomArray(6)
    atoms.chain_id = ["DNA"] * 3 + ["RNA"] * 3
    atoms.res_id = [1] * 6
    atoms.res_name = ["DA"] * 3 + ["A"] * 3
    atoms.atom_name = ["P", "O5'", "C5'"] * 2
    atoms.element = ["P", "O", "C"] * 2
    atoms.coord = np.array(
        [
            [0, 0, 0],
            [1.5, 0, 0],
            [2.9, 0, 0],
            [0, 10, 0],
            [1.5, 10, 0],
            [2.9, 10, 0],
        ]
    )
    atoms.bonds = struc.BondList(
        6, np.array([[0, 1, 1], [1, 2, 1], [3, 4, 1], [4, 5, 1]])
    )
    model = tmp_path / "model.cif"
    metadata = pdbx.CIFBlock()
    metadata["struct_asym"] = pdbx.CIFCategory(
        {"id": ["DNA", "RNA"], "entity_id": ["1", "2"]}
    )
    metadata["entity"] = pdbx.CIFCategory(
        {"id": ["1", "2"], "type": ["polymer", "polymer"]}
    )
    metadata["entity_poly"] = pdbx.CIFCategory(
        {
            "entity_id": ["1", "2"],
            "type": ["polydeoxyribonucleotide", "polyribonucleotide"],
            "pdbx_seq_one_letter_code": ["A", "A"],
            "pdbx_seq_one_letter_code_can": ["A", "A"],
        }
    )
    metadata["entity_poly_seq"] = pdbx.CIFCategory(
        {"entity_id": ["1", "2"], "num": [1, 1], "mon_id": ["DA", "A"]}
    )
    save_cif_file(atoms, "prediction", model, source_block=metadata)
    prepared = prepare_prediction(model, tmp_path / "prepared")
    assert prepared.ligands == {}
    assert {
        a.GetPDBResidueInfo().GetChainId()
        for a in prepared.receptor_molecule.GetAtoms()
    } == {"DNA", "RNA"}


@pytest.mark.parametrize("override", ["missing", "ccd", "both"])
@pytest.mark.parametrize("record_type", ["ATOM", "HETATM"])
def test_prepare_unknown_ligand(complete_model, tmp_path, override, record_type):
    atoms = get_structure_with_altloc(
        read_mmcif_file(complete_model), include_bonds=True
    )
    receptor = atoms[atoms.chain_id == "1.A"]
    # This fixture's deposited NAG is missing one atom. Use the complete CCD
    # molecule to exercise a successful exact chemistry override.
    ligand = _get_ccd_atomarray("NAG").copy()
    ligand = ligand[struc.filter_heavy(ligand)]
    ligand.chain_id[:] = "1.C"
    ligand.hetero[:] = True
    model = tmp_path / "model.cif"
    save_cif_file(receptor + ligand, "prediction", model)
    cif = read_mmcif_file(model)
    site = cif.block["atom_site"]
    names = site["label_comp_id"].as_array(str)
    assert "NAG" in names
    site["label_comp_id"] = np.where(names == "NAG", "XYZ", names)
    site["auth_comp_id"] = np.where(names == "NAG", "XYZ", names)
    site["group_PDB"] = np.where(
        names == "NAG", record_type, site["group_PDB"].as_array(str)
    )
    if "chem_comp_bond" in cif.block:
        del cif.block["chem_comp_bond"]
    cif.write(model)
    options = {}
    if override in {"ccd", "both"}:
        options["ligand_ccd_codes"] = {"XYZ": "NAG"}
    if override == "both":
        options["ligand_smiles"] = {"XYZ": "CCO"}
    if override == "missing":
        with pytest.raises(MissingBondOrderError, match="XYZ"):
            prepare_prediction(model, tmp_path / "prepared", **options)
    elif override == "both":
        with pytest.raises(ValueError, match="either SMILES or CCD"):
            prepare_prediction(model, tmp_path / "prepared", **options)
    elif override == "ccd":
        prepared = prepare_prediction(model, tmp_path / "prepared", **options)
        ligand = next(Chem.SDMolSupplier(str(prepared.ligands["1.C"])))
        assert ligand is not None and ligand.GetNumBonds() > 0


@pytest.mark.parametrize("record_type", ["ATOM", "HETATM"])
def test_prepare_positional_smiles(complete_model, tmp_path, record_type):
    atoms = get_structure_with_altloc(
        read_mmcif_file(complete_model), include_bonds=True
    )
    receptor = atoms[atoms.chain_id == "1.A"]
    ligand = struc.AtomArray(3)
    ligand.chain_id[:] = "L"
    ligand.res_name[:] = "XYZ"
    ligand.res_id[:] = 1
    ligand.hetero[:] = True
    ligand.atom_name = ["C1", "C2", "O1"]
    ligand.element = ["C", "C", "O"]
    ligand.coord = np.array([[100, 0, 0], [101.5, 0, 0], [102.9, 0, 0]])
    ligand.bonds = struc.BondList(3)
    model = tmp_path / "model.cif"
    save_cif_file(receptor + ligand, "prediction", model)
    cif = read_mmcif_file(model)
    site = cif.block["atom_site"]
    site["group_PDB"] = np.where(
        site["label_asym_id"].as_array(str) == "L",
        record_type,
        site["group_PDB"].as_array(str),
    )
    cif.write(model)
    prepared = prepare_prediction(
        model, tmp_path / "prepared", ligand_smiles={"XYZ": "CCO"}
    )
    molecule = next(Chem.SDMolSupplier(str(prepared.ligands["L"])))
    assert Chem.MolToSmiles(molecule) == "CCO"


@pytest.mark.skipif(shutil.which("ost") is None, reason="requires OpenStructure CLI")
def test_prepared_files_work_with_native_ost(complete_model, tmp_path):
    prepared = prepare_prediction(complete_model, tmp_path / "prepared")
    ligand = str(prepared.ligands["1.C"])
    result = run_openstructure(
        prepared.receptor,
        prepared.receptor,
        tmp_path / "comparison.json",
        action="compare-ligand-structures",
        options=["--model-ligands", ligand, "--reference-ligands", ligand],
    )
    assert result["rmsd"]["assigned_scores"][0]["score"] == pytest.approx(0, abs=1e-3)


def test_prepare_groups_covalently_connected_ligand_chains(complete_model, tmp_path):
    atoms = get_structure_with_altloc(
        read_mmcif_file(complete_model), include_bonds=True
    )
    atoms = atoms[atoms.chain_id == "1.A"]
    ligand = struc.AtomArray(4)
    ligand.chain_id = ["E", "E", "F", "F"]
    ligand.res_id[:] = 1
    ligand.res_name[:] = "EDO"
    ligand.atom_name = ["O1", "C1", "C2", "O2"]
    ligand.element = ["O", "C", "C", "O"]
    ligand.hetero[:] = True
    ligand.coord = np.array([[100, 0, 0], [101.4, 0, 0], [102.9, 0, 0], [104.3, 0, 0]])
    ligand.bonds = struc.BondList(4, np.array([[0, 1, 1], [1, 2, 1], [2, 3, 1]]))
    model = tmp_path / "model.cif"
    save_cif_file(atoms + ligand, "prediction", model)
    prepared = prepare_prediction(model, tmp_path / "prepared")
    assert set(prepared.ligands) == {"E"}
    molecule = next(Chem.SDMolSupplier(str(prepared.ligands["E"])))
    assert Chem.MolToSmiles(molecule) == "OCCO"


def test_preparation_does_not_group_ligands_through_a_metal(complete_model, tmp_path):
    atoms = get_structure_with_altloc(
        read_mmcif_file(complete_model), include_bonds=True
    )
    ligand = atoms[atoms.chain_id == "1.C"].copy()
    ligand.chain_id[:] = "L"
    ligand.coord += 10
    metal = struc.AtomArray(1)
    metal.chain_id[:] = "Z"
    metal.res_id[:] = 1
    metal.res_name[:] = "ZN"
    metal.atom_name[:] = "ZN"
    metal.element[:] = "Zn"
    metal.hetero[:] = True
    metal.coord[:] = [100, 100, 100]
    metal.bonds = struc.BondList(1)
    combined = atoms + ligand + metal
    zinc = len(combined) - 1
    original_ligand = int(np.flatnonzero(atoms.chain_id == "1.C")[0])
    combined.bonds.add_bond(original_ligand, zinc, struc.BondType.COORDINATION)
    combined.bonds.add_bond(len(atoms), zinc, struc.BondType.COORDINATION)
    model = tmp_path / "coordinated.cif"
    save_cif_file(combined, "prediction", model)
    reread = get_structure_with_altloc(read_mmcif_file(model), include_bonds=True)
    assert (
        np.count_nonzero(reread.bonds.as_array()[:, 2] == struc.BondType.COORDINATION)
        == 2
    )
    prepared = prepare_prediction(model, tmp_path / "prepared")
    assert set(prepared.ligands) == {"1.C", "L", "Z"}
    molecules = {
        key: next(Chem.SDMolSupplier(str(path)))
        for key, path in prepared.ligands.items()
    }
    assert molecules["Z"].GetNumAtoms() == 1
    assert all(
        atom.GetAtomicNum() != 30
        for key in ("1.C", "L")
        for atom in molecules[key].GetAtoms()
    )


def test_prepare_explicit_polymer_ligand(complete_model, tmp_path):
    atoms = get_structure_with_altloc(
        read_mmcif_file(complete_model), include_bonds=True
    )
    peptide = atoms[(atoms.chain_id == "1.A") & (atoms.res_id <= 5)].copy()
    peptide.chain_id[:] = "P"
    model = tmp_path / "model.cif"
    save_cif_file(atoms + peptide, "prediction", model)
    prepared = prepare_prediction(model, tmp_path / "prepared", ligand_chains=["P"])
    assert set(prepared.ligands) == {"1.C", "P"}
    receptor = get_structure_with_altloc(read_mmcif_file(prepared.receptor))
    assert set(receptor.chain_id) == {"1.A"}
    with pytest.raises(ValueError, match="absent from model"):
        prepare_prediction(model, tmp_path / "absent", ligand_chains=["absent"])


def test_prepare_coordinate_only_model(complete_model, tmp_path):
    cif = read_mmcif_file(complete_model)
    for category in list(cif.block):
        if category != "atom_site":
            del cif.block[category]
    model = tmp_path / "model.cif"
    cif.write(model)
    assert set(prepare_prediction(model, tmp_path / "prepared").ligands) == {"1.C"}
    del cif.block["atom_site"]["Cartn_x"]
    cif.write(model)
    with pytest.raises(ValueError, match="_atom_site.Cartn_x"):
        prepare_prediction(model, tmp_path / "malformed")


@pytest.mark.parametrize("include_all", [False, True])
def test_select_reference_ligands_only_from_annotation(tmp_path, include_all):
    reference = SimpleNamespace(
        system_id="reference",
        system=pd.DataFrame(
            {
                "ligand_id": ["proper", "ion", "artifact"],
                "ligand_instance_chain": ["1.L", "1.I", "1.E"],
                "ligand_instance_chains": [["1.L", "1.M"], ["1.I"], ["1.E"]],
                "ligand_is_proper": [True, False, False],
            }
        ),
        ligand_sdfs={
            "1.L": str(tmp_path / "whole_molecule.sdf"),
            "1.I": str(tmp_path / "ion.sdf"),
            "1.E": str(tmp_path / "artifact.sdf"),
        },
    )
    selected = reference_ligands(reference, include_all_ligands=include_all)
    assert set(selected) == (
        {"proper", "ion", "artifact"} if include_all else {"proper"}
    )
    assert selected["proper"] == tmp_path / "whole_molecule.sdf"
