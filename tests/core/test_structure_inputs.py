"""The shared table keeps each structure associated with its own ligand poses."""

import gzip
import importlib

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from plinder.core.structure.inputs import (
    StructureInput,
    read_structure_table,
    write_search_structure,
)


@pytest.mark.parametrize("suffix", [".csv", ".tsv", ".parquet"])
def test_table_groups_ligands_and_resolves_relative_paths(tmp_path, suffix):
    frame = pd.DataFrame(
        {
            "input_id": ["one", "one", "two"],
            "structure_path": ["one.pdb", "one.pdb", "two.cif"],
            "ligand_path": ["a.sdf", "b.sdf", ""],
            "reference_id": ["1avd", "1avd", "2e31"],
        }
    )
    path = tmp_path / f"inputs{suffix}"
    if suffix == ".parquet":
        frame.to_parquet(path)
    else:
        frame.to_csv(path, sep="\t" if suffix == ".tsv" else ",", index=False)
    inputs = read_structure_table(path, require_reference=True)
    assert inputs["one"].coordinates == tmp_path / "one.pdb"
    assert inputs["one"].ligand_sdfs == (tmp_path / "a.sdf", tmp_path / "b.sdf")
    assert inputs["one"].reference_id == "1avd"
    assert inputs["two"].ligand_sdfs is None


@pytest.mark.parametrize(
    "field,values,error",
    [
        ("input_id", ["../escape", "one"], "Invalid input_id"),
        ("structure_path", ["a.cif", "b.cif"], "Conflicting structure_path"),
        ("reference_id", ["1avd", "2e31"], "Conflicting reference_id"),
        ("ligand_path", ["a.sdf", "a.sdf"], "Repeated ligand"),
    ],
)
def test_table_rejects_ambiguous_inputs(field, values, error):
    frame = pd.DataFrame(
        {
            "input_id": ["one", "one"],
            "structure_path": ["a.cif", "a.cif"],
            "reference_id": ["1avd", "1avd"],
            "ligand_path": ["a.sdf", "b.sdf"],
        }
    )
    frame[field] = values
    with pytest.raises(ValueError, match=error):
        read_structure_table(frame)


def test_reference_column_required_only_for_evaluation():
    frame = pd.DataFrame({"input_id": ["one"], "structure_path": ["one.cif"]})
    assert read_structure_table(frame)["one"].reference_id is None
    with pytest.raises(ValueError, match="reference_id"):
        read_structure_table(frame, require_reference=True)


def test_search_cif_preserves_separate_ligand_chemistry(test_dir, tmp_path):
    from plinder.data.annotations.cif_utils import (
        atoms_to_rdkit_mol,
        check_cif_bond_orders,
        check_custom_mmcif_fields,
        get_structure_with_altloc,
        read_mmcif_file,
    )

    root = test_dir / "reconstructed_systems/1avd__1__1.A__1.C"
    ligand_path = root / "ligand_files/1.C.sdf"
    second = tmp_path / "other.sdf"
    molecule = Chem.MolFromSmiles("C[NH3+]")
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, (100, 0, 0))
    conf.SetAtomPosition(1, (101.4, 0, 0))
    molecule.AddConformer(conf)
    with Chem.SDWriter(str(second)) as writer:
        writer.write(molecule)
    structure = StructureInput(root / "receptor.cif", [ligand_path, second])
    output = write_search_structure(
        structure, tmp_path / "query.cif", include_ligands=True
    )
    cif = read_mmcif_file(output)
    check_custom_mmcif_fields(
        cif.block, source=output, structure_mode="as_is", require_label_ids=True
    )
    check_cif_bond_orders(cif)
    atoms = get_structure_with_altloc(cif, include_bonds=True)
    assert set(atoms.chain_id) == {"1.A", "L0001", "L0002"}
    for chain, source in [("L0001", ligand_path), ("L0002", second)]:
        expected = next(Chem.SDMolSupplier(str(source)))
        selected = atoms[atoms.chain_id == chain]
        actual = atoms_to_rdkit_mol(selected)
        assert Chem.MolToSmiles(Chem.RemoveHs(actual)) == Chem.MolToSmiles(
            Chem.RemoveHs(expected)
        )
        np.testing.assert_allclose(
            selected.coord, expected.GetConformer().GetPositions(), atol=1e-3
        )
    mapping = pd.read_csv(output.with_suffix(".ligands.tsv"), sep="\t")
    assert mapping["ligand_path"].tolist() == [
        str(ligand_path.resolve()),
        str(second.resolve()),
    ]


@pytest.mark.parametrize(
    "smiles",
    ["[B]12[B]3[B]4[B]1[B]5[B]6[B]2[B]3[B]45C6", "[Be](C)(C)(C)C"],
    ids=["boron-cage", "beryllium"],
)
def test_search_accepts_supported_overvalent_sdfs(test_dir, tmp_path, smiles):
    from biotite.interface.rdkit import from_mol

    from plinder.core.utils.sanitize import sanitize
    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )

    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    sanitize(molecule)
    # These are supported by PLINDER but still rejected by strict RDKit checks.
    with pytest.raises(Chem.AtomValenceException):
        Chem.SanitizeMol(Chem.Mol(molecule))
    coordinates = np.arange(molecule.GetNumAtoms() * 3).reshape(-1, 3).astype(float)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, position in enumerate(coordinates):
        conformer.SetAtomPosition(index, position)
    molecule.AddConformer(conformer)
    sdf = tmp_path / "pose.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    original = sdf.read_bytes()
    receptor = test_dir / "reconstructed_systems/1avd__1__1.A__1.C/receptor.cif"
    output = write_search_structure(
        StructureInput(receptor, [sdf]), tmp_path / "query.cif", include_ligands=True
    )
    atoms = get_structure_with_altloc(read_mmcif_file(output), include_bonds=True)
    ligand = atoms[atoms.chain_id == "L0001"]
    expected = from_mol(molecule, conformer_id=0, add_hydrogen=False)
    np.testing.assert_allclose(ligand.coord, coordinates, atol=1e-3)
    np.testing.assert_array_equal(ligand.element, expected.element)
    np.testing.assert_array_equal(ligand.charge, expected.charge)
    np.testing.assert_array_equal(ligand.bonds.as_array(), expected.bonds.as_array())
    assert sdf.read_bytes() == original


def test_search_rejects_unrepairable_ligand_valence(test_dir, tmp_path):
    molecule = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
    molecule.UpdatePropertyCache(strict=False)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, (float(index), 0.0, 0.0))
    molecule.AddConformer(conformer)
    sdf = tmp_path / "invalid.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    receptor = test_dir / "reconstructed_systems/1avd__1__1.A__1.C/receptor.cif"
    with pytest.raises(Chem.AtomValenceException):
        write_search_structure(
            StructureInput(receptor, [sdf]),
            tmp_path / "query.cif",
            include_ligands=True,
        )


@pytest.mark.parametrize("suffix", [".pdb", ".PDB.gz"])
@pytest.mark.parametrize("include_ligands", [False, True])
def test_pdb_search_sequences_exclude_nonpolymers(
    test_dir, tmp_path, suffix, include_ligands
):
    import biotite.structure as struc
    from biotite.structure.io import pdb

    from plinder.core.scores.custom import write_custom_query_files
    from plinder.data.annotations.cif_utils import (
        get_structure_with_altloc,
        read_mmcif_file,
    )

    # Waters/ions share author chains with proteins; MSE is also HETATM.
    atoms = struc.AtomArray(9)
    atoms.chain_id = np.array(["A"] * 6 + ["B"] * 3)
    atoms.res_name = np.array(
        ["ALA", "HOH", "MSE", "ZN", "GLY", "HOH", "SER", "CL", "THR"]
    )
    atoms.res_id = np.array([10, 501, 10, 502, 20, 503, 40, 601, 45])
    atoms.ins_code[2] = "A"
    atoms.hetero = np.array([False, True, True, True, False, True, False, True, False])
    atoms.atom_name = np.array(["CA", "O", "CA", "ZN", "CA", "O", "CA", "CL", "CA"])
    atoms.element = np.array(["C", "O", "C", "Zn", "C", "O", "C", "Cl", "C"])
    atoms.coord = np.arange(27, dtype=float).reshape(9, 3)
    source = tmp_path / f"receptor{suffix}"
    pdb_file = pdb.PDBFile()
    pdb_file.set_structure(atoms)
    with (
        gzip.open(source, "wt") if suffix.endswith(".gz") else source.open("w")
    ) as stream:
        pdb_file.write(stream)
    original = source.read_bytes()
    sdf = test_dir / "reconstructed_systems/1avd__1__1.A__1.C/ligand_files/1.C.sdf"
    output = write_search_structure(
        StructureInput(source, [sdf] if include_ligands else None),
        tmp_path / "query.cif",
        include_ligands=include_ligands,
    )
    cif = read_mmcif_file(output)
    converted = get_structure_with_altloc(cif, include_bonds=True)
    protein = converted[np.isin(converted.chain_id, ["A", "B"])]
    assert protein.res_name.tolist() == ["ALA", "MSE", "GLY", "SER", "THR"]
    assert protein.chain_id.tolist() == ["A", "A", "A", "B", "B"]
    assert protein.res_id.tolist() == [1, 2, 3, 1, 2]
    assert protein.ins_code.tolist() == ["", "A", "", "", ""]
    np.testing.assert_allclose(protein.coord, atoms.coord[[0, 2, 4, 6, 8]])
    sites = cif.block["atom_site"]
    selected = np.isin(sites["label_asym_id"].as_array(str), ["A", "B"])
    assert sites["auth_seq_id"].as_array(int)[selected].tolist() == [10, 10, 20, 40, 45]
    assert cif.block["entity_poly_seq"]["mon_id"].as_array(str).tolist() == [
        "ALA",
        "MSE",
        "GLY",
        "SER",
        "THR",
    ]
    assert ("L0001" in converted.chain_id) == include_ligands
    queries = write_custom_query_files(
        [output], work_dir=tmp_path / "search", min_chain_length=1
    )
    manifest = pd.read_parquet(queries.chain_manifest).set_index("chain_asym_id")
    assert set(manifest.index) == {"A", "B"}
    assert manifest.loc["A", "sequence_length"] == 3
    assert manifest.loc["B", "sequence"] == "ST"
    assert list(manifest.loc["A", "resolved_residue_numbers"]) == [1, 2, 3]
    assert list(manifest.loc["B", "resolved_residue_numbers"]) == [1, 2]
    assert source.read_bytes() == original


def test_search_dispatches_fasta_folder_and_table(test_dir, tmp_path, monkeypatch):
    module = importlib.import_module("plinder.core.scores.search")
    calls = []
    monkeypatch.setattr(
        module,
        "score_custom_sequence_file",
        lambda source, **kw: calls.append((source, kw)),
    )
    monkeypatch.setattr(
        module,
        "score_custom_cif_files",
        lambda sources, **kw: calls.append((sources, kw)),
    )
    module.search(tmp_path / "query.fasta", output_dir=tmp_path / "fasta_search")
    assert calls[-1][1]["backends"] == ("mmseqs",)
    folder = tmp_path / "models"
    folder.mkdir()
    original = test_dir / "reconstructed_systems/1avd__1__1.A__1.C/receptor.cif"
    (folder / "one.cif").write_bytes(original.read_bytes())
    module.search(folder, output_dir=tmp_path / "folder_search")
    assert calls[-1][1]["include_ligands"] is False
    table = pd.DataFrame({"input_id": ["custom_id"], "structure_path": [str(original)]})
    module.search(table, output_dir=tmp_path / "table_search", mode="ligands")
    assert calls[-1][0][0].name == "custom_id.cif"
    assert calls[-1][1]["include_ligands"] is True
    with pytest.raises(ValueError, match="input table"):
        module.search(folder, output_dir=tmp_path / "bad", mode="ligands")
    with pytest.raises(ValueError, match="FASTA inputs"):
        module.search(
            tmp_path / "query.fa", output_dir=tmp_path / "bad", mode="interfaces"
        )


def test_search_annotation_keeps_dotted_chain_ids(test_dir, tmp_path, monkeypatch):
    from plinder.core.scores.custom import annotate_custom_cif_files
    from plinder.data.annotations import ligand_utils

    monkeypatch.setattr(ligand_utils, "BINDING_AFFINITY", {})
    root = test_dir / "reconstructed_systems/1avd__1__1.A__1.C"
    source = StructureInput(root / "receptor.cif", [root / "ligand_files/1.C.sdf"])
    cif = write_search_structure(source, tmp_path / "paired.cif", include_ligands=True)
    annotations = annotate_custom_cif_files(
        [cif], work_dir=tmp_path / "annotated", include_interfaces=False
    )
    entry = annotations.entries_by_structure["paired"]
    rows = pd.read_parquet(annotations.annotation_table)
    assert not rows.empty
    assert "1.A" in entry.chains
    assert entry.chains["1.A"].holo
    assert all(
        "1.1.A" in system.pocket_residue_number_to_index
        for system in entry.systems.values()
    )
    assert not pd.read_parquet(annotations.entry_chains).empty
