# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from biotite.structure.io import pdbx
from plinder.core.scores import custom
from plinder.core.scores.entries import (
    ChainView,
    EntryView,
    InterfaceView,
    LigandView,
    SystemView,
)


def _write_database_prefix(prefix: Path, *, indexed: bool = False) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".dbtype").write_bytes(b"db")
    if indexed:
        Path(f"{prefix}.idx.dbtype").write_bytes(b"index")


def _write_search_bundle(root: Path, backend: str) -> None:
    root.mkdir(parents=True, exist_ok=True)
    contract = {
        "identity": 1.0,
        "coverage": 1.0,
        "coverage_mode": 0,
        "compressed_search_target": False,
    }
    if backend == "foldseek":
        manifest = {
            **contract,
            "alignment_type": backend,
            "portable": True,
            "search_target": "clustered",
            "conversion_target": "clustered",
        }
        _write_database_prefix(root / "clustered", indexed=True)
    else:
        manifest = {
            **contract,
            "alignment_type": backend,
            "portable": True,
            "search_target": "representatives",
            "conversion_target": "holo_mmseqs",
            "cluster_alignments": "cluster_alignments",
        }
        _write_database_prefix(root / "representatives", indexed=True)
        _write_database_prefix(root / "holo_mmseqs")
        _write_database_prefix(root / "cluster_alignments")
    (root / "exact_cluster.json").write_text(json.dumps(manifest))


def _write_index(root: Path) -> None:
    index = root / "index"
    index.mkdir(parents=True, exist_ok=True)
    for filename in [
        "annotation_table.parquet",
        "entry_chains.parquet",
        "interface_annotation_table.parquet",
        "alignment_chain_lookup.parquet",
    ]:
        (index / filename).touch()


def test_resolve_custom_scoring_assets_accepts_ingest_layout(tmp_path):
    _write_index(tmp_path)
    for backend in custom.SEARCH_BACKENDS:
        _write_search_bundle(
            tmp_path / "dbs" / "subdbs" / f"holo_{backend}", backend
        )
    archive = tmp_path / "ligand_archives" / "ab.parquet"
    archive.parent.mkdir()
    archive.touch()

    assets = custom.resolve_custom_scoring_assets(
        data_dir=tmp_path,
        ligand_pdb_ids=["1abc", "2abd__1__1.A__1.L"],
    )

    assert assets.annotation_table == tmp_path / "index/annotation_table.parquet"
    assert assets.search_databases["foldseek"].search_target.name == "clustered"
    assert assets.search_databases["mmseqs"].cluster_alignments is not None
    assert assets.ligand_archives == {"ab": archive}


def test_remote_resolution_requests_only_bounded_assets(tmp_path, monkeypatch):
    _write_index(tmp_path)
    for backend in custom.SEARCH_BACKENDS:
        _write_search_bundle(
            tmp_path / "search_databases" / f"holo_{backend}", backend
        )
    archive = tmp_path / "ligand_archives" / "xy.parquet"
    archive.parent.mkdir()
    archive.touch()
    requested: list[str] = []

    def get_plinder_path(*, rel: str) -> Path:
        requested.append(rel)
        return tmp_path / rel

    monkeypatch.setattr(custom.cpl, "get_plinder_path", get_plinder_path)

    assets = custom.resolve_custom_scoring_assets(ligand_pdb_ids=["3xyz"])

    assert set(requested) == {
        "index/annotation_table.parquet",
        "index/entry_chains.parquet",
        "index/interface_annotation_table.parquet",
        "index/alignment_chain_lookup.parquet",
        "search_databases/holo_foldseek",
        "search_databases/holo_mmseqs",
        "ligand_archives/xy.parquet",
    }
    assert assets.ligand_archives == {"xy": archive}


def test_ligand_coordinates_are_not_resolved_before_targets_are_known(
    tmp_path, monkeypatch
):
    _write_index(tmp_path)
    for backend in custom.SEARCH_BACKENDS:
        _write_search_bundle(
            tmp_path / "search_databases" / f"holo_{backend}", backend
        )
    requested: list[str] = []

    def get_plinder_path(*, rel: str) -> Path:
        requested.append(rel)
        return tmp_path / rel

    monkeypatch.setattr(custom.cpl, "get_plinder_path", get_plinder_path)

    assets = custom.resolve_custom_scoring_assets()

    assert assets.ligand_archives == {}
    assert not any(path.startswith("ligand_archives/") for path in requested)


def test_search_database_requires_portable_manifest(tmp_path):
    root = tmp_path / "search_databases" / "holo_foldseek"
    _write_search_bundle(root, "foldseek")
    manifest = json.loads((root / "exact_cluster.json").read_text())
    manifest["portable"] = False
    (root / "exact_cluster.json").write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="not portable"):
        custom.resolve_search_database("foldseek", data_dir=tmp_path)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("identity", 0.9),
        ("coverage", 0.9),
        ("coverage_mode", 1),
        ("compressed_search_target", True),
    ],
)
def test_search_database_requires_exact_cluster_contract(tmp_path, field, value):
    root = tmp_path / "search_databases" / "holo_foldseek"
    _write_search_bundle(root, "foldseek")
    manifest = json.loads((root / "exact_cluster.json").read_text())
    manifest[field] = value
    (root / "exact_cluster.json").write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="exact-cluster contract"):
        custom.resolve_search_database("foldseek", data_dir=tmp_path)


def test_search_database_rejects_nested_external_link(tmp_path):
    root = tmp_path / "search_databases" / "holo_foldseek"
    _write_search_bundle(root, "foldseek")
    outside = tmp_path / "outside"
    outside.write_text("external")
    nested = root / "nested"
    nested.mkdir()
    (nested / "external").symlink_to(outside)

    with pytest.raises(ValueError, match="non-portable database link"):
        custom.resolve_search_database("foldseek", data_dir=tmp_path)


def test_search_database_rejects_absolute_internal_link(tmp_path):
    root = tmp_path / "search_databases" / "holo_foldseek"
    _write_search_bundle(root, "foldseek")
    (root / "absolute").symlink_to((root / "clustered.dbtype").resolve())

    with pytest.raises(ValueError, match="non-portable database link"):
        custom.resolve_search_database("foldseek", data_dir=tmp_path)


def test_missing_offline_asset_has_actionable_error(tmp_path, monkeypatch):
    monkeypatch.setattr(custom.cpl, "is_offline", lambda: True)

    with pytest.raises(FileNotFoundError, match="offline cache"):
        custom.resolve_custom_scoring_assets(data_dir=tmp_path, backends=())


def test_write_custom_query_files_uses_protein_label_asym_ids(
    test_dir, tmp_path
):
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"

    inputs = custom.write_custom_query_files([cif], work_dir=tmp_path)
    chains = pd.read_parquet(inputs.chain_manifest)

    assert chains["chain_asym_id"].tolist() == ["A", "B"]
    assert chains["query_id"].tolist() == ["cq00000000", "cq00000001"]
    assert chains["query_chain_id"].tolist() == [
        "pdb_00007cma_xyz-enrich__A",
        "pdb_00007cma_xyz-enrich__B",
    ]
    assert chains["sequence_length"].min() >= 12
    assert set(chains["sequence_source"]) == {"polymer"}
    assert len(list(inputs.chain_cif_dir.glob("*.cif"))) == 2
    assert inputs.sequence_fasta.read_text().count(">cq") == 2
    block = list(
        pdbx.CIFFile.read(str(inputs.chain_cif_dir / "cq00000000.cif")).values()
    )[0]
    assert {
        "entry",
        "entity",
        "entity_poly",
        "entity_poly_seq",
        "chem_comp",
        "struct_asym",
        "atom_site",
    }.issubset(block)
    assert set(block["atom_site"]["label_asym_id"].as_array(str)).issubset(
        block["struct_asym"]["id"].as_array(str)
    )
    assert set(block["atom_site"]["label_entity_id"].as_array(str)).issubset(
        block["entity"]["id"].as_array(str)
    )
    assert set(block["entity_poly_seq"]["mon_id"].as_array(str)).issubset(
        block["chem_comp"]["id"].as_array(str)
    )


def test_write_custom_query_files_accepts_coordinate_only_mmcif(
    test_dir, tmp_path
):
    from plinder.data.annotations.cif_utils import read_mmcif_file

    source = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    cif_file = read_mmcif_file(source)
    block = list(cif_file.values())[0]
    for category_name in list(block):
        if category_name != "atom_site":
            del block[category_name]
    del block["atom_site"]["pdbx_PDB_model_num"]
    del block["atom_site"]["pdbx_PDB_ins_code"]
    coordinate_only = tmp_path / "coordinate_only.cif"
    cif_file.write(str(coordinate_only))

    inputs = custom.write_custom_query_files(
        [coordinate_only], work_dir=tmp_path / "work"
    )
    chains = pd.read_parquet(inputs.chain_manifest)

    assert chains["chain_asym_id"].tolist() == ["A", "B"]
    assert set(chains["sequence_source"]) == {"coordinates"}
    output = pdbx.CIFFile.read(str(inputs.chain_cif_dir / "cq00000000.cif"))
    output_block = list(output.values())[0]
    assert set(output_block["entity_poly_seq"]["num"].as_array(int)) == set(
        range(1, int(chains.loc[0, "sequence_length"]) + 1)
    )
    assert np.min(output_block["atom_site"]["label_seq_id"].as_array(int)) == 1


def test_write_custom_query_files_reports_missing_coordinate_fields(
    test_dir, tmp_path
):
    from plinder.data.annotations.cif_utils import read_mmcif_file

    source = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    cif_file = read_mmcif_file(source)
    del list(cif_file.values())[0]["atom_site"]["Cartn_z"]
    invalid = tmp_path / "missing_z.cif"
    cif_file.write(str(invalid))

    with pytest.raises(ValueError, match=r"_atom_site\.Cartn_z"):
        custom.write_custom_query_files([invalid], work_dir=tmp_path / "work")


def test_write_pdb_query_files_reports_missing_assembly_fields(
    test_dir, tmp_path
):
    from plinder.data.annotations.cif_utils import read_mmcif_file

    source = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    cif_file = read_mmcif_file(source)
    block = list(cif_file.values())[0]
    del block["pdbx_struct_assembly_gen"]
    invalid = tmp_path / "missing_assembly.cif"
    cif_file.write(str(invalid))

    with pytest.raises(ValueError, match="_pdbx_struct_assembly_gen"):
        custom.write_custom_query_files(
            [invalid],
            work_dir=tmp_path / "work",
            structure_mode="pdb",
        )


def test_write_pdb_query_files_deduplicates_assembly_copies(test_dir, tmp_path):
    cif = test_dir / "interfaces/cm/pdb_00007cm8/pdb_00007cm8_xyz-enrich.cif.gz"

    inputs = custom.write_custom_query_files(
        [cif],
        work_dir=tmp_path,
        structure_mode="pdb",
        assembly_ids=["1"],
    )
    chains = pd.read_parquet(inputs.chain_manifest)

    assert chains["chain_asym_id"].tolist() == ["A"]
    assert chains["query_id"].tolist() == ["cq00000000"]


def test_write_custom_query_files_derives_coordinate_only_sequences(
    test_dir, tmp_path, monkeypatch
):
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    monkeypatch.setattr(custom, "_protein_asym_sequences", lambda _block: {})

    inputs = custom.write_custom_query_files([cif], work_dir=tmp_path)
    chains = pd.read_parquet(inputs.chain_manifest)

    assert chains["chain_asym_id"].tolist() == ["A", "B"]
    assert chains["sequence_length"].min() >= 12
    assert set(chains["sequence_source"]) == {"coordinates"}


def test_annotate_custom_cif_files_writes_bonded_ligands(
    test_dir, tmp_path, monkeypatch
):
    import yaml
    from plinder.data.annotations import ligand_utils
    from rdkit import Chem

    monkeypatch.setattr(ligand_utils, "BINDING_AFFINITY", {})

    cif = test_dir / "custom_cif/boltz_8c3u_input_model_0.cif"
    config = yaml.safe_load(
        (test_dir / "custom_cif/boltz_8c3u_input.yaml").read_text()
    )
    smiles = next(
        sequence["ligand"]["smiles"]
        for sequence in config["sequences"]
        if "ligand" in sequence
    )

    annotations = custom.annotate_custom_cif_files(
        [cif],
        work_dir=tmp_path,
        ligand_smiles_dict={"LIG": smiles},
        include_interfaces=False,
    )

    assert set(annotations.entries_by_structure) == {"boltz_8c3u_input_model_0"}
    assert not pd.read_parquet(annotations.annotation_table).empty
    assert not pd.read_parquet(annotations.entry_chains).empty
    ligand_files = sorted(annotations.ligand_sdf_root.rglob("*.sdf"))
    assert ligand_files
    molecules = [
        Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
        for path in ligand_files
    ]
    assert all(molecule is not None for molecule in molecules)
    assert any(
        bond.GetBondType() != Chem.BondType.SINGLE
        for molecule in molecules
        if molecule is not None
        for bond in molecule.GetBonds()
    )


def test_annotate_custom_cif_files_keeps_interface_only_entries(
    test_dir, tmp_path
):
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"

    annotations = custom.annotate_custom_cif_files(
        [cif],
        work_dir=tmp_path,
        include_ligands=False,
        include_interfaces=True,
        interface_annotate_prodigy=False,
    )
    entry = annotations.entries_by_structure["pdb_00007cma_xyz-enrich"]

    assert not entry.systems
    assert set(entry.interfaces) == {
        "pdb_00007cma_xyz-enrich__1__1.A--1.B"
    }
    assert not pd.read_parquet(annotations.interface_annotations).empty


def test_annotate_custom_cif_files_keeps_protein_only_entries(test_dir, tmp_path):
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"

    annotations = custom.annotate_custom_cif_files(
        [cif],
        work_dir=tmp_path,
        include_ligands=False,
        include_interfaces=False,
    )
    entry = annotations.entries_by_structure["pdb_00007cma_xyz-enrich"]

    assert not entry.systems
    assert not entry.interfaces
    assert {"A", "B"} <= set(entry.chains)


def test_create_custom_query_databases_records_backend_identifiers(
    test_dir, tmp_path, monkeypatch
):
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    inputs = custom.write_custom_query_files([cif], work_dir=tmp_path)
    commands: list[list[str]] = []

    monkeypatch.setattr(custom.shutil, "which", lambda backend: f"/bin/{backend}")

    def fake_run(command: list[str]) -> None:
        commands.append(command)
        database = Path(command[3])
        identifiers = (
            ["cq00000000_A", "cq00000001_A"]
            if command[0] == "foldseek"
            else ["cq00000000", "cq00000001"]
        )
        database.with_suffix(".lookup").write_text(
            "".join(
                f"{index}\t{identifier}\t0\n"
                for index, identifier in enumerate(identifiers)
            )
        )

    monkeypatch.setattr(custom, "_run_command", fake_run)
    databases = custom.create_custom_query_databases(
        inputs,
        work_dir=tmp_path,
        threads=3,
    )
    identifier_map = pd.read_parquet(databases.identifier_map)

    assert [command[:2] for command in commands] == [
        ["foldseek", "createdb"],
        ["mmseqs", "createdb"],
    ]
    assert commands[0][commands[0].index("--coord-store-mode") + 1] == "2"
    assert commands[1][commands[1].index("--threads") + 1] == "3"
    assert len(identifier_map) == 4
    assert set(identifier_map["query_id"]) == {"cq00000000", "cq00000001"}


def test_create_custom_query_databases_only_requires_selected_backend(
    test_dir, tmp_path, monkeypatch
):
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    inputs = custom.write_custom_query_files([cif], work_dir=tmp_path)
    executable_checks: list[str] = []
    commands: list[list[str]] = []

    def which(backend: str) -> str:
        executable_checks.append(backend)
        return f"/bin/{backend}"

    def fake_run(command: list[str]) -> None:
        commands.append(command)
        database = Path(command[3])
        database.with_suffix(".lookup").write_text(
            "0\tcq00000000\t0\n1\tcq00000001\t0\n"
        )

    monkeypatch.setattr(custom.shutil, "which", which)
    monkeypatch.setattr(custom, "_run_command", fake_run)

    databases = custom.create_custom_query_databases(
        inputs,
        work_dir=tmp_path,
        backends=["mmseqs"],
    )

    assert executable_checks == ["mmseqs"]
    assert [command[0] for command in commands] == ["mmseqs"]
    assert set(databases.databases) == {"mmseqs"}


def test_map_custom_alignment_hits_maps_query_and_target_chains(tmp_path):
    input_root = tmp_path / "query_inputs"
    input_root.mkdir()
    chain_manifest = input_root / "query_chains.parquet"
    pd.DataFrame(
        {
            "query_id": ["cq00000000"],
            "query_chain_id": ["model__A"],
            "structure_id": ["model"],
            "chain_asym_id": ["A"],
            "sequence_source": ["polymer"],
            "resolved_residue_numbers": [[1, 2]],
        }
    ).to_parquet(chain_manifest, index=False)
    inputs = custom.CustomQueryInputs(
        root=input_root,
        chain_cif_dir=input_root / "chains",
        sequence_fasta=input_root / "query.fasta",
        chain_manifest=chain_manifest,
    )
    database_root = tmp_path / "query_databases"
    database_root.mkdir()
    identifier_map = database_root / "query_identifier_map.parquet"
    pd.DataFrame(
        {
            "backend": ["foldseek"],
            "backend_query_id": ["cq00000000_A"],
            "query_id": ["cq00000000"],
        }
    ).to_parquet(identifier_map, index=False)
    databases = custom.CustomQueryDatabases(
        root=database_root,
        inputs=inputs,
        databases={"foldseek": database_root / "foldseek"},
        identifier_map=identifier_map,
    )
    raw = tmp_path / "raw.parquet"
    pd.DataFrame(
        {
            "query": ["cq00000000_A"],
            "target": ["pdb_00001abc_xyz-enrich.cif.gz_X"],
            "qstart": [1],
            "tstart": [2],
            "qcov": [0.8],
            "fident": [0.5],
            "qaln": ["AC"],
            "taln": ["AC"],
            "lddt": [0.7],
        }
    ).to_parquet(raw, index=False)
    lookup = tmp_path / "alignment_chain_lookup.parquet"
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["B"],
            "chain_auth_id": ["X"],
            "selected_residue_numbers": [[2, 4]],
            "selected_residue_indices": [[1, 3]],
        }
    ).to_parquet(lookup, index=False)
    output = tmp_path / "mapped.parquet"

    custom.map_custom_alignment_hits(
        raw_alignment=raw,
        backend="foldseek",
        query_databases=databases,
        alignment_chain_lookup=lookup,
        output_path=output,
    )
    result = pd.read_parquet(output)

    assert result.loc[0, "query_chain_id"] == "model__A"
    assert result.loc[0, "target_entry"] == "1abc"
    assert result.loc[0, "target_chain_asym_id"] == "B"
    assert result.loc[0, "source"] == "foldseek"
    assert result.loc[0, "query_sequence_source"] == "polymer"
    assert result.loc[0, "query_resolved_residue_numbers"].tolist() == [1, 2]
    assert result.loc[0, "target_selected_residue_numbers"].tolist() == [2, 4]


def test_custom_search_defaults_keep_coverage_and_disable_identity_filter():
    config = custom.CustomProteinSearchConfig()

    assert config.max_seqs == 10_000
    assert config.coverage == 0.0
    assert config.min_seq_id == 0.0


def test_custom_protein_searches_against_local_release(test_dir, tmp_path):
    data_dir_value = os.environ.get("PLINDER_CUSTOM_SEARCH_SMOKE_DATA_DIR")
    if data_dir_value is None:
        pytest.skip("local V3 search databases were not requested")
    data_dir = Path(data_dir_value)
    search_databases: dict[str, custom.SearchDatabaseBundle] = {}
    for backend in custom.SEARCH_BACKENDS:
        root = data_dir / "dbs/subdbs" / f"holo_{backend}"
        manifest = json.loads((root / "exact_cluster.json").read_text())
        search_databases[backend] = custom.SearchDatabaseBundle(
            backend=backend,
            root=root,
            search_target=root / str(manifest["search_target"]),
            conversion_target=root / str(manifest["conversion_target"]),
            cluster_alignments=(
                root / str(manifest["cluster_alignments"])
                if backend == "mmseqs"
                else None
            ),
            manifest=manifest,
        )
    assets = custom.CustomScoringAssets(
        annotation_table=data_dir / "index/annotation_table.parquet",
        entry_chains=data_dir / "index/entry_chains.parquet",
        interface_annotations=data_dir / "index/interface_annotation_table.parquet",
        alignment_chain_lookup=data_dir / "index/alignment_chain_lookup.parquet",
        search_databases=search_databases,
        ligand_archives={},
    )
    cif = test_dir / "interfaces/cm/pdb_00007cma/pdb_00007cma_xyz-enrich.cif.gz"
    inputs = custom.write_custom_query_files([cif], work_dir=tmp_path)
    databases = custom.create_custom_query_databases(
        inputs,
        work_dir=tmp_path,
        threads=4,
    )

    outputs = custom.run_custom_protein_searches(
        query_databases=databases,
        assets=assets,
        output_dir=tmp_path / "hits",
        scratch_dir=tmp_path / "scratch",
        config=custom.CustomProteinSearchConfig(max_seqs=100),
        threads=4,
    )

    assert set(outputs) == set(custom.SEARCH_BACKENDS)
    for output in outputs.values():
        hits = pd.read_parquet(output)
        assert set(hits["query_chain_asym_id"]) == {"A", "B"}
        assert "7cma" in set(hits["target_entry"])


@pytest.mark.parametrize(
    ("backend", "query_number", "target_numbers", "target_indices"),
    [
        ("foldseek", 10, [20], [1]),
        ("mmseqs", 2, [2], [1]),
    ],
)
def test_prepare_custom_score_alignments_maps_selected_residues(
    tmp_path,
    backend,
    query_number,
    target_numbers,
    target_indices,
):
    query_ligand = LigandView(
        id="model__1__1.L",
        pdb_id="model",
        system_id="model__1__1.A__1.L",
        instance_chain="1.L",
        asym_id="L",
        is_proper=True,
        protein_chains_asym_id=["1.A"],
        num_pocket_residues=1,
        num_interactions=0,
        num_unique_interactions=0,
        pocket_residue_number_to_index={"1.A": {query_number: 1}},
    )
    query_system = SystemView(
        id=query_ligand.system_id,
        pdb_id="model",
        system_type="holo",
        protein_chains_asym_id=["1.A"],
        proper_num_pocket_residues=1,
        proper_num_interactions=0,
        proper_num_unique_interactions=0,
        pocket_residue_number_to_index={"1.A": {query_number: 1}},
        ligands={"1.L": query_ligand},
    )
    query_entry = EntryView(
        pdb_id="model",
        chains={"A": ChainView(asym_id="A", auth_id="A", length=2)},
        systems={query_system.id: query_system},
        author_to_asym={"A": "A"},
    )
    hit_path = tmp_path / f"{backend}_hits.parquet"
    row = {
        "structure_id": "model",
        "query_chain_asym_id": "A",
        "target_entry": "1abc",
        "target_chain_asym_id": "B",
        "target_selected_residue_numbers": target_numbers,
        "target_selected_residue_indices": target_indices,
        "qstart": 1,
        "tstart": 1,
        "qcov": 1.0,
        "fident": 1.0,
        "qaln": "AC",
        "taln": "AC",
    }
    if backend == "foldseek":
        row["lddt"] = 0.9
    else:
        row["query_sequence_source"] = "polymer"
        row["query_resolved_residue_numbers"] = [1, 2]
    pd.DataFrame([row]).to_parquet(hit_path, index=False)

    outputs = custom.prepare_custom_score_alignments(
        {backend: hit_path},
        entries_by_structure={"model": query_entry},
        output_dir=tmp_path / "score_alignments",
    )
    result = pd.read_parquet(outputs[backend])

    assert result.loc[0, "query_entry"] == "model"
    assert result.loc[0, "query_chain_mapped"] == "A"
    assert result.loc[0, "target_chain_mapped"] == "B"
    assert result.loc[0, "query_selected_residue_numbers"].tolist() == [query_number]
    assert result.loc[0, "target_selected_residue_numbers"].tolist() == [
        target_numbers[0]
    ]
    assert result.loc[0, "selected_residue_identity"] == b"\x01"
    assert result.loc[0, "seqsim"] == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("backend", "target_numbers", "target_indices", "custom_number"),
    [
        ("foldseek", [20], [1], 200),
        ("mmseqs", [2], [1], 2),
    ],
)
def test_prepare_custom_protein_score_alignments_reverses_direction(
    tmp_path,
    backend,
    target_numbers,
    target_indices,
    custom_number,
):
    query_entry = EntryView(
        pdb_id="model",
        chains={"A": ChainView(asym_id="A", auth_id="A", length=2)},
        systems={},
        author_to_asym={"A": "A"},
    )
    hit_path = tmp_path / f"{backend}_hits.parquet"
    row = {
        "structure_id": "model",
        "query_chain_asym_id": "A",
        "query_sequence_source": "polymer",
        "query_resolved_residue_numbers": [100, 200],
        "target_entry": "1abc",
        "target_chain_asym_id": "B",
        "target_selected_residue_numbers": target_numbers,
        "target_selected_residue_indices": target_indices,
        "qstart": 1,
        "tstart": 1,
        "qcov": 0.5,
        "tcov": 0.75,
        "fident": 1.0,
        "qaln": "AC",
        "taln": "AC",
    }
    if backend == "foldseek":
        row["lddt"] = 0.9
    pd.DataFrame([row]).to_parquet(hit_path, index=False)

    outputs = custom.prepare_custom_protein_score_alignments(
        {backend: hit_path},
        entries_by_structure={"model": query_entry},
        output_dir=tmp_path / "protein_score_alignments",
    )
    result = pd.read_parquet(outputs[backend])

    assert result.loc[0, "query_entry"] == "1abc"
    assert result.loc[0, "target_entry"] == "model"
    assert result.loc[0, "query_chain_mapped"] == "B"
    assert result.loc[0, "target_chain_mapped"] == "A"
    assert result.loc[0, "qcov"] == pytest.approx(0.75)
    assert result.loc[0, "tcov"] == pytest.approx(0.5)
    assert result.loc[0, "query_selected_residue_numbers"].tolist() == [
        target_numbers[0]
    ]
    assert result.loc[0, "target_selected_residue_numbers"].tolist() == [
        custom_number
    ]
    assert result.loc[0, "selected_residue_identity"] == b"\x01"
    assert result.loc[0, "fident_qcov"] == pytest.approx(0.75)


def test_prepare_custom_score_alignments_maps_coordinate_fasta_positions(tmp_path):
    query_ligand = LigandView(
        id="model__1__1.L",
        pdb_id="model",
        system_id="model__1__1.A__1.L",
        instance_chain="1.L",
        asym_id="L",
        is_proper=True,
        protein_chains_asym_id=["1.A"],
        num_pocket_residues=1,
        num_interactions=0,
        num_unique_interactions=0,
        pocket_residue_number_to_index={"1.A": {40: 2}},
    )
    query_system = SystemView(
        id=query_ligand.system_id,
        pdb_id="model",
        system_type="holo",
        protein_chains_asym_id=["1.A"],
        proper_num_pocket_residues=1,
        proper_num_interactions=0,
        proper_num_unique_interactions=0,
        pocket_residue_number_to_index={"1.A": {40: 2}},
        ligands={"1.L": query_ligand},
    )
    query_entry = EntryView(
        pdb_id="model",
        chains={"A": ChainView(asym_id="A", auth_id="A", length=3)},
        systems={query_system.id: query_system},
        author_to_asym={"A": "A"},
    )
    hit_path = tmp_path / "mmseqs_hits.parquet"
    pd.DataFrame(
        [
            {
                "structure_id": "model",
                "query_chain_asym_id": "A",
                "query_sequence_source": "coordinates",
                "query_resolved_residue_numbers": [10, 20, 40],
                "target_entry": "1abc",
                "target_chain_asym_id": "B",
                "target_selected_residue_numbers": [3],
                "target_selected_residue_indices": [2],
                "qstart": 1,
                "tstart": 1,
                "qcov": 1.0,
                "fident": 1.0,
                "qaln": "ACD",
                "taln": "ACD",
            }
        ]
    ).to_parquet(hit_path, index=False)

    outputs = custom.prepare_custom_score_alignments(
        {"mmseqs": hit_path},
        entries_by_structure={"model": query_entry},
        output_dir=tmp_path / "score_alignments",
    )
    result = pd.read_parquet(outputs["mmseqs"])

    assert result.loc[0, "query_selected_residue_numbers"].tolist() == [40]
    assert result.loc[0, "target_selected_residue_numbers"].tolist() == [3]
    assert result.loc[0, "selected_residue_identity"] == b"\x01"


def _scoring_entry(
    *,
    pdb_id: str,
    chain_id: str,
    ligand_id: str,
    ligand_chain: str,
    pocket_number: int,
) -> EntryView:
    system_id = f"{pdb_id}__1__1.{chain_id}__1.{ligand_chain}"
    ligand = LigandView(
        id=ligand_id,
        pdb_id=pdb_id,
        system_id=system_id,
        instance_chain=f"1.{ligand_chain}",
        asym_id=ligand_chain,
        is_proper=True,
        protein_chains_asym_id=[f"1.{chain_id}"],
        num_pocket_residues=1,
        num_interactions=0,
        num_unique_interactions=0,
        pocket_residue_number_to_index={f"1.{chain_id}": {pocket_number: 1}},
    )
    system = SystemView(
        id=system_id,
        pdb_id=pdb_id,
        system_type="holo",
        protein_chains_asym_id=[f"1.{chain_id}"],
        proper_num_pocket_residues=1,
        proper_num_interactions=0,
        proper_num_unique_interactions=0,
        pocket_residue_number_to_index={f"1.{chain_id}": {pocket_number: 1}},
        ligands={ligand.instance_chain: ligand},
    )
    return EntryView(
        pdb_id=pdb_id,
        chains={
            chain_id: ChainView(asym_id=chain_id, auth_id=chain_id, length=2)
        },
        systems={system_id: system},
        author_to_asym={chain_id: chain_id},
    )


def _custom_annotations(tmp_path, entry):
    root = tmp_path / "custom_annotations"
    ligand_root = root / "ligands"
    ligand_root.mkdir(parents=True)
    return custom.CustomStructureAnnotations(
        root=root,
        ligand_sdf_root=ligand_root,
        annotation_table=root / "annotation.parquet",
        entry_chains=root / "chains.parquet",
        interface_annotations=root / "interfaces.parquet",
        entries_by_structure={"model": entry},
    )


def _custom_assets(tmp_path):
    return custom.CustomScoringAssets(
        annotation_table=tmp_path / "release_annotation.parquet",
        entry_chains=tmp_path / "release_chains.parquet",
        interface_annotations=tmp_path / "release_interfaces.parquet",
        alignment_chain_lookup=tmp_path / "release_lookup.parquet",
        search_databases={},
        ligand_archives={},
    )


def test_calculate_custom_similarity_scores_reuses_release_metrics(
    tmp_path, monkeypatch
):
    query = _scoring_entry(
        pdb_id="model",
        chain_id="A",
        ligand_id="model__1__1.L",
        ligand_chain="L",
        pocket_number=2,
    )
    target = _scoring_entry(
        pdb_id="1abc",
        chain_id="B",
        ligand_id="1abc__1__1.Z",
        ligand_chain="Z",
        pocket_number=20,
    )
    monkeypatch.setattr(
        custom,
        "_load_release_entry_views",
        lambda _assets, *, pdb_ids: {"1abc": target},
    )
    alignment = tmp_path / "foldseek.parquet"
    pd.DataFrame(
        [
            {
                "query_entry": "model",
                "target_entry": "1abc",
                "query_chain_mapped": "A",
                "target_chain_mapped": "B",
                "source": "foldseek",
                "qcov": 1.0,
                "fident": 1.0,
                "seqsim": 1.0,
                "lddt": 0.9,
                "query_selected_residue_numbers": [2],
                "target_selected_residue_numbers": [20],
                "selected_residue_identity": b"\x01",
            }
        ]
    ).to_parquet(alignment, index=False)

    scores = custom.calculate_custom_similarity_scores(
        {"foldseek": alignment},
        annotations=_custom_annotations(tmp_path, query),
        assets=_custom_assets(tmp_path),
        work_dir=tmp_path / "score_work",
        include_shape=False,
    )

    pocket = scores.loc[scores["metric"] == "pocket_qcov"]
    assert len(pocket) == 1
    assert pocket.iloc[0]["similarity"] == 100
    assert pocket.iloc[0]["query_ligand_id"] == "model__1__1.L"
    assert pocket.iloc[0]["target_ligand_id"] == "1abc__1__1.Z"


def test_calculate_custom_protein_scores_uses_plinder_pocket(
    tmp_path, monkeypatch
):
    plinder_entry = _scoring_entry(
        pdb_id="1abc",
        chain_id="B",
        ligand_id="1abc__1__1.Z",
        ligand_chain="Z",
        pocket_number=20,
    )
    monkeypatch.setattr(
        custom,
        "_load_release_entry_views",
        lambda _assets, *, pdb_ids: {"1abc": plinder_entry},
    )
    alignment = tmp_path / "reverse_foldseek.parquet"
    pd.DataFrame(
        [
            {
                "query_entry": "1abc",
                "target_entry": "model",
                "query_chain_mapped": "B",
                "target_chain_mapped": "A",
                "source": "foldseek",
                "qcov": 1.0,
                "fident": 1.0,
                "seqsim": 1.0,
                "lddt": 0.9,
                "query_selected_residue_numbers": [20],
                "target_selected_residue_numbers": [-1],
                "selected_residue_identity": b"\x01",
            }
        ]
    ).to_parquet(alignment, index=False)

    scores = custom.calculate_custom_protein_similarity_scores(
        {"foldseek": alignment},
        assets=_custom_assets(tmp_path),
        work_dir=tmp_path / "protein_score_work",
        custom_chain_ids={"model_A"},
    )

    pocket = scores.loc[scores["metric"].astype(str) == "pocket_fident"]
    assert len(pocket) == 1
    assert pocket.iloc[0]["similarity"] == 100
    assert pocket.iloc[0]["query_system"] == "1abc__1__1.B__1.Z"
    assert pocket.iloc[0]["query_ligand_id"] == "1abc__1__1.Z"
    assert pocket.iloc[0]["target_system"] == "model_A"
    assert pd.isna(pocket.iloc[0]["target_ligand_id"])


def test_calculate_custom_protein_scores_reads_alignments_once(
    tmp_path, monkeypatch
):
    from plinder.data.annotations import get_similarity_scores

    release_entries = {
        pdb_id: _scoring_entry(
            pdb_id=pdb_id,
            chain_id="B",
            ligand_id=f"{pdb_id}__1__1.Z",
            ligand_chain="Z",
            pocket_number=20,
        )
        for pdb_id in ["1abc", "2def"]
    }
    monkeypatch.setattr(
        custom,
        "_load_release_entry_views",
        lambda _assets, *, pdb_ids: {
            pdb_id: release_entries[pdb_id] for pdb_id in pdb_ids
        },
    )
    alignment = tmp_path / "reverse_foldseek.parquet"
    pd.DataFrame(
        [
            {
                "query_entry": pdb_id,
                "target_entry": "model_with_underscore",
                "query_chain_mapped": "B",
                "target_chain_mapped": "A",
                "source": "foldseek",
                "qcov": 1.0,
                "fident": 1.0,
                "seqsim": 1.0,
                "lddt": 0.9,
                "query_selected_residue_numbers": [20],
                "target_selected_residue_numbers": [-1],
                "selected_residue_identity": b"\x01",
            }
            for pdb_id in release_entries
        ]
    ).to_parquet(alignment, index=False)
    original_load = get_similarity_scores.Scorer.load_alignments
    load_calls = []

    def counted_load(self, *args, **kwargs):
        load_calls.append((args, kwargs))
        return original_load(self, *args, **kwargs)

    monkeypatch.setattr(
        get_similarity_scores.Scorer,
        "load_alignments",
        counted_load,
    )

    scores = custom.calculate_custom_protein_similarity_scores(
        {"foldseek": alignment},
        assets=_custom_assets(tmp_path),
        work_dir=tmp_path / "protein_score_work",
        custom_chain_ids={"model_with_underscore_A"},
    )

    assert len(load_calls) == 1
    assert set(scores["query_system"]) == {
        "1abc__1__1.B__1.Z",
        "2def__1__1.B__1.Z",
    }
    assert set(scores["target_system"]) == {"model_with_underscore_A"}


def test_write_custom_aligned_pocket_residues(tmp_path, monkeypatch):
    plinder_entry = _scoring_entry(
        pdb_id="1abc",
        chain_id="B",
        ligand_id="1abc__1__1.Z",
        ligand_chain="Z",
        pocket_number=20,
    )
    monkeypatch.setattr(
        custom,
        "_load_release_entry_views",
        lambda _assets, *, pdb_ids: {"1abc": plinder_entry},
    )
    alignment = tmp_path / "reverse_foldseek.parquet"
    pd.DataFrame(
        [
            {
                "query_entry": "1abc",
                "target_entry": "model_with_underscore",
                "query_chain_mapped": "B",
                "target_chain_mapped": "A",
                "source": "foldseek",
                "query_selected_residue_numbers": [20],
                "target_selected_residue_numbers": [42],
                "selected_residue_identity": b"\x01",
            }
        ]
    ).to_parquet(alignment, index=False)
    protein_scores = tmp_path / "protein_scores.parquet"
    pd.DataFrame(
        [
            {
                "query_system": "1abc__1__1.B__1.Z",
                "query_ligand_id": "1abc__1__1.Z",
                "target_system": "model_with_underscore_A",
                "metric": "pocket_fident",
            }
        ]
    ).to_parquet(protein_scores, index=False)

    output = custom.write_custom_aligned_pocket_residues(
        {"foldseek": alignment},
        protein_scores=protein_scores,
        assets=_custom_assets(tmp_path),
        output_path=tmp_path / "aligned_pocket_residues.parquet",
    )
    result = pd.read_parquet(output)

    assert result.to_dict("records") == [
        {
            "plinder_system_id": "1abc__1__1.B__1.Z",
            "plinder_ligand_id": "1abc__1__1.Z",
            "plinder_entry_id": "1abc",
            "plinder_chain_instance": "1.B",
            "plinder_chain_asym_id": "B",
            "plinder_residue_number": 20,
            "custom_structure_id": "model_with_underscore",
            "custom_chain_asym_id": "A",
            "custom_residue_number": 42,
            "residue_identical": True,
            "source": "foldseek",
        }
    ]


def test_calculate_custom_interface_scores_uses_compact_maps(
    tmp_path, monkeypatch
):
    query_interface = InterfaceView(
        id="model__1__1.A_1.C",
        pdb_id="model",
        biounit_id="1",
        chain_1="1.A",
        chain_2="1.C",
        chain_1_residue_number_to_index={2: 1},
        chain_2_residue_number_to_index={3: 2},
        num_contact_residue_pairs=1,
    )
    target_interface = InterfaceView(
        id="1abc__1__1.B_1.D",
        pdb_id="1abc",
        biounit_id="1",
        chain_1="1.B",
        chain_2="1.D",
        chain_1_residue_number_to_index={20: 1},
        chain_2_residue_number_to_index={30: 2},
        num_contact_residue_pairs=1,
    )
    query = EntryView(
        pdb_id="model",
        chains={
            chain: ChainView(asym_id=chain, auth_id=chain, length=3)
            for chain in ["A", "C"]
        },
        systems={},
        author_to_asym={"A": "A", "C": "C"},
        interfaces={query_interface.id: query_interface},
    )
    target = EntryView(
        pdb_id="1abc",
        chains={
            chain: ChainView(asym_id=chain, auth_id=chain, length=3)
            for chain in ["B", "D"]
        },
        systems={},
        author_to_asym={"B": "B", "D": "D"},
        interfaces={target_interface.id: target_interface},
    )
    monkeypatch.setattr(
        custom,
        "_load_release_entry_views",
        lambda _assets, *, pdb_ids: {"1abc": target},
    )
    alignment = tmp_path / "interface_foldseek.parquet"
    pd.DataFrame(
        [
            {
                "query_entry": "model",
                "target_entry": "1abc",
                "query_chain_mapped": query_chain,
                "target_chain_mapped": target_chain,
                "source": "foldseek",
                "query_selected_residue_numbers": [query_number],
                "target_selected_residue_numbers": [target_number],
            }
            for query_chain, target_chain, query_number, target_number in [
                ("A", "B", 2, 20),
                ("C", "D", 3, 30),
            ]
        ]
    ).to_parquet(alignment, index=False)

    scores = custom.calculate_custom_interface_similarity_scores(
        {"foldseek": alignment},
        annotations=_custom_annotations(tmp_path, query),
        assets=_custom_assets(tmp_path),
    )

    assert len(scores) == 1
    assert scores.loc[0, "similarity"] == 100
    assert scores.loc[0, "query_system"] == query_interface.id
    assert scores.loc[0, "target_system"] == target_interface.id
