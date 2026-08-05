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
