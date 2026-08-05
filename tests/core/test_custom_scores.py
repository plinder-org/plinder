# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
from pathlib import Path

import pytest
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
