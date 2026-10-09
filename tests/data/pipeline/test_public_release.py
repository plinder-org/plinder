from pathlib import Path

import pytest

from plinder.core.release import RELEASE_PATHS
from plinder.data.pipeline.public_release import (
    _FILE_ARTIFACTS,
    _PARQUET_DIRECTORIES,
    _SAMPLING_DIRECTORIES,
    _SEARCH_DATABASES,
    prepare_public_release,
)


def _working_release(root: Path) -> None:
    for name in _FILE_ARTIFACTS:
        path = root / RELEASE_PATHS[name]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"table")
    for name in _PARQUET_DIRECTORIES:
        path = root / RELEASE_PATHS.get(name, name) / "part.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"scores")
        (path.parent / "local.json").write_text("/scicore/private/path")
    for name in _SAMPLING_DIRECTORIES:
        directory = root / RELEASE_PATHS.get(name, name)
        cover = directory / "directed_set_cover/metric=pocket_qcov/threshold=50.parquet"
        cover.parent.mkdir(parents=True, exist_ok=True)
        cover.write_bytes(b"cover")
        cache = directory / "directed_set_cover/reductions/part.parquet"
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(b"cache")
    for name in _SEARCH_DATABASES:
        folder = root / "search_databases" / name
        folder.mkdir(parents=True)
        (folder / "db").write_bytes(b"database")
        (folder / "alias").symlink_to("db")
        (folder / "local.json").write_text("/scicore/private/path")
    (root / "search_databases" / "holo_steam").mkdir()
    (root / "search_databases" / "holo_steam" / "db").write_bytes(b"tea")
    (root / "scores").mkdir()
    (root / "scores" / "temporary.parquet").write_bytes(b"work")


def test_public_release_contains_only_reader_artifacts(tmp_path: Path) -> None:
    working = tmp_path / "working"
    public = tmp_path / "public"
    _working_release(working)

    count = prepare_public_release(working, public)

    assert count > len(_FILE_ARTIFACTS)
    assert (public / RELEASE_PATHS["annotation_table"]).read_bytes() == b"table"
    assert (public / "search_databases/holo_foldseek/alias").read_bytes() == b"database"
    assert (public / "search_databases/holo_foldseek/alias").is_symlink()
    assert (public / "search_databases/holo_foldseek/db").stat().st_ino == (
        working / "search_databases/holo_foldseek/db"
    ).stat().st_ino
    assert not list(public.rglob("*.json"))
    assert not (public / "search_databases/holo_steam").exists()
    assert not (public / "scores").exists()
    for name in _SAMPLING_DIRECTORIES:
        directory = public / RELEASE_PATHS.get(name, name)
        assert (
            directory / "directed_set_cover/metric=pocket_qcov/threshold=50.parquet"
        ).read_bytes() == b"cover"
        assert not (directory / "directed_set_cover/reductions").exists()


def test_public_release_allows_an_empty_monomer_universe(tmp_path: Path) -> None:
    working = tmp_path / "working"
    _working_release(working)
    for name in ("monomer_mmseqs", "monomer_foldseek"):
        for path in sorted((working / "search_databases" / name).iterdir()):
            path.unlink()
        (working / "search_databases" / name).rmdir()
    monomer_scores = working / RELEASE_PATHS["monomer_similarity_scores"]
    for path in monomer_scores.iterdir():
        path.unlink()
    monomer_scores.rmdir()

    prepare_public_release(working, tmp_path / "public")

    assert (tmp_path / "public/search_databases/holo_mmseqs/db").is_file()
    assert not (tmp_path / "public/search_databases/monomer_mmseqs").exists()


def test_public_release_replacement_drops_old_files(tmp_path: Path) -> None:
    working = tmp_path / "working"
    public = tmp_path / "public"
    _working_release(working)
    prepare_public_release(working, public)
    (public / "stale.parquet").write_bytes(b"old")

    with pytest.raises(FileExistsError):
        prepare_public_release(working, public)
    prepare_public_release(working, public, replace=True, copy=True)

    assert not (public / "stale.parquet").exists()
    assert (public / "search_databases/holo_foldseek/db").stat().st_ino != (
        working / "search_databases/holo_foldseek/db"
    ).stat().st_ino


def test_missing_public_file_does_not_replace_release(tmp_path: Path) -> None:
    working = tmp_path / "working"
    public = tmp_path / "public"
    _working_release(working)
    prepare_public_release(working, public)
    (working / RELEASE_PATHS["annotation_table"]).unlink()

    with pytest.raises(FileNotFoundError, match="annotation_table.parquet"):
        prepare_public_release(working, public, replace=True)

    assert (public / RELEASE_PATHS["annotation_table"]).is_file()


def test_public_release_rejects_external_database_links(tmp_path: Path) -> None:
    working = tmp_path / "working"
    public = tmp_path / "public"
    _working_release(working)
    alias = working / "search_databases/holo_foldseek/alias"
    alias.unlink()
    alias.symlink_to("../../scores/temporary.parquet")

    with pytest.raises(ValueError, match="link leaves the public release"):
        prepare_public_release(working, public)

    assert not public.exists()


def test_public_release_rejects_absolute_internal_database_links(
    tmp_path: Path,
) -> None:
    working = tmp_path / "working"
    public = tmp_path / "public"
    _working_release(working)
    alias = working / "search_databases/holo_foldseek/alias"
    alias.unlink()
    alias.symlink_to(working / "search_databases/holo_foldseek/db")

    with pytest.raises(ValueError, match="absolute link"):
        prepare_public_release(working, public, copy=True)

    assert not public.exists()
