# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path

import pytest
from plinder.core.release import RELEASE_PATHS, RELEASE_TABLES, PlinderRelease


def test_release_paths_are_relative_and_unique():
    paths = list(RELEASE_PATHS.values())
    assert len(paths) == len(set(paths))
    for path in paths:
        assert not Path(path).is_absolute()
        assert ".." not in Path(path).parts


def test_release_tables_have_artifact_row_description_and_key():
    assert RELEASE_TABLES["annotation"]["row_description"] == "ligand"
    assert RELEASE_TABLES["annotation"]["primary_key"] == ("ligand_id",)
    for table in RELEASE_TABLES.values():
        assert table["artifact"] in RELEASE_PATHS
        assert table["row_description"]
        assert table["primary_key"]


def test_ligand_mmp_pairs_are_a_release_table():
    assert RELEASE_PATHS["ligand_mmp_pairs"] == "index/ligand_mmp_pairs.parquet"
    assert RELEASE_TABLES["ligand_mmp_pairs"]["primary_key"] == (
        "ligand_smiles_id_1",
        "ligand_smiles_id_2",
        "transformation",
        "shared_core_smiles",
    )


def test_system_validation_is_a_release_table():
    assert RELEASE_PATHS["system_validation"] == "index/system_validation.parquet"
    assert RELEASE_TABLES["system_validation"]["primary_key"] == ("system_id",)


def test_path_does_not_require_artifact_to_exist(tmp_path):
    release = PlinderRelease(tmp_path)
    assert release.path("annotation_table") == (
        tmp_path / "index" / "annotation_table.parquet"
    )


def test_parameterized_artifact_path(tmp_path):
    release = PlinderRelease(tmp_path)
    assert release.path(
        "alignment_shard",
        search_db="holo",
        alignment_type="foldseek",
        shard="ab",
    ) == (
        tmp_path
        / "alignments"
        / "search_db=holo"
        / "alignment_type=foldseek"
        / "shard=ab.parquet"
    )


@pytest.mark.parametrize(
    "parameters",
    [
        {},
        {"search_db": "holo", "alignment_type": "foldseek"},
        {
            "search_db": "../holo",
            "alignment_type": "foldseek",
            "shard": "ab",
        },
    ],
)
def test_parameterized_artifact_rejects_missing_or_unsafe_values(tmp_path, parameters):
    with pytest.raises(ValueError):
        PlinderRelease(tmp_path).path("alignment_shard", **parameters)


def test_fetch_checks_explicit_local_artifact(tmp_path):
    annotation = tmp_path / "index" / "annotation_table.parquet"
    annotation.parent.mkdir()
    annotation.touch()
    release = PlinderRelease(tmp_path)
    assert release.fetch("annotation_table") == annotation

    with pytest.raises(FileNotFoundError, match="entry_chains"):
        release.fetch("entry_chains")


def test_fetch_accepts_directory(tmp_path):
    ligand_scores = tmp_path / "ligand_scores"
    ligand_scores.mkdir()
    assert PlinderRelease(tmp_path).fetch("ligand_scores") == ligand_scores


def test_unknown_artifact_and_table_fail_clearly(tmp_path):
    release = PlinderRelease(tmp_path)
    with pytest.raises(KeyError, match="unknown release artifact"):
        release.path("systems")
    with pytest.raises(KeyError, match="unknown release artifact"):
        release.path("splits")
    with pytest.raises(KeyError, match="unknown release table"):
        release.table("systems")


def test_external_source_mmcif_cache_is_not_a_release_artifact(tmp_path):
    with pytest.raises(KeyError, match="unknown release artifact"):
        PlinderRelease(tmp_path).path(
            "source_mmcif", shard="ab", filename="1abc_v1-0.cif.gz"
        )
