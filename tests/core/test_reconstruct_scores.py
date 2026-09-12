# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from pathlib import Path

import pandas as pd
import pytest
from plinder.core import release
from plinder.core.scores import reconstruct
from plinder.core.scores.entries import LigandView


def test_prefetch_resolves_only_requested_alignment_shards(tmp_path, monkeypatch):
    requested: list[str] = []

    def get_plinder_path(*, rel: str) -> Path:
        requested.append(rel)
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return path

    monkeypatch.setattr(release.cpl, "get_plinder_path", get_plinder_path)

    paths = reconstruct.prefetch_similarity_alignments(
        [
            "1abc__1__1.A__1.L",
            "2abd__1__1.B__1.M",
            "3xyz__1__1.C__1.N",
        ]
    )

    assert requested == [
        "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet",
        "alignments/search_db=holo/alignment_type=mmseqs/shard=ab.parquet",
        "alignments/search_db=holo/alignment_type=foldseek/shard=xy.parquet",
        "alignments/search_db=holo/alignment_type=mmseqs/shard=xy.parquet",
    ]
    assert paths["1abc"] == paths["2abd"]
    assert paths["1abc"] != paths["3xyz"]
    assert all("scores" not in path for path in requested)


def test_prefetch_requires_requested_shard_in_offline_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(release.cpl, "is_offline", lambda: True)

    with pytest.raises(FileNotFoundError, match="offline cache"):
        reconstruct.prefetch_similarity_alignments(
            ["1abc__1__1.A__1.L"], data_dir=tmp_path
        )


def test_prefetch_accepts_one_available_alignment_backend(tmp_path):
    foldseek = (
        tmp_path / "alignments/search_db=holo/alignment_type=foldseek/shard=ab.parquet"
    )
    foldseek.parent.mkdir(parents=True)
    foldseek.touch()

    paths = reconstruct.prefetch_similarity_alignments(
        ["1abc__1__1.A__1.L"], data_dir=tmp_path
    )

    assert paths == {"1abc": {"foldseek": foldseek}}


def test_canonical_ligand_resolver_materializes_only_requested_member(tmp_path):
    archive = tmp_path / "ligand_archives" / "ab.parquet"
    archive.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "pdb_id": ["1abc", "2abc"],
            "ligand_asym_id": ["L", "M"],
            "sdf": [b"requested", b"other"],
        }
    ).to_parquet(archive, index=False)

    ligand = LigandView(
        id="1abc__1__1.A__1.L__1.L",
        pdb_id="1abc",
        system_id="1abc__1__1.A__1.L",
        instance_chain="1.L",
        asym_id="L",
        is_proper=True,
        protein_chains_asym_id=["1.A"],
        num_pocket_residues=1,
        num_interactions=1,
        num_unique_interactions=1,
    )
    resolver = reconstruct._canonical_ligand_resolver(
        release_root=tmp_path,
        data_dir=tmp_path,
    )

    extracted = resolver(ligand)

    assert extracted == tmp_path / "ligand_archives/1abc/ligand_files/L.sdf"
    assert extracted.read_text() == "requested"
    assert not (tmp_path / "ligand_archives/2abc").exists()
