from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

from plinder.core.utils.schemas import PROTEIN_SIMILARITY_EXPORT_SCHEMA
from plinder.data import monomer_similarity
from plinder.data.monomer_similarity import (
    _chain_ids,
    _score_alignment_batches,
    finalize_monomer_similarity_scores,
)


@pytest.mark.parametrize("backend", ["mmseqs", "foldseek"])
def test_monomer_pair_scores_map_release_chains(tmp_path: Path, backend: str) -> None:
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "2def"],
            "chain_auth_id": ["R", "S"],
            "chain_asym_id": ["A", "C"],
        }
    ).to_parquet(index / "alignment_chain_lookup.parquet", index=False)
    query = "1abc_R" if backend == "mmseqs" else "pdb_00001abc_xyz-enrich_R"
    target = "2def_S" if backend == "mmseqs" else "pdb_00002def_xyz-enrich_MODEL_1_S"
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.DataFrame(
        {
            "query": [query],
            "target": [target],
            "qcov": [0.8],
            "tcov": [0.6],
            "fident": [0.75],
            "qaln": ["ACDE"],
            "taln": ["ACDF"],
            **({"lddt": [0.7]} if backend == "foldseek" else {}),
        }
    ).to_parquet(raw / "batch.parquet", index=False)
    output = tmp_path / "scores.parquet"
    assert (
        _score_alignment_batches(
            raw, output, backend=backend, chain_ids=_chain_ids(tmp_path, backend)
        )
        == 1
    )
    assert pq.read_schema(output).equals(PROTEIN_SIMILARITY_EXPORT_SCHEMA)
    row = pd.read_parquet(output).iloc[0]
    assert (row.query_entry, row.query_chain_mapped) == ("1abc", "A")
    assert (row.target_entry, row.target_chain_mapped) == ("2def", "C")
    assert (row.qcov, row.tcov, row.fident) == (80, 60, 75)
    assert row.lddt == 70 if backend == "foldseek" else pd.isna(row.lddt)


def test_repeated_aligned_sequences_are_scored_once(
    tmp_path: Path, monkeypatch
) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    alignments = pd.DataFrame(
        {
            "query": ["1abc_R", "1abc_R"],
            "target": ["2def_S", "2def_S"],
            "qcov": [0.8, 0.8],
            "tcov": [0.6, 0.6],
            "fident": [0.75, 0.75],
            "qaln": ["ACDE", "ACDE"],
            "taln": ["ACDF", "ACDF"],
        }
    )
    for index in range(len(alignments)):
        alignments.iloc[[index]].to_parquet(raw / f"batch-{index}.parquet", index=False)
    calls = []
    monkeypatch.setattr(
        monomer_similarity,
        "get_sequence_similarity_helper",
        lambda query, target: calls.append((query, target)) or 0.5,
    )
    assert (
        _score_alignment_batches(
            raw,
            tmp_path / "scores.parquet",
            backend="mmseqs",
            chain_ids={"1abc_R": ("1abc", "A"), "2def_S": ("2def", "C")},
        )
        == 2
    )
    assert calls == [("ACDE", "ACDF")]
    assert pq.ParquetFile(tmp_path / "scores.parquet").num_row_groups == 1


def test_monomer_score_writer_accepts_empty_alignments(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.DataFrame(
        {
            column: pd.Series(dtype=dtype)
            for column, dtype in {
                "query": "string",
                "target": "string",
                "qcov": "float64",
                "tcov": "float64",
                "fident": "float64",
                "qaln": "string",
                "taln": "string",
            }.items()
        }
    ).to_parquet(raw / "empty.parquet", index=False)
    output = tmp_path / "scores.parquet"
    assert _score_alignment_batches(raw, output, backend="mmseqs", chain_ids={}) == 0
    assert pq.read_schema(output).equals(PROTEIN_SIMILARITY_EXPORT_SCHEMA)


def test_foldseek_missing_lddt_remains_null(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.DataFrame(
        {
            "query": ["pdb_00001abc_xyz-enrich_R"],
            "target": ["pdb_00002def_xyz-enrich_S"],
            "qcov": [0.8],
            "tcov": [0.6],
            "fident": [0.75],
            "qaln": ["ACDE"],
            "taln": ["ACDF"],
            "lddt": [None],
        }
    ).to_parquet(raw / "batch.parquet", index=False)
    output = tmp_path / "scores.parquet"
    assert (
        _score_alignment_batches(
            raw,
            output,
            backend="foldseek",
            chain_ids={
                "pdb_00001abc_xyz-enrich_R": ("1abc", "A"),
                "pdb_00002def_xyz-enrich_S": ("2def", "C"),
            },
        )
        == 1
    )
    assert pd.isna(pd.read_parquet(output).iloc[0].lddt)


def test_finalize_monomer_scores_has_only_alignment_partition(tmp_path: Path) -> None:
    staging = tmp_path / "exports/.monomer_similarity_scores_staging"
    for backend in ("mmseqs", "foldseek"):
        for target_kind in ("holo", "monomer"):
            source = (
                staging / f"alignment_type={backend}" / f"target_kind={target_kind}"
            )
            source.mkdir(parents=True)
            pq.write_table(
                PROTEIN_SIMILARITY_EXPORT_SCHEMA.empty_table(),
                source / "part-0.parquet",
            )

    assert finalize_monomer_similarity_scores(tmp_path, query_batch_count=1) == 0
    installed = tmp_path / "exports/monomer_similarity_scores"
    assert len(list(installed.rglob("*.parquet"))) == 4
    assert not list(installed.rglob("target_kind=*"))
    assert set(pd.read_parquet(installed).columns) == {
        *PROTEIN_SIMILARITY_EXPORT_SCHEMA.names,
        "alignment_type",
    }


def test_finalize_monomer_scores_requires_every_batch(tmp_path: Path) -> None:
    staging = tmp_path / "exports/.monomer_similarity_scores_staging"
    staging.mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        finalize_monomer_similarity_scores(tmp_path, query_batch_count=1)
    assert staging.is_dir()
    assert not (tmp_path / "exports/monomer_similarity_scores").exists()


def test_finalize_monomer_scores_rejects_extra_or_zero_batches(tmp_path: Path) -> None:
    staging = tmp_path / "exports/.monomer_similarity_scores_staging"
    for backend in ("mmseqs", "foldseek"):
        for target_kind in ("holo", "monomer"):
            source = (
                staging / f"alignment_type={backend}" / f"target_kind={target_kind}"
            )
            source.mkdir(parents=True)
            for index in range(2):
                pq.write_table(
                    PROTEIN_SIMILARITY_EXPORT_SCHEMA.empty_table(),
                    source / f"part-{index}.parquet",
                )

    with pytest.raises(ValueError, match="positive"):
        finalize_monomer_similarity_scores(tmp_path, query_batch_count=0)
    with pytest.raises(ValueError, match="unexpected monomer score batches"):
        finalize_monomer_similarity_scores(tmp_path, query_batch_count=1)

    assert staging.is_dir()
    assert not (tmp_path / "exports/monomer_similarity_scores").exists()
