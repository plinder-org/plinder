# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import gzip
import os
import shutil
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

GOLDEN_PATH = (
    Path(__file__).resolve().parent.parent
    / "test_data"
    / "scoring_fixture"
    / "golden_scores.parquet"
)
PDB_IDS = ["2gdo", "4qyf", "1atp"]
SORT_KEYS = [
    "query_system",
    "query_ligand_id",
    "target_system",
    "target_ligand_id",
    "metric",
]


@pytest.fixture
def scoring_fixture(
    mock_alternative_datasets, cif_2gdo, cif_4qyf, cif_1atp, test_env, tmp_path
):
    from plinder.data.databases import make_db
    from plinder.data.get_system_annotations import GetPlinderAnnotation

    data_dir = test_env  # = tmp_path / bucket/test/v0 with dbs/ already populated

    cifs = {"2gdo": cif_2gdo, "4qyf": cif_4qyf, "1atp": cif_1atp}
    entries: dict[str, "object"] = {}
    for pdb_id in PDB_IDS:
        entry_save = mock_alternative_datasets(pdb_id)
        anno = GetPlinderAnnotation(cifs[pdb_id], "", save_folder=entry_save)
        anno.annotate()
        entries[pdb_id] = anno.entry

    # Write the published-format index parquet so core.scores.entries.load_entry_views
    # (via core.scores.index.query_index) can read it back. Mirrors the layout
    # the scoring pipeline expects the ligand table, entry metadata, and a
    # splits parquet stub (load logic merges splits even though we do not filter
    # on them).
    rows = pd.concat([entry.to_df() for entry in entries.values()])
    (data_dir / "index").mkdir()
    rows.to_parquet(data_dir / "index" / "annotation_table.parquet", index=False)
    pd.concat([entry.metadata_to_df() for entry in entries.values()]).to_parquet(
        data_dir / "index" / "entry_metadata.parquet", index=False
    )
    pd.concat([entry.chains_to_df() for entry in entries.values()]).to_parquet(
        data_dir / "index" / "entry_chains.parquet", index=False
    )
    from plinder.data.annotations.interface_utils import protein_interfaces_to_table

    pq.write_table(
        pa.concat_tables(
            [
                protein_interfaces_to_table(entry.interfaces)
                for entry in entries.values()
            ]
        ),
        data_dir / "index" / "interface_annotation_table.parquet",
    )
    (data_dir / "splits").mkdir()
    pd.DataFrame(
        {"system_id": rows["system_id"].unique(), "split": "train"}
    ).to_parquet(data_dir / "splits" / "split.parquet", index=False)

    # Build a tiny seqres FASTA covering every chain in the two entries (not
    # just the ones chains_for_alignment returns) so that make_sub_db actually
    # has to filter the lookup — mirroring production's real pdb_seqres.txt.gz
    # which carries all chains across all PDB entries.
    from plinder.data.annotations.cif_utils import read_mmcif_file
    from plinder.data.annotations.protein_utils import get_seqres_from_cif

    seqres_lines: list[str] = []
    for pdb_id, cif in cifs.items():
        cif_obj = read_mmcif_file(Path(cif))
        block = list(cif_obj.values())[0]
        for auth_id, seq in get_seqres_from_cif(block).items():
            seqres_lines.append(
                f">{pdb_id}_{auth_id} mol:protein length:{len(seq)}\n{seq}\n"
            )
    seqres_path = data_dir / "dbs" / "seqres" / "pdb_seqres.txt.gz"
    seqres_path.unlink()
    with gzip.open(seqres_path, "wt") as f:
        f.write("".join(seqres_lines))

    make_db(
        input_dir=seqres_path,
        output_dir=data_dir / "dbs" / "mmseqs",
        db="mmseqs",
        tmp_dir=tmp_path / "mmseqs_full_tmp",
    )
    foldseek_in = data_dir / "raw_cifs"
    foldseek_in.mkdir()
    for pdb_id, cif in cifs.items():
        cif_path = Path(cif)
        shutil.copyfile(cif_path, foldseek_in / cif_path.name)
    make_db(
        input_dir=foldseek_in,
        output_dir=data_dir / "dbs" / "foldseek",
        db="foldseek",
        tmp_dir=tmp_path / "foldseek_full_tmp",
    )

    return data_dir, rows


def test_scoring_regression(scoring_fixture, tmp_path):
    from plinder.core.scores.entries import entry_views_from_df
    from plinder.core.scores.reconstruct import reconstruct_similarity_scores
    from plinder.data.annotations.get_similarity_scores import Scorer
    from plinder.data.pipeline import tasks
    from plinder.data.pipeline.utils import get_db_sources

    data_dir, annotation_rows = scoring_fixture
    entries = entry_views_from_df(
        annotation_rows,
        entry_chains=pd.read_parquet(data_dir / "index" / "entry_chains.parquet"),
        interface_annotations=pd.read_parquet(
            data_dir / "index" / "interface_annotation_table.parquet"
        ),
    )

    db_sources = get_db_sources(data_dir=data_dir, sub_databases=["holo"])
    scorer = Scorer(
        entries=entries,
        source_to_full_db_file=db_sources,
        db_dir=data_dir / "dbs" / "subdbs",
        scores_dir=data_dir / "scores",
    )
    scorer.make_dbs()
    output_folder = data_dir / "scores_output"
    output_folder.mkdir()
    scorer.run_alignments(
        entry_ids=PDB_IDS, search_db="holo", output_folder=output_folder
    )

    dfs = []
    for pdb_id in PDB_IDS:
        path = scorer.get_score_df(data_dir, pdb_id, "holo", overwrite=True)
        if path.exists():
            dfs.append(pd.read_parquet(path))
    assert dfs, "Scorer produced no score parquets — check fixture wiring"
    df = (
        pd.concat(dfs)
        .reset_index(drop=True)
        .sort_values(SORT_KEYS)
        .reset_index(drop=True)
    )
    assert not df[["query_ligand_id", "target_ligand_id"]].isna().any().any()
    assert {
        "shape",
        "color",
        "sucos_shape",
        "sucos_shape_pocket_qcov",
    }.issubset(set(df["metric"]))

    if os.environ.get("PLINDER_REGEN_SCORING") or not GOLDEN_PATH.exists():
        GOLDEN_PATH.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(GOLDEN_PATH, index=False)
        pytest.skip(f"regenerated reference scores at {GOLDEN_PATH}")

    expected = pd.read_parquet(GOLDEN_PATH)
    pd.testing.assert_frame_equal(df, expected, check_like=True)

    for shard in tasks.scatter_collate_alignments(data_dir=data_dir):
        tasks.collate_alignments(data_dir=data_dir, partition=shard)
    for release_shard in (data_dir / "alignments").rglob("*.parquet"):
        release_columns = set(pd.read_parquet(release_shard).columns)
        assert {
            "query_selected_residue_numbers",
            "target_selected_residue_numbers",
            "selected_residue_identity",
        }.issubset(release_columns)
        assert {
            "qrnum",
            "trnum",
            "qaa",
            "taa",
            "qaln",
            "taln",
            "query",
            "target",
        }.isdisjoint(release_columns)
    system_ids = set(annotation_rows["system_id"])
    reconstructed = (
        reconstruct_similarity_scores(
            system_ids,
            system_ids,
            data_dir=data_dir,
        )
        .sort_values(SORT_KEYS)
        .reset_index(drop=True)
    )

    def normalize_nulls(frame: pd.DataFrame) -> pd.DataFrame:
        return frame.astype(object).where(frame.notna(), None)

    similarity_mismatch = reconstructed["similarity"].ne(df["similarity"])
    assert not similarity_mismatch.any(), pd.concat(
        {
            "reconstructed": reconstructed.loc[similarity_mismatch, df.columns],
            "direct": df.loc[similarity_mismatch],
        },
        names=["score_path"],
    ).to_string()
    pd.testing.assert_frame_equal(
        normalize_nulls(reconstructed[df.columns]),
        normalize_nulls(df),
        check_dtype=False,
        check_like=True,
    )

    selected = df.iloc[0]
    bounded = (
        reconstruct_similarity_scores(
            [selected["query_system"]],
            [selected["target_system"]],
            query_ligand_ids=[selected["query_ligand_id"]],
            target_ligand_ids=[selected["target_ligand_id"]],
            data_dir=data_dir,
        )
        .sort_values(SORT_KEYS)
        .reset_index(drop=True)
    )
    expected_bounded = (
        df[
            (df["query_system"] == selected["query_system"])
            & (df["target_system"] == selected["target_system"])
            & (df["query_ligand_id"] == selected["query_ligand_id"])
            & (df["target_ligand_id"] == selected["target_ligand_id"])
        ]
        .sort_values(SORT_KEYS)
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(
        normalize_nulls(bounded[df.columns]),
        normalize_nulls(expected_bounded),
        check_dtype=False,
        check_like=True,
    )
