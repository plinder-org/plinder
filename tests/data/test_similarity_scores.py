# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem

from plinder.core.scores.entries import (
    ChainView,
    EntryView,
    InterfaceView,
    LigandView,
    SystemView,
    entry_views_from_df,
    load_entry_views,
)
from plinder.core.scores.reconstruct import (
    calculate_interface_similarity_scores,
    reconstruct_interface_similarity_scores,
)
from plinder.core.utils.schemas import PROTEIN_SIMILARITY_SCHEMA
from plinder.data.annotations import get_similarity_scores as scoring_module
from plinder.data.annotations.get_similarity_scores import (
    Scorer,
    _alignment_search_command,
    _stream_alignment_tsv_to_dataset,
    annotate_cofactor_similarity,
    annotate_ligand_3d_score_ability,
    annotate_ligand_similarity,
    build_ligand_similarity_annotations,
    compute_ligand_fingerprints,
    get_feature_map_score,
    ligand_scores,
    load_ligands_from_index,
    prepare_ligand_for_3d_scoring,
    run_alignment,
    write_ecfp4_fingerprint_table,
)
from plinder.data.pipeline.config import FoldseekConfig, MMSeqsConfig

SDF_FILE = (
    Path(__file__).resolve().parents[1]
    / "test_data"
    / "mini_system_files_new"
    / "1fbz__1__1.A__1.C"
    / "ligand_files"
    / "1.C.sdf"
)
HEM_SDF_FILE = (
    Path(__file__).resolve().parents[1]
    / "test_data"
    / "reconstructed_systems"
    / "19hc__1__1.B__1.T"
    / "ligand_files"
    / "1.T.sdf"
)


@pytest.mark.parametrize(
    ("query", "target", "expected_identity", "expected_similarity"),
    [
        ("ARND", "ARNE", 0.75, 1.0),
        ("ACDE", "AGDE", 0.75, 0.75),
        ("A-CD", "AC-D", 1.0, 1.0),
        ("AJUO", "ALCK", 1.0, 1.0),
    ],
)
def test_sequence_similarity_uses_cached_blosum_lookup(
    query: str,
    target: str,
    expected_identity: float,
    expected_similarity: float,
) -> None:
    identity, similarity = scoring_module.get_sequence_similarity(query, target)

    assert identity == pytest.approx(expected_identity)
    assert similarity == pytest.approx(expected_similarity)


def test_sequence_similarity_helper_returns_zero_for_incompatible_alignment() -> None:
    assert scoring_module.get_sequence_similarity_helper("ACD", "AC") == 0


def test_atomic_copy_file_allows_concurrent_writers(tmp_path: Path) -> None:
    sources = [tmp_path / "first", tmp_path / "second"]
    sources[0].write_bytes(b"first")
    sources[1].write_bytes(b"second")
    target = tmp_path / "score.parquet"
    old_shared_staging = tmp_path / "score.parquet.tmp"
    old_shared_staging.write_bytes(b"unrelated")

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(scoring_module._atomic_copy_file, source, target)
            for source in sources
        ]
        for future in futures:
            future.result()

    assert target.read_bytes() in {b"first", b"second"}
    assert old_shared_staging.read_bytes() == b"unrelated"
    assert not list(tmp_path.glob(".score.parquet.*.tmp"))


def test_protein_pair_scores_choose_best_duplicate_backend_hit(tmp_path: Path) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path,
        scores_dir=tmp_path,
    )
    alignments = pd.DataFrame(
        {
            "qcov": [0.95, 0.94, 0.8],
            "fident": [0.5, 0.5, 0.9],
            "seqsim": [0.6, 0.6, 0.9],
            "fident_qcov": [0.475, 0.47, 0.72],
            "seqsim_qcov": [0.57, 0.564, 0.72],
            "lddt": [0.67, 0.68, np.nan],
            "lddt_qcov": [0.6365, 0.6392, np.nan],
        },
        index=["foldseek", "foldseek", "mmseqs"],
    )

    scores = scorer.get_protein_scores_pair(alignments)

    assert scores["protein_qcov_foldseek"] == 0.94
    assert scores["protein_lddt_foldseek"] == 0.68
    assert scores["protein_fident_qcov_mmseqs"] == 0.72


def _entry_with_pocket_map(
    pdb_id: str,
    asym_id: str,
    number_to_index: dict[int, int],
) -> EntryView:
    system = SystemView(
        id=f"{pdb_id}-system",
        pdb_id=pdb_id,
        system_type="holo",
        protein_chains_asym_id=[f"1.{asym_id}"],
        proper_num_pocket_residues=len(number_to_index),
        proper_num_interactions=0,
        proper_num_unique_interactions=0,
        pocket_residue_number_to_index={f"1.{asym_id}": number_to_index},
    )
    return EntryView(
        pdb_id=pdb_id,
        chains={},
        systems={system.id: system},
        author_to_asym={asym_id: asym_id},
    )


@pytest.mark.parametrize(
    (
        "alignment_type",
        "qstart",
        "tstart",
        "query_numbers",
        "target_numbers",
        "identity",
    ),
    [
        (
            "foldseek",
            1,
            1,
            [30, 50],
            [400, 500],
            bytes([1, 1]),
        ),
        (
            "mmseqs",
            8,
            98,
            [10],
            [100],
            bytes([1]),
        ),
    ],
)
def test_map_row_vectorizes_sparse_pocket_positions(
    tmp_path: Path,
    alignment_type: str,
    qstart: int,
    tstart: int,
    query_numbers: list[int],
    target_numbers: list[int],
    identity: bytes,
) -> None:
    scorer = Scorer(
        entries={
            "query": _entry_with_pocket_map("query", "A", {10: 1, 30: 3, 50: 5}),
            "target": _entry_with_pocket_map("target", "B", {100: 0, 400: 3, 500: 4}),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    row = pd.Series(
        {
            "query_entry": "query",
            "target_entry": "target",
            "query_chain_mapped": "A",
            "target_chain_mapped": "B",
            "qstart": qstart,
            "tstart": tstart,
            "qaln": "A-BCDEF",
            "taln": "AB-CD-F",
        }
    )

    mapped = scorer.map_row(row, aln_type=alignment_type, search_db="holo")

    assert mapped["query_selected_residue_numbers"] == query_numbers
    assert mapped["target_selected_residue_numbers"] == target_numbers
    assert mapped["selected_residue_identity"] == identity


def test_search_uses_configured_zero_minimum_sequence_identity(tmp_path) -> None:
    common = {
        "query_db": tmp_path / "query",
        "target_db": tmp_path / "target",
        "search_db": tmp_path / "search",
        "tmp_dir": tmp_path / "scratch",
        "threads": 4,
    }
    foldseek = _alignment_search_command(
        aln_type="foldseek",
        alignment_config=FoldseekConfig(min_seq_id=0.0),
        **common,
    )
    mmseqs = _alignment_search_command(
        aln_type="mmseqs",
        alignment_config=MMSeqsConfig(min_seq_id=0.0),
        **common,
    )

    assert foldseek[foldseek.index("--max-seqs") + 1] == "10000"
    assert "-e" in foldseek
    assert foldseek[foldseek.index("--min-seq-id") + 1] == "0.0"
    assert "--cov-mode" in foldseek
    assert mmseqs[mmseqs.index("--min-seq-id") + 1] == "0.0"
    assert "--cov-mode" in mmseqs

    clustered_foldseek = _alignment_search_command(
        aln_type="foldseek",
        alignment_config=FoldseekConfig(min_seq_id=0.0),
        expand_exact_clusters=True,
        **common,
    )
    assert clustered_foldseek[clustered_foldseek.index("--cluster-search") + 1] == "1"


def test_mmseqs_search_expands_and_realigns_exact_cluster_members(
    tmp_path, monkeypatch
) -> None:
    commands: list[list[str]] = []
    monkeypatch.setattr(
        scoring_module.subprocess,
        "check_call",
        lambda command, **_kwargs: commands.append(command),
    )
    monkeypatch.setattr(
        scoring_module, "_stream_alignment_tsv_to_dataset", lambda *args, **kwargs: None
    )

    run_alignment(
        aln_type="mmseqs",
        query_db=tmp_path / "query",
        target_db=tmp_path / "full",
        search_target_db=tmp_path / "representatives",
        cluster_alignment_db=tmp_path / "cluster_alignments",
        search_db=tmp_path / "search",
        aln_file=tmp_path / "alignments.tsv",
        alignment_config=MMSeqsConfig(min_seq_id=0.0),
        tmp_dir=tmp_path / "scratch",
        remove_tmp=False,
        threads=3,
    )

    assert [command[1] for command in commands] == [
        "search",
        "expandaln",
        "align",
        "convertalis",
    ]
    assert commands[0][3] == str(tmp_path / "representatives")
    assert commands[1][5] == str(tmp_path / "cluster_alignments")
    assert commands[2][3] == str(tmp_path / "full")
    assert commands[0][commands[0].index("--min-seq-id") + 1] == "0.0"
    assert "--min-seq-id" not in commands[1]
    assert commands[2][commands[2].index("--min-seq-id") + 1] == "0.0"
    assert "--min-seq-id" not in commands[3]


def test_mmseqs_search_accepts_an_unclustered_target(tmp_path, monkeypatch) -> None:
    commands: list[list[str]] = []
    monkeypatch.setattr(
        scoring_module.subprocess,
        "check_call",
        lambda command, **_kwargs: commands.append(command),
    )
    monkeypatch.setattr(
        scoring_module, "_stream_alignment_tsv_to_dataset", lambda *args, **kwargs: None
    )

    run_alignment(
        aln_type="mmseqs",
        query_db=tmp_path / "query",
        target_db=tmp_path / "selected-targets",
        search_target_db=tmp_path / "selected-targets",
        cluster_alignment_db=None,
        search_db=tmp_path / "search",
        aln_file=tmp_path / "alignments.tsv",
        alignment_config=MMSeqsConfig(min_seq_id=0.0),
        tmp_dir=tmp_path / "scratch",
        remove_tmp=False,
        threads=2,
    )

    assert [command[1] for command in commands] == ["search", "convertalis"]
    assert commands[0][3] == str(tmp_path / "selected-targets")
    assert commands[1][3] == str(tmp_path / "selected-targets")


def test_no_hit_search_writes_typed_empty_raw_and_mapped_checkpoints(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={"1abc": object()},
        source_to_full_db_file={"holo_mmseqs": tmp_path / "full"},
        db_dir=tmp_path / "dbs" / "subdbs",
        scores_dir=tmp_path / "scores",
    )
    monkeypatch.setattr(
        scoring_module.databases,
        "get_db_ids",
        lambda *_args, **_kwargs: {"1abc_A"},
    )
    monkeypatch.setattr(
        scoring_module.databases,
        "make_sub_db",
        lambda *_args, **_kwargs: set(),
    )
    monkeypatch.setattr(scoring_module, "run_alignment", lambda **_kwargs: None)

    scorer.run_alignments(
        entry_ids=["1abc"],
        search_db="holo",
        output_folder=tmp_path / "work",
        alignment_types=["mmseqs"],
    )

    raw = scorer.db_dir / "holo_mmseqs" / "aln" / "1abc.parquet"
    assert pd.read_parquet(raw).empty
    assert set(scoring_module.pq.read_schema(raw).names) == {
        "query",
        "target",
        "qlen",
        "fident",
        "alnlen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bits",
        "qcov",
        "tcov",
        "qaln",
        "taln",
        "target_pdb_id",
    }
    mapped = scorer.map_alignment_df(raw, "mmseqs", "holo")
    assert mapped.empty
    assert list(mapped.index.names) == [
        "query_entry",
        "target_entry",
        "query_chain_mapped",
        "target_chain_mapped",
        "source",
    ]
    assert {
        "seqsim",
        "query_selected_residue_numbers",
        "target_selected_residue_numbers",
        "selected_residue_identity",
    } <= set(mapped.columns)
    assert {"evalue", "bits", "tcov"}.isdisjoint(mapped.columns)


def test_search_can_probe_an_alternate_target_without_writing_empty_results(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={"holo_mmseqs": tmp_path / "full"},
        db_dir=tmp_path / "dbs" / "subdbs",
        scores_dir=tmp_path / "scores",
    )
    monkeypatch.setattr(
        scoring_module.databases,
        "make_sub_db",
        lambda *_args, **_kwargs: set(),
    )
    target_roots = []

    def target_paths(root, *_args):
        target_roots.append(root)
        return root / "target", root / "selected", None

    monkeypatch.setattr(
        scoring_module.databases,
        "exact_search_database_paths",
        target_paths,
    )

    def search(**kwargs):
        output = kwargs["aln_file"].with_suffix(".parquet") / "query_pdb_id=1abc"
        output.mkdir(parents=True)
        pd.DataFrame(
            {"query": ["1abc_A"], "target": ["3ghi_A"], "target_pdb_id": ["3ghi"]}
        ).to_parquet(output / "part.parquet", index=False)

    monkeypatch.setattr(scoring_module, "run_alignment", search)
    result_root = tmp_path / "probe-results"
    hits = scorer.run_alignments(
        entry_ids=["1abc", "2def"],
        search_db="holo",
        output_folder=tmp_path / "work",
        alignment_types=["mmseqs"],
        query_chain_auth_ids={"1abc": {"A"}, "2def": {"B"}},
        target_database_dir=tmp_path / "changed-targets",
        result_database_dir=result_root,
        write_empty_results=False,
    )

    assert target_roots == [tmp_path / "changed-targets"]
    assert hits == {"1abc"}
    assert (result_root / "holo_mmseqs/aln/1abc.parquet").is_file()
    assert not (result_root / "holo_mmseqs/aln/2def.parquet").exists()
    assert not (scorer.db_dir / "holo_mmseqs/aln").exists()


def test_unavailable_query_backend_writes_typed_empty_checkpoint(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={"1abc": object()},
        source_to_full_db_file={"holo_foldseek": tmp_path / "full"},
        db_dir=tmp_path / "dbs" / "subdbs",
        scores_dir=tmp_path / "scores",
    )
    monkeypatch.setattr(
        scoring_module.databases,
        "get_db_ids",
        lambda *_args, **_kwargs: {"pdb_00001abc_xyz-enrich_A"},
    )
    monkeypatch.setattr(
        scoring_module.databases,
        "make_sub_db",
        lambda db_ids, *_args, **_kwargs: db_ids,
    )
    monkeypatch.setattr(
        scoring_module,
        "run_alignment",
        lambda **_kwargs: pytest.fail("an unavailable backend must not be searched"),
    )

    scorer.run_alignments(
        entry_ids=["1abc"],
        search_db="apo",
        output_folder=tmp_path / "work",
        alignment_types=["foldseek"],
    )

    raw = scorer.db_dir / "apo_foldseek" / "aln" / "1abc.parquet"
    assert raw.is_file()
    assert pd.read_parquet(raw).empty
    assert pq.read_schema(raw).equals(scoring_module._raw_alignment_schema("foldseek"))


@pytest.mark.parametrize(
    ("alignment_type", "expected_ids"),
    [
        ("foldseek", {"pdb_00001abc_xyz-enrich_R", "pdb_00001abc_xyz-enrich_S"}),
        ("mmseqs", {"1abc_R", "1abc_S"}),
    ],
)
def test_run_alignments_accepts_explicit_query_chains_without_entry_views(
    tmp_path, monkeypatch, alignment_type, expected_ids
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={f"holo_{alignment_type}": tmp_path / "full"},
        db_dir=tmp_path / "dbs" / "subdbs",
        scores_dir=tmp_path / "scores",
    )
    observed = []
    monkeypatch.setattr(
        scoring_module.databases,
        "get_db_ids",
        lambda *_args, **_kwargs: pytest.fail("explicit chains must avoid EntryView"),
    )

    def unavailable(ids, *_args, **_kwargs):
        observed.append(ids)
        return ids

    monkeypatch.setattr(scoring_module.databases, "make_sub_db", unavailable)

    scorer.run_alignments(
        entry_ids=["1abc"],
        search_db="apo",
        output_folder=tmp_path / "work",
        alignment_types=[alignment_type],
        query_chain_auth_ids={"1abc": {"R", "S"}},
    )

    assert observed == [expected_ids]
    assert pd.read_parquet(
        scorer.db_dir / f"apo_{alignment_type}/aln/1abc.parquet"
    ).empty


def test_alignment_mapping_preserves_author_chain_ids_with_underscores(
    tmp_path,
) -> None:
    def mapping_entry(pdb_id: str, author_id: str, asym_id: str) -> EntryView:
        system = SystemView(
            id="alignment",
            pdb_id=pdb_id,
            system_type="holo",
            protein_chains_asym_id=[f"1.{asym_id}"],
            proper_num_pocket_residues=1,
            proper_num_interactions=0,
            proper_num_unique_interactions=0,
            pocket_residue_number_to_index={f"1.{asym_id}": {1: 0}},
        )
        return EntryView(
            pdb_id=pdb_id,
            chains={},
            systems={"alignment": system},
            author_to_asym={author_id: asym_id},
        )

    scorer = Scorer(
        entries={
            "2aaz": mapping_entry("2aaz", "AUTH_WITH_UNDERSCORES", "A"),
            "1abc": mapping_entry("1abc", "TARGET_WITH_UNDERSCORES", "B"),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    raw = tmp_path / "raw.parquet"
    pd.DataFrame(
        {
            "query": ["pdb_00002aaz_xyz-enrich_MODEL_1_AUTH_WITH_UNDERSCORES"],
            "target": ["pdb_00001abc_xyz-enrich_TARGET_WITH_UNDERSCORES"],
            "qlen": [1],
            "fident": [1.0],
            "alnlen": [1],
            "qstart": [1],
            "qend": [1],
            "tstart": [1],
            "tend": [1],
            "evalue": [0.0],
            "bits": [1],
            "qcov": [1.0],
            "tcov": [1.0],
            "qaln": ["A"],
            "taln": ["A"],
            "lddt": [1.0],
            "target_pdb_id": ["1abc"],
        }
    ).to_parquet(raw, index=False)

    mapped = scorer.map_alignment_df(raw, "foldseek", "holo")

    assert mapped.index.to_list() == [("2aaz", "1abc", "A", "B", "foldseek")]


def test_foldseek_identifier_cleanup_does_not_remove_cif_from_pdb_id(
    tmp_path,
) -> None:
    def mapping_entry(pdb_id: str, asym_id: str) -> EntryView:
        system = SystemView(
            id="alignment",
            pdb_id=pdb_id,
            system_type="holo",
            protein_chains_asym_id=[f"1.{asym_id}"],
            proper_num_pocket_residues=1,
            proper_num_interactions=0,
            proper_num_unique_interactions=0,
            pocket_residue_number_to_index={f"1.{asym_id}": {1: 0}},
        )
        return EntryView(
            pdb_id=pdb_id,
            chains={},
            systems={"alignment": system},
            author_to_asym={asym_id: asym_id},
        )

    scorer = Scorer(
        entries={
            "1cif": mapping_entry("1cif", "A"),
            "1abc": mapping_entry("1abc", "B"),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    raw = tmp_path / "raw.parquet"
    pd.DataFrame(
        {
            "query": ["pdb_00001cif_xyz-enrich_A"],
            "target": ["pdb_00001abc_xyz-enrich_B"],
            "qlen": [1],
            "fident": [1.0],
            "alnlen": [1],
            "qstart": [1],
            "qend": [1],
            "tstart": [1],
            "tend": [1],
            "evalue": [0.0],
            "bits": [1],
            "qcov": [1.0],
            "tcov": [1.0],
            "qaln": ["A"],
            "taln": ["A"],
            "lddt": [1.0],
            "target_pdb_id": ["1abc"],
        }
    ).to_parquet(raw, index=False)

    mapped = scorer.map_alignment_df(raw, "foldseek", "holo")

    assert mapped.index.to_list() == [("1cif", "1abc", "A", "B", "foldseek")]


def test_foldseek_cluster_target_is_also_used_for_conversion(
    tmp_path, monkeypatch
) -> None:
    commands: list[list[str]] = []
    stream_options = []
    monkeypatch.setattr(
        scoring_module.subprocess,
        "check_call",
        lambda command, **_kwargs: commands.append(command),
    )
    monkeypatch.setattr(
        scoring_module,
        "_stream_alignment_tsv_to_dataset",
        lambda *_args, **kwargs: stream_options.append(kwargs),
    )

    run_alignment(
        aln_type="foldseek",
        query_db=tmp_path / "query",
        target_db=tmp_path / "full-build-db",
        search_target_db=tmp_path / "clustered",
        search_db=tmp_path / "search",
        aln_file=tmp_path / "alignments.tsv",
        alignment_config=FoldseekConfig(),
        tmp_dir=tmp_path / "scratch",
        remove_tmp=False,
        threads=2,
        include_target_pdb_id=False,
    )

    assert [command[1] for command in commands] == ["search", "convertalis"]
    assert "--cluster-search" in commands[0]
    assert commands[1][3] == str(tmp_path / "clustered")
    assert stream_options == [{"aln_type": "foldseek", "include_target_pdb_id": False}]


def test_alignment_tsv_is_streamed_to_query_partitions(tmp_path) -> None:
    tsv_path = tmp_path / "alignment.tsv"
    tsv_path.write_text(
        "query\ttarget\n"
        "pdb_00001abc_A\tpdb_00002def_B\n"
        "pdb_00003ghi_C\tpdb_00004jkl_D\n"
    )
    dataset_path = tmp_path / "alignment.parquet"

    _stream_alignment_tsv_to_dataset(
        tsv_path,
        dataset_path,
        aln_type="foldseek",
        include_target_pdb_id=True,
    )

    first = pd.read_parquet(dataset_path / "query_pdb_id=1abc")
    assert first["target_pdb_id"].tolist() == ["2def"]
    assert (dataset_path / "query_pdb_id=3ghi").is_dir()


@pytest.mark.parametrize("aln_type", ["foldseek", "mmseqs"])
def test_empty_alignment_tsv_writes_readable_dataset(tmp_path, aln_type) -> None:
    columns = [field.name for field in scoring_module._raw_alignment_schema(aln_type)]
    columns.remove("target_pdb_id")
    tsv_path = tmp_path / "alignment.tsv"
    tsv_path.write_text("\t".join(columns) + "\n")
    dataset_path = tmp_path / "alignment.parquet"

    _stream_alignment_tsv_to_dataset(
        tsv_path,
        dataset_path,
        aln_type=aln_type,
        include_target_pdb_id=True,
    )

    result = pd.read_parquet(dataset_path)
    assert result.empty
    assert result.columns.tolist() == [
        field.name for field in scoring_module._raw_alignment_schema(aln_type)
    ]


def test_ligand_scoring_inputs_include_only_proper_holo_ligands() -> None:
    annotation = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "1abc"],
            "system_id": ["proper", "artifact", "ion-system"],
            "system_type": ["holo", "holo", "ion"],
            "ligand_is_proper": [True, False, True],
            "ligand_smiles": ["CCO", "O", "[Na+]"],
            "ligand_unique_ccd_code": ["LIG", "HOH", "NA"],
            "ligand_id": ["1abc__1__1.L", "1abc__1__1.W", "1abc__1__1.N"],
            "ligand_asym_id": ["L", "W", "N"],
        }
    )

    ligands = load_ligands_from_index(annotation=annotation)

    assert ligands["ligand_id"].tolist() == ["1abc__1__1.L"]


def test_entry_views_accept_annotation_dataframe(
    cif_2gdo, tmp_path, monkeypatch
) -> None:
    from plinder.data.annotations import ligand_utils
    from plinder.data.annotations.aggregate_annotations import Entry

    monkeypatch.setattr(ligand_utils, "BINDING_AFFINITY", {})
    entry = Entry.from_cif_file(cif_2gdo)
    annotation = entry.to_df()
    assert [column for column in annotation.columns if column.startswith("entry_")] == [
        "entry_pdb_id"
    ]
    assert "entry_release_date" in entry.metadata_to_df().columns
    assert not any(column.startswith("entry_chains_") for column in annotation)
    entry_chains = entry.chains_to_df()
    assert len(entry_chains) <= len(entry.chains)
    assert entry_chains["chain_type"].str.lower().str.contains("polypeptide").all()

    # Exercise the Arrow representation used by the ingest pipeline,
    # including the nested UniProt accession lists.
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    annotation_path = index_dir / "annotation_table.parquet"
    chain_path = index_dir / "entry_chains.parquet"
    interface_path = index_dir / "interface_annotation_table.parquet"
    annotation.to_parquet(annotation_path, index=False)
    entry_chains.to_parquet(chain_path, index=False)
    from plinder.data.annotations.interface_utils import protein_interfaces_to_table

    pq.write_table(protein_interfaces_to_table(entry.interfaces), interface_path)
    view = entry_views_from_df(
        pd.read_parquet(annotation_path),
        entry_chains=pd.read_parquet(chain_path),
        interface_annotations=pd.read_parquet(interface_path),
    )[entry.pdb_id]

    assert view.systems
    for system_id, system in entry.systems.items():
        expected = {ligand.instance_chain for ligand in system.ligands}
        assert set(view.systems[system_id].ligands) == expected
    for chain_type in ["holo", "apo", "pred"]:
        for aln_type in ["foldseek", "mmseqs"]:
            assert sorted(view.chains_for_alignment(chain_type, aln_type)) == sorted(
                entry.chains_for_alignment(chain_type, aln_type)
            )

    loaded = load_entry_views(pdb_ids=[entry.pdb_id], data_dir=tmp_path)[entry.pdb_id]
    assert loaded.chains == view.chains
    assert loaded.interfaces == view.interfaces

    chain_path.unlink()
    with pytest.raises(FileNotFoundError, match="missing entry chain index"):
        load_entry_views(pdb_ids=[entry.pdb_id], data_dir=tmp_path)


def test_entry_view_builds_holo_apo_and_pred_alignment_ids() -> None:
    system = SystemView(
        id="1abc__1__1.A__1.L",
        pdb_id="1abc",
        system_type="holo",
        protein_chains_asym_id=["1.A"],
        proper_num_pocket_residues=1,
        proper_num_interactions=1,
        proper_num_unique_interactions=1,
    )
    view = EntryView(
        pdb_id="1abc",
        chains={
            "A": ChainView(
                asym_id="A",
                auth_id="R",
                entity_id="1",
                chain_type="polypeptide(L)",
                length=200,
                holo=True,
                uniprot_ids=("P12345",),
            ),
            "B": ChainView(
                asym_id="B",
                auth_id="B",
                entity_id="2",
                chain_type="polypeptide(L)",
                length=150,
                holo=False,
            ),
            # A non-holo copy of a holo entity must not enter the apo DB.
            "C": ChainView(
                asym_id="C",
                auth_id="C",
                entity_id="1",
                chain_type="polypeptide(L)",
                length=200,
                holo=False,
            ),
            # Search DBs are protein-only.
            "N": ChainView(
                asym_id="N",
                auth_id="N",
                entity_id="3",
                chain_type="polyribonucleotide",
                length=80,
                holo=False,
            ),
        },
        systems={system.id: system},
        author_to_asym={"R": "A", "B": "B", "C": "C"},
    )

    assert view.chains_for_alignment("holo", "foldseek") == [
        "pdb_00001abc_xyz-enrich_R"
    ]
    assert view.chains_for_alignment("holo", "mmseqs") == ["1abc_R"]
    assert view.chains_for_alignment("apo", "foldseek") == ["pdb_00001abc_xyz-enrich_B"]
    assert view.chains_for_alignment("apo", "mmseqs") == ["1abc_B"]
    assert view.chains_for_alignment("pred", "foldseek") == ["AF-P12345-F1-model_v4_A"]
    assert view.chains_for_alignment("pred", "mmseqs") == ["P12345"]


def test_interface_scores_choose_swapped_assignment_and_best_backend() -> None:
    query = InterfaceView(
        id="1abc__1__1.A--1.B",
        pdb_id="1abc",
        biounit_id="1",
        chain_1="1.A",
        chain_2="1.B",
        chain_1_residue_number_to_index={1: 0, 2: 1},
        chain_2_residue_number_to_index={3: 2, 4: 3},
        num_contact_residue_pairs=4,
    )
    target = InterfaceView(
        id="2def__1__1.X--1.Y",
        pdb_id="2def",
        biounit_id="1",
        chain_1="1.X",
        chain_2="1.Y",
        chain_1_residue_number_to_index={10: 9, 20: 19},
        chain_2_residue_number_to_index={30: 29, 40: 39},
        num_contact_residue_pairs=4,
    )
    alignments = pd.DataFrame(
        [
            # Foldseek's swapped assignment covers both interface sides.
            {
                "query_entry": "1abc",
                "target_entry": "2def",
                "query_chain_mapped": "A",
                "target_chain_mapped": "Y",
                "source": "foldseek",
                "query_selected_residue_numbers": [1, 2],
                "target_selected_residue_numbers": [30, 40],
            },
            {
                "query_entry": "1abc",
                "target_entry": "2def",
                "query_chain_mapped": "B",
                "target_chain_mapped": "X",
                "source": "foldseek",
                "query_selected_residue_numbers": [3, 4],
                "target_selected_residue_numbers": [10, 20],
            },
            # MMseqs has a complete direct assignment but only 50% product.
            {
                "query_entry": "1abc",
                "target_entry": "2def",
                "query_chain_mapped": "A",
                "target_chain_mapped": "X",
                "source": "mmseqs",
                "query_selected_residue_numbers": [1, 2],
                "target_selected_residue_numbers": [10, -1],
            },
            {
                "query_entry": "1abc",
                "target_entry": "2def",
                "query_chain_mapped": "B",
                "target_chain_mapped": "Y",
                "source": "mmseqs",
                "query_selected_residue_numbers": [3, 4],
                "target_selected_residue_numbers": [30, 40],
            },
            # Reverse coverage is directional and only 25%.
            {
                "query_entry": "2def",
                "target_entry": "1abc",
                "query_chain_mapped": "X",
                "target_chain_mapped": "A",
                "source": "foldseek",
                "query_selected_residue_numbers": [10, 20],
                "target_selected_residue_numbers": [1, -1],
            },
            {
                "query_entry": "2def",
                "target_entry": "1abc",
                "query_chain_mapped": "Y",
                "target_chain_mapped": "B",
                "source": "foldseek",
                "query_selected_residue_numbers": [30, 40],
                "target_selected_residue_numbers": [3, -1],
            },
        ]
    )

    forward = calculate_interface_similarity_scores(
        alignments,
        query_interfaces={query.id: query},
        target_interfaces={target.id: target},
    )
    assert forward.to_dict("records") == [
        {
            "query_system": query.id,
            "target_system": target.id,
            "mapping": "1.A:1.Y;1.B:1.X",
            "source": "foldseek",
            "metric": "interface_qcov",
            "iface1_qcov": 1.0,
            "iface2_qcov": 1.0,
            "similarity": 100,
        }
    ]

    reverse = calculate_interface_similarity_scores(
        alignments,
        query_interfaces={target.id: target},
        target_interfaces={query.id: query},
    )
    assert reverse.loc[0, "similarity"] == 25
    assert reverse.loc[0, "iface1_qcov"] == 0.5
    assert reverse.loc[0, "iface2_qcov"] == 0.5
    assert reverse.loc[0, "mapping"] == "1.X:1.A;1.Y:1.B"

    incomplete = calculate_interface_similarity_scores(
        alignments.iloc[[0]],
        query_interfaces={query.id: query},
        target_interfaces={target.id: target},
    )
    assert incomplete.empty

    unrelated = pd.DataFrame(
        [
            {
                "query_entry": "1abc",
                "target_entry": "2def",
                "query_chain_mapped": query_chain,
                "target_chain_mapped": target_chain,
                "source": "foldseek",
                "query_selected_residue_numbers": [999],
                "target_selected_residue_numbers": [999],
            }
            for query_chain, target_chain in (("A", "X"), ("B", "Y"))
        ]
    )
    assert calculate_interface_similarity_scores(
        unrelated,
        query_interfaces={query.id: query},
        target_interfaces={target.id: target},
    ).empty


def test_entry_views_support_interface_only_entries() -> None:
    interface_id = "1abc__1__1.A--2.B"
    interfaces = pd.DataFrame(
        [
            {
                "entry_pdb_id": "1abc",
                "system_id": interface_id,
                "system_biounit_id": "1",
                "interface_chain_1": "1.A",
                "interface_chain_2": "2.B",
                "interface_chain_1_residue_numbers": [10, 11, 12],
                "interface_chain_1_residue_indices": [9, 10, 11],
                "interface_chain_2_residue_numbers": [20, 21, 22],
                "interface_chain_2_residue_indices": [19, 20, 21],
                "interface_num_contact_residue_pairs": 5,
            }
        ]
    )
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc"],
            "chain_asym_id": ["A", "B"],
            "chain_auth_id": ["X", "Y"],
            "chain_entity_id": ["1", "2"],
            "chain_type": ["polypeptide(L)", "polypeptide(L)"],
            "chain_length": [100, 80],
            "chain_is_holo": [True, True],
            "chain_uniprot_ids": [[], []],
        }
    )

    entry = entry_views_from_df(
        pd.DataFrame(columns=["entry_pdb_id"]),
        entry_chains=chains,
        interface_annotations=interfaces,
    )["1abc"]

    assert not entry.systems
    assert set(entry.interfaces) == {interface_id}
    assert entry.chains_for_alignment("holo", "mmseqs") == ["1abc_X", "1abc_Y"]
    assert entry.selected_index_to_number_per_chain == {
        "A": {9: 10, 10: 11, 11: 12},
        "B": {19: 20, 20: 21, 21: 22},
    }


def test_entry_views_support_protein_chain_only_entries() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["model"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["X"],
            "chain_entity_id": ["1"],
            "chain_type": ["polypeptide(L)"],
            "chain_length": [100],
            "chain_is_holo": [True],
            "chain_uniprot_ids": [[]],
        }
    )

    entry = entry_views_from_df(
        pd.DataFrame(columns=["entry_pdb_id"]),
        entry_chains=chains,
    )["model"]

    assert not entry.systems
    assert not entry.interfaces
    assert set(entry.chains) == {"A"}
    assert entry.chains["A"].length == 100


def test_reconstruct_interface_scores_from_release_shard(tmp_path: Path) -> None:
    index = tmp_path / "index"
    index.mkdir()
    pd.DataFrame({"entry_pdb_id": pd.Series(dtype="string")}).to_parquet(
        index / "annotation_table.parquet", index=False
    )
    pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "2def", "2def"],
            "chain_asym_id": ["A", "B", "X", "Y"],
            "chain_auth_id": ["A", "B", "X", "Y"],
            "chain_entity_id": ["1", "2", "1", "2"],
            "chain_type": ["polypeptide(L)"] * 4,
            "chain_length": [100] * 4,
            "chain_is_holo": [True] * 4,
            "chain_uniprot_ids": [[] for _ in range(4)],
        }
    ).to_parquet(index / "entry_chains.parquet", index=False)
    query_id = "1abc__1__1.A--1.B"
    target_id = "2def__1__1.X--1.Y"
    interface_rows = [
        {
            "entry_pdb_id": "1abc",
            "system_id": query_id,
            "system_biounit_id": "1",
            "interface_chain_1": "1.A",
            "interface_chain_2": "1.B",
            "interface_chain_1_residue_numbers": [1, 2, 3],
            "interface_chain_1_residue_indices": [0, 1, 2],
            "interface_chain_2_residue_numbers": [4, 5, 6],
            "interface_chain_2_residue_indices": [3, 4, 5],
            "interface_num_contact_residue_pairs": 3,
        },
        {
            "entry_pdb_id": "2def",
            "system_id": target_id,
            "system_biounit_id": "1",
            "interface_chain_1": "1.X",
            "interface_chain_2": "1.Y",
            "interface_chain_1_residue_numbers": [10, 20, 30],
            "interface_chain_1_residue_indices": [9, 19, 29],
            "interface_chain_2_residue_numbers": [40, 50, 60],
            "interface_chain_2_residue_indices": [39, 49, 59],
            "interface_num_contact_residue_pairs": 3,
        },
    ]
    from plinder.data.annotations.interface_utils import INTERFACE_ANNOTATION_SCHEMA

    pq.write_table(
        pa.Table.from_pylist(interface_rows, schema=INTERFACE_ANNOTATION_SCHEMA),
        index / "interface_annotation_table.parquet",
    )
    alignment = (
        tmp_path
        / "alignments"
        / "search_db=holo"
        / "alignment_type=foldseek"
        / "shard=ab.parquet"
    )
    alignment.parent.mkdir(parents=True)
    alignment_rows = [
        {
            "query_entry": "1abc",
            "target_entry": "2def",
            "query_chain_mapped": "A",
            "target_chain_mapped": "X",
            "source": "foldseek",
            "qcov": 1.0,
            "fident": 1.0,
            "seqsim": 1.0,
            "query_selected_residue_numbers": [1, 2, 3],
            "target_selected_residue_numbers": [10, 20, 30],
            "selected_residue_identity": bytes([1, 1, 1]),
            "lddt": 1.0,
        },
        {
            "query_entry": "1abc",
            "target_entry": "2def",
            "query_chain_mapped": "B",
            "target_chain_mapped": "Y",
            "source": "foldseek",
            "qcov": 1.0,
            "fident": 1.0,
            "seqsim": 1.0,
            "query_selected_residue_numbers": [4, 5, 6],
            "target_selected_residue_numbers": [40, 50, 60],
            "selected_residue_identity": bytes([1, 1, 1]),
            "lddt": 1.0,
        },
    ]
    from plinder.core.utils.schemas import mapped_alignment_schema

    pq.write_table(
        pa.Table.from_pylist(
            alignment_rows,
            schema=mapped_alignment_schema(alignment_type="foldseek"),
        ),
        alignment,
    )

    scores = reconstruct_interface_similarity_scores(
        [query_id],
        [target_id],
        data_dir=tmp_path,
    )

    assert scores.to_dict("records") == [
        {
            "query_system": query_id,
            "target_system": target_id,
            "mapping": "1.A:1.X;1.B:1.Y",
            "source": "foldseek",
            "metric": "interface_qcov",
            "iface1_qcov": 1.0,
            "iface2_qcov": 1.0,
            "similarity": 100,
        }
    ]


def test_entry_views_keep_ligand_pockets_separate() -> None:
    common = {
        "entry_pdb_id": "1abc",
        "system_id": "1abc__1__1.A__1.B_1.C",
        "system_type": "holo",
        "system_protein_chains_asym_id": ["1.A"],
        "system_protein_chains_auth_id": ["A"],
        "system_protein_chains_length": [200],
        "system_proper_num_pocket_residues": 3,
        "system_proper_num_interactions": 2,
        "system_proper_num_unique_interactions": 2,
        "ligand_is_proper": True,
        "ligand_protein_chains_asym_id": ["1.A"],
        "ligand_num_interactions": 1,
        "ligand_num_unique_interactions": 1,
    }
    rows = pd.DataFrame(
        [
            {
                **common,
                "ligand_id": "1abc__1__1.B",
                "ligand_instance_chain": "1.B",
                "ligand_asym_id": "B",
                "ligand_num_pocket_residues": 2,
                "ligand_is_shape_comparable": False,
                "ligand_neighboring_residues": ["1.A_10_9_10", "1.A_20_19_20"],
                "ligand_interactions": ["1.A_10_type:hydrogen_bonds"],
            },
            {
                **common,
                "ligand_id": "1abc__1__1.C",
                "ligand_instance_chain": "1.C",
                "ligand_asym_id": "C",
                "ligand_num_pocket_residues": 1,
                "ligand_is_shape_comparable": True,
                "ligand_neighboring_residues": ["1.A_30_29_30"],
                "ligand_interactions": ["1.A_30_type:hydrophobic_contacts"],
            },
        ]
    )

    system = entry_views_from_df(rows)["1abc"].systems[common["system_id"]]

    assert set(system.ligands) == {"1.B", "1.C"}
    assert system.ligands["1.B"].pocket_residue_number_to_index == {
        "1.A": {10: 9, 20: 19}
    }
    assert system.ligands["1.C"].pocket_residue_number_to_index == {"1.A": {30: 29}}
    assert not system.ligands["1.B"].is_shape_comparable
    assert system.ligands["1.C"].is_shape_comparable
    # The existing system view remains the union for backward compatibility.
    assert system.pocket_residue_number_to_index == {"1.A": {10: 9, 20: 19, 30: 29}}


def test_entry_views_exclude_na_receptors_from_scoring(tmp_path) -> None:
    rows = pd.DataFrame(
        [
            {
                "entry_pdb_id": "1abc",
                "system_id": "1abc__1__1.A_1.N__1.L",
                "system_type": "holo",
                "system_receptor_type": "protein+dna",
                "system_protein_chains_asym_id": ["1.A", "1.N"],
                "system_protein_chains_auth_id": ["A", "N"],
                "system_protein_chains_length": [100, 20],
                "system_proper_num_pocket_residues": 2,
                "system_proper_num_interactions": 2,
                "system_proper_num_unique_interactions": 2,
                "ligand_id": "1abc__1__1.L",
                "ligand_instance_chain": "1.L",
                "ligand_asym_id": "L",
                "ligand_is_proper": True,
                "ligand_protein_chains_asym_id": ["1.A", "1.N"],
                "ligand_num_pocket_residues": 2,
                "ligand_num_interactions": 2,
                "ligand_num_unique_interactions": 2,
                "ligand_neighboring_residues": ["1.A_10_9_10", "1.N_2_1_2"],
                "ligand_interactions": [
                    "1.A_10_type:hydrogen_bonds",
                    "1.N_2_type:hydrogen_bonds",
                ],
            }
        ]
    )
    entry_chains = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["A"],
            "chain_entity_id": ["1"],
            "chain_type": ["polypeptide(L)"],
            "chain_receptor_type": ["protein"],
            "chain_length": [100],
            "chain_num_unresolved_residues": [0],
            "chain_is_holo": [True],
            "chain_uniprot_ids": [["P12345"]],
        }
    )

    with pytest.raises(ValueError, match="entry_chains is required for 1abc"):
        entry_views_from_df(rows)

    entry = entry_views_from_df(rows, entry_chains=entry_chains)["1abc"]
    system = next(iter(entry.systems.values()))
    ligand = next(iter(system.ligands.values()))
    scorer = Scorer(
        entries={"1abc": entry},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )

    assert system.receptor_type == "protein+dna"
    assert system.protein_chains_asym_id == ["1.A"]
    assert system.pocket_residue_number_to_index == {"1.A": {10: 9}}
    assert system.proper_num_pocket_residues == 1
    assert ligand.protein_chains_asym_id == ["1.A"]
    assert ligand.pocket_residue_number_to_index == {"1.A": {10: 9}}
    assert ligand.num_interactions == 1
    assert scorer.get_protein_chain_length("1abc", ["1.A", "1.N"]) == 100
    assert entry.chains_for_alignment("holo", "mmseqs") == ["1abc_A"]


def test_na_only_entry_has_no_similarity_chains() -> None:
    rows = pd.DataFrame(
        [
            {
                "entry_pdb_id": "1dna",
                "system_id": "1dna__1__1.N__1.L",
                "system_type": "holo",
                "system_receptor_type": "dna",
                "system_protein_chains_asym_id": ["1.N"],
                "system_protein_chains_auth_id": ["N"],
                "system_protein_chains_length": [20],
                "system_proper_num_pocket_residues": 1,
                "system_proper_num_interactions": 1,
                "system_proper_num_unique_interactions": 1,
                "ligand_id": "1dna__1__1.L",
                "ligand_instance_chain": "1.L",
                "ligand_asym_id": "L",
                "ligand_is_proper": True,
                "ligand_protein_chains_asym_id": ["1.N"],
                "ligand_num_pocket_residues": 1,
                "ligand_num_interactions": 1,
                "ligand_num_unique_interactions": 1,
                "ligand_neighboring_residues": ["1.N_2_1_2"],
                "ligand_interactions": ["1.N_2_type:hydrogen_bonds"],
            }
        ]
    )
    entry_chains = pd.DataFrame(
        {
            "entry_pdb_id": ["1dna"],
            "chain_asym_id": ["N"],
            "chain_auth_id": ["N"],
            "chain_entity_id": ["1"],
            "chain_type": ["polydeoxyribonucleotide"],
            "chain_receptor_type": ["dna"],
            "chain_length": [20],
            "chain_num_unresolved_residues": [0],
            "chain_is_holo": [True],
            "chain_uniprot_ids": [[]],
        }
    )

    entry = entry_views_from_df(rows, entry_chains=entry_chains)["1dna"]
    system = next(iter(entry.systems.values()))

    assert system.receptor_type == "dna"
    assert system.protein_chains_asym_id == []
    assert system.proper_num_pocket_residues == 0
    assert entry.chains["N"].receptor_type == "dna"
    assert entry.chains_for_alignment("holo", "foldseek") == []
    assert entry.chains_for_alignment("holo", "mmseqs") == []


def _ligand(
    ligand_id: str,
    instance_chain: str,
    pocket: dict[str, dict[int, int]],
    interactions: dict[str, dict[int, Counter[str]]],
) -> LigandView:
    return LigandView(
        id=ligand_id,
        pdb_id=ligand_id[:4],
        system_id=f"{ligand_id[:4]}_system",
        instance_chain=instance_chain,
        asym_id=instance_chain.split(".", 1)[1],
        is_proper=True,
        protein_chains_asym_id=list(pocket),
        num_pocket_residues=sum(len(residues) for residues in pocket.values()),
        num_interactions=sum(
            sum(counter.values())
            for chain in interactions.values()
            for counter in chain.values()
        ),
        num_unique_interactions=sum(
            len(counter)
            for chain in interactions.values()
            for counter in chain.values()
        ),
        pocket_residue_number_to_index=pocket,
        interactions_counter=interactions,
    )


def test_ligand_pair_pocket_scores_do_not_use_system_union(tmp_path) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        include_pli_fident=True,
    )
    alignments = {
        ("1.A", "1.X"): pd.DataFrame(
            [
                {
                    "query_selected_residue_numbers": [10, 20],
                    "target_selected_residue_numbers": [110, 120],
                    "selected_residue_identity": bytes([1, 1]),
                }
            ],
            index=["foldseek"],
        )
    }
    query = _ligand(
        "1abc__1__1.B",
        "1.B",
        {"1.A": {10: 9, 20: 19}},
        {
            "1.A": {
                10: Counter({"hydrogen_bond": 1}),
                20: Counter({"hydrophobic": 1}),
            }
        },
    )
    target_full = _ligand(
        "2def__1__1.Y",
        "1.Y",
        {"1.X": {110: 109, 120: 119}},
        {
            "1.X": {
                110: Counter({"hydrogen_bond": 1}),
                120: Counter({"hydrophobic": 1}),
            }
        },
    )
    target_partial = _ligand(
        "2def__1__1.Z",
        "1.Z",
        {"1.X": {110: 109}},
        {"1.X": {110: Counter({"hydrogen_bond": 1})}},
    )

    full_pocket, full_pli = scorer._get_ligand_pair_pocket_pli_scores_for_mapping(
        alignments, query, target_full
    )
    partial_pocket, partial_pli = scorer._get_ligand_pair_pocket_pli_scores_for_mapping(
        alignments, query, target_partial
    )

    assert full_pocket["pocket_qcov_foldseek"] == pytest.approx(1.0)
    assert partial_pocket["pocket_qcov_foldseek"] == pytest.approx(0.5)
    assert full_pli["pli_qcov_foldseek"] == pytest.approx(1.0)
    assert partial_pli["pli_qcov_foldseek"] == pytest.approx(0.5)
    assert full_pli["pli_fident_foldseek"] == pytest.approx(1.0)
    assert partial_pli["pli_fident_foldseek"] == pytest.approx(1.0)


def test_ligand_pair_pocket_scores_ignore_null_compact_maps(tmp_path) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    alignments = {
        ("1.A", "1.X"): pd.DataFrame(
            [
                {
                    "query_selected_residue_numbers": np.nan,
                    "target_selected_residue_numbers": np.nan,
                    "selected_residue_identity": np.nan,
                }
            ],
            index=["foldseek"],
        )
    }
    query = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    target = _ligand("2def__1__1.Y", "1.Y", {"1.X": {110: 109}}, {})

    pocket_scores, pli_scores = scorer._get_ligand_pair_pocket_pli_scores_for_mapping(
        alignments, query, target
    )

    assert pocket_scores == {}
    assert pli_scores == {}


def test_ligand_pocket_scores_report_identity_over_all_pli_residues(tmp_path) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        include_pli_fident=True,
    )
    query = _ligand(
        "1abc__1__1.B",
        "1.B",
        {"1.A": {10: 9, 20: 19, 30: 29}},
        {
            "1.A": {
                10: Counter({"hydrogen_bond": 1}),
                30: Counter({"hydrophobic": 1}),
            }
        },
    )
    alignments = {
        ("1.A", "0.X"): pd.DataFrame(
            [
                {
                    "query_selected_residue_numbers": [10, 20],
                    "target_selected_residue_numbers": [100, 200],
                    "selected_residue_identity": bytes([1, 1]),
                }
            ],
            index=["mmseqs"],
        )
    }

    scores = scorer.get_ligand_pocket_scores(alignments, query)

    assert scores["pocket_fident_mmseqs"] == pytest.approx(2 / 3)
    assert scores["pli_fident_mmseqs"] == pytest.approx(0.5)


def test_ligand_pair_pocket_mapping_maximizes_coverage_before_similarity(
    tmp_path,
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    query = _ligand(
        "1abc__1__1.L",
        "1.L",
        {"1.A": {1: 0}, "1.B": {2: 1}},
        {},
    )
    target = _ligand(
        "2def__1__1.M",
        "1.M",
        {"1.X": {10: 9}, "1.Y": {20: 19}},
        {},
    )
    rows = [
        # Whole-protein similarity favors this direct Foldseek assignment,
        # but neither alignment maps into the corresponding target pocket.
        ("A", "X", "foldseek", 1, 99, 1.0),
        ("B", "Y", "foldseek", 2, 99, 0.9),
        # The lower-scoring swapped assignment covers both pocket residues.
        ("A", "Y", "foldseek", 1, 20, 0.2),
        ("B", "X", "foldseek", 2, 10, 0.1),
        # MMseqs covers only one residue with its preferred assignment.
        ("A", "X", "mmseqs", 1, 10, 0.95),
        ("B", "Y", "mmseqs", 2, 99, 0.9),
    ]
    alignments = pd.DataFrame(
        [
            {
                "query_chain_mapped": query_chain,
                "target_chain_mapped": target_chain,
                "source": source,
                "query_selected_residue_numbers": [query_number],
                "target_selected_residue_numbers": [target_number],
                "selected_residue_identity": bytes([1]),
                "qcov": 1.0,
                "fident": similarity,
                "fident_qcov": similarity,
                "lddt_qcov": similarity,
            }
            for (
                query_chain,
                target_chain,
                source,
                query_number,
                target_number,
                similarity,
            ) in rows
        ]
    ).set_index(["query_chain_mapped", "target_chain_mapped", "source"])
    alignments.sort_index(inplace=True)

    pocket_scores, _, mappings = scorer.get_ligand_pair_pocket_pli_scores(
        alignments,
        query,
        target,
    )

    assert pocket_scores["pocket_qcov_foldseek"] == pytest.approx(1.0)
    assert pocket_scores["pocket_qcov_mmseqs"] == pytest.approx(0.5)
    assert mappings["pocket_qcov_foldseek"] == [
        ("1.A", "1.Y"),
        ("1.B", "1.X"),
    ]


def test_ligand_pair_shape_scores_gate_sdf_access_and_cache(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    query = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    target = _ligand("2def__1__1.Y", "1.Y", {"1.X": {110: 109}}, {})
    resolved: list[str] = []

    def resolve_sdf(_data_dir: Path, ligand: LigandView) -> Path:
        resolved.append(ligand.id)
        return SDF_FILE

    monkeypatch.setattr(scorer, "resolve_ligand_sdf", resolve_sdf)
    shape_calls = 0
    sucos_calls = 0

    def shape_scores(*_args) -> tuple[float, float]:
        nonlocal shape_calls
        shape_calls += 1
        return 0.8, 0.6

    def sucos_score(*_args) -> float:
        nonlocal sucos_calls
        sucos_calls += 1
        return 0.7

    monkeypatch.setattr(scoring_module, "align_molecules", shape_scores)
    monkeypatch.setattr(scoring_module, "get_sucos_score", sucos_score)

    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.0) == {}
    assert (
        scorer.get_ligand_pair_shape_scores(tmp_path, query, target, float("nan")) == {}
    )
    assert resolved == []

    query.is_shape_comparable = False
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == {}
    assert resolved == []
    query.is_shape_comparable = True

    scores = scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5)
    assert set(scores) == {
        "shape",
        "color",
        "sucos_shape",
        "sucos_shape_pocket_qcov",
    }
    assert scores["sucos_shape_pocket_qcov"] == pytest.approx(
        scores["sucos_shape"] * 0.5
    )
    assert all(0 <= value <= 1 for value in scores.values())
    assert resolved == [query.id, target.id]

    # The same canonical ASU ligand pair can occur in several systems. Its
    # expensive base scores are reused while pocket weighting remains local.
    repeated_target = replace(
        target,
        id="2def__2__2.Y",
        system_id="2def__2",
        instance_chain="2.Y",
    )
    repeated_scores = scorer.get_ligand_pair_shape_scores(
        tmp_path, query, repeated_target, 0.25
    )
    assert repeated_scores["shape"] == scores["shape"]
    assert repeated_scores["color"] == scores["color"]
    assert repeated_scores["sucos_shape"] == scores["sucos_shape"]
    assert repeated_scores["sucos_shape_pocket_qcov"] == pytest.approx(0.7 * 0.25)
    assert resolved == [query.id, target.id]
    assert shape_calls == 1
    assert sucos_calls == 1


def test_ligand_pair_shape_scores_drop_nonfinite_results(tmp_path, monkeypatch) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        ligand_sdf_resolver=lambda _ligand: SDF_FILE,
    )
    query = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    target = _ligand("2def__1__1.Y", "1.Y", {"1.X": {110: 109}}, {})

    monkeypatch.setattr(
        scoring_module, "align_molecules", lambda *_args: (1.0, float("nan"))
    )
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == {}

    scorer._ligand_shape_score_cache.clear()
    monkeypatch.setattr(scoring_module, "align_molecules", lambda *_args: (0.8, 0.6))
    monkeypatch.setattr(scoring_module, "get_sucos_score", lambda *_args: float("nan"))
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == {
        "shape": 0.8,
        "color": 0.6,
    }


def test_ligand_sdf_resolver_uses_only_canonical_entry_path(tmp_path) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    ligand = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    canonical = tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files" / "B.sdf"
    legacy = (
        tmp_path / "raw_entries" / "ab" / ligand.system_id / "ligand_files" / "1.B.sdf"
    )
    legacy.parent.mkdir(parents=True)
    legacy.write_text("legacy")

    assert scorer.resolve_ligand_sdf(tmp_path, ligand) is None
    canonical.parent.mkdir(parents=True)
    canonical.write_text("canonical")
    assert scorer.resolve_ligand_sdf(tmp_path, ligand) == canonical


def test_ligand_molecules_bulk_load_from_packed_sdf_parquet(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    ligand = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    archive = tmp_path / "ligand_archives" / "ab.parquet"
    archive.parent.mkdir()
    pd.DataFrame(
        {
            "pdb_id": ["1abc"],
            "ligand_asym_id": ["B"],
            "sdf": [SDF_FILE.read_bytes()],
        }
    ).to_parquet(archive, index=False)
    monkeypatch.setattr(
        scorer,
        "resolve_ligand_sdf",
        lambda *_args: pytest.fail("packed SDF should avoid raw-file resolution"),
    )

    scorer._preload_packed_ligand_sdfs(tmp_path, [ligand])
    molecule = scorer._get_ligand_mol(tmp_path, ligand)

    assert molecule is not None
    assert molecule.GetNumConformers() == 1
    assert (ligand.pdb_id, ligand.asym_id) not in scorer._ligand_sdf_block_cache


def test_get_score_df_loads_entries_from_ingest_data_dir(tmp_path, monkeypatch) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    alignment_path = scorer.db_dir / "holo_foldseek" / "aln" / "1abc.parquet"
    alignment_path.parent.mkdir(parents=True)
    pd.DataFrame({"target_pdb_id": ["2def"]}).to_parquet(alignment_path, index=False)
    calls: list[tuple[set[str], Path]] = []

    def fake_load_entry_views(*, pdb_ids, data_dir):
        calls.append((set(pdb_ids), data_dir))
        return {}

    monkeypatch.setattr(scoring_module, "load_entry_views", fake_load_entry_views)
    monkeypatch.setattr(
        scorer,
        "map_alignment_df",
        lambda *_args: pd.DataFrame({"mapped": [True]}),
    )
    monkeypatch.setattr(
        scorer,
        "aggregate_scores",
        lambda *_args, **_kwargs: pd.DataFrame(),
    )

    scorer.get_score_df(tmp_path, "1abc", "holo")

    assert calls == [({"1abc", "2def"}, tmp_path)]


def test_get_score_df_loads_mapped_targets_when_mapping_is_separate(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    mapped_path = scorer.db_dir / "holo_foldseek" / "mapped_aln" / "1abc.parquet"
    mapped_path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "query_entry": ["1abc"] * 3,
            "target_entry": ["2def", "3ghi", "2def"],
            "query_chain_mapped": ["A"] * 3,
            "target_chain_mapped": ["B"] * 3,
            "source": ["foldseek"] * 3,
            "qcov": [1.0] * 3,
        }
    ).set_index(
        [
            "query_entry",
            "target_entry",
            "query_chain_mapped",
            "target_chain_mapped",
            "source",
        ]
    ).to_parquet(mapped_path, index=True)
    calls: list[tuple[set[str], Path]] = []

    def fake_load_entry_views(*, pdb_ids, data_dir):
        calls.append((set(pdb_ids), data_dir))
        return {pdb_id: object() for pdb_id in pdb_ids}

    monkeypatch.setattr(scoring_module, "load_entry_views", fake_load_entry_views)
    monkeypatch.setattr(
        scorer,
        "aggregate_scores",
        lambda *_args, **_kwargs: pd.DataFrame(),
    )

    scorer.get_score_df(tmp_path, "1abc", "holo", map_alignments=False)

    assert calls == [({"1abc", "2def", "3ghi"}, tmp_path)]


def test_get_score_df_records_and_reraises_aggregation_failure(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={"1abc": object()},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )

    def fail(*_args, **_kwargs):
        raise RuntimeError("broken aggregation")

    monkeypatch.setattr(scorer, "aggregate_scores", fail)

    with pytest.raises(RuntimeError, match="broken aggregation"):
        scorer.get_score_df(tmp_path, "1abc", "holo", map_alignments=False)

    failure = (
        tmp_path / "scratch" / "scores" / "aggregate_scores_failures" / "holo_1abc.txt"
    )
    assert failure.read_text() == (
        "RuntimeError('broken aggregation'): broken aggregation"
    )


def test_get_score_df_atomically_replaces_stale_cache_and_marks_empty_completion(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={"1abc": object()},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    output = scorer.db_dir / "search_db=holo" / "1abc.parquet"
    output.parent.mkdir(parents=True)
    pd.DataFrame({"legacy": [True]}).to_parquet(output, index=False)
    calls = 0

    def no_scores(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        return None

    monkeypatch.setattr(scorer, "aggregate_scores", no_scores)
    scratch = tmp_path / "node-scratch"

    assert (
        scorer.get_score_df(
            tmp_path,
            "1abc",
            "holo",
            overwrite=False,
            map_alignments=False,
            scratch_dir=scratch,
        )
        == output
    )
    assert pd.read_parquet(output).empty
    assert set(scoring_module.pq.read_schema(output).names) == set(
        PROTEIN_SIMILARITY_SCHEMA.names
    )
    assert not list(scratch.glob("*.parquet"))

    scorer.get_score_df(
        tmp_path,
        "1abc",
        "holo",
        overwrite=False,
        map_alignments=False,
        scratch_dir=scratch,
    )
    assert calls == 1


def test_get_score_df_retains_and_tracks_requested_metrics(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={"1abc": object()},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    calls = 0

    def scores(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        rows = []
        for metric in ["pocket_fident", "protein_qcov_weighted_sum"]:
            rows.append(
                {
                    "query_system": "1abc__1__1.A__1.L",
                    "query_ligand_id": "1abc__1__1.L",
                    "target_system": "2def_A",
                    "target_ligand_id": None,
                    "protein_mapping": "1.A:0.A",
                    "mapping": "1.A:0.A",
                    "protein_mapper": "foldseek",
                    "source": "foldseek",
                    "metric": metric,
                    "similarity": 95,
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(scorer, "aggregate_scores", scores)
    output = scorer.get_score_df(
        tmp_path,
        "1abc",
        "apo",
        overwrite=False,
        map_alignments=False,
        score_metrics={"pocket_fident"},
    )
    assert pd.read_parquet(output)["metric"].astype(str).tolist() == ["pocket_fident"]
    assert (pq.read_schema(output).metadata or {})[
        scoring_module.SCORE_METRICS_METADATA_KEY
    ] == scoring_module.score_metrics_metadata({"pocket_fident"})

    scorer.get_score_df(
        tmp_path,
        "1abc",
        "apo",
        overwrite=False,
        map_alignments=False,
        score_metrics={"pocket_fident"},
    )
    assert calls == 1

    scorer.get_score_df(
        tmp_path,
        "1abc",
        "apo",
        overwrite=False,
        map_alignments=False,
        score_metrics={"protein_qcov_weighted_sum"},
    )
    assert calls == 2
    assert pd.read_parquet(output)["metric"].astype(str).tolist() == [
        "protein_qcov_weighted_sum"
    ]


def test_get_score_df_defers_ligand_3d_and_writes_full_precision_candidates(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={"1abc": object()},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    calls = 0

    def candidate_scores(*_args, **kwargs):
        nonlocal calls
        calls += 1
        assert kwargs["data_dir"] is None
        assert kwargs["include_holo_protein_scores"] is False
        kwargs["ligand_3d_candidates"].append(
            {
                "query_system": "1abc_system",
                "query_ligand_id": "1abc__1__1.B",
                "query_entry": "1abc",
                "query_ligand_asym_id": "B",
                "target_system": "2def_system",
                "target_ligand_id": "2def__1__1.Y",
                "target_entry": "2def",
                "target_ligand_asym_id": "Y",
                "protein_mapping": "1.A:1.X",
                "protein_mapper": "foldseek",
                "pocket_qcov": 2 / 3,
            }
        )
        kwargs["ligand_pair_scores"].append(
            {
                "query_system": "1abc_system",
                "query_ligand_id": "1abc__1__1.B",
                "query_entry": "1abc",
                "query_ligand_asym_id": "B",
                "target_system": "2def_system",
                "target_ligand_id": "2def__1__1.Y",
                "target_entry": "2def",
                "target_ligand_asym_id": "Y",
                "pocket_qcov": 27,
                "pocket_fident_qcov": 19,
                "pli_qcov": 11,
            }
        )
        return None

    monkeypatch.setattr(scorer, "aggregate_scores", candidate_scores)

    output = scorer.get_score_df(
        tmp_path,
        "1abc",
        "holo",
        overwrite=False,
        map_alignments=False,
        defer_ligand_3d=True,
    )
    candidates = (
        scorer.scores_dir
        / "ligand_3d_candidates"
        / "search_db=holo"
        / "shard=ab"
        / "1abc.parquet"
    )

    assert output.is_file()
    assert pd.read_parquet(candidates)["pocket_qcov"].item() == pytest.approx(2 / 3)
    ligand_pair_scores = pd.read_parquet(
        scorer.scores_dir
        / "ligand_pair_scores"
        / "search_db=holo"
        / "shard=ab"
        / "1abc.parquet"
    )
    assert ligand_pair_scores[
        ["pocket_qcov", "pocket_fident_qcov", "pli_qcov"]
    ].to_dict("records") == [
        {"pocket_qcov": 27, "pocket_fident_qcov": 19, "pli_qcov": 11}
    ]
    assert (scoring_module.pq.read_schema(output).metadata or {}).get(
        b"plinder.ligand_3d"
    ) == b"deferred"
    assert (scoring_module.pq.read_schema(output).metadata or {}).get(
        scoring_module.HOLO_PROTEIN_SCORES_METADATA_KEY
    ) == b"excluded"
    scorer.get_score_df(
        tmp_path,
        "1abc",
        "holo",
        overwrite=False,
        map_alignments=False,
        defer_ligand_3d=True,
    )
    assert calls == 1

    scorer.minimum_thresholds["protein_lddt_weighted_sum"] = 0.2
    scorer.get_score_df(
        tmp_path,
        "1abc",
        "holo",
        overwrite=False,
        map_alignments=False,
        defer_ligand_3d=True,
    )
    assert calls == 2


def test_repair_score_df_targets_replaces_only_affected_target_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from types import SimpleNamespace

    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "dbs" / "subdbs",
        scores_dir=tmp_path / "scores",
    )
    score_path = scorer.db_dir / "search_db=holo" / "1abc.parquet"
    score_path.parent.mkdir(parents=True)

    def score_row(target_system: str, similarity: int) -> dict[str, object]:
        return {
            "query_system": "1abc__1__1.A__1.B",
            "query_ligand_id": "1abc__1__1.B",
            "target_system": target_system,
            "target_ligand_id": f"{target_system.split('__', 1)[0]}__1__1.Y",
            "protein_mapping": "1.A:1.X",
            "mapping": None,
            "protein_mapper": "foldseek",
            "source": "foldseek",
            "metric": "pocket_qcov",
            "similarity": similarity,
        }

    pd.DataFrame(
        [
            score_row("2def__1__1.X__1.Y", 60),
            score_row("3ghi__1__1.X__1.Y", 70),
            {
                **score_row("3ghi__1__1.X__1.Y", 90),
                "metric": "protein_fident_weighted_sum",
            },
        ]
    ).to_parquet(
        score_path,
        index=False,
        schema=PROTEIN_SIMILARITY_SCHEMA.with_metadata(
            {b"plinder.ligand_3d": b"deferred"}
        ),
    )
    candidate_path = (
        scorer.scores_dir / "ligand_3d_candidates/search_db=holo/shard=ab/1abc.parquet"
    )
    candidate_path.parent.mkdir(parents=True)

    def candidate_row(target_entry: str, target_system: str) -> dict[str, object]:
        return {
            "query_system": "1abc__1__1.A__1.B",
            "query_ligand_id": "1abc__1__1.B",
            "query_entry": "1abc",
            "query_ligand_asym_id": "B",
            "target_system": target_system,
            "target_ligand_id": f"{target_entry}__1__1.Y",
            "target_entry": target_entry,
            "target_ligand_asym_id": "Y",
            "protein_mapping": "1.A:1.X",
            "protein_mapper": "foldseek",
            "pocket_qcov": 0.6,
        }

    pd.DataFrame(
        [
            candidate_row("2def", "2def__1__1.X__1.Y"),
            candidate_row("3ghi", "3ghi__1__1.X__1.Y"),
        ]
    ).to_parquet(
        candidate_path,
        index=False,
        schema=scoring_module.schemas.LIGAND_3D_CANDIDATE_SCHEMA,
    )
    ligand_pair_score_path = (
        scorer.scores_dir / "ligand_pair_scores/search_db=holo/shard=ab/1abc.parquet"
    )
    ligand_pair_score_path.parent.mkdir(parents=True)

    def ligand_pair_row(target_entry: str, target_system: str) -> dict[str, object]:
        row = candidate_row(target_entry, target_system)
        return {
            key: value
            for key, value in row.items()
            if key not in {"protein_mapping", "protein_mapper", "pocket_qcov"}
        } | {
            "pocket_qcov": 60,
            "pocket_fident_qcov": 50,
            "pli_qcov": 40,
        }

    pd.DataFrame(
        [
            ligand_pair_row("2def", "2def__1__1.X__1.Y"),
            ligand_pair_row("3ghi", "3ghi__1__1.X__1.Y"),
        ]
    ).to_parquet(
        ligand_pair_score_path,
        index=False,
        schema=scoring_module.schemas.LIGAND_PAIR_SCORE_SCHEMA,
    )
    entries = {
        "1abc": SimpleNamespace(systems={"1abc__1__1.A__1.B": object()}),
        "2def": SimpleNamespace(systems={"2def__2__1.X__1.Y": object()}),
    }
    monkeypatch.setattr(
        scoring_module,
        "load_entry_views",
        lambda *, pdb_ids, data_dir: {
            pdb_id: entries[pdb_id] for pdb_id in pdb_ids if pdb_id in entries
        },
    )

    def repaired_scores(*_args, **kwargs):
        assert kwargs["target_system_ids"] == {"2def__2__1.X__1.Y"}
        kwargs["ligand_3d_candidates"].append(
            candidate_row("2def", "2def__2__1.X__1.Y")
        )
        kwargs["ligand_pair_scores"].append(
            ligand_pair_row("2def", "2def__2__1.X__1.Y")
        )
        return pd.DataFrame([score_row("2def__2__1.X__1.Y", 80)])

    monkeypatch.setattr(scorer, "aggregate_scores", repaired_scores)

    scorer.repair_score_df_targets(
        tmp_path,
        "1abc",
        affected_target_entries={"2def"},
        scratch_dir=tmp_path / "scratch",
    )

    repaired = pd.read_parquet(score_path)
    assert set(repaired["target_system"]) == {
        "2def__2__1.X__1.Y",
        "3ghi__1__1.X__1.Y",
    }
    assert not repaired["metric"].str.startswith("protein_").any()
    candidates = pd.read_parquet(candidate_path)
    assert set(candidates["target_system"]) == {
        "2def__2__1.X__1.Y",
        "3ghi__1__1.X__1.Y",
    }
    assert set(pd.read_parquet(ligand_pair_score_path)["target_system"]) == {
        "2def__2__1.X__1.Y",
        "3ghi__1__1.X__1.Y",
    }
    assert scoring_module.SCORE_THRESHOLDS_METADATA_KEY not in (
        scoring_module.pq.read_schema(score_path).metadata or {}
    )


def test_repair_score_df_targets_can_create_bounded_query_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from types import SimpleNamespace

    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "dbs" / "subdbs",
        scores_dir=tmp_path / "scores",
    )
    entries = {
        "1abc": SimpleNamespace(systems={"1abc__1__1.A__1.B": object()}),
        "2def": SimpleNamespace(systems={"2def__1__1.X__1.Y": object()}),
    }
    monkeypatch.setattr(
        scoring_module,
        "load_entry_views",
        lambda *, pdb_ids, data_dir: {
            pdb_id: entries[pdb_id] for pdb_id in pdb_ids if pdb_id in entries
        },
    )
    score_row = {
        "query_system": "1abc__1__1.A__1.B",
        "query_ligand_id": "1abc__1__1.B",
        "target_system": "2def__1__1.X__1.Y",
        "target_ligand_id": "2def__1__1.Y",
        "protein_mapping": "1.A:1.X",
        "mapping": None,
        "protein_mapper": "foldseek",
        "source": "foldseek",
        "metric": "pocket_qcov",
        "similarity": 75,
    }
    candidate_row = {
        "query_system": "1abc__1__1.A__1.B",
        "query_ligand_id": "1abc__1__1.B",
        "query_entry": "1abc",
        "query_ligand_asym_id": "B",
        "target_system": "2def__1__1.X__1.Y",
        "target_ligand_id": "2def__1__1.Y",
        "target_entry": "2def",
        "target_ligand_asym_id": "Y",
        "protein_mapping": "1.A:1.X",
        "protein_mapper": "foldseek",
        "pocket_qcov": 0.75,
    }

    def repaired_scores(*_args, **kwargs):
        assert kwargs["query_system_ids"] == {"1abc__1__1.A__1.B"}
        assert kwargs["query_ligand_ids"] == {"1abc__1__1.B"}
        assert kwargs["target_system_ids"] == {"2def__1__1.X__1.Y"}
        assert kwargs["target_ligand_ids"] == {"2def__1__1.Y"}
        kwargs["ligand_3d_candidates"].append(candidate_row)
        kwargs["ligand_pair_scores"].append(
            {
                key: value
                for key, value in candidate_row.items()
                if key not in {"protein_mapping", "protein_mapper", "pocket_qcov"}
            }
            | {
                "pocket_qcov": 75,
                "pocket_fident_qcov": 65,
                "pli_qcov": 55,
            }
        )
        return pd.DataFrame([score_row])

    monkeypatch.setattr(scorer, "aggregate_scores", repaired_scores)
    output = scorer.repair_score_df_targets(
        tmp_path,
        "1abc",
        affected_target_entries={"2def"},
        scratch_dir=tmp_path / "scratch",
        allow_missing=True,
        query_system_ids={"1abc__1__1.A__1.B"},
        query_ligand_ids={"1abc__1__1.B"},
        target_system_ids={"2def__1__1.X__1.Y"},
        target_ligand_ids={"2def__1__1.Y"},
    )

    assert pd.read_parquet(output)["similarity"].tolist() == [75]
    candidate_path = (
        tmp_path / "scores/ligand_3d_candidates/search_db=holo/shard=ab/1abc.parquet"
    )
    assert pd.read_parquet(candidate_path)["pocket_qcov"].tolist() == [0.75]
    ligand_pair_score_path = (
        tmp_path / "scores/ligand_pair_scores/search_db=holo/shard=ab/1abc.parquet"
    )
    assert pd.read_parquet(ligand_pair_score_path)["pli_qcov"].tolist() == [55]


def test_map_alignment_files_replaces_stale_schema_without_force(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    raw_path = scorer.db_dir / "holo_foldseek" / "aln" / "1abc.parquet"
    raw_path.parent.mkdir(parents=True)
    pd.DataFrame({"target_pdb_id": ["2def"]}).to_parquet(raw_path, index=False)
    mapped_path = scorer.db_dir / "holo_foldseek" / "mapped_aln" / "1abc.parquet"
    mapped_path.parent.mkdir(parents=True)
    pd.DataFrame({"qrnum": [[(0, 1)]]}).to_parquet(mapped_path, index=False)
    calls = []

    monkeypatch.setattr(scoring_module, "load_entry_views", lambda **_kwargs: {})
    monkeypatch.setattr(
        scorer,
        "map_alignment_df",
        lambda *args: calls.append(args) or pd.DataFrame({"mapped": [True]}),
    )

    assert scorer.map_alignment_files(
        tmp_path,
        "1abc",
        "holo",
        overwrite=False,
        scratch_dir=tmp_path / "scratch",
    ) == [mapped_path]
    assert len(calls) == 1
    assert pd.read_parquet(mapped_path)["mapped"].tolist() == [True]


def test_feature_map_score_handles_molecules_without_features() -> None:
    helium = Chem.MolFromSmiles("[He]")
    assert helium is not None
    conformer = Chem.Conformer(helium.GetNumAtoms())
    conformer.SetAtomPosition(0, (0.0, 0.0, 0.0))
    helium.AddConformer(conformer)

    assert get_feature_map_score(helium, helium) == 0.0


def test_cached_reference_pharmacophore_features_preserve_score() -> None:
    reference = Chem.MolFromMolFile(str(SDF_FILE))
    mobile = Chem.MolFromMolFile(str(SDF_FILE))
    assert reference is not None
    assert mobile is not None
    conformer = mobile.GetConformer()
    for atom_index in range(mobile.GetNumAtoms()):
        position = conformer.GetAtomPosition(atom_index)
        conformer.SetAtomPosition(
            atom_index,
            (-position.y + 0.7, position.x - 0.4, position.z + 0.2),
        )

    expected = get_feature_map_score(reference, mobile)
    reference_features = scoring_module._pharmacophore_features(reference)
    observed = get_feature_map_score(
        reference,
        mobile,
        scoring_module.FeatMaps.FeatMapScoreMode.All,
        reference_features,
        None,
    )

    assert observed == pytest.approx(expected, abs=1e-12)


def test_align_molecules_passes_unmodified_mobile_to_shapealign(monkeypatch) -> None:
    reference = Chem.MolFromMolFile(str(SDF_FILE))
    mobile = Chem.MolFromMolFile(str(SDF_FILE))
    assert reference is not None
    assert mobile is not None
    conformer = mobile.GetConformer()
    for atom_index in range(mobile.GetNumAtoms()):
        position = conformer.GetAtomPosition(atom_index)
        conformer.SetAtomPosition(
            atom_index,
            (-position.y + 0.7, position.x - 0.4, position.z + 0.2),
        )
    expected_positions = np.asarray(conformer.GetPositions()).copy()
    observed_positions = None

    def observe_shapealign(_reference, observed_mobile, *_args):
        nonlocal observed_positions
        observed_positions = np.asarray(
            observed_mobile.GetConformer().GetPositions()
        ).copy()
        return 0.8, 0.6

    monkeypatch.setattr(scoring_module.rdShapeAlign, "AlignMol", observe_shapealign)

    assert scoring_module.align_molecules(reference, mobile) == (0.8, 0.6)
    assert observed_positions is not None
    assert np.array_equal(observed_positions, expected_positions)


def test_ligand_feature_detection_is_cached_per_canonical_ligand(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        ligand_sdf_resolver=lambda _ligand: SDF_FILE,
    )
    query = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    target = _ligand("2def__1__1.Y", "1.Y", {"1.X": {110: 109}}, {})
    original = scoring_module._pharmacophore_features
    feature_calls = 0

    def observed_features(molecule):
        nonlocal feature_calls
        feature_calls += 1
        return original(molecule)

    monkeypatch.setattr(scoring_module, "_pharmacophore_features", observed_features)
    monkeypatch.setattr(scoring_module, "align_molecules", lambda *_args: (0.8, 0.6))
    monkeypatch.setattr(scoring_module, "get_sucos_score", lambda *_args: 0.7)

    scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5)
    scorer.get_ligand_pair_shape_scores(tmp_path, target, query, 0.5)

    assert feature_calls == 2


def test_canonical_ligand_pair_scoring_deduplicates_pairs(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        ligand_sdf_resolver=lambda _ligand: SDF_FILE,
        shape_score_threads=2,
    )
    pairs = pd.DataFrame(
        {
            "query_entry": ["1abc", "1abc"],
            "query_ligand_asym_id": ["B", "B"],
            "target_entry": ["2def", "2def"],
            "target_ligand_asym_id": ["Y", "Y"],
        }
    )
    align_calls = 0

    def align_once(*_args):
        nonlocal align_calls
        align_calls += 1
        return 0.8, 0.6

    monkeypatch.setattr(scoring_module, "align_molecules", align_once)
    monkeypatch.setattr(scoring_module, "get_sucos_score", lambda *_args: 0.7)

    scores = scorer.score_canonical_ligand_pairs(tmp_path, pairs)

    assert align_calls == 1
    assert scores.to_dict("records") == [
        {
            "query_entry": "1abc",
            "query_ligand_asym_id": "B",
            "target_entry": "2def",
            "target_ligand_asym_id": "Y",
            "shape": 0.8,
            "color": 0.6,
            "sucos_shape": 0.7,
        }
    ]


def test_shape_scoring_uses_custom_radius_and_preserves_hem_iron() -> None:
    molecule = Chem.MolFromMolFile(str(HEM_SDF_FILE))
    assert molecule is not None
    assert any(atom.GetAtomicNum() == 26 for atom in molecule.GetAtoms())

    prepared = prepare_ligand_for_3d_scoring(molecule)

    assert prepared is not None
    assert any(atom.GetAtomicNum() == 26 for atom in prepared.GetAtoms())
    shape, color = scoring_module.align_molecules(
        Chem.Mol(prepared), Chem.Mol(prepared)
    )
    assert shape == pytest.approx(1.0)
    assert color == pytest.approx(1.0)
    assert scoring_module.get_sucos_score(prepared, prepared) == pytest.approx(1.0)


def test_sdf_loaders_fall_back_to_unsanitized_bad_valence(tmp_path) -> None:
    editable = Chem.RWMol()
    nitrogen = editable.AddAtom(Chem.Atom(7))
    carbons = [editable.AddAtom(Chem.Atom(6)) for _ in range(4)]
    for carbon in carbons:
        editable.AddBond(nitrogen, carbon, Chem.BondType.SINGLE)
    molecule = editable.GetMol()
    molecule.UpdatePropertyCache(strict=False)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, (float(index), 0.0, 0.0))
    molecule.AddConformer(conformer)
    block = Chem.MolToMolBlock(molecule)
    sdf_file = tmp_path / "bad-valence.sdf"
    sdf_file.write_text(block + "\n$$$$\n")

    assert Chem.MolFromMolBlock(block) is None
    loaded_block = scoring_module.load_sdf_molecule_block(
        (block + "\n$$$$\n").encode(), label="bad-valence"
    )
    loaded_file = scoring_module.load_sdf_molecule(sdf_file)

    assert loaded_block is not None
    assert loaded_file is not None
    assert loaded_block.GetNumConformers() == 1
    assert loaded_file.GetNumConformers() == 1


def test_ligand_3d_score_ability_uses_canonical_sdf_and_caches_by_path(
    tmp_path, monkeypatch
) -> None:
    valid_sdf = tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files" / "B.sdf"
    valid_sdf.parent.mkdir(parents=True)
    valid_sdf.write_bytes(SDF_FILE.read_bytes())
    calls: list[Path] = []
    original = scoring_module.is_ligand_shape_comparable

    def observed(sdf_file: Path) -> bool:
        calls.append(sdf_file)
        return original(sdf_file)

    monkeypatch.setattr(scoring_module, "is_ligand_shape_comparable", observed)
    ligands = pd.DataFrame(
        {
            "pdb_id": ["1abc", "1abc", "1abc"],
            "ligand_asym_id": ["B", "B", "C"],
            "ligand_id": ["first", "second", "missing"],
        }
    )

    annotated = annotate_ligand_3d_score_ability(ligands, data_dir=tmp_path)

    assert annotated["ligand_is_shape_comparable"].tolist() == [True, True, False]
    assert calls == [valid_sdf, valid_sdf.with_name("C.sdf")]


def test_ligand_3d_score_ability_reuses_success_for_same_molecular_graph(
    tmp_path, monkeypatch
) -> None:
    first_sdf = tmp_path / "first.sdf"
    second_sdf = tmp_path / "second.sdf"
    first_sdf.write_bytes(SDF_FILE.read_bytes())
    second_sdf.write_bytes(SDF_FILE.read_bytes())
    scoring_module._LIGAND_SHAPE_COMPARABILITY_CACHE.clear()
    align_calls = 0
    sucos_calls = 0
    original_align = scoring_module.align_molecules
    original_sucos = scoring_module.get_sucos_score

    def observed_align(*args, **kwargs):
        nonlocal align_calls
        align_calls += 1
        return original_align(*args, **kwargs)

    def observed_sucos(*args, **kwargs):
        nonlocal sucos_calls
        sucos_calls += 1
        return original_sucos(*args, **kwargs)

    monkeypatch.setattr(scoring_module, "align_molecules", observed_align)
    monkeypatch.setattr(scoring_module, "get_sucos_score", observed_sucos)

    try:
        assert scoring_module.is_ligand_shape_comparable(first_sdf)
        assert scoring_module.is_ligand_shape_comparable(second_sdf)
    finally:
        scoring_module._LIGAND_SHAPE_COMPARABILITY_CACHE.clear()

    assert align_calls == 1
    assert sucos_calls == 1


def test_cofactor_similarity_uses_ccd_reference_fingerprints() -> None:
    from rdkit import DataStructs

    from plinder.core.structure.smallmols_similarity import mol2morgan_fp

    smiles = ["CCO", "c1ccccc1"]
    unique_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1],
            "ligand_rdkit_canonical_smiles": smiles,
            "fingerprint": [
                DataStructs.BitVectToBinaryText(mol2morgan_fp(value, nbits=1024))
                for value in smiles
            ],
        }
    )

    annotated = annotate_cofactor_similarity(
        unique_ligands,
        cofactor_smiles={"COF": "CCO"},
        threshold=90,
    ).set_index("ligand_smiles_id")

    assert annotated.loc[0, "ligand_max_cofactor_similarity"] == 100
    assert annotated.loc[0, "ligand_most_similar_cofactor"] == "COF"
    assert bool(annotated.loc[0, "ligand_is_cofactor_like"])
    assert not bool(annotated.loc[1, "ligand_is_cofactor_like"])


def test_ligand_scores_use_bulk_tanimoto_for_unique_smiles(
    tmp_path, monkeypatch
) -> None:
    from rdkit import DataStructs

    from plinder.core.structure.smallmols_similarity import mol2morgan_fp

    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    smiles = ["CCO", "CCN", "c1ccccc1"]
    fingerprint_table = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1, 2],
            "ligand_rdkit_canonical_smiles": smiles,
            "fingerprint": [
                DataStructs.BitVectToBinaryText(mol2morgan_fp(value, nbits=1024))
                for value in smiles
            ],
        }
    )
    write_ecfp4_fingerprint_table(
        fingerprint_table,
        fingerprint_dir / "ligands_per_smiles.parquet",
    )
    calls = 0
    original_bulk = scoring_module.DataStructs.BulkTanimotoSimilarity

    def observed_bulk(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_bulk(*args, **kwargs)

    monkeypatch.setattr(
        scoring_module.DataStructs,
        "BulkTanimotoSimilarity",
        observed_bulk,
    )
    output_path = tmp_path / "scores.parquet"

    ligand_scores(
        ligand_ids=[0, 1, 2],
        data_dir=tmp_path,
        output_path=output_path,
        minimum_similarity=30,
    )
    scores = pd.read_parquet(output_path)

    assert calls == 3
    assert set(
        scores.loc[
            scores["query_ligand_id"] == scores["target_ligand_id"],
            "query_ligand_id",
        ]
    ) == {0, 1, 2}
    assert not (
        (scores["query_ligand_id"] == 0) & (scores["target_ligand_id"] == 2)
    ).any()
    assert scores.columns.tolist() == [
        "query_ligand_id",
        "target_ligand_id",
        "tanimoto_similarity_ecfp4_1024",
    ]
    assert (
        scoring_module.pq.read_schema(output_path).metadata
        == scoring_module.ECFP4_PARQUET_METADATA
    )


def test_mhfp6_scores_use_minhash_jaccard_on_shared_node_universe(tmp_path) -> None:
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    # SMILES 0 and 3 are identical, so their MinHash Jaccard must be exactly 1.0;
    # the alkane/benzene pair must fall below the 30% threshold and be dropped.
    smiles = ["CCO", "CCN", "c1ccccc1", "CCO"]
    unique_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": np.arange(len(smiles), dtype=np.int32),
            "ligand_rdkit_canonical_smiles": smiles,
        }
    )
    scoring_module.write_mhfp6_fingerprints(unique_ligands, fingerprint_dir)

    mhfp6_path = fingerprint_dir / scoring_module.MHFP6_FINGERPRINT_FILE
    assert mhfp6_path.is_file()
    # the fingerprint table carries the MHFP6 provenance keys (alongside pandas
    # metadata), matching how the scorer validates it
    fingerprint_metadata = scoring_module.pq.read_schema(mhfp6_path).metadata
    assert all(
        fingerprint_metadata.get(key) == value
        for key, value in scoring_module.MHFP6_PARQUET_METADATA.items()
    )

    output_path = tmp_path / "mhfp6_scores.parquet"
    scoring_module.mhfp6_ligand_scores(
        ligand_ids=[0, 1, 2, 3],
        data_dir=tmp_path,
        output_path=output_path,
        minimum_similarity=30,
    )
    scores = pd.read_parquet(output_path)

    assert scores.columns.tolist() == [
        "query_ligand_id",
        "target_ligand_id",
        scoring_module.MHFP6_METRIC,
    ]
    # every node is self-identical, and the two identical SMILES match at 100%
    self_edges = scores[scores["query_ligand_id"] == scores["target_ligand_id"]]
    assert set(self_edges["query_ligand_id"]) == {0, 1, 2, 3}
    assert (self_edges[scoring_module.MHFP6_METRIC] == 100.0).all()
    identical = scores[
        (scores["query_ligand_id"] == 0) & (scores["target_ligand_id"] == 3)
    ][scoring_module.MHFP6_METRIC]
    assert identical.tolist() == [100.0]
    # dissimilar ethanol/benzene pair is filtered by the minimum-similarity gate
    assert not (
        (scores["query_ligand_id"] == 0) & (scores["target_ligand_id"] == 2)
    ).any()
    assert (
        scoring_module.pq.read_schema(output_path).metadata
        == scoring_module.MHFP6_PARQUET_METADATA
    )


def test_mhfp6_scores_reject_non_mhfp6_fingerprint_metadata(tmp_path) -> None:
    # A stale/mislabelled fingerprint table must fail loudly rather than silently
    # score the wrong fingerprint.
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    table = pd.DataFrame(
        {
            "ligand_smiles_id": np.array([0], dtype=np.int32),
            "ligand_rdkit_canonical_smiles": ["CCO"],
            "mhfp6": [b"\x00" * 4],
        }
    )
    scoring_module.pq.write_table(
        scoring_module.pa.Table.from_pandas(table, preserve_index=False),
        fingerprint_dir / scoring_module.MHFP6_FINGERPRINT_FILE,
    )
    with pytest.raises(ValueError, match="MHFP6"):
        scoring_module.mhfp6_ligand_scores(
            ligand_ids=[0],
            data_dir=tmp_path,
            output_path=tmp_path / "out.parquet",
        )


def test_ligand_similarity_annotations_exclude_fingerprint_bytes() -> None:
    unique_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1, 2],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCN", "c1ccccc1"],
            "fingerprint": [b"", b"", b""],
            "ligand_is_cofactor_like": [True, False, False],
        }
    )

    annotations = build_ligand_similarity_annotations(
        unique_ligands=unique_ligands,
    )

    assert annotations.columns.tolist() == [
        "ligand_smiles_id",
        "ligand_rdkit_canonical_smiles",
        "ligand_is_cofactor_like",
    ]
    assert annotations["ligand_smiles_id"].tolist() == [0, 1, 2]


def test_annotate_ligand_similarity_requires_complete_mhfp6_shards(tmp_path):
    from plinder.core.utils.schemas import TANIMOTO_SCORE_SCHEMA

    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    smiles = ["CCO", "CCN", "c1ccccc1"]
    unique_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": np.arange(len(smiles), dtype=np.int32),
            "ligand_rdkit_canonical_smiles": smiles,
            "fingerprint": [b"", b"", b""],
            "ligand_is_cofactor_like": [False, False, False],
        }
    )
    unique_ligands.to_parquet(
        fingerprint_dir / "ligands_per_smiles.parquet", index=False
    )
    ecfp4_dir = tmp_path / "ligand_scores"
    ecfp4_dir.mkdir()
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "query_ligand_id": node,
                    "target_ligand_id": node,
                    "tanimoto_similarity_ecfp4_1024": 100.0,
                }
                for node in range(len(smiles))
            ],
            schema=TANIMOTO_SCORE_SCHEMA,
        ),
        ecfp4_dir / "part.parquet",
    )

    # MHFP6 scoring is opt-in: no shards at all is acceptable.
    annotate_ligand_similarity(data_dir=tmp_path)

    scoring_module.write_mhfp6_fingerprints(
        unique_ligands[["ligand_smiles_id", "ligand_rdkit_canonical_smiles"]],
        fingerprint_dir,
    )
    mhfp6_dir = tmp_path / scoring_module.MHFP6_SCORES_DIR
    mhfp6_dir.mkdir()
    scoring_module.mhfp6_ligand_scores(
        ligand_ids=[0],
        data_dir=tmp_path,
        output_path=mhfp6_dir / "part.parquet",
    )
    with pytest.raises(
        ValueError, match="MHFP6 score shards do not cover the fingerprint set"
    ):
        annotate_ligand_similarity(data_dir=tmp_path)

    scoring_module.mhfp6_ligand_scores(
        ligand_ids=[0, 1, 2],
        data_dir=tmp_path,
        output_path=mhfp6_dir / "part.parquet",
    )
    annotation_path = annotate_ligand_similarity(data_dir=tmp_path)

    assert pd.read_parquet(annotation_path)["ligand_smiles_id"].tolist() == [0, 1, 2]
    assert (
        json.loads((ecfp4_dir / scoring_module.LIGAND_SCORE_MANIFEST).read_text())[
            "fingerprint_column"
        ]
        == "fingerprint"
    )
    assert (
        json.loads((mhfp6_dir / scoring_module.LIGAND_SCORE_MANIFEST).read_text())[
            "fingerprint_column"
        ]
        == "mhfp6"
    )


def test_ligand_similarity_pipeline_does_not_write_per_system_mapping(
    tmp_path, monkeypatch
) -> None:
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    pd.DataFrame(
        {
            "entry_pdb_id": ["1aaa", "2bbb", "3ccc"],
            "system_id": ["1aaa_system", "2bbb_system", "3ccc_system"],
            "system_type": ["holo", "holo", "holo"],
            "ligand_is_proper": [True, True, True],
            "ligand_smiles": ["CCO", "CCO", "c1ccccc1"],
            "ligand_unique_ccd_code": ["LIG", "LIG", "BEN"],
            "ligand_id": ["1aaa__1__1.L", "2bbb__1__1.L", "3ccc__1__1.L"],
            "ligand_asym_id": ["L", "L", "L"],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    from plinder.data.annotations import ligand_utils

    # A cofactor's reference SMILES now comes from the CCD via _get_ccd_smiles
    # (the components.parquet lookup was retired); parse_cofactors only supplies
    # the code set. Stub both so COF resolves to ethanol, matching the "CCO"
    # ligand and exercising the exact-structure cofactor-like flag.
    monkeypatch.setattr(ligand_utils, "parse_cofactors", lambda _data_dir: {"COF"})
    monkeypatch.setattr(
        ligand_utils,
        "_get_ccd_smiles",
        lambda code: "CCO" if code == "COF" else None,
    )
    compute_ligand_fingerprints(data_dir=tmp_path)

    unique_ligands = pd.read_parquet(
        tmp_path / "fingerprints" / "ligands_per_smiles.parquet"
    )
    assert len(unique_ligands) == 2
    assert not (tmp_path / "fingerprints" / "ligands_per_system.parquet").exists()

    score_dir = tmp_path / "ligand_scores"
    score_dir.mkdir()
    retained_score = score_dir / "retained.parquet"
    retained_score.write_bytes(b"unchanged fingerprint score basis")
    stale_annotations = (
        tmp_path / "fingerprints" / "ligand_similarity_annotations.parquet"
    )
    stale_annotations.write_bytes(b"stale")
    compute_ligand_fingerprints(data_dir=tmp_path)
    assert retained_score.is_file()
    assert not stale_annotations.exists()
    retained_score.unlink()

    ligand_scores(
        ligand_ids=[int(unique_ligands["ligand_smiles_id"].iloc[0])],
        data_dir=tmp_path,
        output_path=score_dir / "part.parquet",
    )
    with pytest.raises(
        ValueError,
        match="BulkTanimoto score shards do not cover the fingerprint set",
    ):
        annotate_ligand_similarity(data_dir=tmp_path)

    (score_dir / "part.parquet").unlink()
    ligand_scores(
        ligand_ids=unique_ligands["ligand_smiles_id"].tolist(),
        data_dir=tmp_path,
        output_path=score_dir / "part.parquet",
    )
    annotation_path = annotate_ligand_similarity(data_dir=tmp_path)
    annotations = pd.read_parquet(annotation_path).set_index(
        "ligand_rdkit_canonical_smiles"
    )

    assert bool(annotations.loc["CCO", "ligand_is_cofactor_like"])
    assert annotations.loc["CCO", "ligand_smiles_id"] == 0
    assert annotations.loc["c1ccccc1", "ligand_smiles_id"] == 1
    assert not any("cluster" in column for column in annotations.columns)

    index = pd.read_parquet(index_dir / "annotation_table.parquet")
    index.loc[index.index[-1], "ligand_smiles"] = "CCN"
    index.to_parquet(index_dir / "annotation_table.parquet", index=False)
    compute_ligand_fingerprints(data_dir=tmp_path, retain_score_shards=True)
    updated_ids = pd.read_parquet(
        tmp_path / "fingerprints/ligands_per_smiles.parquet"
    ).set_index("ligand_rdkit_canonical_smiles")["ligand_smiles_id"]
    assert updated_ids.to_dict() == {"CCN": 0, "CCO": 1}
    assert (score_dir / "part.parquet").is_file()

    index.loc[index.index[-1], "ligand_smiles"] = "CCC"
    index.to_parquet(index_dir / "annotation_table.parquet", index=False)
    compute_ligand_fingerprints(data_dir=tmp_path)
    assert not list(score_dir.glob("*.parquet"))


def _system(pdb_id: str, ligands: list[LigandView]) -> SystemView:
    protein_chains = sorted(
        {chain for ligand in ligands for chain in ligand.protein_chains_asym_id}
    )
    return SystemView(
        id=f"{pdb_id}_system",
        pdb_id=pdb_id,
        system_type="holo",
        protein_chains_asym_id=protein_chains,
        proper_num_pocket_residues=0,
        proper_num_interactions=0,
        proper_num_unique_interactions=0,
        ligands={ligand.instance_chain: ligand for ligand in ligands},
    )


def _entry(pdb_id: str, system: SystemView) -> EntryView:
    asym_ids = [chain.split(".", 1)[1] for chain in system.protein_chains_asym_id]
    return EntryView(
        pdb_id=pdb_id,
        chains={
            asym_id: ChainView(asym_id=asym_id, auth_id=asym_id, length=100)
            for asym_id in asym_ids
        },
        systems={system.id: system},
        author_to_asym={asym_id: asym_id for asym_id in asym_ids},
    )


def test_scorer_limits_query_systems_by_protein_and_proper_ligand_chains(
    tmp_path,
) -> None:
    ligands = [
        _ligand(
            f"1abc__1__1.{chr(ord('C') + index)}",
            f"1.{chr(ord('C') + index)}",
            {"1.A": {10 + index: 9 + index}},
            {},
        )
        for index in range(6)
    ]
    five_ligands = _system("1abc", ligands[:5])
    six_ligands = _system("1abc", ligands)
    six_ligands_one_improper = _system(
        "1abc", [*ligands[:5], replace(ligands[5], is_proper=False)]
    )
    no_proper_ligands = _system("1abc", [replace(ligands[0], is_proper=False)])
    six_proteins = replace(
        five_ligands,
        protein_chains_asym_id=[f"1.{chain}" for chain in "ABCDEF"],
    )
    no_proteins = replace(five_ligands, protein_chains_asym_id=[])
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        max_query_protein_chains=5,
        max_query_proper_ligand_chains=5,
    )

    assert scorer.system_is_query_scoreable(five_ligands)
    assert not scorer.system_is_query_scoreable(six_ligands)
    assert scorer.system_is_query_scoreable(six_ligands_one_improper)
    assert not scorer.system_is_query_scoreable(no_proper_ligands)
    assert not scorer.system_is_query_scoreable(six_proteins)
    assert not scorer.system_is_query_scoreable(no_proteins)
    assert scorer.system_is_target_scoreable(six_ligands)
    assert scorer.system_is_target_scoreable(six_proteins)
    assert not scorer.system_is_target_scoreable(no_proteins)


def test_holo_scores_are_emitted_per_ligand_pair(tmp_path, monkeypatch) -> None:
    query_ligands = [
        _ligand("1abc__1__1.C", "1.C", {"1.A": {10: 9}}, {}),
        _ligand("1abc__1__1.D", "1.D", {"1.B": {20: 19}}, {}),
    ]
    target_ligands = [
        _ligand("2def__1__1.Z", "1.Z", {"1.X": {110: 109}}, {}),
        _ligand("2def__1__1.W", "1.W", {"1.Y": {120: 119}}, {}),
    ]
    query_system = _system("1abc", query_ligands)
    target_system = _system("2def", target_ligands)
    scorer = Scorer(
        entries={
            "1abc": _entry("1abc", query_system),
            "2def": _entry("2def", target_system),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    alignments = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [("2def", "A", "X"), ("2def", "B", "Y")],
            names=["target_entry", "query_chain_mapped", "target_chain_mapped"],
        )
    )
    protein_calls: list[tuple[tuple[str, ...], tuple[str, ...], int]] = []
    pocket_qcov_value = 1.0

    def protein_scores(
        _alignments: pd.DataFrame,
        _query_system: SystemView,
        target_protein_chains: list[str],
        query_length: int,
        query_protein_chains: list[str] | None = None,
    ) -> tuple[dict, dict, dict, str]:
        assert query_protein_chains is not None
        pair = (query_protein_chains[0], target_protein_chains[0])
        protein_calls.append(
            (
                tuple(query_protein_chains),
                tuple(target_protein_chains),
                query_length,
            )
        )
        return (
            {"protein_qcov_foldseek_weighted_sum": [pair]},
            {"protein_qcov_foldseek_weighted_sum": 0.75},
            {pair: pd.DataFrame()},
            "protein_qcov_foldseek",
        )

    def pocket_scores(
        _alns: dict,
        query_ligand: LigandView,
        target_ligand: LigandView,
    ) -> tuple[dict[str, float], dict[str, float], dict]:
        qcov = pocket_qcov_value * float(
            query_ligand.id == query_ligands[0].id
            and target_ligand.id == target_ligands[0].id
        )
        mapping = [
            (
                query_ligand.protein_chains_asym_id[0],
                target_ligand.protein_chains_asym_id[0],
            )
        ]
        return (
            {
                "pocket_qcov_foldseek": qcov,
                "pocket_fident_qcov_foldseek": 0.19 * qcov,
            },
            {"pli_qcov_foldseek": 0.11 * qcov},
            {"pocket_qcov_foldseek": mapping},
        )

    shape_calls: list[tuple[str, str, float]] = []

    def shape_scores(
        _data_dir: Path,
        query_ligand: LigandView,
        target_ligand: LigandView,
        pocket_qcov: float,
    ) -> dict[str, float]:
        shape_calls.append((query_ligand.id, target_ligand.id, pocket_qcov))
        return {
            "shape": 0.8,
            "color": 0.7,
            "sucos_shape": 0.6,
            "sucos_shape_pocket_qcov": 0.6 * pocket_qcov,
        }

    monkeypatch.setattr(scorer, "get_protein_scores", protein_scores)
    monkeypatch.setattr(scorer, "get_ligand_pair_pocket_pli_scores", pocket_scores)
    monkeypatch.setattr(scorer, "get_ligand_pair_shape_scores", shape_scores)

    scores = list(
        scorer.get_scores_holo(
            query_system,
            alignments,
            data_dir=tmp_path,
            include_protein_scores=True,
        )
    )

    assert {
        (score["query_ligand_id"], score["target_ligand_id"]) for score in scores
    } == {(query.id, target.id) for query in query_ligands for target in target_ligands}
    assert len(protein_calls) == 4
    assert {
        (query_chains, target_chains)
        for query_chains, target_chains, _ in protein_calls
    } == {
        (("1.A",), ("1.X",)),
        (("1.A",), ("1.Y",)),
        (("1.B",), ("1.X",)),
        (("1.B",), ("1.Y",)),
    }
    assert all(query_length == 100 for _, _, query_length in protein_calls)
    assert shape_calls == [(query_ligands[0].id, target_ligands[0].id, 1.0)]
    shape_row = next(score for score in scores if "shape" in score)
    assert shape_row["sucos_shape_pocket_qcov"] == pytest.approx(0.6)

    ligand_scores = list(
        scorer.get_scores_holo(
            query_system,
            alignments,
            include_protein_scores=False,
        )
    )
    assert all("protein_qcov_weighted_sum" not in score for score in ligand_scores)
    assert any("pocket_qcov" in score for score in ligand_scores)
    assert len(protein_calls) == 4

    candidates = []
    ligand_pair_scores = []
    protein_only_scores = list(
        scorer.get_scores_holo(
            query_system,
            alignments,
            ligand_3d_candidates=candidates,
            ligand_pair_scores=ligand_pair_scores,
        )
    )

    assert len(shape_calls) == 1
    assert all("shape" not in score for score in protein_only_scores)
    assert candidates == [
        {
            "query_system": query_system.id,
            "query_ligand_id": query_ligands[0].id,
            "query_entry": "1abc",
            "query_ligand_asym_id": "C",
            "target_system": target_system.id,
            "target_ligand_id": target_ligands[0].id,
            "target_entry": "2def",
            "target_ligand_asym_id": "Z",
            "protein_mapping": "1.A:1.X",
            "protein_mapper": "foldseek",
            "pocket_qcov": 1.0,
        }
    ]
    assert ligand_pair_scores == [
        {
            "query_system": query_system.id,
            "query_ligand_id": query_ligands[0].id,
            "query_entry": "1abc",
            "query_ligand_asym_id": "C",
            "target_system": target_system.id,
            "target_ligand_id": target_ligands[0].id,
            "target_entry": "2def",
            "target_ligand_asym_id": "Z",
            "pocket_qcov": 100,
            "pocket_fident_qcov": 19,
            "pli_qcov": 11,
        }
    ]

    pocket_qcov_value = 0.004
    tiny_ligand_pair_scores = []
    list(
        scorer.get_scores_holo(
            query_system,
            alignments,
            ligand_pair_scores=tiny_ligand_pair_scores,
        )
    )
    assert len(tiny_ligand_pair_scores) == 1
    assert tiny_ligand_pair_scores[0]["pocket_qcov"] == 0
    assert tiny_ligand_pair_scores[0]["pocket_fident_qcov"] == 0
    assert tiny_ligand_pair_scores[0]["pli_qcov"] == 0


def test_holo_threaded_scoring_reuses_canonical_and_receptor_pairs(
    tmp_path, monkeypatch
) -> None:
    query_ligand = _ligand("1abc__1__1.C", "1.C", {"1.A": {10: 9}}, {})
    target_ligand = _ligand("2def__1__1.Z", "1.Z", {"1.X": {110: 109}}, {})
    repeated_target = replace(
        target_ligand,
        id="2def__2__2.Z",
        system_id="2def_alt",
        instance_chain="2.Z",
    )
    query_system = _system("1abc", [query_ligand])
    target_system = _system("2def", [target_ligand])
    repeated_system = replace(_system("2def", [repeated_target]), id="2def_alt")
    target_entry = _entry("2def", target_system)
    target_entry.systems[repeated_system.id] = repeated_system
    scorer = Scorer(
        entries={
            "1abc": _entry("1abc", query_system),
            "2def": target_entry,
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        ligand_sdf_resolver=lambda _ligand: SDF_FILE,
        shape_score_threads=2,
    )
    alignments = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [("2def", "A", "X")],
            names=["target_entry", "query_chain_mapped", "target_chain_mapped"],
        )
    )
    protein_calls = 0
    shape_calls = 0

    def protein_scores(*_args, **_kwargs) -> tuple[dict, dict, dict, str]:
        nonlocal protein_calls
        protein_calls += 1
        pair = ("1.A", "1.X")
        return (
            {"protein_qcov_foldseek_weighted_sum": [pair]},
            {"protein_qcov_foldseek_weighted_sum": 0.75},
            {pair: pd.DataFrame()},
            "protein_qcov_foldseek",
        )

    def pocket_scores(
        _alns: dict, _query: LigandView, target: LigandView
    ) -> tuple[dict[str, float], dict[str, float], dict]:
        return (
            {"pocket_qcov_foldseek": (0.5 if target.id == target_ligand.id else 0.25)},
            {},
            {},
        )

    def align_once(*_args) -> tuple[float, float]:
        nonlocal shape_calls
        shape_calls += 1
        return 0.8, 0.6

    monkeypatch.setattr(scorer, "get_protein_scores", protein_scores)
    monkeypatch.setattr(scorer, "get_ligand_pair_pocket_pli_scores", pocket_scores)
    monkeypatch.setattr(scoring_module, "align_molecules", align_once)
    monkeypatch.setattr(scoring_module, "get_sucos_score", lambda *_args: 0.7)

    scores = list(
        scorer.get_scores_holo(
            query_system,
            alignments,
            data_dir=tmp_path,
            include_protein_scores=True,
        )
    )

    assert len(scores) == 2
    assert protein_calls == 1
    assert shape_calls == 1
    assert sorted(
        float(score["sucos_shape_pocket_qcov"]) for score in scores
    ) == pytest.approx([0.175, 0.35])


def test_holo_weighted_sum_retains_unmatched_query_receptor_length(
    tmp_path, monkeypatch
) -> None:
    query_ligand = _ligand(
        "1abc__1__1.C",
        "1.C",
        {"1.A": {10: 9}, "1.B": {20: 19}},
        {},
    )
    target_ligand = _ligand("2def__1__1.Z", "1.Z", {"1.X": {110: 109}}, {})
    query_system = _system("1abc", [query_ligand])
    target_system = _system("2def", [target_ligand])
    scorer = Scorer(
        entries={
            "1abc": _entry("1abc", query_system),
            "2def": _entry("2def", target_system),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    # Only query chain A has a target hit. Chain B must nevertheless remain
    # in the 200-residue directed query denominator.
    alignments = pd.DataFrame(
        [
            {
                "qcov": 0.5,
                "fident": 0.5,
                "seqsim": 0.5,
                "fident_qcov": 0.25,
                "seqsim_qcov": 0.25,
                "lddt": 0.5,
                "lddt_qcov": 0.5,
            }
        ],
        index=pd.MultiIndex.from_tuples(
            [("2def", "A", "X", "foldseek")],
            names=[
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
                "source",
            ],
        ),
    )
    monkeypatch.setattr(
        scorer,
        "get_ligand_pair_pocket_pli_scores",
        lambda *_args: ({}, {}, {}),
    )

    scores = list(
        scorer.get_scores_holo(
            query_system,
            alignments,
            include_protein_scores=True,
        )
    )

    assert len(scores) == 1
    assert scores[0]["protein_lddt_qcov_weighted_max"] == pytest.approx(0.5)
    assert scores[0]["protein_lddt_qcov_weighted_sum"] == pytest.approx(0.25)


def test_apo_pred_scores_are_emitted_per_query_ligand(tmp_path, monkeypatch) -> None:
    query_ligands = [
        _ligand(
            "1abc__1__1.C",
            "1.C",
            {"1.A": {10: 9}, "1.B": {20: 19}},
            {},
        ),
        _ligand("1abc__1__1.D", "1.D", {"1.B": {20: 19}}, {}),
    ]
    query_system = _system("1abc", query_ligands)
    scorer = Scorer(
        entries={"1abc": _entry("1abc", query_system)},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    target_system = replace(
        query_system,
        id="model_a_system",
        pdb_id="model_a",
        protein_chains_asym_id=["0.X"],
        ligands={},
    )
    scorer.entries["model_a"] = _entry("model_a", target_system)
    alignments = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [("model_a", "A", "X"), ("model_b", "B", "Y")],
            names=["target_entry", "query_chain_mapped", "target_chain_mapped"],
        )
    )
    protein_calls: list[tuple[tuple[str, ...], int]] = []

    def protein_scores(
        _alignments: pd.DataFrame,
        _query_system: SystemView,
        target_protein_chains: list[str],
        query_length: int,
        query_protein_chains: list[str] | None = None,
    ) -> tuple[dict, dict, dict, str]:
        assert query_protein_chains is not None
        protein_calls.append((tuple(query_protein_chains), query_length))
        pair = (query_protein_chains[0], target_protein_chains[0])
        return (
            {"protein_qcov_foldseek_weighted_sum": [pair]},
            {"protein_qcov_foldseek_weighted_sum": 0.75},
            {pair: pd.DataFrame()},
            "protein_qcov_foldseek",
        )

    monkeypatch.setattr(scorer, "get_protein_scores", protein_scores)
    monkeypatch.setattr(
        scorer,
        "get_ligand_pocket_scores",
        lambda _alns, _ligand: {"pocket_fident_foldseek": 0.5},
    )

    scores = list(scorer.get_scores_apo_pred(query_system, alignments))

    assert {score["query_ligand_id"] for score in scores} == {
        ligand.id for ligand in query_ligands
    }
    assert {score["target_ligand_id"] for score in scores} == {None}
    assert {score["target_system"] for score in scores} == {
        "model_a_X",
        "model_b_Y",
    }
    assert (("1.A", "1.B"), 200) in protein_calls


def test_aggregate_scores_keeps_ligand_ids_and_shape_metrics(
    tmp_path, monkeypatch
) -> None:
    ligand = _ligand("1abc__1__1.C", "1.C", {"1.A": {10: 9}}, {})
    system = _system("1abc", [ligand])
    scorer = Scorer(
        entries={"1abc": _entry("1abc", system)},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        minimum_threshold=0.3,
    )
    monkeypatch.setattr(
        scorer,
        "load_alignments",
        lambda **_kwargs: pd.DataFrame(
            {"present": [True]}, index=pd.Index(["1abc"], name="query_entry")
        ),
    )
    monkeypatch.setattr(
        scorer,
        "get_scores",
        lambda *_args, **_kwargs: iter(
            [
                {
                    "query_system": system.id,
                    "query_ligand_id": ligand.id,
                    "target_system": "2def_system",
                    "target_ligand_id": "2def__1__1.Z",
                    "protein_mapping": "1.A:1.X",
                    "protein_mapper": "foldseek",
                    "protein_qcov_weighted_sum": 0.1,
                    "protein_qcov_weighted_sum_source": "foldseek",
                    "protein_qcov_weighted_sum_mapping": "1.A:1.X",
                    "shape": 0.3,
                    "color": 0.31,
                    "sucos_shape": 0.32,
                    "sucos_shape_pocket_qcov": 0.33,
                }
            ]
        ),
    )

    scores = scorer.aggregate_scores("1abc", data_dir=tmp_path)

    assert scores is not None
    assert set(scores["metric"]) >= {
        "shape",
        "color",
        "sucos_shape",
        "sucos_shape_pocket_qcov",
    }
    assert set(scores["query_ligand_id"]) == {ligand.id}
    assert set(scores["target_ligand_id"]) == {"2def__1__1.Z"}
    assert scores.loc[scores["metric"] == "shape", "similarity"].item() == 30
    assert scores.loc[scores["metric"] == "color", "similarity"].item() == 31
    assert scores.loc[scores["metric"] == "sucos_shape", "similarity"].item() == 32
    assert (
        scores.loc[scores["metric"] == "sucos_shape_pocket_qcov", "similarity"].item()
        == 33
    )
    assert "protein_qcov_weighted_sum" not in set(scores["metric"])
    assert {"query_ligand_id", "target_ligand_id"}.issubset(
        PROTEIN_SIMILARITY_SCHEMA.names
    )
    output_file = tmp_path / "scores.parquet"
    scores.to_parquet(output_file, index=False, schema=PROTEIN_SIMILARITY_SCHEMA)
    written = pd.read_parquet(output_file)
    assert set(written["query_ligand_id"]) == {ligand.id}
