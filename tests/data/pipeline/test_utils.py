# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
import pytest
from plinder.data.pipeline import utils


@pytest.mark.parametrize(
    "funcname, run, skip, expect",
    [
        ("a", [], [], True),
        ("a", ["a"], [], True),
        ("a", [], ["b"], True),
        ("a", [], ["a"], False),
        ("scatter_a", ["a", "b"], ["a"], False),
        ("join_a", [], [], True),
    ],
)
def test_should_run_stage(funcname, run, skip, expect):
    assert utils.should_run_stage(funcname, run, skip) == expect


@pytest.mark.parametrize(
    "expect",
    [
        (False),
        (True),
    ],
)
def test_entry_exists(expect, tmp_path):
    path = tmp_path / "aa" / "aaaa.parquet"
    path.parent.mkdir(parents=True)
    if expect:
        path.write_text("")
        chain_path = tmp_path / "aa" / "aaaa" / "entry_chains.parquet"
        chain_path.parent.mkdir()
        chain_path.write_text("")
    assert (
        utils.entry_exists(
            entry_dir=tmp_path,
            pdb_id="aaaa",
        )
        == expect
    )


def test_entry_exists_requires_chain_sidecar(tmp_path):
    annotation = tmp_path / "aa" / "aaaa.parquet"
    annotation.parent.mkdir(parents=True)
    annotation.touch()

    assert not utils.entry_exists(entry_dir=tmp_path, pdb_id="aaaa")


@pytest.mark.parametrize(
    "contents",
    [
        ["aaaa", "bbbb"],
        ["pdb_0000cccc", "pdb_0000dddd"],
    ],
)
def test_hash_contents(contents):
    utils.hash_contents(contents)


@pytest.mark.parametrize(
    "two_char_codes, expect",
    [
        (["aa", "bb"], 2),
        (["aa"], 1),
        ([], 2),
        (None, 2),
    ],
)
def test_get_local_contents(two_char_codes, expect, tmp_path):
    a = tmp_path / "aa" / "aaaa" / "aaaa.cif"
    b = tmp_path / "bb" / "bbbb" / "bbbb.cif"
    a.parent.mkdir(parents=True)
    b.parent.mkdir(parents=True)
    a.touch()
    b.touch()
    contents = utils.get_local_contents(
        data_dir=tmp_path,
        two_char_codes=two_char_codes,
    )
    assert len(contents) == expect


def test_get_local_contents_pdb_ids(tmp_path):
    a = tmp_path / "aa" / "pdb_0000aaaa" / "aaaa.cif"
    b = tmp_path / "bb" / "pdb_0000bbbb" / "bbbb.cif"
    a.parent.mkdir(parents=True)
    b.parent.mkdir(parents=True)
    a.touch()
    b.touch()
    contents = utils.get_local_contents(data_dir=tmp_path, as_four_char_ids=True)
    assert contents == ["aaaa", "bbbb"]


def test_create_index_collates_per_entry_parquets(tmp_path, monkeypatch):
    first = tmp_path / "raw_entries" / "aa" / "1aaa.parquet"
    second = tmp_path / "raw_entries" / "bb" / "2bbb.parquet"
    first.parent.mkdir(parents=True)
    second.parent.mkdir(parents=True)
    pd.DataFrame({"system_id": ["1aaa__1__A"]}).to_parquet(first, index=False)
    pd.DataFrame({"system_id": ["2bbb__1__B"]}).to_parquet(second, index=False)
    chain_columns = {
        "chain_asym_id": ["A"],
        "chain_auth_id": ["A"],
        "chain_entity_id": ["1"],
        "chain_type": ["polypeptide(L)"],
        "chain_length": [100],
        "chain_num_unresolved_residues": [0],
        "chain_is_holo": [True],
        "chain_uniprot_ids": [["P12345"]],
    }
    for pdb_id, code in [("1aaa", "aa"), ("2bbb", "bb")]:
        chain_path = tmp_path / "raw_entries" / code / pdb_id / "entry_chains.parquet"
        chain_path.parent.mkdir(parents=True)
        pd.DataFrame({"entry_pdb_id": [pdb_id], **chain_columns}).to_parquet(
            chain_path, index=False
        )
    monkeypatch.setattr(utils, "add_aggregated_columns", lambda index: index)

    index = utils.create_index(data_dir=tmp_path, force_update=True)

    assert index["system_id"].tolist() == ["1aaa__1__A", "2bbb__1__B"]
    assert (tmp_path / "index" / "annotation_table.parquet").is_file()
    entry_chains = pd.read_parquet(tmp_path / "index" / "entry_chains.parquet")
    assert entry_chains["entry_pdb_id"].tolist() == ["1aaa", "2bbb"]
