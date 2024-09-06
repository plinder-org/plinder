# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pytest
from plinder.data.pipeline import utils

_ENTRY = """
{{
    "pdb_id": "{pdb_id}",
    "release_date": "2002-11-29",
    "oligomeric_state": "monomeric",
    "determination_method": "X-RAY DIFFRACTION",
    "keywords": "MEMBRANE PROTEIN",
    "pH": "8",
    "resolution": 2.7
}}
""".format


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
    path = tmp_path / "aa" / "aaaa.json"
    path.parent.mkdir(parents=True)
    if expect:
        path.write_text("")
    assert (
        utils.entry_exists(
            entry_dir=tmp_path,
            pdb_id="aaaa",
        )
        == expect
    )


def test_load_entries(tmp_path):
    a = tmp_path / "raw_entries" / "aa" / "aaaa.json"
    b = tmp_path / "raw_entries" / "bb" / "bbbb.json"
    a.parent.mkdir(parents=True)
    b.parent.mkdir(parents=True)
    a.write_text(_ENTRY(pdb_id="aaaa"))
    b.write_text(_ENTRY(pdb_id="bbbb"))
    ret = utils.load_entries(
        data_dir=tmp_path,
        pdb_ids=["aaaa", "bbbb"],
    )
    assert len(ret) == 2


@pytest.mark.parametrize(
    "kwargs, expect",
    [
        ({}, 2),
        ({"two_char_codes": ["2g"]}, 2),
        ({"pdb_ids": ["22gs"]}, 1),
    ],
)
def test_load_entries_from_zips(kwargs, expect, tmp_path, entry_zip):
    zip_dir = tmp_path / "entries"
    zip_dir.mkdir(parents=True)
    (zip_dir / entry_zip.name).write_bytes(entry_zip.read_bytes())
    entries = utils.load_entries_from_zips(data_dir=tmp_path, **kwargs)
    assert len(entries) == expect


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
    a = tmp_path / "aa" / "aaaa" / "aaaa.json"
    b = tmp_path / "bb" / "bbbb" / "bbbb.json"
    a.parent.mkdir(parents=True)
    b.parent.mkdir(parents=True)
    a.write_text(_ENTRY(pdb_id="aaaa"))
    b.write_text(_ENTRY(pdb_id="bbbb"))
    contents = utils.get_local_contents(
        data_dir=tmp_path,
        two_char_codes=two_char_codes,
    )
    assert len(contents) == expect


def test_get_local_contents_pdb_ids(tmp_path):
    a = tmp_path / "aa" / "pdb_0000aaaa" / "aaaa.json"
    b = tmp_path / "bb" / "pdb_0000bbbb" / "bbbb.json"
    a.parent.mkdir(parents=True)
    b.parent.mkdir(parents=True)
    a.write_text(_ENTRY(pdb_id="aaaa"))
    b.write_text(_ENTRY(pdb_id="bbbb"))
    contents = utils.get_local_contents(data_dir=tmp_path, as_four_char_ids=True)
    assert contents == ["aaaa", "bbbb"]
