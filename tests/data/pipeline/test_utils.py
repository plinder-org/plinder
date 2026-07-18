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
        pd.DataFrame({"system_receptor_type": ["protein"]}).to_parquet(
            path, index=False
        )
        chain_path = tmp_path / "aa" / "aaaa" / "entry_chains.parquet"
        chain_path.parent.mkdir()
        pd.DataFrame({"chain_receptor_type": ["protein"]}).to_parquet(
            chain_path, index=False
        )
        pd.DataFrame(
            {
                "entry_pdb_id": ["aaaa"],
                "biounit_id": ["1"],
                "chain_instance": ["1.A"],
                "chain_asym_id": ["A"],
                "chain_role": ["receptor"],
            }
        ).to_parquet(chain_path.parent / "entry_biounit_chains.parquet", index=False)
        pd.DataFrame({"entry_pdb_id": ["aaaa"]}).to_parquet(
            chain_path.parent / "entry_source.parquet", index=False
        )
    assert (
        utils.entry_exists(
            entry_dir=tmp_path,
            pdb_id="aaaa",
        )
        == expect
    )


def test_entry_exists_requires_entry_chain_table(tmp_path):
    annotation = tmp_path / "aa" / "aaaa.parquet"
    annotation.parent.mkdir(parents=True)
    annotation.touch()

    assert not utils.entry_exists(entry_dir=tmp_path, pdb_id="aaaa")


def test_entry_exists_requires_entry_source_table(tmp_path):
    annotation = tmp_path / "aa" / "aaaa.parquet"
    annotation.parent.mkdir(parents=True)
    annotation.touch()
    entry_dir = tmp_path / "aa" / "aaaa"
    entry_dir.mkdir()
    (entry_dir / "entry_chains.parquet").touch()

    assert not utils.entry_exists(entry_dir=tmp_path, pdb_id="aaaa")


def test_entry_exists_requires_current_biounit_chain_table(tmp_path):
    annotation = tmp_path / "aa" / "aaaa.parquet"
    annotation.parent.mkdir(parents=True)
    pd.DataFrame({"system_receptor_type": ["protein"]}).to_parquet(
        annotation,
        index=False,
    )
    entry_dir = annotation.parent / "aaaa"
    entry_dir.mkdir()
    pd.DataFrame({"chain_receptor_type": ["protein"]}).to_parquet(
        entry_dir / "entry_chains.parquet",
        index=False,
    )
    pd.DataFrame({"chain_instance": ["1.A"]}).to_parquet(
        entry_dir / "entry_biounit_chains.parquet",
        index=False,
    )
    pd.DataFrame({"entry_pdb_id": ["aaaa"]}).to_parquet(
        entry_dir / "entry_source.parquet",
        index=False,
    )

    assert not utils.entry_exists(entry_dir=tmp_path, pdb_id="aaaa")


def test_entry_exists_invalidates_pre_receptor_type_cache(tmp_path):
    annotation = tmp_path / "aa" / "aaaa.parquet"
    annotation.parent.mkdir(parents=True)
    pd.DataFrame({"system_id": ["aaaa__1__1.A__1.L"]}).to_parquet(
        annotation, index=False
    )
    entry_dir = tmp_path / "aa" / "aaaa"
    entry_dir.mkdir()
    pd.DataFrame({"chain_type": ["polypeptide(L)"]}).to_parquet(
        entry_dir / "entry_chains.parquet", index=False
    )
    pd.DataFrame({"entry_pdb_id": ["aaaa"]}).to_parquet(
        entry_dir / "entry_source.parquet", index=False
    )

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
    pd.DataFrame(
        {
            "system_id": ["1aaa__1__A"],
            "system_pocket_ECOD": ["e1aaaA1"],
            "system_pocket_PANTHER": ["PTHR00001"],
            "system_pocket_kinase_name": ["example kinase"],
        }
    ).to_parquet(first, index=False)
    pd.DataFrame(
        {
            "system_id": ["2bbb__1__B"],
            "ligand_is_kinase_inhibitor": [True],
            "system_has_kinase_inhibitor": [True],
        }
    ).to_parquet(second, index=False)
    chain_columns = {
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
    for pdb_id, code in [("1aaa", "aa"), ("2bbb", "bb")]:
        chain_path = tmp_path / "raw_entries" / code / pdb_id / "entry_chains.parquet"
        chain_path.parent.mkdir(parents=True)
        pd.DataFrame({"entry_pdb_id": [pdb_id], **chain_columns}).to_parquet(
            chain_path, index=False
        )
        pd.DataFrame(
            {
                "entry_pdb_id": [pdb_id],
                "source_mmcif_major_revision": [1],
                "source_mmcif_minor_revision": [0],
            }
        ).to_parquet(chain_path.parent / "entry_source.parquet", index=False)
        pd.DataFrame(
            {
                "entry_pdb_id": [pdb_id],
                "biounit_id": ["1"],
                "chain_instance": ["1.A"],
                "chain_asym_id": ["A"],
                "chain_role": ["receptor"],
            }
        ).to_parquet(chain_path.parent / "entry_biounit_chains.parquet", index=False)
    monkeypatch.setattr(utils, "add_aggregated_columns", lambda index: index)

    index = utils.create_index(data_dir=tmp_path, force_update=True)

    assert index["system_id"].tolist() == ["1aaa__1__A", "2bbb__1__B"]
    assert not {
        column
        for column in index.columns
        if any(marker in column.casefold() for marker in ("ecod", "panther", "kinase"))
    }
    assert (tmp_path / "index" / "annotation_table.parquet").is_file()
    entry_chains = pd.read_parquet(tmp_path / "index" / "entry_chains.parquet")
    assert entry_chains["entry_pdb_id"].tolist() == ["1aaa", "2bbb"]
    biounit_chains = pd.read_parquet(
        tmp_path / "index" / "entry_biounit_chains.parquet"
    )
    assert biounit_chains["entry_pdb_id"].tolist() == ["1aaa", "2bbb"]
    entry_sources = pd.read_parquet(tmp_path / "index" / "entry_sources.parquet")
    assert entry_sources["entry_pdb_id"].tolist() == ["1aaa", "2bbb"]
    assert entry_sources["source_mmcif_major_revision"].tolist() == [1, 1]


def test_create_entry_chain_index_handles_empty_entry_table(tmp_path):
    entry_table = tmp_path / "raw_entries" / "dn" / "1dna" / "entry_chains.parquet"
    entry_table.parent.mkdir(parents=True)
    pd.DataFrame().to_parquet(entry_table, index=False)

    chains = utils.create_entry_chain_index(data_dir=tmp_path, force_update=True)

    assert chains.empty
    assert chains.columns.tolist() == [
        "entry_pdb_id",
        "chain_asym_id",
        "chain_auth_id",
        "chain_entity_id",
        "chain_type",
        "chain_receptor_type",
        "chain_length",
        "chain_num_unresolved_residues",
        "chain_is_holo",
        "chain_uniprot_ids",
    ]


def test_create_entry_biounit_chain_index_handles_empty_input(tmp_path):
    chains = utils.create_entry_biounit_chain_index(
        data_dir=tmp_path,
        force_update=True,
    )

    assert chains.empty
    assert chains.columns.tolist() == [
        "entry_pdb_id",
        "biounit_id",
        "chain_instance",
        "chain_asym_id",
        "chain_role",
    ]


def test_finalize_index_creates_nonredundant_data_from_local_clusters(tmp_path):
    index_dir = tmp_path / "index"
    cluster_file = (
        tmp_path
        / "clusters/cluster=components/directed=True/metric=pli_qcov"
        / "threshold=100.parquet"
    )
    index_dir.mkdir(parents=True)
    cluster_file.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "system_id": ["1aaa__1__1.A__1.X", "1aaa__2__1.A__1.X"],
            "system_id_no_biounit": ["1aaa__1.A__1.X", "1aaa__1.A__1.X"],
            "system_biounit_id": ["1", "2"],
            "system_type": ["holo", "holo"],
            "ligand_id": ["1aaa__1__1.X", "1aaa__2__1.X"],
            "ligand_is_proper": [True, False],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCO"],
        }
    ).to_parquet(index_dir / "annotation_table.parquet", index=False)
    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_rdkit_canonical_smiles": ["CCO"],
            "ligand_smiles_id": [0],
            "ligand_is_cofactor_like": [False],
            "ligand_tanimoto_ecfp4_1024_90_cluster": ["c0"],
            "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids": [1],
        }
    ).to_parquet(fingerprint_dir / "ligand_similarity_annotations.parquet", index=False)
    ligand_dir = tmp_path / "ligands"
    ligand_dir.mkdir()
    pd.DataFrame(
        {
            "ligand_id": ["1aaa__1__1.X"],
            "ligand_is_3d_score_able": [True],
        }
    ).to_parquet(ligand_dir / "part.parquet", index=False)
    pd.DataFrame(
        {
            "system_id": ["1aaa__1__1.A__1.X", "1aaa__2__1.A__1.X"],
            "label": ["c0", "c0"],
            "metric": ["pli_qcov", "pli_qcov"],
            "cluster": ["components", "components"],
            "directed": [True, True],
            "threshold": [100, 100],
        }
    ).to_parquet(cluster_file, index=False)

    utils.finalize_index(data_dir=tmp_path)
    utils.create_nonredundant_dataset(data_dir=tmp_path)

    finalized = pd.read_parquet(index_dir / "annotation_table.parquet")
    assert finalized["pli_qcov__100__strong__component"].tolist() == ["c0", "c0"]
    assert finalized.loc[0, "ligand_smiles_id"] == 0
    assert pd.isna(finalized.loc[1, "ligand_smiles_id"])
    assert finalized["ligand_is_3d_score_able"].tolist() == [True, False]
    assert finalized["uniqueness"].nunique() == 1
    nonredundant = pd.read_parquet(index_dir / "annotation_table_nonredundant.parquet")
    assert len(nonredundant) == 1
