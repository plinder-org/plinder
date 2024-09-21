# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import pandas as pd
import pytest
from plinder.core.utils import io as core_io
from plinder.data.pipeline import io


class _resp:
    def __init__(self, body=None, text=None, content=None):
        self.status_code = 200
        self.called = 0
        self.body = body or {}
        self.text = text or ""
        self.content = content or b""

    def raise_for_status(self):
        pass

    def json(self):
        self.called += 1
        if isinstance(self.body, list):
            return self.body[self.called - 1]
        return self.body


def test_retry_fail(monkeypatch):
    monkeypatch.setattr(
        "plinder.core.utils.io.sleep",
        lambda _: None,
    )

    @core_io.retry
    def f():
        raise Exception("intentional")

    with pytest.raises(Exception):
        f()


def test_download_cofactors_fetch(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.requests.get",
        lambda _: _resp(body={"test": "foo"}),
    )
    obj = io.download_cofactors(data_dir=tmp_path)
    assert "test" in obj


def test_download_cofactors_cached(tmp_path):
    cofactors = tmp_path / "dbs" / "cofactors" / "cofactors.json"
    cofactors.parent.mkdir(parents=True)
    cofactors.write_text('{"test": "foo"}')
    obj = io.download_cofactors(data_dir=tmp_path)
    assert "test" in obj


def test_download_ecod_data_fetch(tmp_path, raw_ecod_data, monkeypatch):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.requests.get", lambda _: _resp(text=raw_ecod_data)
    )
    ecod_path = io.download_ecod_data(data_dir=tmp_path)
    assert len(pd.read_parquet(ecod_path).index)


def test_download_ecod_data_cached(test_env, raw_ecod_path):
    ecod_path = io.download_ecod_data(data_dir=test_env)
    assert len(pd.read_parquet(ecod_path).index)


def test_ecod(test_env, raw_ecod_path):
    ecod_path = io.download_ecod_data(data_dir=test_env)
    df = pd.read_parquet(ecod_path)
    reference = "acid protease"
    assert (
        df.loc[df["domainid"] == "e1udzA1", "domain"].values[0].strip()
        == reference.strip()
    )
    assert (
        df.loc[(df["pdb"] == "1udz") & (df["chain"] == "A"), "domain"].values[0]
        == reference.strip()
    )


def test_download_panther_data_fetch(tmp_path, monkeypatch, raw_panther_data):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.requests.get",
        lambda _: _resp(content=raw_panther_data),
    )
    panther_path = io.download_panther_data(data_dir=tmp_path)
    assert panther_path.is_file()
    assert len(pd.read_parquet(panther_path).index)


def test_download_panther_data_cached(test_env, raw_panther_path):
    panther_path = io.download_panther_data(data_dir=test_env)
    assert len(pd.read_parquet(panther_path).index)


def test_panther(test_env, raw_panther_path):
    panther_path = io.download_panther_data(data_dir=test_env)
    df = pd.read_parquet(panther_path)
    reference = "E5R0L8"
    assert df.loc[df["uniprotac"] == reference, "panther"].values[0] == "NOHIT"
    assert (
        df.loc[df["uniprotac"] == reference, "panther_class"].values[0]
        == "Arthroderma_gypseum_ARTGP"
    )


def test_download_seqres_data_fetch(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.requests.get", lambda _: _resp(content=b"foo")
    )
    raw_seqres_path = io.download_seqres_data(data_dir=tmp_path)
    assert raw_seqres_path.is_file()
    assert raw_seqres_path.read_text() == "foo"


def test_download_seqres_data_cached(tmp_path):
    raw_seqres_path = tmp_path / "dbs" / "seqres" / "pdb_seqres.txt.gz"
    raw_seqres_path.parent.mkdir(parents=True)
    raw_seqres_path.write_bytes(b"foo")
    raw_seqres_path = io.download_seqres_data(data_dir=tmp_path)
    assert raw_seqres_path.read_text() == "foo"


def test_download_kinase_data_fetch(tmp_path, monkeypatch):
    resp = _resp(
        body=[
            [
                {
                    "name": "foo",
                    "HGNC": "foo",
                    "family": "foo",
                    "group": "foo",
                    "kinase_class": "foo",
                    "full_name": "foo",
                    "uniprot": "foo",
                    "kinase_ID": "2",
                }
            ],
            [{"kinase_ID": "2"}],
            [
                {
                    "pdb": "foob",
                    "chain": "A",
                    "structure_ID": "1",
                    "kinase_ID": "2",
                    "missing_atoms": "5",
                    "missing_residues": "6",
                    "rmsd1": "1.00",
                    "rmsd2": "2.00",
                    "resolution": "3.0",
                    "quality_score": "9.2",
                    "Grich_distance": "1.0",
                    "Grich_rotation": "5.6",
                    "Grich_angle": "23.4",
                }
            ],
        ]
    )
    monkeypatch.setattr("plinder.data.pipeline.io.requests.get", lambda *_, **__: resp)
    kinase_path = io.download_kinase_data(data_dir=tmp_path)
    assert kinase_path.is_file()


def test_kinase(test_env, all_kinase_paths):
    kinase_path = io.download_kinase_data(data_dir=test_env)
    df = pd.read_parquet(kinase_path)
    assert df.loc[df["pdbid_chainid"] == "3mvh_A", "kinase"].values[0] == "AKT1"
    assert df.loc[df["pdbid_chainid"] == "3mvh_A", "uniprot"].values[0] == "P31749"


def test_kinase_ligand(test_env, all_kinase_paths):
    kinase_path = io.download_kinase_data(data_dir=test_env)
    df = pd.read_parquet(kinase_path.parent / "kinase_ligand_ccd_codes.parquet")
    assert (
        df.loc[df["PDB-code"] == "IHZ", "SMILES"].values[0]
        == "FC(F)(F)c1cc(NC(=O)c2cc(Nc3cncc(c3)C(=O)N)c(cc2)C)ccc1"
    )
    assert (
        df.loc[df["PDB-code"] == "IHZ", "InChIKey"].values[0]
        == "SAAYRHKJHDIDPH-UHFFFAOYSA-N"
    )


def test_rsync_rcsb(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.check_output",
        lambda *_, **__: None,
    )
    io.rsync_rcsb(kind="cif", two_char_code="aa", data_dir=tmp_path)


def test_list_rcsb(monkeypatch):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.check_output",
        lambda *_, **__: "d aaaa \nd baan",
    )
    pdbs = io.list_rcsb(kind="cif", two_char_code="aa")
    assert pdbs == ["aaaa", "baan"]


def test_downloads_components_cif_fetch(tmp_path, monkeypatch, mini_component_cif_gz):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.requests.get",
        lambda _: _resp(content=mini_component_cif_gz.read_bytes()),
    )
    components_path = io.download_components_cif(data_dir=tmp_path)
    assert components_path.is_file()
    assert (components_path.parent / "components.parquet").is_file()
    assert len(pd.read_parquet(components_path.parent / "components.parquet").index)


def test_download_components_cif_cached(
    tmp_path, mini_component_cif, mini_components_pqt
):
    components_path = tmp_path / "dbs" / "components" / "components.cif"
    components_pqt = tmp_path / "dbs" / "components" / "components.parquet"
    components_path.parent.mkdir(parents=True)
    components_path.write_bytes(mini_component_cif.read_bytes())
    components_pqt.write_bytes(mini_components_pqt.read_bytes())
    components_path = io.download_components_cif(data_dir=tmp_path)
    assert components_path.is_file()
    assert components_pqt.is_file()
    assert len(pd.read_parquet(components_path.parent / "components.parquet").index)


def test_ccd(test_env, mock_alternative_datasets):
    mock_alternative_datasets("test")
    reference = "C[C@H]1[C@H]([C@H]([C@@H]([C@@H](O1)OC)O)O)O"
    path = io.download_components_cif(data_dir=test_env)
    df = pd.read_parquet(path.parent / "components.parquet")
    canonical_smile = df[df["binder_id"] == "MFU"].squeeze()["canonical_smiles"]
    assert canonical_smile == reference


def test_download_uniprot_fasta_data_fetch(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "plinder.data.pipeline.io.requests.get", lambda _: _resp(content=b"foo")
    )
    raw_path = io.download_uniprot_fasta_data(data_dir=tmp_path)
    assert raw_path.is_file()
    assert raw_path.read_text() == "foo"


def test_download_uniprot_fasta_cached(tmp_path):
    raw_path = tmp_path / "dbs" / "uniprot" / "pdb_uniprot.txt.gz"
    raw_path.parent.mkdir(parents=True)
    raw_path.write_bytes(b"foo")
    raw_path = io.download_uniprot_fasta_data(data_dir=tmp_path)
    assert raw_path.read_text() == "foo"


def test_download_alphafold_cif_file(tmp_path):
    cif_file = core_io.download_alphafold_cif_file("P05067", tmp_path)
    assert cif_file is not None and cif_file.exists()
