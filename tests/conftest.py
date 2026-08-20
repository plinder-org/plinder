# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import json
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest

test_asset_fp = Path(__file__).absolute().parent / "test_data"
test_output_fp = Path(__file__).absolute().parent / "xx/output"


def _write_test_entry_metadata(release_dir: Path) -> None:
    annotation = pd.read_parquet(release_dir / "index" / "annotation_table.parquet")
    entry_columns = [column for column in annotation if column.startswith("entry_")]
    metadata = annotation.loc[:, entry_columns].drop_duplicates()
    if metadata["entry_pdb_id"].duplicated().any():
        raise ValueError("test annotation has inconsistent entry metadata")
    dates = pd.read_csv(
        Path(__file__).resolve().parents[1]
        / "src/plinder/data/annotations/static_files/dates.csv"
    ).loc[:, ["entry_pdb_id", "entry_release_date"]]
    metadata = metadata.drop(columns="entry_release_date", errors="ignore").merge(
        dates,
        on="entry_pdb_id",
        how="left",
        validate="one_to_one",
    )
    metadata.to_parquet(release_dir / "index" / "entry_metadata.parquet", index=False)


@pytest.fixture(scope="session")
def plinder_src():
    plinder_root = Path(__file__).absolute().parent.parent
    return plinder_root


@pytest.fixture(scope="session")
def test_dir():
    return test_asset_fp


@pytest.fixture(scope="session")
def out_dir():
    return test_output_fp


@pytest.fixture(scope="session")
def cif_1qz5():
    return test_asset_fp / "xx/pdb_00001qz5/pdb_00001qz5_xyz-enrich.cif.gz"


@pytest.fixture(scope="session")
def validation_1qz5():
    return test_asset_fp / "validation/1qz5_validation.xml.gz"


@pytest.fixture(scope="session")
def cif_1qz5_unzipped():
    return test_asset_fp / "xx/pdb_00001qz5/pdb_00001qz5_xyz-enrich.cif"


@pytest.fixture(scope="session")
def cif_assembly_1qz5():
    return test_asset_fp / "xx/pdb_00001qz5/1qz5-assembly.cif"


@pytest.fixture(scope="session")
def cif_assembly_4ci1():
    return test_asset_fp / "xx/pdb_00004ci1/4ci1-assembly.cif"


@pytest.fixture(scope="session")
def cif_4ci1():
    return test_asset_fp / "xx/pdb_00004ci1/pdb_00004ci1_xyz-enrich.cif.gz"


@pytest.fixture(scope="session")
def cif_5a7w():
    return test_asset_fp / "xx/pdb_00005a7w/pdb_00005a7w_xyz-enrich.cif.gz"


@pytest.fixture(scope="session")
def cif_assembly_5a7w():
    return test_asset_fp / "xx/pdb_00005a7w/5a7w-assembly.cif"


@pytest.fixture(scope="session")
def fingerprint_prop_5a7w():
    return test_asset_fp / "xx/pdb_00005a7w/5a7w_interactions.txt"


# Hydrogenated 5a7w lig(35M)
@pytest.fixture(scope="session")
def sdf_5a7w_lig():
    return test_asset_fp / "xx/pdb_00005a7w/5a7w_A_35M.sdf"


# To test fragmented ligands, carbs
@pytest.fixture(scope="session")
def cif_6fx1():
    return test_asset_fp / "xx/pdb_00006fx1/pdb_00006fx1_xyz-enrich.cif.gz"


# To test fragmented ligands
@pytest.fixture(scope="session")
def cif_assembly_6fx1():
    return test_asset_fp / "xx/pdb_00006fx1/6fx1-assembly.cif"


# To test covalent bond annotation
@pytest.fixture(scope="session")
def cif_6f6r():
    return test_asset_fp / "xx/pdb_00006f6r/pdb_00006f6r_xyz-enrich.cif.gz"


# To test peptide - 9 a.a. with UniProt (maybe ligand?)
@pytest.fixture(scope="session")
def cif_assembly_6i41():
    return test_asset_fp / "xx/pdb_00006i41/6i41-assembly.cif"


# To test peptide - 13 a.a. with UniProt (not ligand)
@pytest.fixture(scope="session")
def cif_2p1q():
    return test_asset_fp / "xx/pdb_00002p1q/pdb_00002p1q_xyz-enrich.cif.gz"


# To test peptide - 13 a.a. no UniProt (synthetic -> ligand)
@pytest.fixture(scope="session")
def cif_6u6k():
    return test_asset_fp / "xx/pdb_00006u6k/pdb_00006u6k_xyz-enrich.cif.gz"


# To test peptide - 7 a.a. with BIRD + covalent (definitely ligand)
@pytest.fixture(scope="session")
def cif_6lu7():
    return test_asset_fp / "xx/pdb_00006lu7/pdb_00006lu7_xyz-enrich.cif.gz"


# To covalent ligand
@pytest.fixture(scope="session")
def cif_7gj7():
    return test_asset_fp / "xx/pdb_00007gj7/pdb_00007gj7_xyz-enrich.cif.gz"


# To noncovalent ligand
@pytest.fixture(scope="session")
def cif_7gl9():
    return test_asset_fp / "xx/pdb_00007gl9/pdb_00007gl9_xyz-enrich.cif.gz"


# To test covalent linkages + 5-char CCD codes
# (PCSK9 + enlicitide/MK-0616: ligand is presented as chain + covalent/modified residues)
@pytest.fixture(scope="session")
def cif_10sb():
    return test_asset_fp / "xx/pdb_000010sb/pdb_000010sb_xyz-enrich.cif.gz"


# To test peptide - 9 a.a. with UniProt (maybe ligand?)
@pytest.fixture(scope="session")
def cif_6i41():
    return test_asset_fp / "xx/pdb_00006i41/pdb_00006i41_xyz-enrich.cif.gz"


# Ligand (STI) that needs fixing - otherwise invalid
@pytest.fixture(scope="session")
def cif_2hyy():
    return test_asset_fp / "xx/pdb_00002hyy/pdb_00002hyy_xyz-enrich.cif.gz"


# Ligand (EF2) that needs to have good valency read from cif structure
@pytest.fixture(scope="session")
def cif_7bqu():
    return test_asset_fp / "xx/pdb_00007bqu/pdb_00007bqu_xyz-enrich.cif.gz"


# Ligand (JEF) that is missing atoms (unresolved)
@pytest.fixture(scope="session")
def cif_1ngx():
    return test_asset_fp / "xx/pdb_00001ngx/pdb_00001ngx_xyz-enrich.cif.gz"


# Ligand (FAD) with resolved but distorted geometry
@pytest.fixture(scope="session")
def cif_3grt():
    return test_asset_fp / "xx/pdb_00003grt/pdb_00003grt_xyz-enrich.cif.gz"


# Test inding affinity
@pytest.fixture(scope="session")
def cif_4jvn():
    return test_asset_fp / "xx/pdb_00004jvn/pdb_00004jvn_xyz-enrich.cif.gz"


# Ligand (peptidic) that is disconnected by bonds
@pytest.fixture(scope="session")
def cif_4nhc():
    return test_asset_fp / "xx/pdb_00004nhc/pdb_00004nhc_xyz-enrich.cif.gz"


# To test dna
@pytest.fixture(scope="session")
def cif_assembly_5fkw():
    return test_asset_fp / "xx/pdb_00005fkw/5fkw-assembly.cif"


# To test dna
@pytest.fixture(scope="session")
def cif_5fkw():
    return test_asset_fp / "xx/pdb_00005fkw/5fkw-assembly.cif"


# To test covalent bond annotation
@pytest.fixture(scope="session")
def cif_metadata():
    return test_asset_fp / "xx/pdb_00006f6r/6f6r-metadata.cif"


# To weird connectivity upon replacing missing atoms
@pytest.fixture(scope="session")
def cif_6ue5():
    return test_asset_fp / "xx/pdb_00006ue5/pdb_00006ue5_xyz-enrich.cif.gz"


# To test missing protein chains
@pytest.fixture(scope="session")
def cif_1utr():
    return test_asset_fp / "xx/pdb_00001utr/pdb_00001utr_xyz-enrich.cif.gz"


# To test rna annotation
@pytest.fixture(scope="session")
def cif_2leb():
    return test_asset_fp / "xx/pdb_00002leb/pdb_00002leb_xyz-enrich.cif.gz"


# To test binding site dna, water, and ptm
@pytest.fixture(scope="session")
def cif_5btg():
    return test_asset_fp / "xx/pdb_00005btg/pdb_00005btg_xyz-enrich.cif.gz"


# To test binding site water
@pytest.fixture(scope="session")
def cif_6jjn():
    return test_asset_fp / "xx/pdb_00006jjn/pdb_00006jjn_xyz-enrich.cif.gz"


# To test EC number
@pytest.fixture(scope="session")
def cif_3g32():
    return test_asset_fp / "xx/pdb_00003g32/pdb_00003g32_xyz-enrich.cif.gz"


@pytest.fixture(scope="session")
def cif_2y4i_system():
    return test_asset_fp / "xx/pdb_00002y4i/pdb_00002y4i_xyz-enrich.cif.gz"


# TODO: PLIP is no longer used — these fixtures test interaction detection (now via peppr)
# CHK1 inhib 1
@pytest.fixture(scope="session")
def cif_2gdo():
    return test_asset_fp / "xx/pdb_00002gdo/pdb_00002gdo_xyz-enrich.cif.gz"


# CHK1 inhib 2
@pytest.fixture(scope="session")
def cif_4qyf():
    return test_asset_fp / "xx/pdb_00004qyf/pdb_00004qyf_xyz-enrich.cif.gz"


@pytest.fixture(scope="session")
def rcsb_ccd_reference_csv():
    return test_asset_fp / "rcsb_ccd_smiles_reference.csv"


@pytest.fixture(scope="session")
def resolved_smiles_csv():
    return test_asset_fp / "resolved_smiles_reference.csv"


@pytest.fixture(scope="session")
def cif_2y4i():
    return test_asset_fp / "xx/pdb_00002y4i/pdb_00002y4i_xyz-enrich.cif.gz"


# To test for hydrogen removal before saving
@pytest.fixture(scope="session")
def cif_7az3():
    return test_asset_fp / "xx/pdb_00007az3/pdb_00007az3_xyz-enrich.cif.gz"


# To test for too many hydrogens
@pytest.fixture(scope="session")
def cif_6ntj():
    return test_asset_fp / "xx/pdb_00006ntj/pdb_00006ntj_xyz-enrich.cif.gz"


# To test nucleic acid receptor detection (issue #61)
@pytest.fixture(scope="session")
def cif_8ufz():
    return test_asset_fp / "xx/pdb_00008ufz/pdb_00008ufz_xyz-enrich.cif.gz"


# To test multi-ligand system grouping (GPCR with adjacent binding sites)
@pytest.fixture(scope="session")
def cif_7fee():
    return test_asset_fp / "xx/pdb_00007fee/pdb_00007fee_xyz-enrich.cif.gz"


# To test cofactor-only system classification (HEM in hemoglobin)
@pytest.fixture(scope="session")
def cif_19hc():
    return test_asset_fp / "xx/pdb_000019hc/pdb_000019hc_xyz-enrich.cif.gz"


# To test ATP+metal cofactor system grouping (PKA)
@pytest.fixture(scope="session")
def cif_1atp():
    return test_asset_fp / "xx/pdb_00001atp/pdb_00001atp_xyz-enrich.cif.gz"


@pytest.fixture(scope="session")
def mini_component_cif():
    return test_asset_fp / "components.cif"


@pytest.fixture(scope="session")
def mini_component_cif_gz():
    return test_asset_fp / "components.cif.gz"


@pytest.fixture(scope="session")
def mini_components_pqt():
    return test_asset_fp / "components.parquet"


@pytest.fixture
def test_env(tmp_path, monkeypatch):
    monkeypatch.setenv("PLINDER_MOUNT", tmp_path.as_posix())
    monkeypatch.setenv("PLINDER_BUCKET", "bucket")
    monkeypatch.setenv("PLINDER_RELEASE", "test")
    monkeypatch.setenv("PLINDER_RELEASE_NUMBER", "")
    from plinder.core.utils import config

    config._config._clear()
    return tmp_path / "bucket" / "test"


@pytest.fixture
def components_path(test_env, mini_component_cif, mini_components_pqt):
    components_path = test_env / "dbs" / "components" / "components.cif"
    components_pqt = test_env / "dbs" / "components" / "components.parquet"
    components_path.parent.mkdir(parents=True)
    components_path.write_text(mini_component_cif.read_text())
    components_pqt.write_bytes(mini_components_pqt.read_bytes())
    return components_path


@pytest.fixture
def cofactors_path(test_env):
    cofactors_path = test_env / "dbs" / "cofactors" / "cofactors.json"
    cofactors_path.parent.mkdir(parents=True)
    mini_cofactors = """\
{
    "Coenzyme A": [
        {
            "cofactors": [
                "01A",
                "EKY"
            ],
            "EC": [
                "1.1.1.34"
            ]
        }
    ],
    "Orthoquinone residues (LTQ, TTQ, CTQ)": [
        {
            "cofactors": [
                "0AF",
                "TOQ",
                "TQQ",
                "TRQ"
            ],
            "EC": [
                "1.4.99.3"
            ]
        }
    ],
    "Heme": [
        {
            "cofactors": [
                "HEM",
                "HEC",
                "HEB",
                "HEA"
            ],
            "EC": [
                "1.11.1.5",
                "1.11.2.2"
            ]
        }
    ],
    "Adenosine nucleotides": [
        {
            "cofactors": [
                "ATP",
                "ADP",
                "AMP"
            ],
            "EC": [
                "2.7.1.1",
                "2.7.11.1"
            ]
        }
    ]
}"""
    cofactors_path.write_text(mini_cofactors)
    with cofactors_path.open("r") as f:
        cofactors = json.load(f)
        assert len(cofactors) == 4  # CoA, TTQ, Heme (19hc), ATP (1atp)
    return cofactors_path


@pytest.fixture
def affinity_path(test_env):
    affinity_path = test_env / "dbs" / "affinity" / "affinity.json"
    affinity_path.parent.mkdir(parents=True)
    affinity = """\
{
"pchembl": {
    "4JVM_XDI": 5.3979400087,
    "4JVN_YUG": 7.638272164,
    "4JVO_A5A": 5.7706810755,
    "4JVP_SO4": 1.2732414543,
    "4JVQ_1ML": 6.2006594505,
    "4JVR_1MT": 8.0268721464}
}"""
    affinity = json.loads(affinity)
    with open(affinity_path, "w") as f:
        json.dump(affinity, f)
    return affinity_path


@pytest.fixture
def seqres_path(test_env):
    seqres_path = test_env / "dbs" / "seqres" / "pdb_seqres.txt.gz"
    fixture_path = test_asset_fp / "pdb_seqres.txt.gz"
    seqres_path.parent.mkdir(parents=True)
    seqres_path.write_bytes(fixture_path.read_bytes())
    return seqres_path


@pytest.fixture
def mock_alternative_datasets(
    test_env,
    seqres_path,
    components_path,
    cofactors_path,
    affinity_path,
):
    def inner(pdb_id: str):
        entry_dir = test_env / "raw_entries" / pdb_id[-3:-1]
        entry_dir.mkdir(parents=True)
        return entry_dir

    return inner


@pytest.fixture
def read_plinder_mount(monkeypatch, tmp_path):
    source = test_asset_fp / "plinder" / "mount"
    adir = tmp_path / "plinder" / "mount"
    shutil.copytree(source, adir, copy_function=os.symlink)
    _write_test_entry_metadata(adir)

    monkeypatch.setenv("PLINDER_MOUNT", tmp_path.as_posix())
    monkeypatch.setenv("PLINDER_RELEASE", "mount")
    monkeypatch.setenv("PLINDER_RELEASE_NUMBER", "")
    monkeypatch.setenv("PLINDER_BUCKET", "plinder")
    monkeypatch.setenv("PLINDER_OFFLINE", "true")
    from plinder.core.utils import config, cpl

    config._config._clear()
    monkeypatch.setattr(cpl, "_CLIENTS", {})
    cfg = config.get_config()
    assert Path(cfg.data.plinder_dir) == adir

    for path in adir.rglob("*_done"):
        path.unlink()

    return adir


@pytest.fixture
def read_plinder_eval_mount(monkeypatch, tmp_path):
    plinder_mount = tmp_path / "plinder_mount"
    adir = plinder_mount / "eval"
    shutil.copytree(test_asset_fp / "eval", adir)
    annotation_path = adir / "index" / "annotation_table.parquet"
    annotation = pd.read_parquet(annotation_path)
    annotation = annotation.rename(
        columns={"ligand_rdkit_canonical_smiles": "ligand_smiles"}
    )
    instance_chains = annotation["ligand_id"].str.rsplit("__", n=1).str[-1]
    annotation["ligand_instance_chain"] = instance_chains
    annotation["ligand_instance"] = instance_chains.str.split(".").str[0].astype(int)
    annotation["ligand_asym_id"] = instance_chains.str.rsplit(".", n=1).str[-1]
    annotation.to_parquet(annotation_path, index=False)
    for archive in ("a3.zip", "ai.zip"):
        shutil.unpack_archive(
            adir / "systems" / archive,
            adir / "reconstructed_systems",
        )
    _write_test_entry_metadata(adir)
    ligand_archive_dir = adir / "ligand_archives"
    ligand_archive_dir.mkdir()
    for system_id in ("1a3b__1__1.B__1.D", "1ai5__1__1.A_1.B__1.D"):
        pdb_id = system_id[:4]
        ligand_file = (
            adir / "reconstructed_systems" / system_id / "ligand_files" / "1.D.sdf"
        )
        pd.DataFrame(
            {
                "pdb_id": [pdb_id],
                "ligand_asym_id": ["D"],
                "sdf": [ligand_file.read_bytes()],
            }
        ).to_parquet(ligand_archive_dir / f"{pdb_id[1:3]}.parquet", index=False)
    monkeypatch.setenv("PLINDER_MOUNT", plinder_mount.as_posix())
    monkeypatch.setenv("PLINDER_RELEASE", "")
    monkeypatch.setenv("PLINDER_RELEASE_NUMBER", "")
    monkeypatch.setenv("PLINDER_BUCKET", "eval")
    monkeypatch.setenv("PLINDER_OFFLINE", True)
    from plinder.core.utils import config, cpl

    config._config._clear()
    monkeypatch.setattr(cpl, "_CLIENTS", {})
    cfg = config.get_config()
    assert Path(cfg.data.plinder_dir) == adir

    return adir


@pytest.fixture
def write_plinder_mount(monkeypatch, tmp_path):
    read_plinder_mount = test_asset_fp / "plinder" / "mount"
    write_plinder_mount = tmp_path / "plinder" / "mount"
    write_plinder_mount.mkdir(parents=True)
    monkeypatch.setenv("PLINDER_MOUNT", tmp_path.as_posix())
    monkeypatch.setenv("PLINDER_RELEASE", "mount")
    monkeypatch.setenv("PLINDER_RELEASE_NUMBER", "")
    monkeypatch.setenv("PLINDER_BUCKET", "plinder")
    from plinder.core.utils import config, cpl

    config._config._clear()
    monkeypatch.setattr(cpl, "_CLIENTS", {})
    for path in read_plinder_mount.rglob("*"):
        if path.is_dir() or path.name.endswith("_done"):
            continue
        write_path = write_plinder_mount / path.relative_to(read_plinder_mount)
        write_path.parent.mkdir(exist_ok=True, parents=True)
        write_path.write_bytes(path.read_bytes())
    _write_test_entry_metadata(write_plinder_mount)
    return write_plinder_mount


@pytest.fixture
def cached_plinder_system(read_plinder_mount, tmp_path):
    """Build a system using explicit current-release cache locations."""
    from plinder.core import PlinderSystem

    def build(system_id: str) -> PlinderSystem:
        source = test_asset_fp / "reconstructed_systems" / system_id
        reconstruction_dir = tmp_path / "reconstructed_systems" / system_id
        shutil.copytree(source, reconstruction_dir)

        canonical_ligand_dir = tmp_path / "canonical_ligands" / system_id
        canonical_ligand_dir.mkdir(parents=True)
        for ligand_file in (source / "ligand_files").glob("*.sdf"):
            asym_id = ligand_file.stem.rsplit(".", maxsplit=1)[-1]
            target = canonical_ligand_dir / f"{asym_id}.sdf"
            if not target.exists():
                shutil.copyfile(ligand_file, target)

        return PlinderSystem(
            system_id=system_id,
            reconstruction_dir=reconstruction_dir,
            canonical_ligand_dir=canonical_ligand_dir,
        )

    return build


@pytest.fixture(autouse=True)
def mock_ccd_lookups(monkeypatch):
    from plinder.data.annotations.ligand_utils import sort_ccd_codes

    data = json.loads((test_asset_fp / "ccd_lookups.json").read_text())
    synonyms = [set(s) for s in data["ccd_synonyms"]]
    monkeypatch.setattr(
        "plinder.data.annotations.ligand_utils.CCD_SYNONYMS_DICT",
        {code: sort_ccd_codes(list(s))[0] for s in synonyms for code in s},
    )
    monkeypatch.setattr(
        "plinder.data.annotations.ligand_utils.COFACTORS",
        set(data["cofactors"]),
    )
    monkeypatch.setattr(
        "plinder.data.annotations.ligand_utils.ARTIFACTS",
        set(data["artifacts"]),
    )
    monkeypatch.setattr(
        "plinder.data.annotations.ligand_utils.BINDING_AFFINITY",
        None,
    )


@pytest.fixture(scope="session")
def system_1a3b():
    return "1a3b__1__1.B__1.D"


@pytest.fixture(scope="session")
def predicted_pose_1a3b():
    return test_asset_fp / "eval/predicted_poses/1a3b__1__1.B__1.D/rank1.sdf"


@pytest.fixture(scope="session")
def predicted_named_pose_1a3b():
    return test_asset_fp / "eval/predicted_poses/1a3b__1__1.B__1.D/rank1_named.sdf"


@pytest.fixture(scope="session")
def system_1ai5():
    return "1ai5__1__1.A_1.B__1.D"


@pytest.fixture(scope="session")
def predicted_pose_1ai5():
    return test_asset_fp / "eval/predicted_poses/1ai5__1__1.A_1.B__1.D/rank1.sdf"


@pytest.fixture(autouse=True)
def failfast(monkeypatch):
    monkeypatch.setattr("plinder.core.utils.io.sleep", lambda x: None)

    def f(*args, **kwargs):
        class obj:
            status_code = 200
            text = ""

            def raise_for_status(self):
                pass

        return obj()

    monkeypatch.setattr("plinder.data.pipeline.io.requests.get", f)


@pytest.fixture(scope="session")
def cif_atom_array(cif_1qz5_unzipped):
    from plinder.core.structure import atoms

    print(cif_1qz5_unzipped)
    atom_array = atoms.atom_array_from_cif_file(
        cif_1qz5_unzipped, use_author_fields=False
    )
    print(atom_array)
    return atom_array
