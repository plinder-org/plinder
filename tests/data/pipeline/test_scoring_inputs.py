import json
import sys
from shutil import copytree

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from plinder.core.utils.files import file_sha256
from plinder.data.annotations.interface_utils import INTERFACE_ANNOTATION_SCHEMA
from plinder.data.pipeline import collate, score, tasks


def _write_inputs(root, *, revised=False, ligands=True, interfaces=True):
    index = root / "index"
    index.mkdir(parents=True, exist_ok=True)
    schema = pa.schema(
        [
            ("entry_pdb_id", pa.string()),
            ("system_id", pa.string()),
            ("system_type", pa.string()),
            ("ligand_id", pa.string()),
            ("ligand_asym_id", pa.string()),
            ("ligand_is_proper", pa.bool_()),
            ("ligand_is_shape_comparable", pa.bool_()),
            ("ligand_protein_chains_asym_id", pa.list_(pa.string())),
            ("ligand_neighboring_residues", pa.list_(pa.string())),
            ("ligand_interacting_residues", pa.list_(pa.string())),
            ("ligand_interactions", pa.list_(pa.string())),
        ]
    )
    asym, residue = ("M", 20) if revised else ("L", 10)
    pocket = f"1.A_{residue}_{residue - 1}_{residue + 100}"
    rows = (
        [
            {
                "entry_pdb_id": "1abc",
                "system_id": f"1abc__1__1.A__1.{asym}",
                "system_type": "holo",
                "ligand_id": f"1abc__1__1.{asym}",
                "ligand_asym_id": asym,
                "ligand_is_proper": True,
                "ligand_is_shape_comparable": True,
                "ligand_protein_chains_asym_id": ["1.A"],
                "ligand_neighboring_residues": [pocket],
                "ligand_interacting_residues": [pocket],
                "ligand_interactions": [f"1.A_{residue}_hydrogen_bonds"],
            }
        ]
        if ligands
        else []
    )
    pq.write_table(
        pa.Table.from_pylist(rows, schema=schema), index / "annotation_table.parquet"
    )
    other = "3ghi" if revised else "2def"
    pq.write_table(
        pa.table(
            {
                "entry_pdb_id": ["1abc", other, other],
                "chain_asym_id": ["A", "A", "B"],
                "chain_auth_id": ["X", "Y", "Z"],
                "chain_receptor_type": ["protein"] * 3,
            }
        ),
        index / "entry_chains.parquet",
    )
    rows = (
        [
            {
                "entry_pdb_id": other,
                "system_id": f"{other}__1__1.A--1.B",
                "system_biounit_id": "1",
                "interface_chain_1": "1.A",
                "interface_chain_2": "1.B",
                "interface_chain_1_residue_numbers": [30],
                "interface_chain_1_residue_indices": [29],
                "interface_chain_2_residue_numbers": [40],
                "interface_chain_2_residue_indices": [39],
                "interface_num_contact_residue_pairs": 1,
            }
        ]
        if interfaces
        else []
    )
    pq.write_table(
        pa.Table.from_pylist(rows, schema=INTERFACE_ANNOTATION_SCHEMA),
        index / "interface_annotation_table.parquet",
    )


@pytest.mark.parametrize(
    "ligands,interfaces", [(True, True), (False, True), (True, False), (False, False)]
)
def test_prepare_scoring_inputs_refreshes_and_resumes(
    tmp_path, monkeypatch, capsys, ligands, interfaces
):
    workspace = tmp_path / "updated"
    _write_inputs(workspace)
    marker = workspace / "index/collation.json"
    marker.write_text(json.dumps({"status": collate.REPAIR_REQUIRED_STATUS}))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "score",
            "prepare-scoring-inputs",
            str(workspace),
            "--threads",
            "2",
            "--memory-limit",
            "1GB",
        ],
    )
    configurations = []
    configure = collate._configure_duckdb

    def record_configuration(connection, **kwargs):
        configure(connection, **kwargs)
        configurations.append(kwargs)

    monkeypatch.setattr(collate, "_configure_duckdb", record_configuration)
    score.main()
    assert len(configurations) == 3
    assert all(c["threads"] == 2 and c["memory_limit"] == "1GB" for c in configurations)
    assert json.loads(capsys.readouterr().out)["status"] == "complete"

    # Replace a ligand, remove an entry and add an interface-only entry.
    _write_inputs(workspace, revised=True, ligands=ligands, interfaces=interfaces)
    fresh = tmp_path / "fresh"
    copytree(workspace / "index", fresh / "index")
    for path in (fresh / "index").glob("*.parquet"):
        if path.name not in {
            "annotation_table.parquet",
            "entry_chains.parquet",
            "interface_annotation_table.parquet",
        }:
            path.unlink()
    inputs = [
        workspace / "index" / name
        for name in (
            "annotation_table.parquet",
            "entry_chains.parquet",
            "interface_annotation_table.parquet",
            "collation.json",
        )
    ]
    before = {p: file_sha256(p) for p in inputs}
    score.main()
    tasks.make_alignment_chain_lookup(
        data_dir=fresh,
        scratch_dir=None,
        threads=1,
        memory_limit="1GB",
        force_update=True,
    )
    outputs = [
        tasks.ALIGNMENT_CHAIN_LOOKUP_RELATIVE,
        tasks.LIGAND_POCKET_REPRESENTATIVES_RELATIVE,
        tasks.LIGAND_POCKET_MEMBERSHIP_RELATIVE,
        tasks.LIGAND_POCKET_RESIDUES_RELATIVE,
        tasks.INTERFACE_REPRESENTATIVES_RELATIVE,
        tasks.INTERFACE_HALF_REPRESENTATIVES_RELATIVE,
        tasks.INTERFACE_MEMBERSHIP_RELATIVE,
    ]
    for relative in outputs:
        assert pq.read_table(workspace / relative).equals(
            pq.read_table(fresh / relative)
        )
    lookup = pq.read_table(
        workspace / tasks.ALIGNMENT_CHAIN_LOOKUP_RELATIVE
    ).to_pylist()
    assert {row["entry_pdb_id"] for row in lookup} == {"1abc", "3ghi"}
    assert lookup[0]["selected_residue_numbers"] == ([20] if ligands else [])
    assert lookup[1]["selected_residue_numbers"] == ([30] if interfaces else [])
    mtimes = {p: (workspace / p).stat().st_mtime_ns for p in outputs}
    score.main()
    assert mtimes == {p: (workspace / p).stat().st_mtime_ns for p in outputs}
    (workspace / tasks.LIGAND_POCKET_RESIDUES_RELATIVE).unlink()
    score.main()
    assert tasks._completed_alignment_chain_lookup(workspace) is not None
    assert before == {p: file_sha256(p) for p in inputs}
