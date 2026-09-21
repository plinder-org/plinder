from __future__ import annotations

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from plinder.core.utils.schemas import (
    INTERFACE_APO_LINK_SCHEMA,
    mapped_alignment_schema,
)
from plinder.data.interface_apo import (
    InterfaceApoSelectionConfig,
    build_interface_apo_candidate_manifest,
    build_interface_apo_query_manifest,
    write_interface_apo_structure_table,
)


def test_interface_query_manifest_has_one_row_per_side() -> None:
    queries = build_interface_apo_query_manifest(
        pd.DataFrame(
            {
                "entry_pdb_id": ["1abc"],
                "system_id": ["1abc__1__1.A--2.B"],
                "interface_chain_1": ["1.A"],
                "interface_chain_2": ["2.B"],
            }
        )
    )

    assert queries.to_dict("records") == [
        {
            "query_entry": "1abc",
            "reference_system_id": "1abc__1__1.A--2.B",
            "reference_side": 1,
            "reference_chain_instance": "1.A",
            "reference_chain_asym_id": "A",
        },
        {
            "query_entry": "1abc",
            "reference_system_id": "1abc__1__1.A--2.B",
            "reference_side": 2,
            "reference_chain_instance": "2.B",
            "reference_chain_asym_id": "B",
        },
    ]


def test_interface_apo_candidates_require_a_protein_contact_free_instance() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["2def", "2def", "2def"],
            "chain_asym_id": ["A", "B", "C"],
            "chain_auth_id": ["R", "S", "T"],
            "chain_receptor_type": ["protein", "protein", "protein"],
            "chain_is_ligand_like": [False, True, False],
        }
    )
    membership = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"] * 5,
            "biounit_id": ["1", "2", "1", "1", "1"],
            "chain_instance": ["1.A", "1.A", "1.B", "1.C", "2.C"],
            "chain_asym_id": ["A", "A", "B", "C", "C"],
            "chain_role": ["receptor", "receptor", "ligand", "receptor", "receptor"],
            "chain_num_contacting_proteins": [1, 0, 0, 1, 0],
            "chain_num_contacting_ions": [0, 1, 0, 0, 0],
            "chain_num_contacting_artifacts": [0, 0, 0, 0, 1],
            "chain_num_contacting_other_ligands": [0, 0, 0, 0, 0],
        }
    )
    candidates = build_interface_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=pd.DataFrame(
            {"entry_pdb_id": ["2def"], "entry_resolution": [1.8]}
        ),
    )

    assert candidates["target_system"].tolist() == ["2def_A", "2def_C"]
    assert candidates["source_chain_instance"].tolist() == ["1.A", "2.C"]
    assert candidates["source_num_contacting_proteins"].tolist() == [0, 0]

    ligand_contact_columns = [
        "chain_num_contacting_ions",
        "chain_num_contacting_artifacts",
        "chain_num_contacting_other_ligands",
    ]
    membership[ligand_contact_columns] = pd.NA
    candidates = build_interface_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=pd.DataFrame(
            {"entry_pdb_id": ["2def"], "entry_resolution": [1.8]}
        ),
    )
    assert candidates["target_system"].tolist() == ["2def_A", "2def_C"]
    assert (
        candidates[
            [
                "source_num_contacting_ions",
                "source_num_contacting_artifacts",
                "source_num_contacting_other_ligands",
            ]
        ]
        .isna()
        .all(axis=None)
    )

    membership["chain_num_contacting_proteins"] = pd.Series(
        [pd.NA, 0, 0, 1, 0], dtype="Int64"
    )
    with pytest.raises(
        ValueError,
        match="biological-assembly membership has invalid chain_num_contacting_proteins",
    ):
        build_interface_apo_candidate_manifest(
            chains,
            biounit_chains=membership,
            entry_metadata=pd.DataFrame(
                {"entry_pdb_id": ["2def"], "entry_resolution": [1.8]}
            ),
        )


def _alignment_row(
    query_chain: str,
    target_entry: str,
    target_chain: str,
    *,
    qcov: float = 1.0,
    tcov: float = 1.0,
    fident: float = 1.0,
    lddt: float | None = None,
) -> dict[str, object]:
    row: dict[str, object] = {
        "query_entry": "1abc",
        "target_entry": target_entry,
        "query_chain_mapped": query_chain,
        "target_chain_mapped": target_chain,
        "source": "foldseek" if lddt is not None else "mmseqs",
        "qcov": qcov,
        "tcov": tcov,
        "fident": fident,
        "seqsim": fident,
        "query_selected_residue_numbers": [],
        "target_selected_residue_numbers": [],
        "selected_residue_identity": b"",
    }
    if lddt is not None:
        row["lddt"] = lddt
    return row


def test_interface_apo_links_apply_two_sided_coverage_and_rank_clean_chains(
    tmp_path,
) -> None:
    workspace = tmp_path / "release's files"
    workspace.mkdir()
    queries = build_interface_apo_query_manifest(
        pd.DataFrame(
            {
                "entry_pdb_id": ["1abc"],
                "system_id": ["1abc__1__1.A--1.B"],
                "interface_chain_1": ["1.A"],
                "interface_chain_2": ["1.B"],
            }
        )
    )
    query_path = workspace / "queries.parquet"
    queries.to_parquet(query_path, index=False)
    candidates = pd.DataFrame(
        {
            "target_system": ["2def_X", "3ghi_Y", "1abc_C"],
            "source_entry_id": ["2def", "3ghi", "1abc"],
            "source_chain_asym_id": ["X", "Y", "C"],
            "source_chain_auth_id": ["R", "S", "C"],
            "source_biounit_id": ["1", "1", "1"],
            "source_chain_instance": ["1.X", "1.Y", "1.C"],
            "source_num_contacting_proteins": [0, 0, 0],
            "source_num_contacting_ions": [1, 0, 0],
            "source_num_contacting_artifacts": [0, 0, 0],
            "source_num_contacting_other_ligands": [0, 0, 0],
            "source_resolution": [1.5, 2.5, 1.0],
        }
    )
    candidate_path = workspace / "candidates.parquet"
    candidates.to_parquet(candidate_path, index=False)
    mmseqs_rows = [
        _alignment_row("A", "2def", "X", fident=0.99),
        _alignment_row("A", "3ghi", "Y", fident=0.99, qcov=0.8, tcov=1.0),
        _alignment_row("A", "3ghi", "Y", fident=0.98, qcov=1.0, tcov=0.8),
        _alignment_row("A", "1abc", "C"),
        _alignment_row("B", "2def", "X", tcov=0.79),
    ]
    mmseqs_path = workspace / "mmseqs" / "shard=ab.parquet"
    mmseqs_path.parent.mkdir()
    pq.write_table(
        pa.Table.from_pylist(
            mmseqs_rows,
            schema=mapped_alignment_schema(alignment_type="mmseqs"),
        ),
        mmseqs_path,
    )
    foldseek_path = workspace / "foldseek" / "shard=ab.parquet"
    foldseek_path.parent.mkdir()
    pq.write_table(
        pa.Table.from_pylist(
            [_alignment_row("A", "3ghi", "Y", lddt=0.72)],
            schema=mapped_alignment_schema(alignment_type="foldseek"),
        ),
        foldseek_path,
    )

    output = workspace / "interface_apo_structures.parquet"
    write_interface_apo_structure_table(
        mmseqs_path.parent,
        foldseek_alignments=foldseek_path.parent,
        queries=query_path,
        candidates=candidate_path,
        output_path=output,
        config=InterfaceApoSelectionConfig(max_per_side=2),
        scratch_dir=workspace / "scratch's",
        threads=1,
        memory_limit="1GB",
    )

    links = pd.read_parquet(output)
    assert pq.read_schema(output) == INTERFACE_APO_LINK_SCHEMA
    assert links["reference_side"].tolist() == [1, 1]
    assert links["linked_structure_id"].tolist() == ["3ghi_Y", "2def_X"]
    assert links["rank"].tolist() == [1, 2]
    assert links["foldseek_lddt"].iloc[0] == pytest.approx(0.72)
    assert links["mmseqs_fident"].iloc[0] == pytest.approx(0.99)
    assert links["mmseqs_query_coverage"].iloc[0] == pytest.approx(0.8)
    assert links["mmseqs_target_coverage"].iloc[0] == pytest.approx(1.0)


def test_interface_apo_links_require_mmseqs_results_for_nonempty_inputs(
    tmp_path,
) -> None:
    queries = pd.DataFrame(
        {
            "query_entry": ["1abc"],
            "reference_system_id": ["1abc__1__1.A--1.B"],
            "reference_side": [1],
            "reference_chain_instance": ["1.A"],
            "reference_chain_asym_id": ["A"],
        }
    )
    candidates = pd.DataFrame(
        {
            "target_system": ["2def_X"],
            "source_entry_id": ["2def"],
            "source_chain_asym_id": ["X"],
            "source_chain_auth_id": ["R"],
            "source_biounit_id": ["1"],
            "source_chain_instance": ["1.X"],
            "source_num_contacting_proteins": [0],
            "source_num_contacting_ions": [0],
            "source_num_contacting_artifacts": [0],
            "source_num_contacting_other_ligands": [0],
            "source_resolution": [2.0],
        }
    )

    with pytest.raises(FileNotFoundError, match="no MMseqs alignment shards"):
        write_interface_apo_structure_table(
            tmp_path / "missing-mmseqs",
            foldseek_alignments=tmp_path / "missing-foldseek",
            queries=queries,
            candidates=candidates,
            output_path=tmp_path / "interface_apo_structures.parquet",
        )
