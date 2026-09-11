from __future__ import annotations

import pandas as pd
import pyarrow.parquet as pq
import pytest
from plinder.core.utils.schemas import STRUCTURE_LINK_SCHEMA
from plinder.data.linked_apo import (
    REQUIRED_SCORE_METRICS,
    LinkedApoSelectionConfig,
    build_apo_candidate_manifest,
    select_linked_apo_structures,
    write_linked_apo_structure_table,
)

SYSTEM = "1abc__1__1.A__1.L_1.M"


def _membership(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    for column in [
        "chain_num_contacting_ions",
        "chain_num_contacting_artifacts",
        "chain_num_contacting_other_ligands",
    ]:
        if column not in frame:
            frame[column] = 0
    return frame


def _annotation(*ligand_ids: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "system_id": [SYSTEM] * len(ligand_ids),
            "ligand_id": list(ligand_ids),
            "ligand_is_proper": [True] * len(ligand_ids),
            "entry_pdb_id": ["1abc"] * len(ligand_ids),
            "system_biounit_id": ["1"] * len(ligand_ids),
            "ligand_is_ion": [False] * len(ligand_ids),
            "ligand_is_artifact": [False] * len(ligand_ids),
            "ligand_neighboring_residues": [[] for _ in ligand_ids],
            "ligand_protein_chains_asym_id": [[] for _ in ligand_ids],
        }
    )


def _holo_annotation(*chain_ids: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "entry_pdb_id": ["2def"] * len(chain_ids),
            "system_biounit_id": ["1"] * len(chain_ids),
            "ligand_id": [f"ligand-{chain_id}" for chain_id in chain_ids],
            "ligand_is_proper": [True] * len(chain_ids),
            "ligand_is_ion": [False] * len(chain_ids),
            "ligand_is_artifact": [False] * len(chain_ids),
            "ligand_neighboring_residues": [[] for _ in chain_ids],
            "ligand_protein_chains_asym_id": [
                [f"1.{chain_id}"] for chain_id in chain_ids
            ],
        }
    )


def _candidate(
    target_system: str,
    *,
    ions: int = 0,
    artifacts: int = 0,
    other_ligands: int = 0,
    resolution: float | None = 2.0,
) -> dict[str, object]:
    entry_id, chain_id = target_system.rsplit("_", maxsplit=1)
    return {
        "target_system": target_system,
        "source_entry_id": entry_id,
        "source_chain_asym_id": chain_id,
        "source_chain_auth_id": chain_id,
        "source_biounit_id": "1",
        "source_chain_instance": f"1.{chain_id}",
        "source_num_contacting_ions": ions,
        "source_num_contacting_artifacts": artifacts,
        "source_num_contacting_other_ligands": other_ligands,
        "source_resolution": resolution,
    }


def _scores(
    target_system: str,
    ligand_values: dict[str, dict[str, int]],
) -> pd.DataFrame:
    rows = []
    defaults = {
        "pocket_fident": 100,
        "protein_fident_weighted_sum": 100,
        "protein_fident_qcov_weighted_sum": 100,
        "protein_lddt_weighted_sum": 100,
    }
    for ligand_id, overrides in ligand_values.items():
        for metric, similarity in (defaults | overrides).items():
            rows.append(
                {
                    "query_system": SYSTEM,
                    "query_ligand_id": ligand_id,
                    "target_system": target_system,
                    "metric": metric,
                    "similarity": similarity,
                }
            )
    return pd.DataFrame(rows)


def test_candidate_manifest_matches_apo_database_definition() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"] * 7,
            "chain_asym_id": list("ABCDEFG"),
            "chain_auth_id": list("RSTUVWX"),
            "chain_entity_id": ["1", "2", "1", "3", "4", "5", "6"],
            "chain_receptor_type": [
                "protein",
                "protein",
                "protein",
                "dna",
                "protein",
                "protein",
                "protein",
            ],
            "chain_is_holo": [True, False, False, False, True, False, False],
            "chain_is_ligand_like": [False, False, False, False, False, True, False],
        }
    )
    membership = _membership(
        pd.DataFrame(
            {
                "entry_pdb_id": ["2def"] * 6,
                "biounit_id": ["1"] * 6,
                "chain_instance": ["1.A", "1.B", "1.C", "1.D", "1.E", "1.F"],
                "chain_asym_id": list("ABCDEF"),
                "chain_role": ["receptor"] * 5 + ["ligand"],
            }
        )
    )
    metadata = pd.DataFrame({"entry_pdb_id": ["2def"], "entry_resolution": [1.8]})

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=metadata,
        annotation=_holo_annotation("A", "E"),
    )

    # B is apo. C copies holo entity 1; D is DNA; F is a peptide ligand; G is
    # absent from the deposited biological assemblies.
    assert candidates["target_system"].tolist() == ["2def_B"]
    assert candidates.loc[0, "source_chain_auth_id"] == "S"


def test_candidate_manifest_prefers_ligand_free_assembly() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["R"],
            "chain_entity_id": ["1"],
            "chain_receptor_type": ["protein"],
            "chain_is_holo": [False],
            "chain_is_ligand_like": [False],
        }
    )
    membership = _membership(
        pd.DataFrame(
            {
                "entry_pdb_id": ["2def"] * 3,
                "biounit_id": ["1", "1", "2"],
                "chain_instance": ["1.A", "1.I", "1.A"],
                "chain_asym_id": ["A", "I", "A"],
                "chain_role": ["receptor", "ligand", "receptor"],
            }
        )
    )
    membership.loc[
        membership["biounit_id"].eq("1") & membership["chain_instance"].eq("1.A"),
        "chain_num_contacting_other_ligands",
    ] = 1
    metadata = pd.DataFrame({"entry_pdb_id": ["2def"], "entry_resolution": [2.0]})

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=metadata,
        annotation=_holo_annotation(),
    )

    assert candidates.loc[0, "source_biounit_id"] == "2"
    assert candidates.loc[0, "source_num_contacting_other_ligands"] == 0


def test_candidate_manifest_uses_all_ligand_contact_counts() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["R"],
            "chain_entity_id": ["1"],
            "chain_receptor_type": ["protein"],
            "chain_is_holo": [False],
            "chain_is_ligand_like": [False],
        }
    )
    membership = _membership(
        pd.DataFrame(
            {
                "entry_pdb_id": ["2def"],
                "biounit_id": ["1"],
                "chain_instance": ["1.A"],
                "chain_asym_id": ["A"],
                "chain_role": ["receptor"],
            }
        )
    )
    membership.loc[0, "chain_num_contacting_ions"] = 1
    membership.loc[0, "chain_num_contacting_artifacts"] = 1
    membership.loc[0, "chain_num_contacting_other_ligands"] = 1

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=pd.DataFrame(
            {
                "entry_pdb_id": ["2def"],
                "entry_resolution": [1.8],
            }
        ),
        annotation=_holo_annotation(),
    )

    assert candidates.loc[0, "source_num_contacting_ions"] == 1
    assert candidates.loc[0, "source_num_contacting_artifacts"] == 1
    assert candidates.loc[0, "source_num_contacting_other_ligands"] == 1


def test_candidate_manifest_ranks_repeated_chain_instances_before_selection() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"],
            "chain_asym_id": ["A"],
            "chain_auth_id": ["R"],
            "chain_entity_id": ["1"],
            "chain_receptor_type": ["protein"],
            "chain_is_ligand_like": [False],
        }
    )
    membership = _membership(
        pd.DataFrame(
            {
                "entry_pdb_id": ["2def", "2def"],
                "biounit_id": ["1", "1"],
                "chain_instance": ["1.A", "2.A"],
                "chain_asym_id": ["A", "A"],
                "chain_role": ["receptor", "receptor"],
                "chain_num_contacting_ions": [1, 0],
            }
        )
    )

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=pd.DataFrame(
            {
                "entry_pdb_id": ["2def"],
                "entry_resolution": [1.8],
            }
        ),
        annotation=_holo_annotation(),
    )

    assert candidates.loc[0, "source_chain_instance"] == "2.A"
    assert candidates.loc[0, "source_num_contacting_ions"] == 0


def test_candidate_manifest_uses_ligand_holo_rows_not_chain_flag() -> None:
    chains = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"] * 3,
            "chain_asym_id": ["A", "B", "C"],
            "chain_auth_id": ["A", "B", "C"],
            "chain_entity_id": ["1", "1", "2"],
            "chain_receptor_type": ["protein"] * 3,
            "chain_is_holo": [True] * 3,
            "chain_is_ligand_like": [False] * 3,
        }
    )
    membership = _membership(
        pd.DataFrame(
            {
                "entry_pdb_id": ["2def"] * 3,
                "biounit_id": ["1"] * 3,
                "chain_instance": ["1.A", "1.B", "1.C"],
                "chain_asym_id": ["A", "B", "C"],
                "chain_role": ["receptor"] * 3,
            }
        )
    )

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=pd.DataFrame(
            {
                "entry_pdb_id": ["2def"],
                "entry_resolution": [1.8],
            }
        ),
        annotation=_holo_annotation("A"),
    )

    assert candidates["target_system"].tolist() == ["2def_C"]


def test_requires_every_proper_ligand_pocket() -> None:
    scores = _scores("2def_A", {"ligand-1": {}, "ligand-2": {"pocket_fident": 94}})

    links = select_linked_apo_structures(
        scores,
        annotation=_annotation("ligand-1", "ligand-2"),
        candidates=pd.DataFrame([_candidate("2def_A")]),
    )

    assert links.empty


def test_rejects_candidate_when_a_proper_ligand_is_absent() -> None:
    links = select_linked_apo_structures(
        _scores("2def_A", {"ligand-1": {}}),
        annotation=_annotation("ligand-1", "ligand-2"),
        candidates=pd.DataFrame([_candidate("2def_A")]),
    )

    assert links.empty


def test_ignores_nonproper_ligands() -> None:
    annotation = pd.concat(
        [
            _annotation("ligand-1"),
            pd.DataFrame(
                {
                    "system_id": [SYSTEM],
                    "ligand_id": ["artifact"],
                    "ligand_is_proper": [False],
                }
            ),
        ],
        ignore_index=True,
    )

    links = select_linked_apo_structures(
        _scores("2def_A", {"ligand-1": {}}),
        annotation=annotation,
        candidates=pd.DataFrame([_candidate("2def_A")]),
    )

    assert links["linked_structure_id"].tolist() == ["2def_A"]
    assert links.loc[0, "num_ligand_pockets"] == 1


def test_applies_every_protein_threshold() -> None:
    config = LinkedApoSelectionConfig(
        min_pocket_fident=0,
        min_protein_fident_weighted_sum=0,
        min_protein_fident_qcov_weighted_sum=81,
        min_protein_lddt_weighted_sum=0,
    )
    links = select_linked_apo_structures(
        _scores(
            "2def_A",
            {"ligand-1": {"protein_fident_qcov_weighted_sum": 80}},
        ),
        annotation=_annotation("ligand-1"),
        candidates=pd.DataFrame([_candidate("2def_A")]),
        config=config,
    )

    assert links.empty


def test_prefers_chain_local_ligand_class_before_resolution() -> None:
    candidates = pd.DataFrame(
        [
            _candidate("2def_A", other_ligands=1, resolution=1.0),
            _candidate("3ghi_B", artifacts=1, resolution=1.0),
            _candidate("4jkl_C", ions=1, resolution=1.0),
            _candidate("5mno_D", resolution=3.0),
        ]
    )
    scores = pd.concat(
        [_scores(target, {"ligand-1": {}}) for target in candidates["target_system"]],
        ignore_index=True,
    )

    links = select_linked_apo_structures(
        scores,
        annotation=_annotation("ligand-1"),
        candidates=candidates,
    )

    assert links["linked_structure_id"].tolist() == [
        "5mno_D",
        "4jkl_C",
        "3ghi_B",
        "2def_A",
    ]
    assert links["rank"].tolist() == [1, 2, 3, 4]


def test_resolutionless_apo_is_a_fallback() -> None:
    candidates = pd.DataFrame(
        [
            _candidate("2def_A", resolution=None),
            _candidate("3ghi_B", resolution=3.0),
        ]
    )
    scores = pd.concat(
        [_scores(target, {"ligand-1": {}}) for target in candidates["target_system"]],
        ignore_index=True,
    )

    links = select_linked_apo_structures(
        scores,
        annotation=_annotation("ligand-1"),
        candidates=candidates,
    )

    assert links["linked_structure_id"].tolist() == ["3ghi_B", "2def_A"]


def test_excludes_apo_chain_from_same_entry_as_holo_system() -> None:
    links = select_linked_apo_structures(
        _scores("1abc_B", {"ligand-1": {}}),
        annotation=_annotation("ligand-1"),
        candidates=pd.DataFrame([_candidate("1abc_B")]),
    )

    assert links.empty


def test_limits_links_per_holo_system() -> None:
    candidates = pd.DataFrame(
        [
            _candidate("2def_A", resolution=1.0),
            _candidate("3ghi_A", resolution=2.0),
        ]
    )
    scores = pd.concat(
        [_scores(target, {"ligand-1": {}}) for target in candidates["target_system"]],
        ignore_index=True,
    )

    links = select_linked_apo_structures(
        scores,
        annotation=_annotation("ligand-1"),
        candidates=candidates,
        config=LinkedApoSelectionConfig(max_per_system=1),
    )

    assert links["linked_structure_id"].tolist() == ["2def_A"]


def test_rejects_duplicate_required_score_rows() -> None:
    scores = _scores("2def_A", {"ligand-1": {}})
    scores = pd.concat([scores, scores.iloc[[0]]], ignore_index=True)

    with pytest.raises(ValueError, match="duplicate required metrics"):
        select_linked_apo_structures(
            scores,
            annotation=_annotation("ligand-1"),
            candidates=pd.DataFrame([_candidate("2def_A")]),
        )


def test_score_path_reads_only_required_metrics(tmp_path, monkeypatch) -> None:
    score_path = tmp_path / "apo_scores.parquet"
    _scores("2def_A", {"ligand-1": {}}).to_parquet(score_path, index=False)
    original_read_parquet = pd.read_parquet
    score_filters = []

    def read_parquet(path, *args, **kwargs):
        if path == score_path:
            score_filters.append(kwargs.get("filters"))
        return original_read_parquet(path, *args, **kwargs)

    monkeypatch.setattr(pd, "read_parquet", read_parquet)

    links = select_linked_apo_structures(
        score_path,
        annotation=_annotation("ligand-1"),
        candidates=pd.DataFrame([_candidate("2def_A")]),
    )

    assert links["linked_structure_id"].tolist() == ["2def_A"]
    assert score_filters == [[("metric", "in", list(REQUIRED_SCORE_METRICS))]]


def test_score_dataset_directory_is_supported(tmp_path) -> None:
    score_dir = tmp_path / "scores" / "shard=de"
    score_dir.mkdir(parents=True)
    _scores("2def_A", {"ligand-1": {}}).to_parquet(
        score_dir / "part.parquet", index=False
    )

    links = select_linked_apo_structures(
        score_dir.parent,
        annotation=_annotation("ligand-1"),
        candidates=pd.DataFrame([_candidate("2def_A")]),
    )

    assert links["linked_structure_id"].tolist() == ["2def_A"]


def test_writes_declared_release_schema(tmp_path) -> None:
    output = tmp_path / "linked_apo_structures.parquet"

    result = write_linked_apo_structure_table(
        _scores("2def_A", {"ligand-1": {}}),
        annotation=_annotation("ligand-1"),
        candidates=pd.DataFrame([_candidate("2def_A")]),
        output_path=output,
    )

    assert result == output
    assert pq.read_schema(output) == STRUCTURE_LINK_SCHEMA
    links = pd.read_parquet(output)
    assert links.loc[0, "reference_system_id"] == SYSTEM
    assert links.loc[0, "min_pocket_fident"] == 100


def test_parquet_writer_matches_in_memory_selection(tmp_path) -> None:
    annotation = _annotation("ligand-1", "ligand-2")
    candidates = pd.DataFrame(
        [
            _candidate("2def_A", other_ligands=1, resolution=1.0),
            _candidate("3ghi_B", resolution=3.0),
        ]
    )
    scores = pd.concat(
        [
            _scores(target, {"ligand-1": {}, "ligand-2": {}})
            for target in candidates["target_system"]
        ],
        ignore_index=True,
    )
    score_path = tmp_path / "scores" / "shard=de"
    annotation_path = tmp_path / "annotation.parquet"
    candidate_path = tmp_path / "candidates.parquet"
    output = tmp_path / "linked_apo_structures.parquet"
    score_path.mkdir(parents=True)
    scores.to_parquet(score_path / "part.parquet", index=False)
    annotation.to_parquet(annotation_path, index=False)
    candidates.to_parquet(candidate_path, index=False)

    write_linked_apo_structure_table(
        score_path.parent,
        annotation=annotation_path,
        candidates=candidate_path,
        output_path=output,
        scratch_dir=tmp_path / "scratch",
        threads=2,
        memory_limit="1GB",
    )

    assert pq.read_schema(output) == STRUCTURE_LINK_SCHEMA
    links = pd.read_parquet(output)
    assert links["linked_structure_id"].tolist() == ["3ghi_B", "2def_A"]
    assert links["rank"].tolist() == [1, 2]
    assert links["num_ligand_pockets"].tolist() == [2, 2]


def test_parquet_writer_rejects_duplicate_required_score_rows(tmp_path) -> None:
    scores = _scores("2def_A", {"ligand-1": {}})
    scores = pd.concat([scores, scores.iloc[[0]]], ignore_index=True)
    score_path = tmp_path / "scores.parquet"
    annotation_path = tmp_path / "annotation.parquet"
    candidate_path = tmp_path / "candidates.parquet"
    scores.to_parquet(score_path, index=False)
    _annotation("ligand-1").to_parquet(annotation_path, index=False)
    pd.DataFrame([_candidate("2def_A")]).to_parquet(candidate_path, index=False)

    with pytest.raises(ValueError, match="duplicate required metrics"):
        write_linked_apo_structure_table(
            score_path,
            annotation=annotation_path,
            candidates=candidate_path,
            output_path=tmp_path / "links.parquet",
            scratch_dir=tmp_path / "scratch",
            memory_limit="1GB",
        )


def test_parquet_writer_preserves_schema_without_apo_candidates(tmp_path) -> None:
    score_path = tmp_path / "scores.parquet"
    annotation_path = tmp_path / "annotation.parquet"
    candidate_path = tmp_path / "candidates.parquet"
    output = tmp_path / "linked_apo_structures.parquet"
    _scores("2def_A", {"ligand-1": {}}).to_parquet(score_path, index=False)
    _annotation("ligand-1").to_parquet(annotation_path, index=False)
    pd.DataFrame([_candidate("2def_A")]).iloc[:0].to_parquet(
        candidate_path, index=False
    )

    write_linked_apo_structure_table(
        score_path,
        annotation=annotation_path,
        candidates=candidate_path,
        output_path=output,
        scratch_dir=tmp_path / "scratch",
        memory_limit="1GB",
    )

    assert pq.read_schema(output) == STRUCTURE_LINK_SCHEMA
    assert pd.read_parquet(output).empty
