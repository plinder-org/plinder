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


def _annotation(*ligand_ids: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "system_id": [SYSTEM] * len(ligand_ids),
            "ligand_id": list(ligand_ids),
            "ligand_is_proper": [True] * len(ligand_ids),
        }
    )


def _candidate(
    target_system: str,
    *,
    ligand_chains: int = 0,
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
        "source_num_ligand_chains": ligand_chains,
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
            "entry_pdb_id": ["2def"] * 5,
            "chain_asym_id": list("ABCDE"),
            "chain_auth_id": list("RSTUV"),
            "chain_entity_id": ["1", "2", "1", "3", "4"],
            "chain_receptor_type": ["protein", "protein", "protein", "dna", "protein"],
            "chain_is_holo": [True, False, False, False, True],
        }
    )
    membership = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"] * 5,
            "biounit_id": ["1"] * 5,
            "chain_instance": ["1.A", "1.B", "1.C", "1.D", "1.E"],
            "chain_asym_id": list("ABCDE"),
            "chain_role": ["receptor"] * 5,
        }
    )
    metadata = pd.DataFrame(
        {"entry_pdb_id": ["2def"], "entry_resolution": [1.8]}
    )

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=metadata,
    )

    # B is apo. C is rejected because it copies holo entity 1; D is DNA.
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
        }
    )
    membership = pd.DataFrame(
        {
            "entry_pdb_id": ["2def"] * 3,
            "biounit_id": ["1", "1", "2"],
            "chain_instance": ["1.A", "1.I", "1.A"],
            "chain_asym_id": ["A", "I", "A"],
            "chain_role": ["receptor", "ligand", "receptor"],
        }
    )
    metadata = pd.DataFrame(
        {"entry_pdb_id": ["2def"], "entry_resolution": [2.0]}
    )

    candidates = build_apo_candidate_manifest(
        chains,
        biounit_chains=membership,
        entry_metadata=metadata,
    )

    assert candidates.loc[0, "source_biounit_id"] == "2"
    assert candidates.loc[0, "source_num_ligand_chains"] == 0


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


def test_prefers_empty_apo_assembly_before_resolution() -> None:
    candidates = pd.DataFrame(
        [
            _candidate("2def_A", ligand_chains=1, resolution=1.0),
            _candidate("3ghi_B", ligand_chains=0, resolution=3.0),
            _candidate("4jkl_C", ligand_chains=0, resolution=2.0),
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

    assert links["linked_structure_id"].tolist() == ["4jkl_C", "3ghi_B", "2def_A"]
    assert links["rank"].tolist() == [1, 2, 3]


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
        [_candidate("2def_A", resolution=1.0), _candidate("3ghi_A", resolution=2.0)]
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
    assert score_filters == [
        [("metric", "in", list(REQUIRED_SCORE_METRICS))]
    ]


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
