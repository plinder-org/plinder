# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from collections import Counter

import pandas as pd
import pytest
from plinder.core.scores.entries import LigandView, entry_views_from_df
from plinder.data.utils.annotations.get_similarity_scores import Scorer


def test_entry_views_accept_annotation_dataframe(cif_2gdo) -> None:
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    entry = Entry.from_cif_file(cif_2gdo)
    view = entry_views_from_df(entry.to_df())[entry.pdb_id]

    assert view.systems
    for system_id, system in entry.systems.items():
        expected = {ligand.instance_chain for ligand in system.ligands}
        assert set(view.systems[system_id].ligands) == expected


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
                "ligand_neighboring_residues": ["1.A_10_9_10", "1.A_20_19_20"],
                "ligand_interactions": ["1.A_10_type:hydrogen_bonds"],
            },
            {
                **common,
                "ligand_id": "1abc__1__1.C",
                "ligand_instance_chain": "1.C",
                "ligand_asym_id": "C",
                "ligand_num_pocket_residues": 1,
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
    # The existing system view remains the union for backward compatibility.
    assert system.pocket_residue_number_to_index == {"1.A": {10: 9, 20: 19, 30: 29}}


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
    )
    alignments = {
        ("1.A", "1.X"): pd.DataFrame(
            [
                {
                    "qrnum": {0: 10, 1: 20},
                    "trnum": {0: 110, 1: 120},
                    "qaln": "AG",
                    "taln": "AG",
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

    full_pocket, full_pli = scorer.get_ligand_pair_pocket_pli_scores(
        alignments, query, target_full
    )
    partial_pocket, partial_pli = scorer.get_ligand_pair_pocket_pli_scores(
        alignments, query, target_partial
    )

    assert full_pocket["pocket_qcov_foldseek"] == pytest.approx(1.0)
    assert partial_pocket["pocket_qcov_foldseek"] == pytest.approx(0.5)
    assert full_pli["pli_qcov_foldseek"] == pytest.approx(1.0)
    assert partial_pli["pli_qcov_foldseek"] == pytest.approx(0.5)
