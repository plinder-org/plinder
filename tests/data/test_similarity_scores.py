# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
import pytest
from plinder.core.scores.entries import (
    ChainView,
    EntryView,
    LigandView,
    SystemView,
    entry_views_from_df,
    load_entry_views,
)
from plinder.core.utils.schemas import PROTEIN_SIMILARITY_SCHEMA
from plinder.data.utils.annotations import get_similarity_scores as scoring_module
from plinder.data.utils.annotations.get_similarity_scores import (
    Scorer,
    annotate_cofactor_similarity,
    annotate_ligand_3d_score_ability,
    annotate_ligand_similarity,
    build_ligand_similarity_annotations,
    compute_ligand_fingerprints,
    get_feature_map_score,
    ligand_scores,
    load_ligands_from_index,
    prepare_ligand_for_3d_scoring,
    write_ecfp4_fingerprint_table,
)
from rdkit import Chem

SDF_FILE = (
    Path(__file__).resolve().parents[1]
    / "test_data"
    / "mini_system_files_new"
    / "1fbz__1__1.A__1.C"
    / "ligand_files"
    / "1.C.sdf"
)
HEM_SDF_FILE = (
    Path(__file__).resolve().parents[1]
    / "test_data"
    / "plinder"
    / "mount"
    / "systems"
    / "19hc__1__1.B__1.T"
    / "ligand_files"
    / "1.T.sdf"
)


def test_ligand_scoring_inputs_include_only_proper_holo_ligands() -> None:
    annotation = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc", "1abc", "1abc"],
            "system_id": ["proper", "artifact", "ion-system"],
            "system_type": ["holo", "holo", "ion"],
            "ligand_is_proper": [True, False, True],
            "ligand_rdkit_canonical_smiles": ["CCO", "O", "[Na+]"],
            "ligand_unique_ccd_code": ["LIG", "HOH", "NA"],
            "ligand_id": ["1abc__1__1.L", "1abc__1__1.W", "1abc__1__1.N"],
            "ligand_asym_id": ["L", "W", "N"],
        }
    )

    ligands = load_ligands_from_index(annotation=annotation)

    assert ligands["ligand_id"].tolist() == ["1abc__1__1.L"]


def test_entry_views_accept_annotation_dataframe(cif_2gdo, tmp_path) -> None:
    from plinder.data.utils.annotations.aggregate_annotations import Entry

    entry = Entry.from_cif_file(cif_2gdo)
    annotation = entry.to_df()
    assert not any(column.startswith("entry_chains_") for column in annotation)
    entry_chains = entry.chains_to_df()
    assert len(entry_chains) <= len(entry.chains)
    assert entry_chains["chain_type"].str.lower().str.contains("polypeptide").all()

    # Exercise the Arrow representation used by the ingest pipeline,
    # including the normalized nested UniProt accession lists.
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    annotation_path = index_dir / "annotation_table.parquet"
    chain_path = index_dir / "entry_chains.parquet"
    annotation.to_parquet(annotation_path, index=False)
    entry_chains.to_parquet(chain_path, index=False)
    view = entry_views_from_df(
        pd.read_parquet(annotation_path),
        entry_chains=pd.read_parquet(chain_path),
    )[entry.pdb_id]

    assert view.systems
    for system_id, system in entry.systems.items():
        expected = {ligand.instance_chain for ligand in system.ligands}
        assert set(view.systems[system_id].ligands) == expected
    for chain_type in ["holo", "apo", "pred"]:
        for aln_type in ["foldseek", "mmseqs"]:
            assert sorted(view.chains_for_alignment(chain_type, aln_type)) == sorted(
                entry.chains_for_alignment(chain_type, aln_type)
            )

    loaded = load_entry_views(pdb_ids=[entry.pdb_id], data_dir=tmp_path)[entry.pdb_id]
    assert loaded.chains == view.chains

    chain_path.unlink()
    with pytest.raises(FileNotFoundError, match="normalized entry chain index"):
        load_entry_views(pdb_ids=[entry.pdb_id], data_dir=tmp_path)


def test_entry_view_builds_holo_apo_and_pred_alignment_ids() -> None:
    system = SystemView(
        id="1abc__1__1.A__1.L",
        pdb_id="1abc",
        system_type="holo",
        protein_chains_asym_id=["1.A"],
        proper_num_pocket_residues=1,
        proper_num_interactions=1,
        proper_num_unique_interactions=1,
    )
    view = EntryView(
        pdb_id="1abc",
        chains={
            "A": ChainView(
                asym_id="A",
                auth_id="R",
                entity_id="1",
                chain_type="polypeptide(L)",
                length=200,
                holo=True,
                uniprot_ids=("P12345",),
            ),
            "B": ChainView(
                asym_id="B",
                auth_id="B",
                entity_id="2",
                chain_type="polypeptide(L)",
                length=150,
                holo=False,
            ),
            # A non-holo copy of a holo entity must not enter the apo DB.
            "C": ChainView(
                asym_id="C",
                auth_id="C",
                entity_id="1",
                chain_type="polypeptide(L)",
                length=200,
                holo=False,
            ),
            # Search DBs are protein-only.
            "N": ChainView(
                asym_id="N",
                auth_id="N",
                entity_id="3",
                chain_type="polyribonucleotide",
                length=80,
                holo=False,
            ),
        },
        systems={system.id: system},
        author_to_asym={"R": "A", "B": "B", "C": "C"},
    )

    assert view.chains_for_alignment("holo", "foldseek") == [
        "pdb_00001abc_xyz-enrich_R"
    ]
    assert view.chains_for_alignment("holo", "mmseqs") == ["1abc_R"]
    assert view.chains_for_alignment("apo", "foldseek") == ["pdb_00001abc_xyz-enrich_B"]
    assert view.chains_for_alignment("apo", "mmseqs") == ["1abc_B"]
    assert view.chains_for_alignment("pred", "foldseek") == ["AF-P12345-F1-model_v4_A"]
    assert view.chains_for_alignment("pred", "mmseqs") == ["P12345"]


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
                "ligand_is_3d_score_able": False,
                "ligand_neighboring_residues": ["1.A_10_9_10", "1.A_20_19_20"],
                "ligand_interactions": ["1.A_10_type:hydrogen_bonds"],
            },
            {
                **common,
                "ligand_id": "1abc__1__1.C",
                "ligand_instance_chain": "1.C",
                "ligand_asym_id": "C",
                "ligand_num_pocket_residues": 1,
                "ligand_is_3d_score_able": True,
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
    assert not system.ligands["1.B"].is_3d_score_able
    assert system.ligands["1.C"].is_3d_score_able
    # The existing system view remains the union for backward compatibility.
    assert system.pocket_residue_number_to_index == {"1.A": {10: 9, 20: 19, 30: 29}}


def test_entry_views_exclude_na_receptors_from_scoring(tmp_path) -> None:
    rows = pd.DataFrame(
        [
            {
                "entry_pdb_id": "1abc",
                "system_id": "1abc__1__1.A_1.N__1.L",
                "system_type": "holo",
                "system_receptor_type": "protein+dna",
                "system_protein_chains_asym_id": ["1.A", "1.N"],
                "system_protein_chains_auth_id": ["A", "N"],
                "system_protein_chains_length": [100, 20],
                "system_proper_num_pocket_residues": 2,
                "system_proper_num_interactions": 2,
                "system_proper_num_unique_interactions": 2,
                "ligand_id": "1abc__1__1.L",
                "ligand_instance_chain": "1.L",
                "ligand_asym_id": "L",
                "ligand_is_proper": True,
                "ligand_protein_chains_asym_id": ["1.A", "1.N"],
                "ligand_num_pocket_residues": 2,
                "ligand_num_interactions": 2,
                "ligand_num_unique_interactions": 2,
                "ligand_neighboring_residues": ["1.A_10_9_10", "1.N_2_1_2"],
                "ligand_interactions": [
                    "1.A_10_type:hydrogen_bonds",
                    "1.N_2_type:hydrogen_bonds",
                ],
            }
        ]
    )
    entry_chains = pd.DataFrame(
        {
            "entry_pdb_id": ["1abc"],
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
    )

    with pytest.raises(ValueError, match="entry_chains is required for 1abc"):
        entry_views_from_df(rows)

    entry = entry_views_from_df(rows, entry_chains=entry_chains)["1abc"]
    system = next(iter(entry.systems.values()))
    ligand = next(iter(system.ligands.values()))
    scorer = Scorer(
        entries={"1abc": entry},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )

    assert system.receptor_type == "protein+dna"
    assert system.protein_chains_asym_id == ["1.A"]
    assert system.pocket_residue_number_to_index == {"1.A": {10: 9}}
    assert system.proper_num_pocket_residues == 1
    assert ligand.protein_chains_asym_id == ["1.A"]
    assert ligand.pocket_residue_number_to_index == {"1.A": {10: 9}}
    assert ligand.num_interactions == 1
    assert scorer.get_protein_chain_length("1abc", ["1.A", "1.N"]) == 100
    assert entry.chains_for_alignment("holo", "mmseqs") == ["1abc_A"]


def test_na_only_entry_has_no_similarity_chains() -> None:
    rows = pd.DataFrame(
        [
            {
                "entry_pdb_id": "1dna",
                "system_id": "1dna__1__1.N__1.L",
                "system_type": "holo",
                "system_receptor_type": "dna",
                "system_protein_chains_asym_id": ["1.N"],
                "system_protein_chains_auth_id": ["N"],
                "system_protein_chains_length": [20],
                "system_proper_num_pocket_residues": 1,
                "system_proper_num_interactions": 1,
                "system_proper_num_unique_interactions": 1,
                "ligand_id": "1dna__1__1.L",
                "ligand_instance_chain": "1.L",
                "ligand_asym_id": "L",
                "ligand_is_proper": True,
                "ligand_protein_chains_asym_id": ["1.N"],
                "ligand_num_pocket_residues": 1,
                "ligand_num_interactions": 1,
                "ligand_num_unique_interactions": 1,
                "ligand_neighboring_residues": ["1.N_2_1_2"],
                "ligand_interactions": ["1.N_2_type:hydrogen_bonds"],
            }
        ]
    )
    entry_chains = pd.DataFrame(
        {
            "entry_pdb_id": ["1dna"],
            "chain_asym_id": ["N"],
            "chain_auth_id": ["N"],
            "chain_entity_id": ["1"],
            "chain_type": ["polydeoxyribonucleotide"],
            "chain_receptor_type": ["dna"],
            "chain_length": [20],
            "chain_num_unresolved_residues": [0],
            "chain_is_holo": [True],
            "chain_uniprot_ids": [[]],
        }
    )

    entry = entry_views_from_df(rows, entry_chains=entry_chains)["1dna"]
    system = next(iter(entry.systems.values()))

    assert system.receptor_type == "dna"
    assert system.protein_chains_asym_id == []
    assert system.proper_num_pocket_residues == 0
    assert entry.chains["N"].receptor_type == "dna"
    assert entry.chains_for_alignment("holo", "foldseek") == []
    assert entry.chains_for_alignment("holo", "mmseqs") == []


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


def test_ligand_pair_shape_scores_gate_sdf_access_and_cache(
    tmp_path, monkeypatch
) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    query = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    target = _ligand("2def__1__1.Y", "1.Y", {"1.X": {110: 109}}, {})
    resolved: list[str] = []

    def resolve_sdf(_data_dir: Path, ligand: LigandView) -> Path:
        resolved.append(ligand.id)
        return SDF_FILE

    monkeypatch.setattr(scorer, "resolve_ligand_sdf", resolve_sdf)

    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.0) == {}
    assert (
        scorer.get_ligand_pair_shape_scores(tmp_path, query, target, float("nan")) == {}
    )
    assert resolved == []

    query.is_3d_score_able = False
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == {}
    assert resolved == []
    query.is_3d_score_able = True

    scores = scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5)
    assert set(scores) == {
        "shape",
        "color",
        "sucos_shape",
        "sucos_shape_pocket_qcov",
    }
    assert scores["sucos_shape_pocket_qcov"] == pytest.approx(
        scores["sucos_shape"] * 0.5
    )
    assert all(0 <= value <= 1 for value in scores.values())
    assert resolved == [query.id, target.id]

    # Re-scoring uses pristine clones of the two cached base molecules.
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == scores
    assert resolved == [query.id, target.id]


def test_ligand_pair_shape_scores_drop_nonfinite_results(tmp_path, monkeypatch) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
        ligand_sdf_resolver=lambda _ligand: SDF_FILE,
    )
    query = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    target = _ligand("2def__1__1.Y", "1.Y", {"1.X": {110: 109}}, {})

    monkeypatch.setattr(
        scoring_module, "align_molecules", lambda *_args: (1.0, float("nan"))
    )
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == {}

    monkeypatch.setattr(scoring_module, "align_molecules", lambda *_args: (0.8, 0.6))
    monkeypatch.setattr(scoring_module, "get_sucos_score", lambda *_args: float("nan"))
    assert scorer.get_ligand_pair_shape_scores(tmp_path, query, target, 0.5) == {
        "shape": 0.8,
        "color": 0.6,
    }


def test_ligand_sdf_resolver_uses_only_canonical_entry_path(tmp_path) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    ligand = _ligand("1abc__1__1.B", "1.B", {"1.A": {10: 9}}, {})
    canonical = tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files" / "B.sdf"
    legacy = (
        tmp_path / "raw_entries" / "ab" / ligand.system_id / "ligand_files" / "1.B.sdf"
    )
    legacy.parent.mkdir(parents=True)
    legacy.write_text("legacy")

    assert scorer.resolve_ligand_sdf(tmp_path, ligand) is None
    canonical.parent.mkdir(parents=True)
    canonical.write_text("canonical")
    assert scorer.resolve_ligand_sdf(tmp_path, ligand) == canonical


def test_get_score_df_loads_entries_from_ingest_data_dir(tmp_path, monkeypatch) -> None:
    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    alignment_path = scorer.db_dir / "holo_foldseek" / "aln" / "1abc.parquet"
    alignment_path.parent.mkdir(parents=True)
    pd.DataFrame({"target_pdb_id": ["2def"]}).to_parquet(alignment_path, index=False)
    calls: list[tuple[set[str], Path]] = []

    def fake_load_entry_views(*, pdb_ids, data_dir):
        calls.append((set(pdb_ids), data_dir))
        return {}

    monkeypatch.setattr(scoring_module, "load_entry_views", fake_load_entry_views)
    monkeypatch.setattr(
        scorer,
        "map_alignment_df",
        lambda *_args: pd.DataFrame({"mapped": [True]}),
    )
    monkeypatch.setattr(
        scorer,
        "aggregate_scores",
        lambda *_args, **_kwargs: pd.DataFrame(),
    )

    scorer.get_score_df(tmp_path, "1abc", "holo")

    assert calls == [({"1abc", "2def"}, tmp_path)]


def test_feature_map_score_handles_molecules_without_features() -> None:
    helium = Chem.MolFromSmiles("[He]")
    assert helium is not None
    conformer = Chem.Conformer(helium.GetNumAtoms())
    conformer.SetAtomPosition(0, (0.0, 0.0, 0.0))
    helium.AddConformer(conformer)

    assert get_feature_map_score(helium, helium) == 0.0


def test_shape_scoring_uses_custom_radius_and_preserves_hem_iron() -> None:
    molecule = Chem.MolFromMolFile(str(HEM_SDF_FILE))
    assert molecule is not None
    assert any(atom.GetAtomicNum() == 26 for atom in molecule.GetAtoms())

    prepared = prepare_ligand_for_3d_scoring(molecule)

    assert prepared is not None
    assert any(atom.GetAtomicNum() == 26 for atom in prepared.GetAtoms())
    shape, color = scoring_module.align_molecules(
        Chem.Mol(prepared), Chem.Mol(prepared)
    )
    assert shape == pytest.approx(1.0)
    assert color == pytest.approx(1.0)
    assert scoring_module.get_sucos_score(prepared, prepared) == pytest.approx(1.0)


def test_ligand_3d_score_ability_uses_canonical_sdf_and_caches_by_path(
    tmp_path, monkeypatch
) -> None:
    valid_sdf = tmp_path / "raw_entries" / "ab" / "1abc" / "ligand_files" / "B.sdf"
    valid_sdf.parent.mkdir(parents=True)
    valid_sdf.write_bytes(SDF_FILE.read_bytes())
    calls: list[Path] = []
    original = scoring_module.is_ligand_3d_score_able

    def observed(sdf_file: Path) -> bool:
        calls.append(sdf_file)
        return original(sdf_file)

    monkeypatch.setattr(scoring_module, "is_ligand_3d_score_able", observed)
    ligands = pd.DataFrame(
        {
            "pdb_id": ["1abc", "1abc", "1abc"],
            "ligand_asym_id": ["B", "B", "C"],
            "ligand_id": ["first", "second", "missing"],
        }
    )

    annotated = annotate_ligand_3d_score_ability(ligands, data_dir=tmp_path)

    assert annotated["ligand_is_3d_score_able"].tolist() == [True, True, False]
    assert calls == [valid_sdf, valid_sdf.with_name("C.sdf")]


def test_ligand_3d_score_ability_reuses_success_for_same_molecular_graph(
    tmp_path, monkeypatch
) -> None:
    first_sdf = tmp_path / "first.sdf"
    second_sdf = tmp_path / "second.sdf"
    first_sdf.write_bytes(SDF_FILE.read_bytes())
    second_sdf.write_bytes(SDF_FILE.read_bytes())
    scoring_module._LIGAND_3D_SCORE_ABILITY_CACHE.clear()
    align_calls = 0
    sucos_calls = 0
    original_align = scoring_module.align_molecules
    original_sucos = scoring_module.get_sucos_score

    def observed_align(*args, **kwargs):
        nonlocal align_calls
        align_calls += 1
        return original_align(*args, **kwargs)

    def observed_sucos(*args, **kwargs):
        nonlocal sucos_calls
        sucos_calls += 1
        return original_sucos(*args, **kwargs)

    monkeypatch.setattr(scoring_module, "align_molecules", observed_align)
    monkeypatch.setattr(scoring_module, "get_sucos_score", observed_sucos)

    try:
        assert scoring_module.is_ligand_3d_score_able(first_sdf)
        assert scoring_module.is_ligand_3d_score_able(second_sdf)
    finally:
        scoring_module._LIGAND_3D_SCORE_ABILITY_CACHE.clear()

    assert align_calls == 1
    assert sucos_calls == 1


def test_cofactor_similarity_uses_ccd_reference_fingerprints() -> None:
    from plinder.core.structure.smallmols_similarity import mol2morgan_fp
    from rdkit import DataStructs

    smiles = ["CCO", "c1ccccc1"]
    unique_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1],
            "ligand_rdkit_canonical_smiles": smiles,
            "fingerprint": [
                DataStructs.BitVectToBinaryText(mol2morgan_fp(value, nbits=1024))
                for value in smiles
            ],
        }
    )

    annotated = annotate_cofactor_similarity(
        unique_ligands,
        cofactor_smiles={"COF": "CCO"},
        threshold=90,
    ).set_index("ligand_smiles_id")

    assert annotated.loc[0, "ligand_max_cofactor_similarity"] == 100
    assert annotated.loc[0, "ligand_most_similar_cofactor"] == "COF"
    assert bool(annotated.loc[0, "ligand_is_cofactor_like"])
    assert not bool(annotated.loc[1, "ligand_is_cofactor_like"])


def test_ligand_scores_use_bulk_tanimoto_for_unique_smiles(
    tmp_path, monkeypatch
) -> None:
    from plinder.core.structure.smallmols_similarity import mol2morgan_fp
    from rdkit import DataStructs

    fingerprint_dir = tmp_path / "fingerprints"
    fingerprint_dir.mkdir()
    smiles = ["CCO", "CCN", "c1ccccc1"]
    fingerprint_table = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1, 2],
            "ligand_rdkit_canonical_smiles": smiles,
            "fingerprint": [
                DataStructs.BitVectToBinaryText(mol2morgan_fp(value, nbits=1024))
                for value in smiles
            ],
        }
    )
    write_ecfp4_fingerprint_table(
        fingerprint_table,
        fingerprint_dir / "ligands_per_smiles.parquet",
    )
    calls = 0
    original_bulk = scoring_module.DataStructs.BulkTanimotoSimilarity

    def observed_bulk(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_bulk(*args, **kwargs)

    monkeypatch.setattr(
        scoring_module.DataStructs,
        "BulkTanimotoSimilarity",
        observed_bulk,
    )
    output_path = tmp_path / "scores.parquet"

    ligand_scores(
        ligand_ids=[0, 1, 2],
        data_dir=tmp_path,
        output_path=output_path,
        minimum_similarity=30,
    )
    scores = pd.read_parquet(output_path)

    assert calls == 3
    assert set(
        scores.loc[
            scores["query_ligand_id"] == scores["target_ligand_id"],
            "query_ligand_id",
        ]
    ) == {0, 1, 2}
    assert not (
        (scores["query_ligand_id"] == 0) & (scores["target_ligand_id"] == 2)
    ).any()
    assert scores.columns.tolist() == [
        "query_ligand_id",
        "target_ligand_id",
        "tanimoto_similarity_ecfp4_1024",
    ]
    assert (
        scoring_module.pq.read_schema(output_path).metadata
        == scoring_module.ECFP4_PARQUET_METADATA
    )


def test_tanimoto_90_cluster_counts_distinct_pdb_ids() -> None:
    unique_ligands = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 1, 2],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCN", "c1ccccc1"],
            "fingerprint": [b"", b"", b""],
            "ligand_is_cofactor_like": [True, False, False],
        }
    )
    occurrences = pd.DataFrame(
        {
            "ligand_smiles_id": [0, 0, 1, 2],
            "pdb_id": ["1aaa", "1aaa", "2bbb", "3ccc"],
        }
    )
    edges = pd.DataFrame(
        {
            "query_ligand_id": [0, 0, 1, 2],
            "target_ligand_id": [0, 1, 1, 2],
            "tanimoto_similarity_ecfp4_1024": [100.0, 91.0, 100.0, 100.0],
        }
    )

    annotations = build_ligand_similarity_annotations(
        unique_ligands=unique_ligands,
        ligand_occurrences=occurrences,
        edges=edges,
        cluster_threshold=90,
    ).set_index("ligand_smiles_id")

    assert (
        annotations.loc[0, "ligand_tanimoto_ecfp4_1024_90_cluster"]
        == annotations.loc[1, "ligand_tanimoto_ecfp4_1024_90_cluster"]
    )
    assert (
        annotations.loc[2, "ligand_tanimoto_ecfp4_1024_90_cluster"]
        != annotations.loc[0, "ligand_tanimoto_ecfp4_1024_90_cluster"]
    )
    assert annotations.loc[0, "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids"] == 2
    assert annotations.loc[1, "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids"] == 2
    assert annotations.loc[2, "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids"] == 1


def test_ligand_similarity_pipeline_does_not_write_per_system_mapping(
    tmp_path, monkeypatch
) -> None:
    ligand_dir = tmp_path / "ligands"
    ligand_dir.mkdir()
    pd.DataFrame(
        {
            "pdb_id": ["1aaa", "2bbb", "3ccc"],
            "ligand_rdkit_canonical_smiles": ["CCO", "CCO", "c1ccccc1"],
        }
    ).to_parquet(ligand_dir / "part.parquet", index=False)
    component_dir = tmp_path / "dbs" / "components"
    component_dir.mkdir(parents=True)
    pd.DataFrame({"binder_id": ["COF"], "canonical_smiles": ["CCO"]}).to_parquet(
        component_dir / "components.parquet", index=False
    )

    from plinder.data.utils.annotations import ligand_utils

    monkeypatch.setattr(ligand_utils, "parse_cofactors", lambda _data_dir: {"COF"})
    compute_ligand_fingerprints(data_dir=tmp_path)

    unique_ligands = pd.read_parquet(
        tmp_path / "fingerprints" / "ligands_per_smiles.parquet"
    )
    assert len(unique_ligands) == 2
    assert not (tmp_path / "fingerprints" / "ligands_per_system.parquet").exists()

    score_dir = tmp_path / "ligand_scores"
    score_dir.mkdir()
    ligand_scores(
        ligand_ids=unique_ligands["ligand_smiles_id"].tolist(),
        data_dir=tmp_path,
        output_path=score_dir / "part.parquet",
    )
    annotation_path = annotate_ligand_similarity(data_dir=tmp_path)
    annotations = pd.read_parquet(annotation_path).set_index(
        "ligand_rdkit_canonical_smiles"
    )

    assert bool(annotations.loc["CCO", "ligand_is_cofactor_like"])
    assert (
        annotations.loc["CCO", "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids"] == 2
    )
    assert (
        annotations.loc["c1ccccc1", "ligand_tanimoto_ecfp4_1024_90_cluster_num_pdb_ids"]
        == 1
    )


def _system(pdb_id: str, ligands: list[LigandView]) -> SystemView:
    protein_chains = sorted(
        {chain for ligand in ligands for chain in ligand.protein_chains_asym_id}
    )
    return SystemView(
        id=f"{pdb_id}_system",
        pdb_id=pdb_id,
        system_type="holo",
        protein_chains_asym_id=protein_chains,
        proper_num_pocket_residues=0,
        proper_num_interactions=0,
        proper_num_unique_interactions=0,
        ligands={ligand.instance_chain: ligand for ligand in ligands},
    )


def _entry(pdb_id: str, system: SystemView) -> EntryView:
    asym_ids = [chain.split(".", 1)[1] for chain in system.protein_chains_asym_id]
    return EntryView(
        pdb_id=pdb_id,
        chains={
            asym_id: ChainView(asym_id=asym_id, auth_id=asym_id, length=100)
            for asym_id in asym_ids
        },
        systems={system.id: system},
        author_to_asym={asym_id: asym_id for asym_id in asym_ids},
    )


def test_holo_scores_are_emitted_per_ligand_pair(tmp_path, monkeypatch) -> None:
    query_ligands = [
        _ligand("1abc__1__1.C", "1.C", {"1.A": {10: 9}}, {}),
        _ligand("1abc__1__1.D", "1.D", {"1.B": {20: 19}}, {}),
    ]
    target_ligands = [
        _ligand("2def__1__1.Z", "1.Z", {"1.X": {110: 109}}, {}),
        _ligand("2def__1__1.W", "1.W", {"1.Y": {120: 119}}, {}),
    ]
    query_system = _system("1abc", query_ligands)
    target_system = _system("2def", target_ligands)
    scorer = Scorer(
        entries={
            "1abc": _entry("1abc", query_system),
            "2def": _entry("2def", target_system),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    alignments = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [("2def", "A", "X"), ("2def", "B", "Y")],
            names=["target_entry", "query_chain_mapped", "target_chain_mapped"],
        )
    )
    protein_calls: list[tuple[tuple[str, ...], tuple[str, ...], int]] = []

    def protein_scores(
        _alignments: pd.DataFrame,
        _query_system: SystemView,
        target_protein_chains: list[str],
        query_length: int,
        query_protein_chains: list[str] | None = None,
    ) -> tuple[dict, dict, dict, str]:
        assert query_protein_chains is not None
        pair = (query_protein_chains[0], target_protein_chains[0])
        protein_calls.append(
            (tuple(query_protein_chains), tuple(target_protein_chains), query_length)
        )
        return (
            {"protein_qcov_foldseek_weighted_sum": [pair]},
            {"protein_qcov_foldseek_weighted_sum": 0.75},
            {pair: pd.DataFrame()},
            "protein_qcov_foldseek",
        )

    def pocket_scores(
        _alns: dict,
        query_ligand: LigandView,
        target_ligand: LigandView,
    ) -> tuple[dict[str, float], dict[str, float]]:
        qcov = float(
            query_ligand.id == query_ligands[0].id
            and target_ligand.id == target_ligands[0].id
        )
        return {"pocket_qcov_foldseek": qcov}, {}

    shape_calls: list[tuple[str, str, float]] = []

    def shape_scores(
        _data_dir: Path,
        query_ligand: LigandView,
        target_ligand: LigandView,
        pocket_qcov: float,
    ) -> dict[str, float]:
        shape_calls.append((query_ligand.id, target_ligand.id, pocket_qcov))
        return {
            "shape": 0.8,
            "color": 0.7,
            "sucos_shape": 0.6,
            "sucos_shape_pocket_qcov": 0.6 * pocket_qcov,
        }

    monkeypatch.setattr(scorer, "get_protein_scores", protein_scores)
    monkeypatch.setattr(scorer, "get_ligand_pair_pocket_pli_scores", pocket_scores)
    monkeypatch.setattr(scorer, "get_ligand_pair_shape_scores", shape_scores)

    scores = list(scorer.get_scores_holo(query_system, alignments, data_dir=tmp_path))

    assert {
        (score["query_ligand_id"], score["target_ligand_id"]) for score in scores
    } == {(query.id, target.id) for query in query_ligands for target in target_ligands}
    assert len(protein_calls) == 4
    assert {
        (query_chains, target_chains)
        for query_chains, target_chains, _ in protein_calls
    } == {
        (("1.A",), ("1.X",)),
        (("1.A",), ("1.Y",)),
        (("1.B",), ("1.X",)),
        (("1.B",), ("1.Y",)),
    }
    assert all(query_length == 100 for _, _, query_length in protein_calls)
    assert shape_calls == [(query_ligands[0].id, target_ligands[0].id, 1.0)]
    shape_row = next(score for score in scores if "shape" in score)
    assert shape_row["sucos_shape_pocket_qcov"] == pytest.approx(0.6)


def test_holo_weighted_sum_retains_unmatched_query_receptor_length(
    tmp_path, monkeypatch
) -> None:
    query_ligand = _ligand(
        "1abc__1__1.C",
        "1.C",
        {"1.A": {10: 9}, "1.B": {20: 19}},
        {},
    )
    target_ligand = _ligand("2def__1__1.Z", "1.Z", {"1.X": {110: 109}}, {})
    query_system = _system("1abc", [query_ligand])
    target_system = _system("2def", [target_ligand])
    scorer = Scorer(
        entries={
            "1abc": _entry("1abc", query_system),
            "2def": _entry("2def", target_system),
        },
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    # Only query chain A has a target hit. Chain B must nevertheless remain
    # in the 200-residue directed query denominator.
    alignments = pd.DataFrame(
        [
            {
                "qcov": 0.5,
                "fident": 0.5,
                "seqsim": 0.5,
                "fident_qcov": 0.25,
                "seqsim_qcov": 0.25,
                "lddt": 0.5,
                "lddt_qcov": 0.5,
            }
        ],
        index=pd.MultiIndex.from_tuples(
            [("2def", "A", "X", "foldseek")],
            names=[
                "target_entry",
                "query_chain_mapped",
                "target_chain_mapped",
                "source",
            ],
        ),
    )
    monkeypatch.setattr(
        scorer,
        "get_ligand_pair_pocket_pli_scores",
        lambda *_args: ({}, {}),
    )

    scores = list(scorer.get_scores_holo(query_system, alignments))

    assert len(scores) == 1
    assert scores[0]["protein_lddt_qcov_weighted_max"] == pytest.approx(0.5)
    assert scores[0]["protein_lddt_qcov_weighted_sum"] == pytest.approx(0.25)


def test_apo_pred_scores_are_emitted_per_query_ligand(tmp_path, monkeypatch) -> None:
    query_ligands = [
        _ligand(
            "1abc__1__1.C",
            "1.C",
            {"1.A": {10: 9}, "1.B": {20: 19}},
            {},
        ),
        _ligand("1abc__1__1.D", "1.D", {"1.B": {20: 19}}, {}),
    ]
    query_system = _system("1abc", query_ligands)
    scorer = Scorer(
        entries={"1abc": _entry("1abc", query_system)},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    alignments = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [("model_a", "A", "X"), ("model_b", "B", "Y")],
            names=["target_entry", "query_chain_mapped", "target_chain_mapped"],
        )
    )
    protein_calls: list[tuple[tuple[str, ...], int]] = []

    def protein_scores(
        _alignments: pd.DataFrame,
        _query_system: SystemView,
        target_protein_chains: list[str],
        query_length: int,
        query_protein_chains: list[str] | None = None,
    ) -> tuple[dict, dict, dict, str]:
        assert query_protein_chains is not None
        protein_calls.append((tuple(query_protein_chains), query_length))
        pair = (query_protein_chains[0], target_protein_chains[0])
        return (
            {"protein_qcov_foldseek_weighted_sum": [pair]},
            {"protein_qcov_foldseek_weighted_sum": 0.75},
            {pair: pd.DataFrame()},
            "protein_qcov_foldseek",
        )

    monkeypatch.setattr(scorer, "get_protein_scores", protein_scores)
    monkeypatch.setattr(
        scorer,
        "get_ligand_pocket_scores",
        lambda _alns, _ligand: {"pocket_fident_foldseek": 0.5},
    )

    scores = list(scorer.get_scores_apo_pred(query_system, alignments))

    assert {score["query_ligand_id"] for score in scores} == {
        ligand.id for ligand in query_ligands
    }
    assert {score["target_ligand_id"] for score in scores} == {None}
    assert {score["target_system"] for score in scores} == {
        "model_a_X",
        "model_b_Y",
    }
    assert (("1.A", "1.B"), 200) in protein_calls


def test_aggregate_scores_keeps_ligand_ids_and_shape_metrics(
    tmp_path, monkeypatch
) -> None:
    ligand = _ligand("1abc__1__1.C", "1.C", {"1.A": {10: 9}}, {})
    system = _system("1abc", [ligand])
    scorer = Scorer(
        entries={"1abc": _entry("1abc", system)},
        source_to_full_db_file={},
        db_dir=tmp_path / "db",
        scores_dir=tmp_path / "scores",
    )
    monkeypatch.setattr(
        scorer,
        "load_alignments",
        lambda **_kwargs: pd.DataFrame(
            {"present": [True]}, index=pd.Index(["1abc"], name="query_entry")
        ),
    )
    monkeypatch.setattr(
        scorer,
        "get_scores",
        lambda *_args, **_kwargs: iter(
            [
                {
                    "query_system": system.id,
                    "query_ligand_id": ligand.id,
                    "target_system": "2def_system",
                    "target_ligand_id": "2def__1__1.Z",
                    "protein_mapping": "1.A:1.X",
                    "protein_mapper": "foldseek",
                    "protein_qcov_weighted_sum": 0.75,
                    "protein_qcov_weighted_sum_source": "foldseek",
                    "protein_qcov_weighted_sum_mapping": "1.A:1.X",
                    "shape": 0.8,
                    "color": 0.7,
                    "sucos_shape": 0.6,
                    "sucos_shape_pocket_qcov": 0.3,
                }
            ]
        ),
    )

    scores = scorer.aggregate_scores("1abc", data_dir=tmp_path)

    assert scores is not None
    assert set(scores["metric"]) >= {
        "shape",
        "color",
        "sucos_shape",
        "sucos_shape_pocket_qcov",
    }
    assert set(scores["query_ligand_id"]) == {ligand.id}
    assert set(scores["target_ligand_id"]) == {"2def__1__1.Z"}
    assert scores.loc[scores["metric"] == "shape", "similarity"].item() == 80
    assert {"query_ligand_id", "target_ligand_id"}.issubset(
        PROTEIN_SIMILARITY_SCHEMA.names
    )
    output_file = tmp_path / "scores.parquet"
    scores.to_parquet(output_file, index=False, schema=PROTEIN_SIMILARITY_SCHEMA)
    written = pd.read_parquet(output_file)
    assert set(written["query_ligand_id"]) == {ligand.id}
