# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Tests for the CCD-anchored ligand relationship databases.

Everything here runs against the real bundled CCD; no synthetic component ids.
"""

import gzip
import json
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from biotite.structure.io.pdbx import CIFFile, get_structure
from plinder.core.structure.smallmols_similarity import smiles2nonstereo
from plinder.data.annotations.ccd_ligand_dbs import (
    LOG,
    CcdParityTable,
    ResidueGraph,
    _ccd_chem_comp_table,
    _ccd_heavy_atoms,
    _ccd_named_mol,
    _ccd_smiles,
    build_ccd_ecfp_db,
    build_ccd_parity_scores,
    build_ccd_tanimoto_scores,
    build_ligand_ccd_match,
    ccd_code_residue_graph,
    ccd_code_sequence,
    ccd_component_table,
    ccd_dbs_dir,
    ccd_fingerprint_path,
    ccd_node_similarity,
    ccd_parity_path,
    ccd_universe_signature,
    component_sequence_identity,
    composite_parity,
    composite_similarity,
    ligand_features,
    ligand_parity,
    ligand_similarity,
    make_ccd_ligand_dbs,
    make_ligand_ccd_match,
    match_composite_by_sequence,
    query_ccd_mmp_pairs,
    query_ccd_parity,
    query_ccd_tanimoto,
    residue_graph_features,
    residue_graph_from_atoms,
)

needs_mmpdb = pytest.mark.skipif(
    shutil.which("mmpdb") is None, reason="mmpdb executable not installed"
)
LIMIT = 1500


@pytest.fixture(scope="module")
def components() -> pd.DataFrame:
    return ccd_component_table(limit=LIMIT)


def test_chem_comp_table_pins_the_biotite_bundle_contract():
    """Fails loudly if a bundle change removes a column this module needs."""
    table = _ccd_chem_comp_table()
    assert len(table) > 40_000
    # Exactly the columns this module consumes - nothing fetched and unused.
    assert set(table.columns) == {"id", "type", "pdbx_release_status"}
    assert {"REL", "OBS"} <= set(table["pdbx_release_status"])


def test_component_table_filters_and_numbers_nodes(components):
    released = set(
        _ccd_chem_comp_table()
        .loc[lambda df: df["pdbx_release_status"] == "REL", "id"]
        .astype(str)
    )
    assert set(components["ccd_id"]) <= released
    assert not {"UNL", "UNX", "DUM"} & set(components["ccd_id"])
    # is_excluded_mol enforces the >=5 heavy-atom floor; no invented cap.
    assert components["num_heavy_atoms"].min() >= 5
    assert all(value == value.lower() for value in components["ccd_type"])
    # 0-based and contiguous: ligand_scores indexes the fingerprint list by node
    # id and rejects any other numbering.
    assert components["ligand_smiles_id"].tolist() == list(range(len(components)))


def test_universe_signature_is_deterministic_and_content_sensitive(components):
    baseline = ccd_universe_signature(components)
    assert baseline == ccd_universe_signature(components.copy())
    changed = components.copy()
    changed.loc[0, "ligand_rdkit_canonical_smiles"] = "CCO"
    assert ccd_universe_signature(changed) != baseline


def test_non_druglike_components_are_dropped(components):
    """Named real components, so this checks the wiring, not is_excluded_mol itself."""
    broad = ccd_component_table(limit=LIMIT, exclude_non_druglike=False)
    dropped = set(broad["ccd_id"]) - set(components["ccd_id"])
    # 06C = CI (1 carbon), 03S = CS(=O)(=O)O (1 carbon)
    assert {"06C", "03S"} <= dropped
    assert "001" in set(components["ccd_id"])


def test_polymer_linking_exclusion_is_opt_in(components):
    ligand_only = ccd_component_table(limit=LIMIT, exclude_polymer_linking=True)
    assert len(ligand_only) < len(components)
    assert not ligand_only["ccd_type"].str.contains("linking").any()
    assert components["ccd_type"].str.contains("linking").any()


def _match_all(ligands: pd.DataFrame, components: pd.DataFrame) -> pd.DataFrame:
    return build_ligand_ccd_match(ligands, components).set_index("ligand_id")


def test_match_by_code_then_by_chemistry_against_the_real_universe(components):
    known, other = components["ccd_id"].iloc[0], components["ccd_id"].iloc[1]
    smiles = dict(
        zip(components["ccd_id"], components["ligand_rdkit_canonical_smiles"])
    )
    nodes = dict(zip(components["ccd_id"], components["ligand_smiles_id"]))
    ligands = pd.DataFrame(
        {
            "ligand_id": ["by_code", "by_chemistry", "composite", "ion", "novel"],
            "ligand_ccd_code": [known, "ZZZZZ", f"{known}-{other}", "NA", "ZZZZZ"],
            "ligand_smiles": [smiles[known], smiles[known], "CC", "[Na+]", "CCOCC"],
        }
    )
    matches = _match_all(ligands, components)
    assert matches.loc["by_code", "match_kind"] == "exact"
    assert matches.loc["by_code", "ccd_node_ids"] == [nodes[known]]
    # an unknown code with a component's chemistry is that component
    assert matches.loc["by_chemistry", "match_kind"] == "exact"
    assert known in matches.loc["by_chemistry", "matched_ccd_ids"]
    # a composite is not its parts: annotated with them, matched as a molecule
    assert matches.loc["composite", "match_kind"] == "novel"
    assert matches.loc["composite", "component_ccd_ids"] == [known, other]
    assert matches.loc["composite", "ccd_node_ids"] == []
    # sodium is a released CCD component the universe filters out, not novel
    assert "NA" not in set(components["ccd_id"])
    assert matches.loc["ion", "match_kind"] == "excluded"
    assert matches.loc["novel", "match_kind"] == "novel"


def test_one_molecule_matches_however_it_was_deposited():
    """Lactose as one code (LAT) or as its linked sugars (GAL-BGC).

    plinder's identity key is stereo-insensitive, so lactose, cellobiose and
    maltose share it; the stereo-aware resolved SMILES picks the right one.
    """
    from rdkit import Chem

    codes = ["LAT", "CBI", "MAL", "GAL", "BGC"]
    components = pd.DataFrame({"ccd_id": codes})
    components["ligand_rdkit_canonical_smiles"] = [_ccd_smiles(c) for c in codes]
    components["ligand_identity"] = components["ligand_rdkit_canonical_smiles"].map(
        smiles2nonstereo
    )
    components["ligand_smiles_id"] = range(len(codes))
    assert len(set(components["ligand_identity"].iloc[:3])) == 1
    # the composite's whole-molecule SMILES, written from a different atom order
    lactose = Chem.MolFromSmiles(_ccd_smiles("LAT"))
    reordered = Chem.MolToSmiles(
        Chem.RenumberAtoms(lactose, list(reversed(range(lactose.GetNumAtoms())))),
        canonical=False,
    )
    assert reordered != _ccd_smiles("LAT")
    ligands = pd.DataFrame(
        {
            "ligand_id": ["one_code", "two_codes", "two_codes_no_stereo", "galactose"],
            "ligand_ccd_code": ["LAT", "GAL-BGC", "GAL-BGC", "GAL"],
            "ligand_smiles": [
                _ccd_smiles("LAT"),
                reordered,
                reordered,
                _ccd_smiles("GAL"),
            ],
            "ligand_resolved_smiles": [_ccd_smiles("LAT"), reordered, None, None],
        }
    )
    matches = _match_all(ligands, components)
    assert matches.loc["two_codes", "match_kind"] == "exact"
    assert matches.loc["two_codes", "matched_ccd_ids"] == ["LAT"]
    assert matches.loc["two_codes", "component_ccd_ids"] == ["GAL", "BGC"]
    # no resolved stereo: every stereoisomer sharing the key is listed
    assert matches.loc["two_codes_no_stereo", "matched_ccd_ids"] == [
        "CBI",
        "LAT",
        "MAL",
    ]
    assert matches.loc["one_code", "matched_ccd_ids"] == ["LAT"]
    # a code match is never overridden by the stereo-insensitive identity
    assert matches.loc["galactose", "matched_ccd_ids"] == ["GAL"]


def test_make_ligand_ccd_match_joins_the_collated_annotation_table(
    tmp_path, components
):
    with pytest.raises(FileNotFoundError):
        make_ligand_ccd_match(data_dir=tmp_path)

    known = components["ccd_id"].iloc[0]
    (tmp_path / "index").mkdir()
    pd.DataFrame(
        {
            "ligand_id": ["a__1__1.A", "a__1__1.A", "b__1__1.B"],
            "ligand_ccd_code": [known, known, "ZZZZZ"],
            "ligand_smiles": ["C", "C", "CCOCC"],
        }
    ).to_parquet(tmp_path / "index" / "annotation_table.parquet", index=False)
    ccd_dbs_dir(tmp_path).mkdir()
    components.to_parquet(ccd_dbs_dir(tmp_path) / "ccd_components.parquet", index=False)

    output = make_ligand_ccd_match(data_dir=tmp_path)

    assert output == ccd_dbs_dir(tmp_path) / "ligand_ccd_match.parquet"
    matches = pd.read_parquet(output).set_index("ligand_id")
    assert matches["match_kind"].to_dict() == {
        "a__1__1.A": "exact",
        "b__1__1.B": "novel",
    }


@pytest.mark.parametrize(
    "query, target, expected",
    [
        ("SER-GLY-ALA", "SER-GLY-ALA", 1.0),
        ("DA-DT-DG", "DA-DC-DG", pytest.approx(2 / 3)),
        # Order matters - a sequence match, not a set intersection.
        ("NAG-BMA", "BMA-NAG", 0.5),
        ("SER-GLY", "TRP-PHE", 0.0),
        ("SER-GLY", "SER-GLY-ALA-LYS", 0.5),
        ("", "ATP", 0.0),
        ("ATP", "ATP", 1.0),
    ],
    ids=[
        "identical",
        "one-substitution",
        "order-matters",
        "disjoint",
        "prefix-normalised-by-longer",
        "empty-query",
        "single-component",
    ],
)
def test_component_sequence_identity(query, target, expected):
    """Each case pins a distinct property; the comparison is class-agnostic."""
    assert (
        component_sequence_identity(ccd_code_sequence(query), ccd_code_sequence(target))
        == expected
    )


def test_match_composite_by_sequence_ranks_and_thresholds():
    """Stereoisomeric sugars that achiral ECFP4 cannot separate, resolved by code."""
    references = {"same": "NAG-NAG-BMA", "one_off": "NAG-NAG-MAN", "other": "DA-DT-DG"}
    matches = match_composite_by_sequence("NAG-NAG-BMA", references, min_identity=0.6)
    assert [name for name, _ in matches] == ["same", "one_off"]
    assert matches[0][1] == 1.0


@needs_mmpdb
def test_tanimoto_scores_match_native_edge_schema(tmp_path, components):
    """The CCD table drives the *native* scorer unchanged."""
    build_ccd_ecfp_db(components, ccd_fingerprint_path(tmp_path))
    scores_dir = build_ccd_tanimoto_scores(
        components, tmp_path, minimum_similarity=30.0, batch_size=750
    )
    shards = sorted(scores_dir.glob("*.parquet"))
    assert len(shards) > 1
    edges = pd.concat([pd.read_parquet(shard) for shard in shards])
    assert list(edges.columns) == [
        "query_ligand_id",
        "target_ligand_id",
        "tanimoto_similarity_ecfp4_1024",
    ]
    nodes = set(components["ligand_smiles_id"])
    assert set(edges["query_ligand_id"]) <= nodes
    assert set(edges["target_ligand_id"]) <= nodes
    self_edges = edges[edges["query_ligand_id"] == edges["target_ligand_id"]]
    assert (self_edges["tanimoto_similarity_ecfp4_1024"] == 100.0).all()
    assert edges["tanimoto_similarity_ecfp4_1024"].min() >= 30.0


def test_parity_scores_every_tanimoto_edge_once(tmp_path):
    """The ECFP4 edges are the prefilter; each unordered pair is scored once."""
    small = ccd_component_table(limit=300)
    build_ccd_ecfp_db(small, ccd_fingerprint_path(tmp_path))
    scores_dir = build_ccd_tanimoto_scores(small, tmp_path, batch_size=100)
    edges = pd.concat(
        [pd.read_parquet(shard) for shard in scores_dir.glob("*.parquet")]
    )
    unordered = edges[edges["query_ligand_id"] < edges["target_ligand_id"]]
    parity = pd.read_parquet(build_ccd_parity_scores(small, tmp_path, threads=2))
    assert list(parity.columns) == [
        "ligand_smiles_id_1",
        "ligand_smiles_id_2",
        "tanimoto_similarity_ecfp4_1024",
        "parity_similarity",
        "fragment_atoms_1",
        "fragment_atoms_2",
    ]
    assert len(parity) == len(unordered) > 0
    # the stored fragment is in CCD atom names, one name per mapped atom on each side
    codes = small.set_index("ligand_smiles_id")["ccd_id"]
    scored = parity[parity["parity_similarity"] > 0]
    assert len(scored) > 0
    for row in scored.head(20).itertuples(index=False):
        assert len(row.fragment_atoms_1) == len(row.fragment_atoms_2) > 0
        for one, other in zip(row.fragment_atoms_1, row.fragment_atoms_2):
            assert len(one) == len(other) > 0
            assert set(one) <= set(_ccd_heavy_atoms(codes[row.ligand_smiles_id_1])[0])
            assert set(other) <= set(_ccd_heavy_atoms(codes[row.ligand_smiles_id_2])[0])
    assert (parity["ligand_smiles_id_1"] < parity["ligand_smiles_id_2"]).all()
    assert not parity.duplicated(["ligand_smiles_id_1", "ligand_smiles_id_2"]).any()
    assert parity["parity_similarity"].between(0.0, 100.0).all()
    # the graph score orders neighbours differently from the fingerprint
    assert (
        parity["parity_similarity"].corr(parity["tanimoto_similarity_ecfp4_1024"])
        < 0.95
    )


@pytest.fixture(scope="module")
def ccd_dbs(tmp_path_factory) -> Path:
    if shutil.which("mmpdb") is None:
        pytest.skip("mmpdb executable not installed")
    data_dir = tmp_path_factory.mktemp("ccd_dbs")
    make_ccd_ligand_dbs(
        data_dir=data_dir,
        scratch_dir=data_dir / "scratch",
        threads=2,
        limit=400,
        build_tanimoto_scores=False,
    )
    return data_dir


@needs_mmpdb
def test_fragment_dictionary_reproduces_every_mmpdb_pair(ccd_dbs):
    pairs = pd.read_parquet(ccd_dbs / "ccd_dbs" / "ccd_mmp_pairs.parquet")
    assert not pairs.empty
    queried = pairs.drop_duplicates("ligand_smiles_id_1").head(12)
    found = query_ccd_mmp_pairs(
        dict(
            zip(queried["ligand_smiles_id_1"].astype(str), queried["ligand_smiles_1"])
        ),
        ccd_dbs,
    )
    expected = pairs[pairs["ligand_smiles_id_1"].isin(queried["ligand_smiles_id_1"])]
    for pair in expected.itertuples(index=False):
        hits = found[
            (found["query_id"] == str(pair.ligand_smiles_id_1))
            & (found["ligand_smiles_id"] == pair.ligand_smiles_id_2)
            & (found["shared_core_smiles"] == pair.shared_core_smiles)
        ]
        # mmpdb writes the pair in canonical order; ours reads query >> component
        reversed_smirks = ">>".join(reversed(pair.transformation.split(">>")))
        assert {pair.transformation, reversed_smirks} & set(
            hits["transformation"]
        ), pair
    # a component never pairs with itself
    assert not (found["query_id"].astype(int) == found["ligand_smiles_id"]).any()


@needs_mmpdb
def test_novel_molecule_pairs_against_the_universe(ccd_dbs):
    components = pd.read_parquet(ccd_dbs / "ccd_dbs" / "ccd_components.parquet")
    # diphenyl ether is not one of the first 400 components
    assert "c1ccc(Oc2ccccc2)cc1" not in set(components["ligand_rdkit_canonical_smiles"])
    novel = query_ccd_mmp_pairs({"dpe": "c1ccccc1Oc1ccccc1"}, ccd_dbs)
    assert not novel.empty
    assert list(novel.columns) == [
        "query_id",
        "query_smiles",
        "ligand_smiles_id",
        "transformation",
        "shared_core_smiles",
        "num_cuts",
        "shared_core_num_heavy_atoms",
        "ccd_id",
        "ligand_smiles",
    ]
    assert set(novel["ligand_smiles_id"]) <= set(components["ligand_smiles_id"])
    assert novel["ccd_id"].notna().all()
    assert (novel["transformation"].str.count(">>") == 1).all()
    # an unparsable query is skipped, not fatal
    assert query_ccd_mmp_pairs({"bad": "C1CC"}, ccd_dbs).empty


@needs_mmpdb
def test_novel_molecule_tanimoto_neighbours_use_the_native_measure(ccd_dbs):
    components = pd.read_parquet(ccd_dbs / "ccd_dbs" / "ccd_components.parquet")
    own = components.iloc[7]
    hits = query_ccd_tanimoto(
        {"self": own["ligand_rdkit_canonical_smiles"], "dpe": "c1ccccc1Oc1ccccc1"},
        ccd_dbs,
        minimum_similarity=30.0,
    )
    assert list(hits.columns) == [
        "query_id",
        "query_smiles",
        "ligand_smiles_id",
        "tanimoto_similarity_ecfp4_1024",
        "ccd_id",
    ]
    best = hits[hits["query_id"] == "self"].iloc[0]
    assert best["ccd_id"] == own["ccd_id"]
    assert best["tanimoto_similarity_ecfp4_1024"] == 100.0
    assert (hits["tanimoto_similarity_ecfp4_1024"] >= 30.0).all()
    assert set(hits["ligand_smiles_id"]) <= set(components["ligand_smiles_id"])
    assert hits["ccd_id"].notna().all()


@needs_mmpdb
def test_novel_molecule_parity_extends_the_precompiled_table(ccd_dbs):
    """Ligands that only exist at dataset time are scored the same way, later."""
    components = pd.read_parquet(ccd_dbs / "ccd_dbs" / "ccd_components.parquet")
    own = components.iloc[7]
    hits = query_ccd_parity(
        {"self": own["ligand_rdkit_canonical_smiles"], "dpe": "c1ccccc1Oc1ccccc1"},
        ccd_dbs,
        threads=2,
    )
    assert list(hits.columns) == [
        "query_id",
        "query_smiles",
        "ligand_smiles_id",
        "tanimoto_similarity_ecfp4_1024",
        "ccd_id",
        "parity_similarity",
    ]
    best = hits[hits["query_id"] == "self"].iloc[0]
    assert best["ccd_id"] == own["ccd_id"]
    assert best["parity_similarity"] == 100.0
    assert (hits["tanimoto_similarity_ecfp4_1024"] >= 30.0).all()
    assert hits["parity_similarity"].between(0.0, 100.0).all()
    dpe = hits[hits["query_id"] == "dpe"]["parity_similarity"]
    assert dpe.is_monotonic_decreasing and len(dpe) > 0


@needs_mmpdb
def test_make_ccd_ligand_dbs_builds_artifacts_and_caches(tmp_path):
    output_dir = make_ccd_ligand_dbs(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch",
        threads=2,
        limit=250,
        build_tanimoto_scores=False,
    )
    paths = [
        output_dir / "ccd_components.parquet",
        output_dir / "fingerprints" / "ligands_per_smiles.parquet",
        output_dir / "ccd_mmp_pairs.parquet",
        output_dir / "ccd_mmp_fragments.parquet",
        output_dir / "ccd_dbs.manifest.json",
    ]
    for path in paths:
        assert path.is_file(), path

    components = pd.read_parquet(paths[0])
    ecfp = pd.read_parquet(paths[1])
    mmp = pd.read_parquet(paths[2])
    assert set(ecfp["ligand_smiles_id"]) == set(components["ligand_smiles_id"])
    assert "fingerprint" in ecfp.columns
    nodes = set(components["ligand_smiles_id"])
    assert set(mmp["ligand_smiles_id_1"]) <= nodes
    assert set(mmp["ligand_smiles_id_2"]) <= nodes

    before = {path.name: path.stat().st_mtime_ns for path in paths}
    make_ccd_ligand_dbs(
        data_dir=tmp_path,
        scratch_dir=tmp_path / "scratch",
        threads=2,
        limit=250,
        build_tanimoto_scores=False,
    )
    assert {path.name: path.stat().st_mtime_ns for path in paths} == before


@needs_mmpdb
def test_ccd_score_cache_tracks_similarity_cutoff(tmp_path):
    def build(cutoff, *, scores=True):
        return make_ccd_ligand_dbs(
            data_dir=tmp_path,
            scratch_dir=tmp_path / "scratch",
            threads=2,
            limit=12,
            minimum_similarity=cutoff,
            build_tanimoto_scores=scores,
        )

    def read_scores(output_dir):
        edges = pd.concat(
            [
                pd.read_parquet(path)
                for path in sorted((output_dir / "ligand_scores").glob("*.parquet"))
            ],
            ignore_index=True,
        )
        parity_pairs = pd.read_parquet(
            ccd_parity_path(tmp_path),
            columns=[
                "ligand_smiles_id_1",
                "ligand_smiles_id_2",
                "tanimoto_similarity_ecfp4_1024",
            ],
        )
        return edges, parity_pairs

    output_dir = build(0.0)
    all_edges, all_parity = read_scores(output_dir)
    similarity = "tanimoto_similarity_ecfp4_1024"
    assert all_parity[similarity].lt(100.0).any()

    # A score-free run must not certify the older score files at a new cutoff.
    build(100.0, scores=False)
    build(100.0)
    high_edges, high_parity = read_scores(output_dir)
    for actual, original in [(high_edges, all_edges), (high_parity, all_parity)]:
        pd.testing.assert_frame_equal(
            actual,
            original.loc[original[similarity].ge(100.0)].reset_index(drop=True),
        )

    build(0.0)
    low_edges, low_parity = read_scores(output_dir)
    pd.testing.assert_frame_equal(low_edges, all_edges)
    pd.testing.assert_frame_equal(low_parity, all_parity)
    manifest = json.loads((output_dir / "ccd_dbs.manifest.json").read_text())
    assert manifest["minimum_similarity"] == 0.0
    assert manifest["build_tanimoto_scores"] is True

    paths = list((output_dir / "ligand_scores").glob("*.parquet")) + [
        ccd_parity_path(tmp_path)
    ]
    modified = {path: path.stat().st_mtime_ns for path in paths}
    build(0.0)
    assert {path: path.stat().st_mtime_ns for path in paths} == modified


@needs_mmpdb
def test_tanimoto_rebuild_does_not_leave_stale_shards(tmp_path):
    """A smaller rebuild must not leave edges referencing dead node ids."""
    big = ccd_component_table(limit=300)
    build_ccd_ecfp_db(big, ccd_fingerprint_path(tmp_path))
    build_ccd_tanimoto_scores(big, tmp_path, batch_size=100)

    small = ccd_component_table(limit=120)
    build_ccd_ecfp_db(small, ccd_fingerprint_path(tmp_path))
    scores_dir = build_ccd_tanimoto_scores(small, tmp_path, batch_size=100)

    edges = pd.concat(
        [pd.read_parquet(shard) for shard in scores_dir.glob("*.parquet")]
    )
    assert edges["query_ligand_id"].max() < len(small)


@pytest.fixture(scope="module")
def atoms_2dty(cif_2dty):
    """2dty as deposited: chains E and F each carry a branched NAG-FUC-NAG."""
    with gzip.open(cif_2dty, "rt") as handle:
        return get_structure(
            CIFFile.read(handle), model=1, include_bonds=True, altloc="first"
        )


def glycan(atoms, *chains):
    return atoms[
        np.isin(atoms.chain_id, list(chains)) & np.isin(atoms.res_name, ["NAG", "FUC"])
    ]


def test_residue_graph_roots_at_the_reducing_end(atoms_2dty):
    graph = residue_graph_from_atoms(glycan(atoms_2dty, "E"))
    assert graph.labels == ("NAG", "FUC", "NAG")
    # The core NAG carries both the fucose and the second NAG: branched, not a
    # chain, and the only residue donating no anomeric carbon.
    assert graph.adjacency == (frozenset({1, 2}), frozenset({0}), frozenset({0}))
    assert graph.roots == (0,)
    assert graph.is_tree and graph.is_connected


def test_connectivity_separates_what_the_ccd_code_cannot(atoms_2dty):
    """The deposited glycan is branched; its own code string reads as linear."""
    deposited = residue_graph_from_atoms(glycan(atoms_2dty, "E"))
    code = "-".join(deposited.labels)

    assert (
        component_sequence_identity(ccd_code_sequence(code), ccd_code_sequence(code))
        == 1.0
    ), "sequence identity cannot see the branch"
    assert composite_similarity(
        deposited, ccd_code_residue_graph(code)
    ) == pytest.approx(4 / 9)
    assert composite_similarity(deposited, deposited) == 1.0


def test_composite_similarity_is_independent_of_residue_order(atoms_2dty):
    """Relabelling the traversal must not move the score; a code string would."""
    graph = residue_graph_from_atoms(glycan(atoms_2dty, "E"))
    order = [2, 0, 1]
    position = {old: new for new, old in enumerate(order)}
    shuffled = ResidueGraph(
        tuple(graph.labels[old] for old in order),
        tuple(frozenset(position[n] for n in graph.adjacency[old]) for old in order),
        tuple(position[root] for root in graph.roots),
    )
    assert shuffled.labels != graph.labels
    assert composite_similarity(graph, shuffled) == 1.0


def test_residue_graph_records_present_atoms_and_link_atoms(atoms_2dty):
    """Linked residues have lost their leaving oxygen; each linkage names its atoms."""
    graph = residue_graph_from_atoms(glycan(atoms_2dty, "E"))
    assert [len(names) for names in graph.atoms] == [14, 10, 14]
    assert all("O1" not in names for names in graph.atoms)
    # fucose alpha1-3 and the second NAG beta1-4, both donating their anomeric C1
    assert graph.links == ((0, 1, "O3", "C1"), (0, 2, "O4", "C1"))


@pytest.fixture(scope="module")
def sugar_table() -> CcdParityTable:
    return CcdParityTable.compute(["NAG", "FUC", "GLC", "BMA"], threads=2)


def test_parity_table_is_symmetric_and_identity_aware(sugar_table):
    forward, backward = (
        sugar_table.mapping("GLC", "NAG"),
        sugar_table.mapping("NAG", "GLC"),
    )
    assert forward and backward == {q: p for p, q in forward.items()}
    assert (
        sugar_table.similarity("GLC", "NAG") == sugar_table.similarity("NAG", "GLC") > 0
    )
    assert sugar_table.similarity("NAG", "NAG") == 1.0
    assert sugar_table.mapping("NAG", "NAG") == {
        n: n for n in _ccd_heavy_atoms("NAG")[0]
    }
    assert sugar_table.mappings("GLC", "NAG")[0] == forward
    # alternatives are stored for residue pairs; the scored fragment leads
    alternatives = CcdParityTable.compute(["FUC", "VAL"]).mappings("FUC", "VAL")
    assert len(alternatives) > 1
    assert len({tuple(sorted(m.items())) for m in alternatives}) == len(alternatives)
    assert all(len(m) <= len(alternatives[0]) for m in alternatives)
    assert (
        sugar_table.similarity("NAG", "XYZ") == 0.0
        and sugar_table.mapping("NAG", "XYZ") == {}
    )


def test_composite_parity_reassembles_the_whole_molecule_score(atoms_2dty, sugar_table):
    """Fragments joined through matched linkages give the exact atom-level score."""
    from plinder.core.structure.smallmols_similarity import rascal_parity_score
    from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

    chain_e, chain_f = glycan(atoms_2dty, "E"), glycan(atoms_2dty, "F")
    core = chain_e[chain_e.res_name != "FUC"]  # NAG-NAG, the fucose branch removed
    graphs = [residue_graph_from_atoms(a) for a in (chain_e, chain_f, core)]
    mols = [atoms_to_rdkit_mol(a) for a in (chain_e, chain_f, core)]
    exact = lambda i, j: rascal_parity_score(mols[i], mols[j], stereo=False)  # noqa: E731
    assert composite_parity(graphs[0], graphs[1], sugar_table) == exact(0, 1) == 1.0
    assert composite_parity(graphs[0], graphs[2], sugar_table) == pytest.approx(
        exact(0, 2)
    )
    assert composite_parity(graphs[2], graphs[0], sugar_table) == pytest.approx(
        exact(0, 2)
    )
    assert 0.7 < exact(0, 2) < 0.75


def test_composite_parity_scores_mono_against_composite(atoms_2dty, sugar_table):
    """A one-residue graph is the mono case; the score tracks the unpruned engine."""
    from plinder.core.structure.smallmols_similarity import rascal_parity_score
    from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

    chain_e = glycan(atoms_2dty, "E")
    graph, whole = residue_graph_from_atoms(chain_e), atoms_to_rdkit_mol(chain_e)
    free = lambda code: ResidueGraph((code,), (frozenset(),), (0,))  # noqa: E731
    for code in ("NAG", "FUC", "GLC"):
        assembled = composite_parity(free(code), graph, sugar_table)
        exact = rascal_parity_score(
            _ccd_named_mol(code), whole, target=0.0, stereo=False
        )
        assert composite_parity(graph, free(code), sugar_table) == assembled
        # the engine may also carry the mono's leaving oxygen onto the linkage
        assert 0 < assembled <= exact + 1e-9
        assert assembled == pytest.approx(exact, abs=0.05)
    assert composite_parity(free("NAG"), graph, sugar_table) > composite_parity(
        free("GLC"), graph, sugar_table
    )


def _linked_chain(codes, *, donor="C1", acceptors=None, leaving="O1"):
    """CCD residues linked donor(i) -> acceptor(i+1), each donor's leaving atom dropped."""
    import biotite.structure as struc
    from plinder.data.annotations.cif_utils import _get_ccd_atomarray

    units = []
    for index, code in enumerate(codes):
        unit = _get_ccd_atomarray(code)
        unit = unit[unit.element != "H"]
        if index < len(codes) - 1:
            unit = unit[unit.atom_name != leaving]
        unit.res_id[:] = index + 1
        unit.chain_id[:] = "A"
        units.append(unit)
    atoms = struc.concatenate(units)
    offsets = np.cumsum([0] + [len(unit) for unit in units[:-1]])
    for index, acceptor in enumerate(acceptors or ["O4"] * (len(codes) - 1)):
        first = offsets[index] + int(np.where(units[index].atom_name == donor)[0][0])
        second = offsets[index + 1] + int(
            np.where(units[index + 1].atom_name == acceptor)[0][0]
        )
        atoms.bonds.add_bond(first, second, struc.BondType.SINGLE)
    return atoms


def test_a_mono_spanning_several_residues_goes_through_the_whole_molecule_engine():
    """Chitotriose (CTO, one code) against NAG-NAG-NAG built from CCD residues.

    The table holds one fragment of CTO per NAG unit, but residues align
    one-to-one, so the assembly credits a single unit. The routed score runs
    the engine on the whole molecules and recovers the identity.
    """
    from plinder.core.structure.smallmols_similarity import rascal_parity_score
    from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

    chain = _linked_chain(["NAG", "NAG", "NAG"])
    graph = residue_graph_from_atoms(chain)
    assert graph.links == ((0, 1, "C1", "O4"), (1, 2, "C1", "O4"))
    assert [len(names) for names in graph.atoms] == [14, 14, 15]
    whole = atoms_to_rdkit_mol(chain, assign_stereo=False)
    cto = _ccd_named_mol("CTO")
    assert (whole.GetNumAtoms(), whole.GetNumBonds()) == (
        cto.GetNumAtoms(),
        cto.GetNumBonds(),
    )
    table = CcdParityTable.compute(["NAG", "CTO"])
    assert len(table.mappings("CTO", "NAG")) == 3  # one fragment per unit
    free = ResidueGraph(("CTO",), (frozenset(),), (0,))
    assembled = composite_parity(free, graph, table)
    assert 0.15 < assembled < 0.25  # one unit of three, on each side's union
    exact = rascal_parity_score(cto, whole, stereo=False)
    assert exact == 1.0
    assert ligand_parity(free, graph, table, cto, whole) == exact
    assert ligand_parity(graph, free, table, whole, cto) == exact
    # two composites still assemble
    assert ligand_parity(graph, graph, table, whole, whole) == composite_parity(
        graph, graph, table
    )


def test_a_linkage_isomer_scores_high_but_not_identical():
    """NAG-NAG-NAG with the last linkage beta1-6 against the beta1-4 chain and CTO.

    The mono table maps NAG onto NAG as the identity, which cannot carry a
    linkage that lands on another oxygen; rewiring the acceptor keeps the chain
    connected at the cost of the two displaced hydroxyl bonds. The engine does
    one better by swapping C4 with C6 along with their oxygens, losing a single
    ring bond, so the assembly sits just below it.
    """
    from plinder.core.structure.smallmols_similarity import (
        parity_similarity,
        rascal_parity_score,
    )
    from plinder.data.annotations.cif_utils import atoms_to_rdkit_mol

    regular = _linked_chain(["NAG", "NAG", "NAG"])
    isomer = _linked_chain(["NAG", "NAG", "NAG"], acceptors=["O4", "O6"])
    graphs = [residue_graph_from_atoms(a) for a in (regular, isomer)]
    mols = [atoms_to_rdkit_mol(a, assign_stereo=False) for a in (regular, isomer)]
    assert graphs[1].links == ((0, 1, "C1", "O4"), (1, 2, "C1", "O6"))
    table = CcdParityTable.compute(["NAG", "CTO"])
    # all 43 atoms stay one fragment; the engine matches 44 of 45 bonds, the
    # rewired assembly 43 (two hydroxyl bonds displaced, one linkage gained)
    exact = rascal_parity_score(mols[0], mols[1], stereo=False)
    assert exact == pytest.approx(parity_similarity(43, 44, 86, 90))
    assembled = composite_parity(graphs[0], graphs[1], table)
    assert assembled == pytest.approx(parity_similarity(43, 43, 86, 90))
    assert composite_parity(graphs[1], graphs[0], table) == assembled
    assert 0.95 < assembled < exact < 1.0
    free = ResidueGraph(("CTO",), (frozenset(),), (0,))
    cto = _ccd_named_mol("CTO")
    assert ligand_parity(free, graphs[0], table, cto, mols[0]) == 1.0
    assert ligand_parity(free, graphs[1], table, cto, mols[1]) == pytest.approx(exact)


def test_residue_graph_requires_a_bond_graph(atoms_2dty):
    unbonded = glycan(atoms_2dty, "E").copy()
    unbonded.bonds = None
    with pytest.raises(ValueError, match="requires bonded atoms"):
        residue_graph_from_atoms(unbonded)


def test_disconnected_residues_are_flagged_not_scored_silently(atoms_2dty, monkeypatch):
    """Two chains' glycans share no bond; that must not pass as one composite."""
    warnings: list[str] = []
    monkeypatch.setattr(LOG, "warning", lambda msg, *a: warnings.append(msg % a))
    graph = residue_graph_from_atoms(glycan(atoms_2dty, "E", "F"))
    assert not graph.is_connected
    assert not graph.is_tree
    assert any("disconnected" in message for message in warnings)


def build_peptide(sequence: str, *, cyclic: bool):
    """Residue graph and whole-molecule SMILES of an RDKit-built peptide."""
    import io

    import biotite.structure.io.pdb as pdb
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom

    mol = Chem.RWMol(Chem.MolFromSequence(sequence))
    named = lambda atom, name: atom.GetPDBResidueInfo().GetName().strip() == name  # noqa: E731
    if cyclic:
        n_term = next(
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetPDBResidueInfo().GetResidueNumber() == 1 and named(a, "N")
        )
        c_term = next(
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetPDBResidueInfo().GetResidueNumber() == len(sequence)
            and named(a, "C")
        )
        for idx in sorted(
            (a.GetIdx() for a in mol.GetAtoms() if named(a, "OXT")), reverse=True
        ):
            mol.RemoveAtom(idx)
        mol.AddBond(n_term, c_term, Chem.BondType.SINGLE)
    peptide = mol.GetMol()
    Chem.SanitizeMol(peptide)
    peptide = Chem.AddHs(peptide)
    rdDistGeom.EmbedMolecule(peptide, randomSeed=0xF00D)
    block = Chem.MolToPDBBlock(Chem.RemoveHs(peptide))
    stripped = Chem.RemoveHs(peptide)
    atoms = pdb.PDBFile.read(io.StringIO(block)).get_structure(
        model=1, include_bonds=True
    )
    return residue_graph_from_atoms(atoms), Chem.MolToSmiles(stripped)


def peptide_graph(sequence: str, *, cyclic: bool):
    """Just the residue graph, for tests that do not need the molecule."""
    return build_peptide(sequence, cyclic=cyclic)[0]


@pytest.fixture(scope="module")
def cyclic_peptide():
    """Head-to-tail cyclo(Ala-Gly-Ser-Phe-Leu)."""
    return peptide_graph("AGSFL", cyclic=True)


def rotated(graph, shift):
    """Same macrocycle, residue numbering rotated by *shift*."""
    size = len(graph.labels)
    position = {old: (old + shift) % size for old in range(size)}
    order = sorted(range(size), key=lambda old: position[old])
    return ResidueGraph(
        tuple(graph.labels[old] for old in order),
        tuple(frozenset(position[x] for x in graph.adjacency[old]) for old in order),
        tuple(sorted(position[root] for root in graph.roots)),
        tuple(frozenset(position[x] for x in graph.donors[old]) for old in order),
    )


def test_macrocycle_has_no_terminus_to_root_at(cyclic_peptide):
    assert cyclic_peptide.labels == ("ALA", "GLY", "SER", "PHE", "LEU")
    assert [len(bonded) for bonded in cyclic_peptide.adjacency] == [2] * 5
    assert cyclic_peptide.is_connected and not cyclic_peptide.is_tree
    # Every residue donates its carbonyl and accepts at its nitrogen, so none is
    # a terminus and all are candidate roots.
    assert cyclic_peptide.roots == (0, 1, 2, 3, 4)


@pytest.mark.parametrize("shift", [0, 1, 2, 3, 4])
def test_macrocycle_score_is_invariant_to_rotation(cyclic_peptide, shift):
    """Rotating where the ring is numbered must not move the score."""
    assert composite_similarity(cyclic_peptide, rotated(cyclic_peptide, shift)) == 1.0


def test_macrocycle_is_distinguished_from_the_real_linear_peptide(cyclic_peptide):
    """A ring and the open chain of the same residues are different molecules.

    The spanning tree drops the ring-closing edge, so this separation rests
    entirely on the two termini having one linkage where the ring has two.
    """
    linear = peptide_graph("AGSFL", cyclic=False)
    assert linear.labels == cyclic_peptide.labels
    assert [len(bonded) for bonded in linear.adjacency] == [1, 2, 2, 2, 1]
    # Three internal residues match fully (degree 2 both); the two termini
    # carry degree 1 against the ring's 2, weighted (1+1)/(1+2).
    assert composite_similarity(cyclic_peptide, linear) == pytest.approx(13 / 15)


@pytest.mark.parametrize(
    "sequence, note",
    [("AGSLF", "last two residues swapped"), ("ASGFL", "middle two swapped")],
)
def test_scrambled_macrocycle_matches_partially(cyclic_peptide, sequence, note):
    """A different ring order over the same residues is neither same nor unrelated."""
    assert sequence not in {
        "AGSFL"[shift:] + "AGSFL"[:shift] for shift in range(5)
    }, "the scramble must not be a rotation"
    scrambled = peptide_graph(sequence, cyclic=True)
    assert sorted(scrambled.labels) == sorted(cyclic_peptide.labels)
    assert 0.0 < composite_similarity(cyclic_peptide, scrambled) < 1.0


def test_ccd_node_similarity_grades_conservative_substitutions():
    """Residue chemistry, not string equality: Tyr is Phe plus a hydroxyl."""
    assert ccd_node_similarity("TYR", "TYR") == 1.0
    assert (
        ccd_node_similarity("TYR", "PHE")
        > ccd_node_similarity("TYR", "TRP")
        > ccd_node_similarity("TYR", "ALA")
    )
    # An underivable component must not invent similarity.
    assert ccd_node_similarity("TYR", "ZZZZZ") == 0.0


def test_graded_kernel_ranks_substitutions_within_the_same_topology():
    """A-X-B vs A-Y-B: identical topology, so the node kernel decides."""
    reference = peptide_graph("AYG", cyclic=False)
    scores = {
        substitution: composite_similarity(
            reference,
            peptide_graph(f"A{substitution}G", cyclic=False),
            node_similarity=ccd_node_similarity,
        )
        for substitution in ("Y", "F", "W", "A")
    }
    assert scores["Y"] == 1.0
    assert scores["F"] > scores["W"] > scores["A"]
    # Exact matching cannot see any of this: all three are one mismatch of three.
    assert all(
        composite_similarity(reference, peptide_graph(f"A{s}G", cyclic=False))
        == pytest.approx(1 / 3 * 2)
        for s in ("F", "W", "A")
    )


def test_single_residue_ligand_still_scores_against_a_composite(cyclic_peptide):
    """A free residue must not be unmatchable against composites containing it.

    The linkage weighting compares degrees, and a free residue has degree 0. An
    unsmoothed ratio makes that pair impossible rather than merely poor, zeroing
    every mono-vs-multi comparison - worse than the whole-molecule fingerprint
    this is meant to improve on.
    """
    free = lambda code: ResidueGraph((code,), (frozenset(),), (0,))  # noqa: E731

    present = composite_similarity(free("PHE"), cyclic_peptide)
    absent = composite_similarity(free("TRP"), cyclic_peptide)
    assert present > 0.0
    # PHE is one of the five ring residues; TRP is not.
    assert present > absent
    # Scaled by how much of the composite one residue can account for.
    assert present == pytest.approx(1 / 3 / 5)
    assert composite_similarity(free("PHE"), free("PHE")) == 1.0


def test_residue_features_are_invariant_to_numbering_but_see_topology(cyclic_peptide):
    """Relabelling must not change features; a different ring order must."""
    assert residue_graph_features(cyclic_peptide) == residue_graph_features(
        rotated(cyclic_peptide, 3)
    )
    assert residue_graph_features(cyclic_peptide) != residue_graph_features(
        peptide_graph("AGSLF", cyclic=True)
    )


@pytest.mark.parametrize(
    "single, parts",
    [("LAT", ("GAL", "BGC")), ("SUC", ("GLC", "FRU"))],
    ids=["lactose", "sucrose"],
)
def test_same_molecule_survives_a_different_decomposition(single, parts):
    """One CCD code or several linked ones can denote the same molecule.

    Heavy atoms confirm the identity: a glycosidic bond loses one oxygen, so
    the single code carries the components' atoms minus one.
    """
    from rdkit import Chem

    heavy = lambda code: Chem.MolFromSmiles(_ccd_smiles(code)).GetNumHeavyAtoms()  # noqa: E731
    assert heavy(single) == sum(heavy(part) for part in parts) - 1

    smiles = _ccd_smiles(single)
    as_one = ResidueGraph((single,), (frozenset(),), (0,))
    as_several = ccd_code_residue_graph("-".join(parts))

    # Residue labels share nothing across the regrouping.
    assert composite_similarity(as_one, as_several) == 0.0
    # Sharing one feature space does not privilege either decomposition.
    combined = ligand_similarity(
        ligand_features(smiles, as_one), ligand_features(smiles, as_several)
    )
    assert combined > 0.7


def test_ligand_similarity_separates_what_fingerprints_alone_cannot(cyclic_peptide):
    """The scramble ECFP4 cannot see, still separated once residues count."""
    features = lambda seq: ligand_features(  # noqa: E731
        *reversed(build_peptide(seq, cyclic=True))
    )
    reference, identical, scrambled = (
        features("AGSFL"),
        features("AGSFL"),
        features("AGSLF"),
    )

    assert ligand_similarity(reference, identical) == 1.0
    assert ligand_similarity(reference, scrambled) < 1.0
    assert ligand_similarity(reference, Counter()) == 0.0


def test_residue_counts_separate_a_repeat_from_its_unit():
    """Presence-only features cannot see multiplicity; counted ones can.

    A pentapeptide and a head-to-tail decapeptide of the same residues have the
    identical set of circular environments and the identical set of residue
    labels, so every presence-only comparison scores them 1.000.
    """
    short_graph, short_smiles = build_peptide("AGSFL", cyclic=True)
    long_graph, long_smiles = build_peptide("AGSFLAGSFL", cyclic=True)
    assert len(long_graph.labels) == 2 * len(short_graph.labels)

    short, long = (
        ligand_features(short_smiles, short_graph),
        ligand_features(long_smiles, long_graph),
    )
    # Same feature *keys* either way - only the counts differ.
    assert set(short) == set(long)
    assert ligand_similarity(short, long) < 1.0

    presence_only = Counter(dict.fromkeys(short, 1)), Counter(dict.fromkeys(long, 1))
    assert ligand_similarity(*presence_only) == 1.0
