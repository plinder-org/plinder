# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from shutil import copytree

from plinder.data import docs


def test_ligand_cluster_column_descriptions():
    import pandas as pd

    columns = [
        "shape__50__strong__component",
        "shape__50__ligand__strong__component",
        "color__70__community",
        "color__70__ligand__community",
    ]

    rows = docs.get_cluster_column_descriptions(pd.DataFrame(columns=columns))

    assert [row[0] for row in rows] == columns
    descriptions = {name: description for name, _, description in rows}
    assert "ligand-level strong component" in descriptions[columns[1]]
    assert "ligand-level community" in descriptions[columns[3]]


def test_make_column_descriptions(read_plinder_mount, tmp_path, monkeypatch):
    from plinder.core.scores import query_index

    generated_tsv_dir = tmp_path / "column_descriptions"
    copytree(docs.TSV_DIR, generated_tsv_dir)
    monkeypatch.setattr(docs, "TSV_DIR", generated_tsv_dir)

    df = query_index(columns=["*"], splits=["*"]).drop(columns=["split"])
    legacy_posebusters = [
        column for column in df.columns if column.startswith("ligand_posebusters_")
    ]
    removed_enrichment_columns = [
        "ligand_is_kinase_inhibitor",
        "system_has_kinase_inhibitor",
        "system_pocket_ECOD",
        "system_pocket_ECOD_t_name",
        "system_pocket_PANTHER",
        "system_pocket_kinase_name",
        "ligand_num_neighboring_ppi_atoms_within_4A_of_gap",
        "ligand_num_neighboring_ppi_atoms_within_8A_of_gap",
        "ligand_num_missing_ppi_interface_residues",
        "ligand_num_pli_atoms_within_4A_of_gap",
        "ligand_num_pli_atoms_within_8A_of_gap",
        "ligand_num_missing_pli_interface_residues",
    ]
    df = df.drop(columns=legacy_posebusters + removed_enrichment_columns)

    schema = docs.get_all_column_descriptions(plindex=df)
    columns = schema["Name"].to_list()
    assert not len(df.columns.difference(columns))
    assert {row[0] for row in docs.DERIVED_LIGAND_COLUMNS}.issubset(columns)
