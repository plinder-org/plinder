# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from plinder.data import docs
from plinder.data.annotations.get_ligand_validation import EntryValidation


def test_entry_descriptions_follow_the_published_table_split() -> None:
    annotation_names = set(docs.get_column_descriptions("annotation")["Name"])
    metadata_names = set(docs.get_column_descriptions("entry_metadata")["Name"])
    moved_names = {
        "entry_release_date",
        "entry_oligomeric_state",
        "entry_determination_method",
        "entry_keywords",
        "entry_pH",
        "entry_resolution",
        "entry_validation_r_minus_rfree",
        "entry_pass_validation_criteria",
        *(
            f"entry_validation_{field_name}"
            for field_name in EntryValidation.model_fields
        ),
    }

    assert moved_names.isdisjoint(annotation_names)
    assert moved_names <= metadata_names
    assert "entry_pdb_id" in annotation_names
    assert "entry_pdb_id" in metadata_names
