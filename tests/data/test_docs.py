# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from plinder.data import docs


def test_make_column_descriptions(read_plinder_mount):
    from plinder.core.scores import query_index

    df = query_index(columns=["*"], splits=["*"]).drop(columns=["split"])

    schema = docs.get_all_column_descriptions(plindex=df)
    columns = schema["Name"].to_list()
    assert not len(df.columns.difference(columns))
