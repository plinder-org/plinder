# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from plinder.data.annotations.chain_metadata import (
    annotate_chain_metadata,
    write_annotated_chain_metadata,
)


def test_unresolved_fraction_handles_missing_and_invalid_counts() -> None:
    table = pa.table(
        {
            "chain_length": [100, 100, 100, 0, -1, None, 100, 100, 100],
            "chain_num_unresolved_residues": [0, 25, 100, 0, 0, 0, None, -1, 101],
        }
    )
    annotated = annotate_chain_metadata(table)
    assert annotated["chain_unresolved_fraction"].type == pa.float64()
    assert annotated["chain_unresolved_fraction"].to_pylist() == [
        0.0,
        0.25,
        1.0,
        None,
        None,
        None,
        None,
        None,
        None,
    ]
    assert annotated.select(table.column_names).equals(table)
    assert annotate_chain_metadata(annotated).equals(annotated)
    stale = annotated.set_column(
        2, "chain_unresolved_fraction", pa.array([0.5] * len(table))
    )
    assert annotate_chain_metadata(stale).equals(annotated)


@pytest.mark.parametrize("rows", [0, 5])
def test_write_annotated_copy_is_bounded_and_preserves_source(
    tmp_path: Path, rows: int
) -> None:
    table = pa.table(
        {
            "chain_length": pa.array([100] * rows, type=pa.int64()),
            "chain_num_unresolved_residues": pa.array([25] * rows, type=pa.int64()),
            "chain_asym_id": pa.array(["A"] * rows, type=pa.string()),
        }
    ).replace_schema_metadata({b"provenance": b"test"})
    source = tmp_path / "source.parquet"
    output = tmp_path / "new" / "entry_chains.parquet"
    pq.write_table(table, source)
    original = source.read_bytes()
    assert write_annotated_chain_metadata(source, output, batch_size=2) == rows
    assert source.read_bytes() == original
    result = pq.read_table(output)
    assert result.equals(annotate_chain_metadata(table), check_metadata=True)
    assert pq.ParquetFile(output).metadata.num_row_groups == (rows + 1) // 2
    with pytest.raises(FileExistsError):
        write_annotated_chain_metadata(source, output)
    with pytest.raises(FileExistsError):
        write_annotated_chain_metadata(source, source)
    with pytest.raises(ValueError, match="positive"):
        write_annotated_chain_metadata(source, tmp_path / "bad.parquet", batch_size=0)
    assert not list(output.parent.glob("*.tmp"))


def test_write_annotated_copy_does_not_publish_on_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = pa.table({"chain_length": [100], "chain_num_unresolved_residues": [25]})
    source = tmp_path / "source.parquet"
    output = tmp_path / "entry_chains.parquet"
    pq.write_table(table, source)

    def fail(*args: object, **kwargs: object) -> None:
        raise OSError("write failed")

    monkeypatch.setattr(pq.ParquetWriter, "write_table", fail)
    with pytest.raises(OSError, match="write failed"):
        write_annotated_chain_metadata(source, output)
    assert not output.exists()
    assert not list(tmp_path.glob("*.tmp"))
