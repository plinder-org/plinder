# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Derive coordinate-completeness annotations without re-ingesting structures."""

from __future__ import annotations

import argparse
from pathlib import Path
from uuid import uuid4

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

UNRESOLVED_FRACTION_COLUMN = "chain_unresolved_fraction"


def annotate_chain_metadata(table: pa.Table) -> pa.Table:
    """Add the unresolved/SEQRES length fraction; invalid counts yield null."""
    length = table["chain_length"]
    unresolved = table["chain_num_unresolved_residues"]
    valid = pc.fill_null(
        pc.and_(
            pc.greater(length, 0),
            pc.and_(pc.greater_equal(unresolved, 0), pc.less_equal(unresolved, length)),
        ),
        False,
    )
    # Use a nonzero denominator even for invalid rows: Arrow evaluates both
    # branches before selecting the valid values.
    denominator = pc.if_else(valid, length, 1).cast(pa.float64())
    fraction = pc.if_else(
        valid,
        pc.divide(unresolved.cast(pa.float64()), denominator),
        pa.scalar(None, type=pa.float64()),
    )
    index = table.schema.get_field_index(UNRESOLVED_FRACTION_COLUMN)
    if index >= 0:
        return table.set_column(index, UNRESOLVED_FRACTION_COLUMN, fraction)
    return table.append_column(UNRESOLVED_FRACTION_COLUMN, fraction)


def write_annotated_chain_metadata(
    source: Path, output: Path, *, batch_size: int = 100_000
) -> int:
    """Write a bounded-memory annotated copy, never overwriting a release table."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if source.resolve() == output.resolve() or output.exists():
        raise FileExistsError(f"choose a new output path: {output}")
    parquet = pq.ParquetFile(source)
    empty = pa.Table.from_batches([], schema=parquet.schema_arrow)
    schema = annotate_chain_metadata(empty).schema
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{uuid4().hex}.tmp")
    rows = 0
    try:
        with pq.ParquetWriter(temporary, schema, compression="zstd") as writer:
            for batch in parquet.iter_batches(batch_size=batch_size):
                table = annotate_chain_metadata(pa.Table.from_batches([batch]))
                writer.write_table(table, row_group_size=batch_size)
                rows += len(table)
        # Atomic, no-clobber publication of the new file.
        output.hardlink_to(temporary)
    finally:
        temporary.unlink(missing_ok=True)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--batch-size", type=int, default=100_000)
    args = parser.parse_args()
    rows = write_annotated_chain_metadata(
        args.source, args.output, batch_size=args.batch_size
    )
    print(f"Annotated {rows:,} chains: {args.output}")


if __name__ == "__main__":
    main()
