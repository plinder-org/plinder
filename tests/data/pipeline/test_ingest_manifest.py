# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0

from pathlib import Path

from plinder.data.pipeline.ingest_manifest import (
    EntryInput,
    balance_entries,
    discover_entries,
    write_manifest,
)


def test_discover_entries_tracks_optional_validation(tmp_path: Path) -> None:
    cif_root = tmp_path / "cif"
    validation_root = tmp_path / "validation"
    for pdb_id, content in (("1abc", b"one"), ("2def", b"second")):
        cif_path = (
            cif_root
            / pdb_id[1:3]
            / f"pdb_0000{pdb_id}"
            / f"pdb_0000{pdb_id}_xyz-enrich.cif.gz"
        )
        cif_path.parent.mkdir(parents=True)
        cif_path.write_bytes(content)
    validation_path = validation_root / "ab" / "1abc" / "1abc_validation.xml.gz"
    validation_path.parent.mkdir(parents=True)
    validation_path.touch()

    entries = discover_entries(cif_root, validation_root)

    assert entries == [
        EntryInput("1abc", 3, True),
        EntryInput("2def", 6, False),
    ]

    unchecked = discover_entries(
        cif_root,
        validation_root,
        check_validation=False,
    )
    assert [entry.validation_exists for entry in unchecked] == [None, None]


def test_balanced_manifest_has_fixed_slices_and_summary(tmp_path: Path) -> None:
    entries = [
        EntryInput(f"{index}abc", size, True)
        for index, size in enumerate((100, 90, 80, 20, 10), start=1)
    ]

    batches = balance_entries(entries, batch_size=2)

    assert [len(batch) for batch in batches] == [2, 2, 1]
    totals = [sum(entry.cif_size for entry in batch) for batch in batches]

    manifest = tmp_path / "pdb_ids.txt"
    summary = write_manifest(
        entries=entries,
        output_path=manifest,
        batch_size=2,
    )
    assert len(manifest.read_text().splitlines()) == len(entries)
    assert summary["batch_count"] == 3
    assert summary["maximum_batch_cif_bytes"] == max(totals)
    assert summary["cif_size_percentiles"]["maximum"] == 100
    assert summary["largest_entries"][0]["cif_size"] == 100
    assert manifest.with_suffix(".txt.json").is_file()
    inventory = manifest.with_suffix(".txt.entries.tsv")
    assert inventory.is_file()
    assert inventory.read_text().splitlines()[0] == (
        "pdb_id\tcif_size\tvalidation_exists"
    )

    manifest_ids = manifest.read_text().splitlines()
    assert [manifest_ids[index : index + 2] for index in range(0, len(entries), 2)] == [
        [entry.pdb_id for entry in batch] for batch in batches
    ]


def test_balance_entries_distributes_large_inputs_across_full_batches() -> None:
    entries = [
        EntryInput(f"{index}abc", size, True)
        for index, size in enumerate((100, 90, 80, 20, 10, 5), start=1)
    ]

    batches = balance_entries(entries, batch_size=2)
    totals = [sum(entry.cif_size for entry in batch) for batch in batches]

    assert [len(batch) for batch in batches] == [2, 2, 2]
    assert max(totals) - min(totals) <= 15
