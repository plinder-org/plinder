"""Keep the bounded documentation input complete and reproducible."""

import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq


def test_docs_score_release_manifest():
    root = Path(__file__).resolve().parents[1] / "test_data/docs_score_release"
    manifest = json.loads((root / "manifest.json").read_text())
    assert manifest["release"] == "2026-09"
    assert {"4agi", "4agc", "7eek", "7cm8"} <= set(manifest["entries"])
    paths = [row["path"] for row in manifest["files"]]
    assert len(paths) == len(set(paths))
    assert "exports/ligand_similarity_scores.parquet" in paths
    assert "exports/interface_similarity_scores.parquet" in paths
    assert "exports/interface_half_similarity_scores.parquet" in paths
    assert (
        "alignment_cigars/search_db=holo/alignment_type=mmseqs/shard=cm.parquet"
        in paths
    )
    for row in manifest["files"]:
        path = root / row["path"]
        assert path.is_relative_to(root)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        assert pq.ParquetFile(path).metadata.num_rows == row["rows"]
    assert sum((root / path).stat().st_size for path in paths) < 5_000_000
