"""Extract real, bounded release rows used by the two score/cluster guides."""

import argparse
import hashlib
import json
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("source", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
con = duckdb.connect()
con.execute("SET threads=2")
con.execute("SET memory_limit='4GB'")
records = []


def read(relative, predicate, parameters):
    path = args.source / relative
    glob = str(path / "**/*.parquet") if path.is_dir() else str(path)
    if relative == "exports/protein_similarity_scores":
        glob = sorted(
            {
                str(file)
                for entry in entries
                for file in path.glob(f"alignment_type=*/shard={entry[1:3]}.parquet")
            }
        )
    return con.execute(
        "SELECT * FROM read_parquet(?, union_by_name=true, hive_partitioning=true) WHERE "
        + predicate,
        [glob, *parameters],
    ).to_arrow_table()


def save(relative, table, selection):
    output = args.output / relative
    output.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, output, compression="zstd")
    records.append(
        {
            "path": relative,
            "rows": table.num_rows,
            "selection": selection,
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        }
    )
    print(relative, table.num_rows, output.stat().st_size, flush=True)


entries = {"4agi", "4agc", "7eek", "7cm8"}
cover_path = (
    "ligand_sampling/directed_set_cover/metric=pocket_qcov/threshold=70.parquet"
)
cover = read(cover_path, "ligand_id = ?", ["4agi__1__1.R"])
assert cover.num_rows
centroids = cover.column("centroid_ligand_id").to_pylist()
entries.update(value.split("__")[0] for value in centroids)
save(cover_path, cover, "4agi__1__1.R assignment at pocket_qcov=70")
annotation = read(
    "index/annotation_table.parquet", "entry_pdb_id IN ?", [sorted(entries)]
)
save(
    "index/annotation_table.parquet",
    annotation,
    "example entries and cover representative entries",
)
ligands = annotation.column("ligand_id").to_pylist()
systems = annotation.column("system_id").to_pylist()
smiles = sorted(
    set(v for v in annotation.column("ligand_smiles_id").to_pylist() if v is not None)
)
interfaces = read(
    "index/interface_annotation_table.parquet", "entry_pdb_id IN ?", [sorted(entries)]
)
save("index/interface_annotation_table.parquet", interfaces, "example entries")
systems.extend(interfaces.column("system_id").to_pylist())
for name, column, values in [
    ("entry_metadata", "entry_pdb_id", sorted(entries)),
    ("entry_sources", "entry_pdb_id", sorted(entries)),
    ("entry_chains", "entry_pdb_id", sorted(entries)),
    ("entry_biounit_chains", "entry_pdb_id", sorted(entries)),
    ("system_validation", "system_id", systems),
    ("ligand_clusters", "ligand_id", ligands),
    ("ligand_affinity", "ligand_id", ligands),
    ("ligand_pocket_membership", "ligand_id", ligands),
    ("interface_clusters", "system_id", systems),
    ("alignment_chain_lookup", "entry_pdb_id", sorted(entries)),
    ("ligand_pocket_residues", "entry_pdb_id", sorted(entries)),
]:
    relative = f"index/{name}.parquet"
    save(
        relative,
        read(relative, f'"{column}" IN ?', [values]),
        f"{column} in example subset",
    )
for name in ("sequence", "structure"):
    relative = f"protein_clusters/{name}.parquet"
    save(
        relative,
        read(relative, "entry_pdb_id IN ?", [sorted(entries)]),
        "example entries",
    )
relative = "index/ligand_mmp_pairs.parquet"
save(
    relative,
    read(
        relative, "ligand_smiles_id_1 IN ? OR ligand_smiles_id_2 IN ?", [smiles, smiles]
    ),
    "either member in example SMILES IDs",
)
for name, column, values in [
    ("ligand_similarity_scores", "ligand_id", ligands),
    ("interface_similarity_scores", "system", systems),
    (
        "interface_half_similarity_scores",
        "half_interface_id",
        [f"{s}::side={side}" for s in systems for side in (1, 2)],
    ),
]:
    relative = f"exports/{name}.parquet"
    save(
        relative,
        read(
            relative,
            f'"query_{column}" IN ? AND "target_{column}" IN ?',
            [values, values],
        ),
        "both endpoints in example subset",
    )
save(
    "ligand_scores/retained.parquet",
    read(
        "ligand_scores",
        "query_ligand_id IN ? AND target_ligand_id IN ?",
        [smiles, smiles],
    ),
    "both endpoints in example SMILES IDs",
)
for name in ("protein_similarity_scores", "monomer_similarity_scores"):
    save(
        f"exports/{name}/examples.parquet",
        read(
            f"exports/{name}",
            "query_entry IN ? AND target_entry IN ?",
            [sorted(entries), sorted(entries)],
        ),
        "both endpoints in example entries",
    )
for entry in ("4agi", "7cm8"):
    relative = (
        "alignment_cigars/search_db=holo/alignment_type=mmseqs/"
        f"shard={entry[1:3]}.parquet"
    )
    save(
        relative,
        read(relative, "query_entry = ? AND target_entry = ?", [entry, entry]),
        f"{entry} self alignments",
    )
metadata = {
    "release": "2026-09",
    "description": "Real release rows, not synthetic scores. Bounded integration-test input for the score and representative notebooks, not a standalone release or benchmark.",
    "entries": sorted(entries),
    "files": records,
}
(args.output / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
print(
    "TOTAL BYTES",
    sum((args.output / r["path"]).stat().st_size for r in records),
    flush=True,
)
