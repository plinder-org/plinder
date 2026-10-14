"""Regenerate the marked README counts table from a local release."""

import argparse
from pathlib import Path

import duckdb

START = "<!-- release-counts:start -->"
END = "<!-- release-counts:end -->"


def release_counts(release_dir: Path) -> str:
    """Count ligand records and receptor-polymer composition per assembly."""
    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        for name in (
            "annotation_table",
            "interface_annotation_table",
            "entry_chains",
            "entry_biounit_chains",
        ):
            con.read_parquet(
                str(release_dir / "index" / f"{name}.parquet")
            ).create_view(name)
        systems, ligands = con.sql("""
            SELECT count(DISTINCT system_id),
                count(*) FILTER (WHERE ligand_is_proper)
            FROM annotation_table
        """).fetchone()
        interfaces = con.sql(
            "SELECT count(*) FROM interface_annotation_table"
        ).fetchone()[0]
        protein, dna, rna, hybrid, complexes = con.sql("""
            WITH assemblies AS (
                SELECT b.entry_pdb_id, b.biounit_id,
                    count(*) FILTER (WHERE c.chain_receptor_type = 'protein') AS proteins,
                    count(*) FILTER (WHERE c.chain_receptor_type = 'dna') AS dna,
                    count(*) FILTER (WHERE c.chain_receptor_type = 'rna') AS rna,
                    count(*) FILTER (WHERE c.chain_receptor_type = 'dna+rna') AS hybrid
                FROM entry_biounit_chains b
                JOIN entry_chains c USING (entry_pdb_id, chain_asym_id)
                WHERE b.chain_role = 'receptor'
                GROUP BY 1, 2
            )
            SELECT
                count(*) FILTER (WHERE proteins = 1 AND dna + rna + hybrid = 0),
                count(*) FILTER (WHERE proteins = 0 AND dna = 1 AND rna + hybrid = 0),
                count(*) FILTER (WHERE proteins = 0 AND rna = 1 AND dna + hybrid = 0),
                count(*) FILTER (WHERE proteins = 0 AND hybrid = 1 AND dna + rna = 0),
                count(*) FILTER (WHERE proteins > 0 AND dna + rna + hybrid > 0)
            FROM assemblies
        """).fetchone()
    rows = [
        ("Ligand systems", systems),
        ("Proper ligands", ligands),
        ("Protein–protein interfaces", interfaces),
        ("Protein-monomer assemblies", protein),
        ("Nucleic-acid-monomer assemblies", dna + rna + hybrid),
        ("Protein–nucleic-acid assemblies", complexes),
    ]
    return "\n".join(
        [
            START,
            "| Release contents | Count |",
            "| --- | ---: |",
            *(f"| {label} | {count:,} |" for label, count in rows),
            "",
            f"Nucleic-acid monomers: {dna:,} DNA, {rna:,} RNA, and {hybrid:,} DNA/RNA hybrids.",
            END,
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release_dir", type=Path)
    parser.add_argument(
        "--readme", type=Path, help="Update this README; otherwise print the table."
    )
    args = parser.parse_args()
    block = release_counts(args.release_dir)
    if args.readme:
        text = args.readme.read_text()
        start = text.index(START)
        end = text.index(END, start) + len(END)
        args.readme.write_text(text[:start] + block + text[end:])
    else:
        print(block)


if __name__ == "__main__":
    main()
