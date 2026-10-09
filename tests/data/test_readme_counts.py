"""Check README counting units and assembly composition."""

from pathlib import Path
from runpy import run_path

import duckdb


def test_release_counts(tmp_path: Path) -> None:
    index = tmp_path / "index"
    index.mkdir()
    tables = {
        "annotation_table": "SELECT * FROM (VALUES ('s1', true), ('s1', false), ('s2', true)) t(system_id, ligand_is_proper)",
        "interface_annotation_table": "SELECT 'i1' AS system_id",
        "entry_chains": """SELECT * FROM (VALUES
            ('p', 'A', 'protein'), ('p', 'L', 'rna'),
            ('d', 'A', 'dna'), ('r', 'A', 'rna'), ('h', 'A', 'dna+rna'),
            ('c', 'A', 'protein'), ('c', 'B', 'dna'),
            ('m', 'A', 'protein'), ('m', 'B', 'protein')
        ) t(entry_pdb_id, chain_asym_id, chain_receptor_type)""",
        "entry_biounit_chains": """SELECT * FROM (VALUES
            ('p', '1', 'A', 'receptor'), ('p', '1', 'L', 'ligand'),
            ('d', '1', 'A', 'receptor'), ('d', '2', 'A', 'receptor'),
            ('r', '1', 'A', 'receptor'), ('h', '1', 'A', 'receptor'),
            ('c', '1', 'A', 'receptor'), ('c', '1', 'B', 'receptor'),
            ('m', '1', 'A', 'receptor'), ('m', '1', 'B', 'receptor')
        ) t(entry_pdb_id, biounit_id, chain_asym_id, chain_role)""",
    }
    with duckdb.connect() as con:
        for name, sql in tables.items():
            con.sql(sql).write_parquet(str(index / f"{name}.parquet"))
    script = Path(__file__).resolve().parents[2] / "scripts" / "update_readme_counts.py"
    counts = run_path(str(script))["release_counts"](tmp_path)
    assert (
        counts
        == """<!-- release-counts:start -->
| Release contents | Count |
| --- | ---: |
| Ligand systems | 2 |
| Proper ligands | 2 |
| Protein–protein interfaces | 1 |
| Protein-monomer assemblies | 1 |
| Nucleic-acid-monomer assemblies | 4 |
| Protein–nucleic-acid assemblies | 1 |

Nucleic-acid monomers: 2 DNA, 1 RNA, and 1 DNA/RNA hybrids.
<!-- release-counts:end -->"""
    )
