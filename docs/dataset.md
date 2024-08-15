---
sd_hide_title: true
---

# Dataset

## Dataset reference

### Directory structure

:::{todo}
- give short description for each file/subdirectory
:::

```bash
2024-04/
|-- v1
    |-- clusters
    |   |-- metric={metric}
    |       |-- directed={directed}
    |           |-- threshold={threshold}.pqt
    |-- dbs
    |   |-- subdbs
    |       |-- {search_db}_{aln_type}
    |           |-- ids.tsv
    |-- entries
    |   |-- {two_char_code}.zip
    |-- fingerprints
    |   |-- ligands_per_inchikey.parquet
    |   |-- ligands_per_inchikey_ecfp4.npy
    |   |-- ligands_per_system.parquet
    |-- index
    |   |-- annotation_table.parquet
    |   |-- annotation_table_nonredundant.parquet
    |-- leakage
    |   |-- {split}_{metric}_{compare_pair}.parquet
    |-- ligand_scores
    |   |-- {fragment}.parquet
    |-- ligands
    |   |-- {fragment}.parquet
    |-- mmp
    |   |-- plinder_mmp_series.parquet
    |   |-- plinder_mms.csv.gz
    |-- scores
    |   |-- search_db=apo
    |   |   |-- {fragment}.parquet
    |   |-- search_db=holo
    |       |-- {fragment}.parquet
    |   |-- search_db=pred
    |       |-- {fragment}.parquet
    |-- splits
    |   |-- {category}
    |       |-- {split_config}.parquet
    |       |-- [{split_config}.yaml]
    |       |-- [{split_flavor}.csv]
    |-- systems
        |-- {two_char_code}.zip
```

# Main table

`annotation_table.parquet` harbours the main table listing all systems and their annotations.
The descriptions for each column in this table are listed here.

:::{raw} html
<div>
:::

```{include} table.html
```

:::{raw} html
</div>
:::