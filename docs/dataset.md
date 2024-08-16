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

(annotation-table-target)=
### Annotation tables (`index/`)

Tables that lists all systems along with their annotations.

- `annotation_table.parquet`: Lists all systems and their annotations.
- `annotation_table_nonredundant.parquet`: Subset of systems without redundant systems.

:::{todo}
- Add `Mandatory` column
- Add `Example` column
- Sort by mandatory columns
- Show `system_id` as first column
- Plugin [datatables](https://datatables.net/) to enable sorting, filtering and pagination
- Then customize table style to fit all columns into the page
:::

:::{include} table.html
:::

### Clusters (`clusters/`)

This directory is nested by

- `metric`: the similarity metrics used for generating the clusters,
- `directed`: whether clustering is done with directed or undirected graph and
- `threshold`: similarity threshold in percent.

:::{list-table} `ligands_per_system.parquet` columns
:widths: 10 5 30
:header-rows: 1

*   - Name
    - Type
    - Description
*   - system_id
    - str
    - The PLINDER system ID
*   - component
    - str
    - The cluster ID of the cluster the system belongs to
:::

### Splits (`splits/`)

This directory contains:

- `plinder-pl50.parquet`: listing the split category for each system
- ...

:::{todo}
- Add missing descriptions for yaml and csv files
:::

:::{list-table} `plinder-pl50.parquet` columns
:widths: 10 5 30
:header-rows: 1

*   - Name
    - Type
    - Description
*   - system_id
    - str
    - The PLINDER system ID
*   - split
    - str
    - Split category: either `train` (training set), `test` (test set),
      `val` (training set) or `removed`
*   - cluster
    - str
    - Cluster used in de-leaking the training vs test dataset
*   - cluster_for_val_split
    - str
    - Cluster used in sampling validation set from training set
:::

### Leakage between split categories (`leakage/`)

Maximum similarity for each system in a evaluation set (i.e. validation or test) to any
system the learning set (i.e. training or validation).
The parquet files are formatted as

```
`plinder_<split_method>__<similarity_metric>__<learning_set>__<evaluation_set>.parquet`
```

where

- `split_method`: the base method the split is based on:
  - `ECOD`: split by ECOD domain
  - ...
- `similarity metric`: the similarity metric reported in the table:
  - `pli_qcov`: ...
  - `pocket_lddt`: ...
  - ...
- `learning_set`: the first split set of systems to be compared:
  - `train`: training set
  - `val`: validation set
- `evaluation set`: the second split set of systems to be compared:
  - `val`: validation set
  - `test`: test set

:::{todo}
- Fill missing methods etc. in list above
:::

:::{list-table} `plinder_<parameters>.parquet` columns
:widths: 10 5 30
:header-rows: 1

*   - Name
    - Type
    - Description
*   - system_id
    - str
    - The PLINDER system ID
*   - similarity
    - int
    - Maximum percentage of similarity between the split sets, based on the given metric
:::

### Small molecule fingerprints (`fingerprints/`)

Tables that contains all the ligand fingerprints used in calculating ligand similarity
stored in `ligand_scores`.

- `ligands_per_inchikey_ecfp4.npy`: `numpy` array of all-vs-all ECFP4 similarity.
- `ligands_per_system.parquet`: TODO.
- `ligands_per_inchikey.parquet`: subset of `ligands_per_system.parquet` with reduced number of columns.

:::{list-table} `ligands_per_system.parquet` columns
:widths: 10 5 30
:header-rows: 1

*   - Name
    - Type
    - Description
*   - system_id
    - str
    - The PLINDER system ID
*   - pdb_id
    -
    -
*   - ligand_rdkit_canonical_smiles
    -
    -
*   - ligand_ccd_code
    -
    -
*   - ligand_id
    -
    -
*   - inchikeys
    -
    -
*   - number_id_by_inchikeys
    -
    - The index of `ligands_per_inchikey_ecfp4.npy` corresponds to
      `number_id_by_inchikeys`, so this can be used to map scores in `ligand_score`
      sub-directory into human readable form.
*   - number_id_by_inchikeys_withl
    -
    -
*   - multi_number_id_by_inchikeys
    -
    -
*   - multi_number_id_by_inchikeys_withl
    -
    -
:::

:::{todo}
- Add sections for missing tables
- Fill missing elements into tables
:::
