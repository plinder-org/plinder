# Dataset tutorial

## Getting the data

The PLINDER data is accessible from a _Google Cloud Platform_
[bucket](https://cloud.google.com/storage/docs/buckets), a container for cloud storage
of data.
The bucket URL of PLINDER is `gs://plinder`.
To interact with the data, we
[install `gsutil`](https://cloud.google.com/storage/docs/gsutil_install) first.

The PLINDER dataset is versioned via two parameters:

- `PLINDER_RELEASE`: the time stamp of the last RCSB sync
- `PLINDER_ITERATION`: iterative development within a release

For the purpose of this tutorial we set `PLINDER_ITERATION` to `tutorial`, to download
only a small manageable excerpt of the entries.

```console
$ export PLINDER_RELEASE=2024-04
$ export PLINDER_ITERATION=tutorial
$ gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/* ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```

The full dataset (`PLINDER_ITERATION=v1`) has a size of hundreds of GB, so you are
advised to have sufficient space for usage of the production dataset.

:::{note}
The version used for the preprint is `gs://plinder/2024-04/v1`, however, we plan to
release a newer version with updated annotations to be used for the
[MLSB challenge](https://www.mlsb.io/) (`gs://plinder/2024-06/v2`) in the future.
:::

## Understanding the directory structure

The directory downloaded from the bucket has the following structure:

```bash
2024-04/                     # the PLINDER release
|-- tutorial                 # the PLINDER iteration
|   |-- clusters             # pre-calculated cluster labels derived from the protein similarity dataset
|   |-- dbs                  # TSVs containing the raw files and IDs in the foldseek and mmseqs sub-databases
|   |-- entries              # raw annotations prior to consolidation (split by `two_char_code` and zipped)
|   |-- fingerprints         # index mapping files for the ligand similarity dataset
|   |-- index                # consolidated tabular annotations
|   |-- leakage              # leakage results
|   |-- ligand_scores        # ligand similarity parquet dataset
|   |-- ligands              # ligand data expanded from entries for computing similarity
|   |-- mmp                  # MMP and MMS data
|   |-- scores               # protein similarity parquet dataset
|   |-- splits               # split files and the configs used to generate them (if available)
|   |-- systems              # structure files for all systems (split by `two_char_code` and zipped)
```

The `systems`, `index`, `clusters`, `splits` and `leakage` directories are most the
important ones for PLINDER utilization and will be covered in the tutorial, while the
rest are for more curios users.

To download specific directories of interest, for example `splits`, run:

```bash
$ gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/splits ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```

## Unpacking the structure files

Similar to the
[PDB NextGen Archive](https://www.wwpdb.org/ftp/pdb-nextgen-archive-site), we split the
structures into subdirectories of chunks to make loading and querying speed palatable.

The structure files can be found in the subfolder
`plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems`.
To unpack the structures run

```bash
cd plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems; for i in ls *zip; do unzip $i; done
```

This will yield directories such as `6lpf__2__1.B__1.D`, which is what we call a PLINDER
system ID in the form
`<PDB ID>__<biological assembly>__<receptor chain ID>__<ligand chain ID>`.
Each system represent a complex between one or multiple proteins and a small molecules,
derived from a biological assembly in the PDB.
The directory contains _mmCIF_, _PDB_ and _SDF_ file formats as well as some additional
metadata files, for e.g. chain mapping and sequences.

## Exploring the annotation table

All systems are listed and annotated in the table contained in the
`index/annotation_table.parquet` file.
The [_Parquet_ format](https://parquet.apache.org/) is an efficient binary data format
for storing table data.
There is a multitude of tools that support reading `.parquet` files.
Here we will use the Python package `pandas` to inspect `annotation_table.parquet`.

```python
>>> df = pd.read_parquet("index/annotation_table.parquet")
>>> df.columns
Index(['entry_pdb_id', 'entry_release_date', 'entry_oligomeric_state',
       'entry_determination_method', 'entry_keywords', 'entry_pH',
       'entry_resolution', 'entry_rfree', 'entry_r', 'entry_clashscore',
       ...
       'ligand_interacting_ligand_chains_Pfam',
       'ligand_neighboring_ligand_chains_Pfam',
       'ligand_interacting_ligand_chains_PANTHER',
       'ligand_neighboring_ligand_chains_PANTHER',
       'system_ligand_chains_SCOP2', 'system_ligand_chains_SCOP2B',
       'protein_lddt_qcov_weighted_sum__100__strong__component',
       'pli_qcov__100__strong__component', 'system_ccd_codes', 'uniqueness'],
      dtype='object', length=488)
```

We see that the table contains hundreds of columns.
Each one is described in more detail in the
[Dataset Reference](#annotation-table-target).
The most important column is the `system_id`, which references the PLINDER systems
in the `systems` directory, we have already seen, but also in the other directories, we
are going to explore.

While `index/annotation_table.parquet` contains annotation for all PLINDER systems,
`index/annotation_table_nonredundant.parquet` contains a smaller set after
ligand-protein redundancy removal.

(cluster-target)=
## Inspecting the clusters

This directory is organized by the similarity metrics used for generating the clusters
and further nested by whether clustering is done with directed or undirected graph and
by the threshold for clustering.

Show nested structure

```console
$ tree clusters

clusters
├── metric=pli_qcov
│   ├── directed=False
│   │   ├── threshold=100.parquet
│   │   ├── threshold=50.parquet
│   │   ├── threshold=70.parquet
│   │   └── threshold=95.parquet
│   └── directed=True
│       ├── threshold=100.parquet
│       ├── threshold=50.parquet
│       ├── threshold=70.parquet
│       └── threshold=95.parquet
└── metric=pocket_fident
    └── directed=False
        ├── threshold=100.parquet
        ├── threshold=50.parquet
        ├── threshold=70.parquet
        └── threshold=95.parquet
```

As example, we will load the clusters based on pocket sequence similarity from an
undirected graph at a similarity threshold of 70 %.

```python
>>> import pandas as pd

>>> clus_file = "clusters/metric=pocket_fident/directed=False/threshold=70.parquet"
>>> df = pd.read_parquet(clus_file)
>>> df
                        system_id component
0               4wk0__1__1.A__1.F        c0
1               7q6e__1__1.B__1.R        c0
2               4dvl__1__1.A__1.I      c314
3           4zqc__1__1.B__1.H_1.I        c0
4       2gtl__1__1.E_1.H__1.X_1.Y        c0
...                           ...       ...
449378      1mqn__2__1.A_1.C__1.K        c0
449379          8sb4__1__1.A__1.O        c0
449380      4lys__1__1.A__1.B_1.D        c0
449381      3h2n__2__2.B__2.H_2.I        c0
449382          5thk__5__1.H__1.S     c1673

[449383 rows x 2 columns]
```

The table assigns a cluster to each system, depicted by the cluster ID in the
`component` column.
This means, all systems with the same cluster ID belong to the same cluster.

## Accessing the splits

The `splits` directory contains an index for _training-validation-test_ splits contained
in a single parquet file.
The _PL50_ split described in the [article](https://doi.org/10.1101/2024.07.17.603955)
can be found in `gs://plinder/2024-04/v1/splits/plinder-pl50.parquet`.

```python
>>> import pandas as pd
>>> df = pd.read_parquet("splits/plinder-pl50.parquet")
>>> df.head()
               system_id  split  cluster cluster_for_val_split
0      3grt__1__1.A__1.B  train  c100718                    c0
1  3grt__1__1.A_2.A__1.C  train  c196491                    c0
2      3grt__1__2.A__2.B  train  c100727                    c0
3  3grt__1__1.A_2.A__2.C  train  c234445                    c0
4      1grx__1__1.A__1.B  train  c186691                  c154
```

The columns are:

- `system_id`: the PLINDER system ID
- `split`: split category, either `train` (training set), `test` (test set),
  `val` (training set) or `removed`.
- `cluster`: cluster used in de-leaking the training vs test dataset.
- `cluster_for_val_split`: cluster used in sampling validation set from training set.

## Checking leakage

For each of the cluster metric described [above](cluster-target), we calculated the
maximum similarity for each system in a evaluation set (i.e. validation or test) to any
system the learning set (i.e. training or validation).
The parquet files are formatted like
`plinder_<split method>__<similarity metric>__<learning set>__<evaluation set>.parquet`.
For example `plinder-ECOD__pli_qcov__train_posebusters.parquet` means the file contains
a maximum `pli_qcov` similarity of posebusters (evaluation set) to any system in the
train set from a split based on the ECOD domain.

```console
$ ls leakage
plinder-ECOD__pli_qcov__train_posebusters.parquet
plinder-ECOD__pli_qcov__train_test.parquet
plinder-ECOD__pli_qcov__train_val.parquet
plinder-ECOD__pli_qcov__val_test.parquet
plinder-ECOD__pocket_lddt__train_posebusters.parquet
plinder-ECOD__pocket_lddt__train_test.parquet
plinder-ECOD__pocket_lddt__train_val.parquet
plinder-ECOD__pocket_lddt__val_test.parquet
plinder-ECOD__pocket_qcov__train_posebusters.parquet
plinder-ECOD__pocket_qcov__train_test.parquet
plinder-ECOD__pocket_qcov__train_val.parquet
plinder-ECOD__pocket_qcov__val_test.parquet
```

```python
>>> import pandas as pd

>>> leak_file = "leakage/plinder-ECOD__pli_qcov__train_posebusters.parquet"
>>> df = pd.read_parquet(leak_file)
>>> df
                      query_system  similarity
0                5s8i__1__1.A__1.B         100
1        5sak__1__1.A__1.B_1.C_1.E         100
2                5sb2__1__1.A__1.G         100
3    5sd5__1__1.A__1.B_1.C_1.D_1.E          89
4                5sis__1__1.A__1.D         100
..                             ...         ...
506      8h0m__1__1.A__1.B_1.C_1.F          75
507          8hfn__1__1.A__1.C_1.D          22
508          8ho0__1__1.A__1.B_1.C          75
509              8slg__1__1.A__1.C          89
510              8slg__1__1.B__1.I          75

[511 rows x 2 columns]
```
