# Dataset tutorial

## Getting the data

The PLINDER data is accessible from a _Google Cloud Platform_
[bucket](https://cloud.google.com/storage/docs/buckets), a container for cloud storage
of data.
The bucket URL of PLINDER is `gs://plinder`.

The PLINDER dataset is versioned via two parameters:

- `PLINDER_RELEASE`: the time stamp of the last RCSB sync
- `PLINDER_ITERATION`: iterative development within a release

There are two ways to obtain the data:

1. Use the `plinder` python package and corresponding API
    - `pip install plinder`
2. Use the `gsutil` command line tool directly
   - [installing `gsutil`](https://cloud.google.com/storage/docs/gsutil_install)

For the purpose of this tutorial we set `PLINDER_ITERATION` to `tutorial`, to download
only a small manageable excerpt of the entries.

Using the `plinder` package:
```bash
# adding --yes will skip all confirmation prompts
plinder_download --release 2024-06 --iteration tutorial --yes
```

Using `gsutil`:
```console
$ export PLINDER_RELEASE=2024-06
$ export PLINDER_ITERATION=tutorial
$ mkdir -p ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
$ gsutil -m cp -r "gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/*" ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```

The full dataset (`PLINDER_ITERATION=v2`) has a size of hundreds of GB, so you are
advised to have sufficient space for usage of the production dataset.

:::{note}
The versions used for the preprint are `gs://plinder/2024-04/v1` (full dataset) and `gs://plinder/2024-04/v0` (non-redundant set used to train DffDock). However, the current version with updated annotations to be used for the
[MLSB challenge](https://www.mlsb.io/) is `gs://plinder/2024-06/v2`.
:::

## Understanding the directory structure

The directory downloaded from the bucket has the following structure:

```bash
2024-06/                     # The PLINDER release
|-- tutorial                 # The PLINDER iteration
|   |-- clusters             # Pre-calculated cluster labels derived from the protein similarity dataset
|   |-- dbs                  # TSVs containing the raw files and IDs in the foldseek and mmseqs sub-databases
|   |-- entries              # Raw annotations prior to consolidation (split by `two_char_code` and zipped)
|   |-- fingerprints         # Index mapping files for the ligand similarity dataset
|   |-- index                # Consolidated tabular annotations
|   |-- ligand_scores        # Ligand similarity parquet dataset
|   |-- ligands              # Ligand data expanded from entries for computing similarity
|   |-- linked_structures    # Apo and predicted structures linked to their holo systems
|   |-- links                # Apo and predicted structures similarity to their holo structures
|   |-- mmp                  # Ligand matched molecular pairs (MMP) and series (MMS) data
|   |-- scores               # Protein similarity parquet dataset
|   |-- splits               # Split files and the configs used to generate them (if available)
|   |-- systems              # Structure files for all systems (split by `two_char_code` and zipped)
```

The `systems`, `index`, `clusters` and `splits` directories are most the
important ones for PLINDER utilization and will be covered in the tutorial, while the
rest are for more curious users.

To download specific directories of interest, for example `splits`, run:

```bash
$ gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/splits ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```

## Unpacking the structure files

If you used the `plinder_download` command, you can skip this section.

Similar to the
[PDB NextGen Archive](https://www.wwpdb.org/ftp/pdb-nextgen-archive-site), we split the
structures into subdirectories of chunks (using two penultimate characters of PDB code) to make loading and querying speed palatable.

The structure files can be found in the subfolder
`~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems`.
To unpack the structures run

```bash
cd ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems; for i in `ls *zip`; do unzip $i; touch ${i//.zip/}_done; done
```

This will yield directories such as `7eek__1__1.A__1.I`, which is what we call a PLINDER
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
Here we will use the Python package [`pandas`](https://pandas.pydata.org) to inspect `annotation_table.parquet`.

```python
>>> df = pd.read_parquet("index/annotation_table.parquet")
>>> df.columns
Index(['entry_pdb_id', 'entry_release_date', 'entry_oligomeric_state',
       'entry_determination_method', 'entry_keywords', 'entry_pH',
       'entry_resolution', 'entry_rfree', 'entry_r', 'entry_clashscore',
       ...
       'ligand_interacting_ligand_chains_UniProt',
       'system_ligand_chains_PANTHER', 'ligand_interacting_ligand_chains_Pfam',
       'ligand_neighboring_ligand_chains_Pfam',
       'ligand_interacting_ligand_chains_PANTHER',
       'ligand_neighboring_ligand_chains_PANTHER',
       'system_ligand_chains_SCOP2', 'system_ligand_chains_SCOP2B',
       'pli_qcov__100__strong__component',
       'protein_lddt_qcov_weighted_sum__100__strong__component'],
      dtype='object', length=500)
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

clusters/
├── cluster=communities
│   └── directed=False
│       ├── metric=pli_qcov
│       │   ├── threshold=100
│       │   │   └── data.parquet
│       │   ├── threshold=50
│       │   │   └── data.parquet
│       │   ├── threshold=70
│       │   │   └── data.parquet
│       │   └── threshold=95
│       │       └── data.parquet
│       ├── metric=pli_unique_qcov
│       │   ├── threshold=100
│       │   │   └── data.parquet
│       │   ├── threshold=50
│       │   │   └── data.parquet
│       │   ├── threshold=70
│       │   │   └── data.parquet
│       │   └── threshold=95
│       │       └── data.parquet
```

As example, we will load the clusters based on pocket sequence similarity from an
undirected graph at a similarity threshold of 70 %.

```python
>>> import pandas as pd

>>> clus_file = "clusters/cluster=communities/directed=False/metric=pli_qcov/threshold=70/data.parquet"
>>> df = pd.read_parquet(clus_file)
>>> df
                            system_id   label    metric      cluster  directed  threshold
0                   3mj2__1__1.A__1.B      c0  pli_qcov  communities     False         70
1       4dh8__1__1.A_1.B__1.C_1.D_1.E      c0  pli_qcov  communities     False         70
2                   7akb__1__1.A__1.C      c0  pli_qcov  communities     False         70
3                   7mgj__2__1.B__1.F      c0  pli_qcov  communities     False         70
4                   4fr4__6__1.F__1.S      c0  pli_qcov  communities     False         70
...                               ...     ...       ...          ...       ...        ...
479806          7xpv__1__2.A__2.C_2.D  c77190  pli_qcov  communities     False         70
479807              4ret__1__1.A__1.I  c77191  pli_qcov  communities     False         70
479808              7ks9__1__1.C__1.R  c77192  pli_qcov  communities     False         70
479809              7s6n__1__1.A__1.H  c77193  pli_qcov  communities     False         70
479810              7sc5__1__1.A__1.H  c77194  pli_qcov  communities     False         70

[479811 rows x 6 columns]
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
>>> df = pd.read_parquet("splits/split.parquet")
>>> df.head()
               system_id            uniqueness  split cluster  ... system_proper_num_interactions  system_proper_ligand_max_molecular_weight  system_has_binding_affinity  system_has_apo_or_pred
0  101m__1__1.A__1.C_1.D  101m__A__C_D_c188899  train     c14  ...                             20                                 616.177293                        False                   False
1      102m__1__1.A__1.C    102m__A__C_c237197  train     c14  ...                             20                                 616.177293                        False                    True
2  103m__1__1.A__1.C_1.D  103m__A__C_D_c252759  train     c14  ...                             16                                 616.177293                        False                   False
3  104m__1__1.A__1.C_1.D  104m__A__C_D_c274687  train     c14  ...                             21                                 616.177293                        False                   False
4  105m__1__1.A__1.C_1.D  105m__A__C_D_c221688  train     c14  ...                             20                                 616.177293                        False                   False

[5 rows x 13 columns]
```

The columns are:

- `system_id`: The PLINDER system ID
- `uniqueness`: An id tag that captures system redundancy based on ligand and pocket similarity
- `split`: Split category, either `train` (training set), `test` (test set)
- `cluster`: Cluster metric used in sampling test dataset.
- `cluster_for_val_split`: Cluster metric used in sampling validation set from training set.
- `system_pass_validation_criteria`: Boolean indicating whether a system pass all the quality criteria
- `system_pass_statistics_criteria`: Boolean indicating whether a system pass the desired statistics criteria
- `system_proper_num_ligand_chains`: Number of chains ligands that are not ions or artifacts
- `system_proper_pocket_num_residues`: Number of pocket residues around ligands that are not ions or artifacts
- `system_proper_num_interactions`: Number of interactions based on ligands that are not ions or artifacts
- `system_proper_ligand_max_molecular_weight`: Maximum molecular weight of ligands that are not ions or artifacts
- `system_has_binding_affinity`: Boolean indicator of whether a system has binding affinity or not
- `system_has_apo_or_pred`: Boolean indicator of whether a system apo or predicted structures linked
