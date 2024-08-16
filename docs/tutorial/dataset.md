# Dataset tutorial

## Dataset versioning

We version the plinder dataset with two controls:

`PLINDER_RELEASE`: the month stamp of the last RCSB sync
`PLINDER_ITERATION`: value that enables iterative development within a release

Changelog:

2024-06/v2 (Upcoming):

- Improved SDF saving to handle some bond order issues
- Updated system definition to be more stable and independent of PLIP
- Added binding affinities from BindingDB and added "has_affinity" as priority for test split
- Annotated all crystal contacts
- Improved covalency detection
  2024-04/v1 (Current): Version with redundancy removal by protein pocket and ligand similarity.

2024-04/v0: Version used to re-train DiffDock in the paper, with redundancy removal based on \<pdbid\>\_\<ligand ccd codes\>

## Getting the data

We provide read access to plinder data via GCP bucket path `gs://plinder`. To interact with the data, users are advised to install `gsutil` via [gsutil installation guide](https://cloud.google.com/storage/docs/gsutil_install).

To download the full dataset, do:

```bash
export PLINDER_RELEASE=2024-04 # Current release
export PLINDER_ITERATION=tutorial # tutorial iteration
gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/* ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```

This full dataset is XGB in size; to download run `export PLINDER_RELEASE=2024-04 && export PLINDER_ITERATION=v1`, followed by the `gsutil` command above.

NOTE!: the version used for the preprint is gs://plinder/2024-04/v1, however, we plan to release a newer version with updated annotations to be used for the MLSB challenge (`gs://plinder/2024-06/v2`) in the future.

Downloading with the commands above will yield the following folder structure shown in **Understanding the directory structure** section, with the `systems`, `splits`, `clusters` and `index` being most important for direct usage while the rest are for more curios users.

## Understanding the directory structure

The directory downloaded from the bucket has the following structure:

```bash
2024-04/                     # the "`plinder` release" (`PLINDER_RELEASE`)
|-- tutorial                 # the "`plinder` iteration" (`PLINDER_ITERATION`); remember to change this to current iteration for production run
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

To download specific directories of interest, for example `splits`, run:

```bash
gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/splits ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/
```

### Unpacking the structure files:

Identical to the PDB NextGen Archive, we split the structures into subdirectories of chunks to make loading and querying speed palatable.

The structure files can be found in the subfolder `plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems`. To unpack the structures do the following

```bash
cd plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems; for i in ls *zip; do unzip $i; done
```

This will yield directories such as `6lpf__2__1.B__1.D`, which is what we call a plinder system id of the form `<pdbid>__<biounit>__<chain_ids_of_receptor>__<chain_ids_of_ligand>`. The directory contains mmcif, pdb and sdf file formats as well as some additional metadata (like chain mapping and fastas).

The naming of the directory is by system - since a given PDB entry can contain multiple protein-ligand complexes, we assign a plinder system id to every system. These system ids are also used for retrieving the splits, as shown in the next step

With the `system_id` contained in these split files, you can load the respective train, val & test splits after unzipping the systems directory. E.g. `~/.local/share/plinder/2024-04/v1/systems/6lpf__2__1.B__1.D/system.cif` will contain the full mmcif of the system. We also provide cif files of seperated receptors (\*/receptor.cif) and ligands (\*/ligand_files/\*.sdf) as well as pdb files (\*/receptor.pdb) but strongly encourage cif, pdb is considered a legacy file format.

Note: a non-redundant and single-ligand smaller subset of this version was used to train diffdock in the paper and is available at `2024-04/v0`.

The folder also contains a .yaml which is the config used to generate the split and can be ignored unless you want to reproduce the splits.

## Accessing the splits files:

`~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/splits` contains an index for splits contained in a single parquet file. The most current split is `gs://plinder/2024-04/v1/splits/plinder-pl50.parquet` containing the pl50 split from the preprint. (Again: note this will be shortly updated to `v2` and the `v2` split will be used for the MLSB leaderboard).

```python
>>> import pandas as pd
>>> df = pd.read_parquet("~/.local/share/plinder/2024-04/tutorial/splits/plinder-pl50.parquet")
>>> df.head()
               system_id  split  cluster cluster_for_val_split
0      3grt__1__1.A__1.B  train  c100718                    c0
1  3grt__1__1.A_2.A__1.C  train  c196491                    c0
2      3grt__1__2.A__2.B  train  c100727                    c0
3  3grt__1__1.A_2.A__2.C  train  c234445                    c0
4      1grx__1__1.A__1.B  train  c186691                  c154
```

The columns are:

- `system_id`: plinder id assigned to every protein-ligand system in a given PDB entry
- `split`: split category, either `train` (training set), `test` (test set), `val` (training set)or `removed`.
- `cluster`: cluster used in de-leaking the training vs test dataset.
- `cluster`: cluster used in sampling validation set from training set

## Inspecting the annotation files

For ease of access to the annotations reported in our publication, we saved selected (488) annotations in a parquet format stored in `index` sub-directory. To list available annotation files, do:

```bash
❯ ls ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/index/
annotation_table.parquet
annotation_table_nonredundant.parquet
```

`annotation_table.parquet` contains annotation for all the systems in RSCB PDB, while `annotation_table_nonredundant.parquet` contains a smaller set after ligand protein redundancy removal. To inspect the data, do

```python
>>> df = pd.read_parquet(anno_file)
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

## Inspecting the entries

```bash
❯ tree  ~/.local/share/plinder/share/${PLINDER_RELEASE}/${PLINDER_ITERATION}

└── 2024-04
    └── v1
        └── entries
            ├── da.zip
            ├── lp.zip
            └── to.zip
```

`plinder.data.utils.annotations.aggregate_annotations.Entry` (to be referred to as `Entry`) provides full access to all the annotations. `entries` sub-directory contains archived serialized version of `Entry`, containing enough information to re-constitute the data class if needed. `Entry` is a higher-level abstraction of the content of a given RCSB PDB. It contains all the information needed to derive all `plinder.data.utils.ligand_utils.System` data objects in one RCSB entry. When unpacked, each of the `.zip` file yields `.json` files of all PDBID represented the the pdb middle two character in the `.zip` file name.

You can re-constitute the `Entry` object directly by running:

```python
>>> import json
>>> from pathlib import Path

>>> import pandas as pd

>>> json_file = "~/.local/share/plinder/2024-04/tutorial/entries/1da0.json"
>>> with Path(json_file).open() as f:
...    entry = json.load(f)
>>> df = pd.json_normalize(entry)
>>> df.T.head()
                                      0
pdb_id                             1da0
release_date                 1992-10-17
oligomeric_state                dimeric
determination_method  X-RAY DIFFRACTION
keywords                            DNA
```

## Inspecting the clusters

clusters
This directory is organized by the similarity metrics used in generating the clusters and further nested by whether clustering is done with directed or undirected graph and the threshold for clustering.

```bash
❯ ls  ~/.local/share/plinder/2024-04/v1/clusters/
metric=pli_qcov/ #  Tracks how many aligned residues from two pockets share similar protein-ligand interaction.
metric=pocket_fident/ # Tracks how similar is the sequence of a given pocket region to a (possibly non-pocket) region of another protein
metric=pocket_fident_qcov/ # Tracks how similar is the sequence of a given pocket region to another pocket-containing region of another protein
metric=pocket_lddt/ # Tracks how well a given pocket structurally resembles a (possibly non-pocket) region of another protein
metric=pocket_lddt_qcov/ # Tracks how well a given pocket structurally resembles a ligand-binding pocket region of another protein.
metric=pocket_qcov/  # Degree of sequence match to query pocket residues
metric=protein_fident_max/ # Maximum protein sequence identity when all-vs-all similarity is computed between multiple proteins
metric=protein_fident_qcov_max/ # Maximum protein sequence identity with reference to query protein, when all-vs-all similarity is computed between multiple proteins
metric=protein_fident_qcov_weighted_max/ # Weighted maximum (by protein length) protein sequence identity with reference to query protein, when all-vs-all similarity is computed between multiple proteins
metric=protein_fident_qcov_weighted_sum/ # Weighted sum (by protein length) protein sequence identity with reference to query protein, when all-vs-all similarity is computed between multiple proteins
metric=protein_fident_weighted_max/ # Weighted max (by protein length) protein sequence identity with reference to query protein, when all-vs-all similarity is computed between multiple proteins
metric=protein_fident_weighted_sum/
metric=protein_lddt_max/
metric=protein_lddt_qcov_max/
metric=protein_lddt_qcov_weighted_max/
metric=protein_lddt_qcov_weighted_sum/
metric=protein_lddt_weighted_max/
metric=protein_lddt_weighted_sum/
metric=protein_qcov_max/
metric=protein_qcov_weighted_max/
metric=protein_qcov_weighted_sum/
metric=protein_seqsim_max/
metric=protein_seqsim_qcov_max/
metric=protein_seqsim_qcov_weighted_max/
metric=protein_seqsim_qcov_weighted_sum/
metric=protein_seqsim_weighted_max/
metric=protein_seqsim_weighted_sum/
```

Show nested structure

```bash
❯ tree  ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/clusters

plinder
└── 2024-04
    └── v1
        └── clusters
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
            ├── metric=pocket_fident
            │   ├── directed=False
            │   │   ├── threshold=100.parquet
            │   │   ├── threshold=50.parquet
            │   │   ├── threshold=70.parquet
            │   │   └── threshold=95.parquet
```

Each of the sub-directories is further nested by the whether clustering is done with directed or undirected graph and the threshold for clustering.

To load and show the content of one of the cluster files, for example, loading cluster based on pocket sequence similarity from an undirected graph at `threshold=70`, do:

```python
>>> import pandas as pd

>>> clus_file = "~/.local/share/plinder/2024-04/tutorial/clusters/metric=pocket_fident/directed=False/threshold=70.parquet"
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

`component` represents the cluster id of each of the `system_id`'s.

## Inspecting leakage files:

For each of the metric described in **Inspecting the clusters** section above, we calculated the maximum similarity for each system in a evaluation (could be validation or test) set to any system the learning (could be training or validation) set. The parquet files are formatted like `<split_method>__<similarity_metric>__<learning_set>__<evaluation_set>.parquet`. For example `plinder-ECOD**pli_qcov__train_posebusters.parquet`mean the file contains maximum`pli_qcov` similarity of posebusters (evaluation set) to any system in the train set from a split based on ECOD domain.

```bash
❯ ls  ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/leakage/
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

>>> leak_file = "/Users/yusuf/.local/share/plinder/2024-04/tutorial/leakage/plinder-ECOD__pli_qcov__train_posebusters.parquet"
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

## Inspecting fingerprint:

```bash
ls  ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/fingerprints
ligands_per_inchikey.parquet
ligands_per_inchikey_ecfp4.npy
ligands_per_system.parquet
```

This contains all the ligand fingerprints used in calculating ligand similarity we stored in `ligand_scores`. `ligands_per_inchikey_ecfp4.npy` is a `numpy` array of all-vs-all ECFP4 similarity. Each row of `ligands_per_system.parquet` contains annotations like `pdb_id`, `system_id`, `ligand_rdkit_canonical_smiles`, `ligand_ccd_code`, `ligand_id`, `inchikeys`, `number_id_by_inchikeys` `number_id_by_inchikeys_withl`, `multi_number_id_by_inchikeys`, `multi_number_id_by_inchikeys_withl`. The index of the `ligands_per_inchikey_ecfp4.npy` correspoonds to `number_id_by_inchikeys`; so, this can be used to map scores in `ligand_score` sub-directory into human readable form. `ligands_per_inchikey.parquet` is just a subset of `ligands_per_system.parquet` with only columns `inchikeys`, `ligand_rdkit_canonical_smiles`, `number_id_by_inchikeys` and `number_id_by_inchikeys_withl`.

```python
import numpy as np
import pandas as pd

>>> fingerprints_arr = np.load("~/.local/share/plinder/2024-04/tutorial/fingerprints/ligands_per_inchikey_ecfp4.npy")
>>> data
array([[0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 0, 0, 0],
       ...,
       [0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 0, 0, 0]], dtype=int8)


>>> data = pd.read_parquet("~/.local/share/plinder/2024-04/tutorial/fingerprints/ligands_per_system.parquet")
>>> data
       pdb_id                  system_id  ... multi_number_id_by_inchikeys multi_number_id_by_inchikeys_withl
0        4n14          4n14__1__1.A__1.B  ...                        48862                             l48862
1        4n15          4n15__1__1.A__1.C  ...                          354                               l354
2        4n16  4n16__1__1.A__1.B_1.F_1.H  ...           2674__29542__30710              l2674__l29542__l30710
3        4n16  4n16__1__1.A__1.B_1.F_1.H  ...           2674__29542__30710              l2674__l29542__l30710
4        4n16  4n16__1__1.A__1.B_1.F_1.H  ...           2674__29542__30710              l2674__l29542__l30710
...       ...                        ...  ...                          ...                                ...
644266   5upy  5upy__1__2.A_3.A__3.B_3.C  ...                 11911__31480                     l11911__l31480
644267   5upy  5upy__1__2.A_3.A__3.B_3.C  ...                 11911__31480                     l11911__l31480
644268   5upy  5upy__1__2.A_4.A__2.B_2.C  ...                 11911__31480                     l11911__l31480
644269   5upy  5upy__1__2.A_4.A__2.B_2.C  ...                 11911__31480                     l11911__l31480
644270   5upz      5upz__1__1.A_1.B__1.C  ...                        30533                             l30533

[644271 rows x 10 columns]
>>> data.columns
Index(['pdb_id', 'system_id', 'ligand_rdkit_canonical_smiles',
       'ligand_ccd_code', 'ligand_id', 'inchikeys', 'number_id_by_inchikeys',
       'number_id_by_inchikeys_withl', 'multi_number_id_by_inchikeys',
       'multi_number_id_by_inchikeys_withl'],
      dtype='object')
```

## Inspecting scores:

### Protein similarity

If you are interested in query the underlying similarity dataset used in computing the clusters, check `scores` and `ligand_score` The `scores` sub-directory contains protein similarity

```bash
❯ ls ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/scores/
search_db=apo
search_db=holo
search_db=pred

❯ ls  ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/scores/search_db=holo/|head
search_db=holo/holo_000d09b6902fb696a248c25ff4d25832.parquet
search_db=holo/holo_002e9b846b84b7cf85856abb02ba6b8e.parquet
search_db=holo/holo_004437a10770fd1cd3c7b6cd47f3e3ba.parquet
search_db=holo/holo_00babc3ceb84a4708a68cfc5710667ce.parquet
search_db=holo/holo_00c24c1662fb0001b2a89331460c030f.parquet
search_db=holo/holo_00ed1d6cf2950410127643989250bfe4.parquet
search_db=holo/holo_0134ec3aeee7f77ea4a6cacaf9bfa97d.parquet
search_db=holo/holo_019a3fbb41352c984514ece5e9d45e4f.parquet
search_db=holo/holo_01c5544cb9a475d47ee22bb75ac52e70.parquet
search_db=holo/holo_01e452953c9f26dc01961dbed2707f5b.parquet
```

```python
>>> df = pd.read_parquet("/Users/yusuf/.local/share/plinder/2024-04/tutorial/scores", filters=[("search_db", "==", "holo"), ("metric", "==", "pli_qcov"), ("similarity", ">=", 70)])
>>> df
                 query_system                  target_system  protein_mapping protein_mapper  ...  source    metric mapping search_db
0       3lpn__1__1.A_1.B__1.C          3lrt__1__1.A_1.B__1.F  1.B:1.B;1.A:1.A       foldseek  ...    both  pli_qcov    None      holo
1       3lpn__1__1.A_1.B__1.C      3mbi__1__1.A_1.C__1.F_1.G  1.A:1.A;1.B:1.C       foldseek  ...    both  pli_qcov    None      holo
2       3lpn__1__1.A_1.B__1.C  6nfe__1__1.A_1.B__1.E_1.J_1.P  1.A:1.A;1.B:1.B       foldseek  ...    both  pli_qcov    None      holo
3       5dap__1__2.A__2.B_2.C          5dap__1__1.A__1.B_1.C          2.A:1.A       foldseek  ...    both  pli_qcov    None      holo
4       5dap__1__2.A__2.B_2.C      5daq__1__1.A__1.B_1.C_1.D          2.A:1.A       foldseek  ...    both  pli_qcov    None      holo
..                        ...                            ...              ...            ...  ...     ...       ...     ...       ...
59          3to9__2__2.A__2.B              3to9__2__1.A__1.B          2.A:1.A       foldseek  ...    both  pli_qcov    None      holo
60          3to9__2__2.A__2.B              3to9__2__3.A__3.B          2.A:3.A       foldseek  ...    both  pli_qcov    None      holo
61          3to9__2__2.A__2.B              5wci__1__1.A__1.C          2.A:1.A       foldseek  ...    both  pli_qcov    None      holo
62          6lpf__2__1.B__1.D              6lr6__2__1.B__1.E          1.B:1.B       foldseek  ...    both  pli_qcov    None      holo
63  4lps__1__1.A__1.C_1.D_1.E      4lps__2__1.B__1.M_1.N_1.O          1.A:1.B       foldseek  ...    both  pli_qcov    None      holo

[64 rows x 9 columns]
>>> df.columns
Index(['query_system', 'target_system', 'protein_mapping', 'protein_mapper',
       'similarity', 'source', 'metric', 'mapping', 'search_db'],
      dtype='object')
```

### Ligand similarity

Ligand similarity is stored in ligand_scores sub-directory and can be queried like:

```python
>>> ligand_similarity_df = pd.read_parquet("/Users/yusuf/.local/share/plinder/2024-04/tutorial/ligand_scores")
>>> ligand_similarity_df
       query_ligand_id  target_ligand_id  tanimoto_similarity_max
0                19987             19987                      100
1                19987             12591                       67
2                19987              8937                       67
3                19987             19079                       64
4                19987             30879                       56
...                ...               ...                      ...
39995             4274             27411                       18
39996             4274             50047                       18
39997             4274              4615                       18
39998             4274             14862                       18
39999             4274               309                       18
```

`query_ligand_id` and `target_ligand_id` maps directly to `number_id_by_inchikeys` from `fingerprints/ligands_per_system.parquet` and we can map the above table to `system_id`'s with

```python
number_id_df = pd.read_parquet("~/.local/share/plinder/2024-04/tutorial/fingerprints/ligands_per_system.parquet")
number_id_to_system_id_map = number_id_df[["system_id", "number_id_by_inchikeys"]].set_index(
    "number_id_by_inchikeys").to_dict()["system_id"]
>>> ligand_similarity_df["query_ligand_id"] = ligand_similarity_df.query_ligand_id.map(number_id_to_system_id_map)
>>> ligand_similarity_df["target_ligand_id"] = ligand_similarity_df.target_ligand_id.map(number_id_to_system_id_map)
>>> ligand_similarity_df
               query_ligand_id               target_ligand_id  tanimoto_similarity_max
0    5uqd__1__1.A__1.B_1.C_1.D      5uqd__1__1.A__1.B_1.C_1.D                      100
1        5uqw__2__1.B__1.E_1.F          5uqw__2__1.B__1.E_1.F                      100
2        5uqw__2__1.B__1.E_1.F          3zjy__3__1.F__1.L_1.M                       98
3        5uqw__2__1.B__1.E_1.F              2p3q__2__2.A__2.H                       93
4        5uqw__2__1.B__1.E_1.F              2bbw__2__1.B__1.D                       93
..                         ...                            ...                      ...
320    4dpg__4__1.H__1.KA_1.LA              3pwk__1__1.B__1.J                       70
321    4dpg__4__1.H__1.KA_1.LA          3vn9__1__1.A__1.B_1.C                       70
322    4dpg__4__1.H__1.KA_1.LA      7z0w__2__1.C_1.D__1.O_1.P                       70
323    4dpg__4__1.H__1.KA_1.LA  6b09__1__1.A_1.B__1.H_1.I_1.J                       70
324    4dpg__4__1.H__1.KA_1.LA              3buz__1__1.A__1.C                       70

[325 rows x 3 columns]
```

Note: The `ligands` is the database of ligands that was used to generate the information in `fingerprints`. We have included it here for completeness sake. You can get the same information from `fingerprints`
