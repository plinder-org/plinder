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
2024-06/
|-- v2
    |-- index
    |   |-- annotation_table.parquet
    |   |-- annotation_table_nonredundant.parquet
    |-- systems
    |   |-- {two_char_code}.zip
    |-- clusters
    |   |-- cluster=communities
    |       |-- ...
    |   |-- cluster=components
    |       |-- ...
    |-- splits
    |   |-- split.parquet
    |   |-- split.yaml
    |-- linked_structures
    |   |-- {two_char_code}.zip
    |-- links
    |   |-- apo_links.parquet
    |   |-- pred_links.parquet

--------------------------------------------------------------------------------
                            miscellaneous data below
--------------------------------------------------------------------------------

    |-- dbs
    |   |-- subdbs
    |       |-- apo.csv
    |       |-- holo.csv
    |       |-- pred.csv
    |-- entries
    |   |-- {two_char_code}.zip
    |-- fingerprints
    |   |-- ligands_per_inchikey.parquet
    |   |-- ligands_per_inchikey_ecfp4.npy
    |   |-- ligands_per_system.parquet
    |-- ligand_scores
    |   |-- {hashid}.parquet
    |-- ligands
    |   |-- {hashid}.parquet
    |-- mmp
    |   |-- plinder_mmp_series.parquet
    |   |-- plinder_mms.csv.gz
    |-- scores
    |   |-- search_db=apo
    |       |-- apo.parquet
    |   |-- search_db=holo
    |       |-- {chunck_id}.parquet
    |   |-- search_db=pred
    |       |-- pred.parquet
```

:::{todo}
- Add sections for missing tables
- Fill missing elements into tables
- Add descriptions for missing files
:::

We will describe the content of the `index`, `systems`, `clusters`, `splits`, `links` and `linked_structures` directories in detail below, the rest are described in the [miscellaneous section](#miscellaneous).

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

`Mandatory`: The column has a non-empty, non-NaN value in for all PLINDER systems.
`Example`: An example non-empty, non-NaN value for the given column in a PLINDER system.


### Systems (`systems/`)
This directory contains all the systems used in the dataset. The systems are grouped into zipped subdirectories by using two penultimate characters of PDB code (`two_char_code`). The purpose of this grouping is to make loading and querying speed palatable.

Each unzipped subdirectory, contains folders named by `system_id` that contain the structure files.
```bash
|-- {two_char_code}
    |-- {system_id}
        |-- chain_mapping.json # mapping between the chains in the receptor and the chains in the system
        |-- ligand_files # mapping between the ligand in the receptor and the ligands in the system
        |-- receptor.cif
        |-- receptor.pdb
        |-- sequences.fasta
        |-- system.cif
        |-- water_mapping.json
```

### Clusters (`clusters/`)

This directory contains pre-calculated cluster labels derived from the protein and pocket similarity dataset.
The nested structure is as follows:

```bash
|-- cluster=communities
    |-- directed=False
        |-- metric={metric}
            |-- threshold={threshold}.parquet
|-- cluster=components
    |-- directed=False
        |-- metric={metric}
            |-- threshold={threshold}.parquet
    |-- directed=True
        |-- metric={metric}
            |-- threshold={threshold}.parquet
```

- `directed`: whether clustering is done with directed or undirected graph
- `metric`: the similarity metrics used for generating the clusters
- `threshold`: similarity threshold in percent done at pre-defined levels

Currently, we provide cluster labels based on undirected commnunities and both directed and undirected components.
This is performed using metrics the metrics below with pre-defined thresholds (eg. 50, 70, 95 and 100 % for `pli_unique_qcov`).

- `pli_qcov`
- `pli_unique_qcov`
- `pocket_fident`
- `pocket_fident_qcov`
- `pocket_lddt`
- `pocket_lddt_qcov`
- `pocket_qcov`
- `protein_fident_max`
- `protein_fident_qcov_max`
- `protein_fident_qcov_weighted_max`
- `protein_fident_qcov_weighted_sum`
- `protein_fident_weighted_max`
- `protein_fident_weighted_sum`
- `protein_lddt_max`
- `protein_lddt_qcov_max`
- `protein_lddt_qcov_weighted_max`
- `protein_lddt_qcov_weighted_sum`
- `protein_lddt_weighted_max`
- `protein_lddt_weighted_sum`
- `protein_qcov_weighted_sum`
- `protein_seqsim_max`
- `protein_seqsim_qcov_max`
- `protein_seqsim_qcov_weighted_max`
- `protein_seqsim_qcov_weighted_sum`
- `protein_seqsim_weighted_max`
- `protein_seqsim_weighted_sum`



### Splits (`splits/`)

This directory contains split files and the configs used to generate them.

- `split.parquet`: listing the split category for each system
- `split.yaml`: the config used to generate the split

:::{todo}
- Add missing descriptions for yaml and csv files
:::

:::{list-table} `split.parquet` columns
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
      `val` (training set) or `removed` (removed for de-leaking purposes)
*   - cluster
    - str
    - Cluster label used in sampling test set
*   - cluster_for_val_split
    - str
    - Cluster label used in sampling validation set from training set
*   - uniqueness
    - str
    - system label used to remove redundant systems from the split
*   - system_pass_validation_criteria
    - bool
    - does as system pass the validation quality criteria for test?
*   - system_pass_statistics_criteria
    - bool
    - does a system fit the statistics criteria for test?
*   - system_proper_num_ligand_chains
    - int
    - number of ligand entries in a system that are not classified as ion or artifact (i.e. "proper" ligands) 
*   - system_proper_pocket_num_residues
    - int
    - total number of pocket residues that are within 6 Å distance to a "proper" ligand(s) in a system
*   - system_proper_num_interactions
    - int
    - total number of PLI interactions to a "proper" ligand(s) in a system
*   - system_proper_ligand_max_molecular_weight
    - float
    - maximum molecular weight of the "proper" ligand(s) in a system
*   - system_has_binding_affinity
    - bool
    - does the system have a ligand with an annotated binding affinity?
*   - system_has_apo_or_pred
    - bool
    - does the system have either `apo` or `pred` structure linked?


:::

```
>>> df.head()
               system_id            uniqueness  split cluster  ... system_proper_num_interactions  system_proper_ligand_max_molecular_weight  system_has_binding_affinity  system_has_apo_or_pred
0  101m__1__1.A__1.C_1.D  101m__A__C_D_c188899  train     c14  ...                             20                                 616.177293                        False                   False
1      102m__1__1.A__1.C    102m__A__C_c237197  train     c14  ...                             20                                 616.177293                        False                    True
2  103m__1__1.A__1.C_1.D  103m__A__C_D_c252759  train     c14  ...                             16                                 616.177293                        False                   False
3  104m__1__1.A__1.C_1.D  104m__A__C_D_c274687  train     c14  ...                             21                                 616.177293                        False                   False
4  105m__1__1.A__1.C_1.D  105m__A__C_D_c221688  train     c14  ...                             20                                 616.177293                        False                   False
```
:::{todo}
Index(['system_id', 'uniqueness', 'split', 'cluster', 'cluster_for_val_split',
       'system_pass_validation_criteria', 'system_pass_statistics_criteria',
       'system_proper_num_ligand_chains', 'system_proper_pocket_num_residues',
       'system_proper_num_interactions',
       'system_proper_ligand_max_molecular_weight',
       'system_has_binding_affinity', 'system_has_apo_or_pred'],
:::

### Linked structures (`linked_structures/`)
This directory contains the linked apo and predicted structures for PLINDER systems. These structures are intended to be used for augmenting the PLINDER dataset, eg. for flexible docking or pocket prediction purposes.
The files are grouped into zipped subdirectories by using `two_char_code` of the system.
Each unzipped subdirectory contains `pred` and `apo` subfolders that in turn contain folders named by `system_id`.
Inside each `apo/{system_id}` and `pred/{system_id}` folder is another directory containing a superposed system: `{source_id}_{chain_id}/superposed.cif`, where `{source_id}` and `{chain_id}` for apo systems is `pdb_id` with a source chain identifier, and for predicted structures, `{source_id}` is `uniprot_id` used in AF2DB with a chain identifier set to `A`.


### Linked systems (`links/`)
This directory contains parquet files linking PLINDER systems to their apo and predicted structures in `linked_structures/`.


```
>>> df.head()
  reference_system_id      id  pocket_fident  pocket_lddt  protein_fident_qcov_weighted_sum  ...  fraction_model_proteins_mapped      lddt   bb_lddt  per_chain_lddt_ave per_chain_bb_lddt_ave
0   6pl9__1__1.A__1.C  2vb1_A          100.0         86.0                             100.0  ...                             1.0  0.903772  0.968844            0.890822              0.959674
1   6ahh__1__1.A__1.G  2vb1_A          100.0         98.0                             100.0  ...                             1.0  0.894349  0.962846            0.883217              0.954721
2   5b59__1__1.A__1.B  2vb1_A          100.0         91.0                             100.0  ...                             1.0  0.903266  0.962318            0.890656              0.955258
3   3ato__1__1.A__1.B  2vb1_A          100.0         99.0                             100.0  ...                             1.0  0.890530  0.954696            0.879496              0.946326
4   6mx9__1__1.A__1.K  2vb1_A          100.0         98.0                             100.0  ...                             1.0  0.904116  0.964309            0.892434              0.955853
```

:::{todo}
Index(['reference_system_id', 'id', 'pocket_fident', 'pocket_lddt',
       'protein_fident_qcov_weighted_sum', 'protein_fident_weighted_sum',
       'protein_lddt_weighted_sum', 'target_id', 'sort_score', 'receptor_file',
       'ligand_files', 'num_reference_ligands', 'num_model_ligands',
       'num_reference_proteins', 'num_model_proteins',
       'fraction_reference_ligands_mapped', 'fraction_model_ligands_mapped',
       'lddt_pli_ave', 'lddt_pli_wave', 'lddt_pli_amd_ave',
       'lddt_pli_amd_wave', 'scrmsd_ave', 'scrmsd_wave', 'lddt_lp_ave',
       'lddt_lp_wave', 'posebusters_mol_pred_loaded',
       'posebusters_mol_cond_loaded', 'posebusters_sanitization',
       'posebusters_all_atoms_connected', 'posebusters_bond_lengths',
       'posebusters_bond_angles', 'posebusters_internal_steric_clash',
       'posebusters_aromatic_ring_flatness',
       'posebusters_double_bond_flatness', 'posebusters_internal_energy',
       'posebusters_protein-ligand_maximum_distance',
       'posebusters_minimum_distance_to_protein',
       'posebusters_minimum_distance_to_organic_cofactors',
       'posebusters_minimum_distance_to_inorganic_cofactors',
       'posebusters_minimum_distance_to_waters',
       'posebusters_volume_overlap_with_protein',
       'posebusters_volume_overlap_with_organic_cofactors',
       'posebusters_volume_overlap_with_inorganic_cofactors',
       'posebusters_volume_overlap_with_waters',
       'fraction_reference_proteins_mapped', 'fraction_model_proteins_mapped',
       'lddt', 'bb_lddt', 'per_chain_lddt_ave', 'per_chain_bb_lddt_ave'],
      dtype='object')
:::

### Miscellaneous

Here we briefly describe subdirectories and their files that are not part of the main dataset but are used in the dataset processing pipeline.
These files should be considered intermediate products and are not intended to be used directly, only for development purposes.

#### Database processed files (`dbs/`)
This directory contains the intermediate files of PDB structures that were successfully processed and scored by foldseek and mmseqs pipeline.
It is used in splitting to make sure that only successfully computed systems are used for splitting.
```bash
|-- subdbs
|   |-- apo.csv
|   |-- holo.csv
|   |-- pred.csv
```
Each file is a CSV with a single column: `pdb_id`.

#### Raw annotations (`entries/`)
This directory contains intermediate raw annotation files prior to consolidation. The files are grouped into zipped subdirectories by using `two_char_code`.
Each subdirectory, contains `{pdb_id}.json` files with raw annotations for every system found in given `pdb_id`.


#### Small molecule fingerprints (`fingerprints/`)

Tables that contains all the ligand fingerprints used in calculating ligand similarity stored in `ligand_scores`.

- `ligands_per_inchikey_ecfp4.npy`: `numpy` array of all-vs-all ECFP4 similarity.
- `ligands_per_system.parquet`: table linking PLINDER systems to their ligands, including ligand ID, SMILES, InChIKey, etc.
- `ligands_per_inchikey.parquet`: subset of `ligands_per_system.parquet` with reduced number of columns.


#### Small molecule data (`ligands/`)
Ligand data expanded from entries for computing similarity, saved in distributed files `{hashid}.parquet`.
Eg.

```
  pdb_id              system_id                      ligand_rdkit_canonical_smiles ligand_ccd_code                   ligand_id                    inchikeys
0   7o00  7o00__1__1.A_1.B__1.D  CC(=O)N[C@H]1CO[C@H](CO)[C@@H](OC2O[C@H](CO)[C...         HSR-HSR  7o00__1__1.A_1.B__1.D__1.D  JHPFQHGUNGJQIZ-BQBDUENHSA-N
1   7o00  7o00__1__1.A_1.B__1.E  CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O             HSR  7o00__1__1.A_1.B__1.E__1.E  OVRNDRQMDRJTHS-FMDGEEDCSA-N
2   7o04      7o04__1__1.A__1.G                        CNCc1cc([N+](=O)[O-])ccc1Cl             4AV      7o04__1__1.A__1.G__1.G  YRTNCUPHKWUHMQ-UHFFFAOYSA-N
3   7o08      7o08__1__1.A__1.C  CC1(C)CCN(Cc2ccc(NCC3(O)CCN(c4cc(NCc5ccccc5)nc...             UXE      7o08__1__1.A__1.C__1.C  GTLDMCHZRAFXCB-UHFFFAOYSA-N
4   7o09      7o09__1__1.A__1.C  CC1(C)CCN(Cc2ccc(N3CCOC4(CCN(c5cc(NCc6ccccc6)n...             UXK      7o09__1__1.A__1.C__1.C  RJEWLHZZXYDBNT-UHFFFAOYSA-N
```

#### Small molecule similarity scores (`ligand_scores/`)

Tables that contains all the ligand similarity scores used in calculating the similarity between two ligands, saved in distributed files `{hashid}.parquet`.
Eg.
```
   query_ligand_id  target_ligand_id  tanimoto_similarity_max
0            35300              6943                      100
1            35300             35300                      100
2            35300             13911                       94
3            35300             44243                       90
4            35300             24003                       90
```

#### Small molecule matched molecular pairs (`mmp/`)

Files that contains all the ligand matched molecular pairs (MMP) and matched molecular series (MMS).
- `plinder_mmp_series.parquet`: matched molecular series (MMS) linked to PLINDER systems,
- `plinder_mms.csv.gz`: compressed [mmpdb](https://github.com/rdkit/mmpdb) index file containing the matched molecular pairs (MMP) of all ligands in PLINDER annotation table.

#### Protein similarity dataset (`scores/`)

TODO

Tables that contains all the protein similarity scores used in calculating the similarity between two proteins.

```bash
|-- search_db=apo
|   |-- apo.parquet
|-- search_db=holo
|   |-- {chunck_id}.parquet
|-- search_db=pred
|   |-- pred.parquet
```
