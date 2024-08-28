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
    |-- index # Consolidated tabular annotations
    |   |-- annotation_table.parquet
    |   |-- annotation_table_nonredundant.parquet
    |-- systems  # Structure files for all systems (split by `two_char_code` and zipped)
    |   |-- {two_char_code}.zip
    |-- clusters # Pre-calculated cluster labels derived from the protein similarity dataset
    |   |-- cluster=communities
    |       |-- ...
    |   |-- cluster=components
    |       |-- ...
    |-- splits # Split files and the configs used to generate them (if available)
    |   |-- split.parquet
    |   |-- split.yaml
    |-- linked_structures # Apo and predicted structures linked to their holo systems
    |   |-- {two_char_code}.zip
    |-- links # Apo and predicted structures similarity to their holo structures
    |   |-- apo_links.parquet
    |   |-- pred_links.parquet

--------------------------------------------------------------------------------
                            miscellaneous data below
--------------------------------------------------------------------------------

    |-- dbs # TSVs containing the raw files and IDs in the foldseek and mmseqs sub-databases
    |   |-- subdbs
    |       |-- apo.csv
    |       |-- holo.csv
    |       |-- pred.csv
    |-- entries # Raw annotations prior to consolidation (split by `two_char_code` and zipped)
    |   |-- {two_char_code}.zip
    |-- fingerprints # Index mapping files for the ligand similarity dataset
    |   |-- ligands_per_inchikey.parquet
    |   |-- ligands_per_inchikey_ecfp4.npy
    |   |-- ligands_per_system.parquet
    |-- ligand_scores # Ligand similarity parquet dataset
    |   |-- {hashid}.parquet
    |-- ligands # Ligand data expanded from entries for computing similarity
    |   |-- {hashid}.parquet
    |-- mmp # Ligand matched molecular pairs (MMP) and series (MMS) data
    |   |-- plinder_mmp_series.parquet
    |   |-- plinder_mms.csv.gz
    |-- scores # Protein similarity parquet dataset
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
        |-- chain_mapping.json # Mapping between the chains in the receptor and the chains in the system
        |-- ligand_files # Mapping between the ligand in the receptor and the ligands in the system
        |-- receptor.cif  # Receptor mmcif file
        |-- receptor.pdb # Receptor pdb file
        |-- sequences.fasta # Receptor sequence fasta
        |-- system.cif # System mmcif file
        |-- water_mapping.json # Receptor binding site water map json file
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

- `cluster`: whether clusters generated from community detection algorithm or via disconnected component
- `directed`: whether clustering is done with directed or undirected graph
- `metric`: the similarity metrics used for generating the clusters
- `threshold`: similarity threshold in percent done at pre-defined levels

Currently, we provide cluster labels based on undirected commnunities and both directed and undirected components.
This is performed using metrics the metrics below with pre-defined thresholds (eg. 50, 70, 95 and 100 %).

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

:::{todo}
```
>>> df.head()
               system_id            uniqueness  split cluster  ... system_proper_num_interactions  system_proper_ligand_max_molecular_weight  system_has_binding_affinity  system_has_apo_or_pred
0  101m__1__1.A__1.C_1.D  101m__A__C_D_c188899  train     c14  ...                             20                                 616.177293                        False                   False
1      102m__1__1.A__1.C    102m__A__C_c237197  train     c14  ...                             20                                 616.177293                        False                    True
2  103m__1__1.A__1.C_1.D  103m__A__C_D_c252759  train     c14  ...                             16                                 616.177293                        False                   False
3  104m__1__1.A__1.C_1.D  104m__A__C_D_c274687  train     c14  ...                             21                                 616.177293                        False                   False
4  105m__1__1.A__1.C_1.D  105m__A__C_D_c221688  train     c14  ...                             20                                 616.177293                        False                   False
```
:::

### Linked structures (`linked_structures/`)
This directory contains the linked apo and predicted structures for PLINDER systems. These structures are intended to be used for augmenting the PLINDER dataset, eg. for flexible docking or pocket prediction purposes.
The files are grouped into zipped subdirectories by using `two_char_code` of the system.
Each unzipped subdirectory contains `pred` and `apo` subfolders that in turn contain folders named by `system_id`.
Inside each `apo/{system_id}` and `pred/{system_id}` folder is another directory containing a superposed system: `{source_id}_{chain_id}/superposed.cif`, where `{source_id}` and `{chain_id}` for apo systems is `pdb_id` with a source chain identifier, and for predicted structures, `{source_id}` is `uniprot_id` used in AF2DB with a chain identifier set to `A`.


### Linked systems (`links/`)
This directory contains parquet files linking PLINDER systems to their apo and predicted structures in `linked_structures/`.


:::{list-table} `{apo|pred}_links.parquet` columns
:widths: 10 5 30
:header-rows: 1

*   - Name
    - Type
    - Description
*   - reference_system_id
    - str
    - The PLINDER system ID
*   - id
    - str
    - The PDB or AF2DB (for `apo` and `pred`, respectively) `{source_id}_{chain_id}` tag
*   - pocket_fident
    - float
    - sequence identity for pocket residues
*   - pocket_lddt
    - float
    - Local Distance Difference Test (lDDT) score for the pocket residues
*   - protein_fident_qcov_weighted_sum
    - float
    - Sum of fident * qcov for all templates, weighted by the number of residues in the template
*   - protein_fident_weighted_sum
    - float
    - Sum of fident for all templates, weighted by the number of residues in the template
*   - protein_lddt_weighted_sum
    - float
    - Sum of lDDT for all residues, weighted by the number of residues in the template
*   - target_id
    - str
    - apo or pred stucture `{source_id}` tag
*   - sort_score
    - float
    - score used to sort linked structures
*   - receptor_file
    - str
    - intermediate aligned linked receptor file path
*   - ligand_files
    - str
    - intermediate file path for ligands used in calculations
*   - num_reference_ligands
    - int
    - number of ligands in reference structure
*   - num_model_ligands
    - int
    - number of ligands in model structure
*   - num_reference_proteins
    - int
    - number of protein chains in reference structure
*   - num_model_proteins
    - int
    - number of protein chains in model structure
*   - fraction_reference_ligands_mapped
    - float
    - Fraction of reference ligands that were successfully mapped to model ligands
*   - fraction_model_ligands_mapped
    - float
    - Fraction of model ligands that were successfully mapped to reference ligands
*   - lddt_pli_ave
    - float
    - Average lDDT score for protein-ligand interactions
*   - lddt_pli_wave
    - float
    - Weighted average lDDT score for protein-ligand interactions
*   - lddt_pli_amd_ave
    - float
    - Average lDDT score for protein-ligand interactions, considering only atoms matched during docking
*   - lddt_pli_amd_wave
    - float
    - Weighted average lDDT score for protein-ligand interactions, considering only atoms matched during docking
*   - scrmsd_ave
    - float
    - Average symmetry-corrected RMSD between reference and model ligands
*   - scrmsd_wave
    - float
    - Weighted average symmetry-corrected RMSD between reference and model ligands
*   - lddt_lp_ave
    - float
    - Average lDDT score for ligand poses
*   - lddt_lp_wave
    - float
    - Weighted average lDDT score for ligand poses
*   - posebusters_mol_pred_loaded
    - bool
    - PoseBusters metric: boolean indicator of whether the predicted ligand could be loaded
*   - posebusters_mol_cond_loaded
    - bool
    - PoseBusters metric: boolean indicator of whether the conditional ligand could be loaded
*   - posebusters_sanitization
    - bool
    - PoseBusters metric: boolean indicator of whether the ligand could be sanitized
*   - posebusters_all_atoms_connected
    - bool
    - PoseBusters metric: boolean indicator of whether all atoms in the ligand are connected
*   - posebusters_bond_lengths
    - bool
    - PoseBusters metric: boolean indicator of whether all bond lengths in the ligand are within 4 standard deviations of the mean
*   - posebusters_bond_angles
    - bool
    - PoseBusters metric: boolean indicator of whether all bond angles in the ligand are within 4 standard deviations of the mean
*   - posebusters_internal_steric_clash
    - bool
    - PoseBusters metric: boolean indicator of whether there are no internal steric clashes in the ligand
*   - posebusters_aromatic_ring_flatness
    - bool
    - PoseBusters metric: boolean indicator of whether all aromatic rings in the ligand are flat
*   - posebusters_double_bond_flatness
    - bool
    - PoseBusters metric: boolean indicator of whether all double bonds in the ligand are flat
*   - posebusters_internal_energy
    - bool
    - PoseBusters metric: boolean indicator of whether the internal energy of the ligand is below 0 kcal/mol
*   - posebusters_protein-ligand_maximum_distance
    - bool
    - PoseBusters metric: boolean indicator of whether the maximum distance between the ligand and the protein is less than 5 Angstrom
*   - posebusters_minimum_distance_to_protein
    - bool
    - PoseBusters metric: boolean indicator of whether the minimum distance between the ligand and the protein is greater than 1.5 Angstrom
*   - posebusters_minimum_distance_to_organic_cofactors
    - float
    - PoseBusters metric: Minimum distance between the ligand and any organic cofactor
*   - posebusters_minimum_distance_to_inorganic_cofactors
    - bool
    - PoseBusters metric: Minimum distance between the ligand and any inorganic cofactor
*   - posebusters_minimum_distance_to_waters
    - float
    - PoseBusters metric: Minimum distance between the ligand and any water molecule
*   - posebusters_volume_overlap_with_protein
    - float
    - PoseBusters metric: Fraction of ligand volume that overlaps with the protein
*   - posebusters_volume_overlap_with_organic_cofactors
    - bool
    - PoseBusters metric: boolean indicator of whether the share of ligand volume that intersects with the organic cofactor is less than 7.5%. The volumes are defined by the van der Waals radii around the heavy atoms scaled by 0.8.
*   - posebusters_volume_overlap_with_inorganic_cofactors
    - bool
    - PoseBusters metric: boolean indicator of whether the share of ligand volume that intersects with the inorganic cofactor is less than 7.5%. The volumes are defined by the van der Waals radii around the heavy atoms scaled by 0.8.
*   - posebusters_volume_overlap_with_waters
    - bool
    - PoseBusters metric: boolean indicator of whether the share of ligand volume that intersects with the linked system waters is less than 7.5%. The volumes are defined by the van der Waals radii around the heavy atoms scaled by 0.8.
*   - fraction_reference_proteins_mapped
    - float
    - Fraction of reference protein chains with corresponding model chains
*   - fraction_model_proteins_mapped
    - float
    - Fraction of model protein chains mapped to corresponding reference chains
*   - lddt
    - float
    - Global lDDT score calculated over all atoms in the structure
*   - bb_lddt
    - float
    - Global lDDT score calculated over backbone atoms (N, CA, C, O) in the structure
*   - per_chain_lddt_ave
    - float
    - Average per-chain lDDT score calculated over all atoms
*   - per_chain_bb_lddt_ave
    - float
    - Average per-chain lDDT score calculated over backbone atoms (N, CA, C, O)
:::


### Miscellaneous

Here we briefly describe subdirectories and their files that are not part of the main dataset but are used in the dataset processing pipeline.
These files should be considered intermediate products and are not intended to be used directly, only for development purposes.

#### Database processed files (`dbs/`)
This directory contains the intermediate files of PDB structures that were successfully processed and scored by Foldseek and MMseqs2 pipeline.
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

Tables that contains all the protein or pocket similarity scores used in calculating the similarity between two systems.

```bash
|-- search_db=apo
|   |-- apo.parquet
|-- search_db=holo
|   |-- {chunck_id}.parquet
|-- search_db=pred
|   |-- pred.parquet
```
All the parquet files have the save columns in the header.
E.g
```
                    query_system target_system protein_mapping protein_mapper  ...    source                            metric  mapping search_db
1070886    1b5d__1__1.A_1.B__1.D        1b49_A         1.A:0.A       foldseek  ...    mmseqs         protein_qcov_weighted_max  1.A:0.A       apo
1070887    1b5d__1__1.A_1.B__1.D        1b49_A         1.A:0.A       foldseek  ...    mmseqs                  protein_qcov_max  1.A:0.A       apo
1070888    1b5d__1__1.A_1.B__1.D        1b49_A         1.A:0.A       foldseek  ...      both       protein_fident_weighted_max  1.A:0.A       apo
1070889    1b5d__1__1.A_1.B__1.D        1b49_A         1.A:0.A       foldseek  ...      both                protein_fident_max  1.A:0.A       apo
1070890    1b5d__1__1.A_1.B__1.D        1b49_A         1.A:0.A       foldseek  ...    mmseqs  protein_fident_qcov_weighted_max  1.A:0.A       apo
...                          ...           ...             ...            ...  ...       ...                               ...      ...       ...
213471528      7eek__1__1.A__1.I        1uor_A         1.A:0.A       foldseek  ...  foldseek    protein_lddt_qcov_weighted_max  1.A:0.A       apo
213471529      7eek__1__1.A__1.I        1uor_A         1.A:0.A       foldseek  ...  foldseek             protein_lddt_qcov_max  1.A:0.A       apo
213471536      7eek__1__1.A__1.I        1uor_A         1.A:0.A       foldseek  ...  foldseek                       pocket_lddt     None       apo
213471540      7eek__1__1.A__1.I        6zl1_A         1.A:0.A       foldseek  ...  foldseek                       pocket_lddt     None       apo
213471541      7eek__1__1.A__1.I        6zl1_B         1.A:0.B       foldseek  ...  foldseek                       pocket_lddt     None       apo
```
:::{list-table} `apo.parquet` columns
:widths: 10 5 30
:header-rows: 1

*   - Name
    - Type
    - Description
*   - query_system
    - str
    - The PLINDER system ID of query system
*   - target_system
    - str
    - The PLINDER system ID of target system
*   - protein_mapping
    - str
    - Chain mapping between query system and target system
*   - protein_mapper
    - str
    - Alignment method used for mapping.
*   - similarity
    - int
    - Similarity metric of interest
*   - source
    - str
    - Source of similarity metric. It could either be `foldseek`, `mmseqs` or `both`
*   - metric
    - str
    - Similarity metric of interest
*   - mapping
    - str
    - ??
*   - search_db
    - str
    - Search database type. Could be `apo`, `holo` or `pred`
