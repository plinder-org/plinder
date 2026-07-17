---
sd_hide_title: true
---

# Dataset

## Dataset reference

### Directory structure

```bash
2024-06/
|-- v2
    |-- index # Consolidated tabular annotations
    |   |-- annotation_table.parquet
    |   |-- annotation_table_nonredundant.parquet
    |   |-- entry_chains.parquet
    |   |-- entry_sources.parquet
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
    |
--------------------------------------------------------------------------------
                            miscellaneous data below
--------------------------------------------------------------------------------
    |
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
    |-- alignments # V3 mapped Foldseek/MMseqs results for score reconstruction
    |   |-- search_db={holo,apo,pred}
    |       |-- alignment_type={foldseek,mmseqs}
    |           |-- shard={two_char_code}.parquet
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

We will describe the content of the `index`, `systems`, `clusters`, `splits`, `links` and `linked_structures` directories in detail below, the rest are described in the [miscellaneous section](#miscellaneous-target).

(annotation-table-target)=

### Annotation tables (`index/`)

Tables that lists all systems along with their annotations.

- `annotation_table.parquet`: Lists all systems and their annotations.
- `annotation_table_nonredundant.parquet`: Subset of systems without redundant systems.
- `entry_chains.parquet`: One row per protein or nucleic-acid receptor chain, including its normalized `chain_receptor_type`. Protein rows also carry the entity, holo partition flag, and UniProt mappings used to construct the Foldseek/MMseqs sub-databases; nucleic-acid rows are never submitted to similarity search.
- `entry_sources.parquet`: One row per PDB entry recording the exact source mmCIF major and minor revision used during ingest. This is normalized entry metadata and is not repeated on ligand rows.

V3 cluster columns are merged only after local scoring and clustering finish.
Receptor-only metrics have system-level labels. Clustered ligand-sensitive metrics,
including pocket, PLI, and `sucos_shape_pocket_qcov`, have both a system projection
used for split safety and a row-specific ligand label. Ligand columns include an
explicit `__ligand__` marker, for example
`sucos_shape_pocket_qcov__50__ligand__strong__component`; the corresponding
projected system column is
`sucos_shape_pocket_qcov__50__strong__component`. The gated `shape`, `color`, and
raw `sucos_shape` values remain diagnostic scores and are not clustered directly.

V3 does not distribute the materialized `scores/` dataset. Instead, it distributes
mapped Foldseek and MMseqs results in `alignments/`, ordered by query and target and
compressed as deterministic two-character Parquet shards. The shard for a PDB ID is
`pdb_id[-3:-1]`, matching the PDB archive convention. This keeps the number of files
manageable while allowing Parquet predicate pushdown to read only row groups relevant
to requested query and target entries. V2 retains the materialized `scores/` layout.

:::{include} table.html
:::

`Mandatory`: The column has a non-empty, non-NaN value in for all PLINDER systems.
`Example`: An example non-empty, non-NaN value for the given column in a PLINDER system.

### Canonical ligand structures

PLINDER stores one canonical asymmetric-unit SDF for each ligand chain. Biological-assembly copies are not stored because they contain the same conformation under different rigid-body transforms.

```bash
|-- {two_char_code}
    |-- {pdb_id}
        |-- ligand_files
            |-- {asym_id}.sdf
```

System and receptor mmCIF files are reconstructed on demand from the deposited PDB
mmCIF and the system selection metadata in the annotation parquet. During V3 ingest,
PLINDER records the exact source structure-model major/minor revision in
`index/entry_sources.parquet`. `PlinderSystem` fetches that exact compressed revision
from the [wwPDB versioned archive](https://www.wwpdb.org/ftp/pdb-versioned-ftp-site)
and caches it under
`<plinder_dir>/source_mmcifs/{two_char_code}/{pdb_id}_v{major}-{minor}.cif.gz`.
The cache is therefore reproducible for a particular PLINDER release rather than a
copy of whichever revision happens to be current when it is requested. An explicitly
supplied `source_mmcif` still takes precedence.

To prepare a subset on an online node before moving the cache to an offline node:

```python
from plinder.core import download_pdb_mmcifs

download_pdb_mmcifs(["2y4i", "1a3b"])
```

Only the unique PDB IDs passed to this function are fetched; the full PDB is never
downloaded. The small release revision manifest is resolved first. Afterward, set
`PLINDER_OFFLINE=true` (or `PLINDER_OFFLINE_MODE=true`). The cached manifest and
source files are reused without network access; requesting an uncached entry raises
an error that includes its expected cache path. Callers explicitly choose whether
reconstructed outputs include interacting or all waters, other biological-assembly
chains, and which output files to write.

V3 similarity scores can likewise be prepared on an online node and reconstructed
later in offline mode:

```python
from plinder.core.scores import (
    prefetch_similarity_alignments,
    reconstruct_similarity_scores,
)

systems = ["2y4i__1__1.B__1.E_1.F", "6cex__1__1.D__1.M"]
prefetch_similarity_alignments(systems)
scores = reconstruct_similarity_scores(systems[:1], systems[1:])
```

Prefetching derives the required two-character shards from the query systems and
downloads each available Foldseek/MMseqs shard once; it never downloads the complete
alignment dataset. If one search backend has no shard, reconstruction uses the other;
it fails only when neither backend is available. Reconstruction then filters those
Parquets by query and target PDB ID, loads only the requested annotation rows, and
emits directed ligand-level scores.
Canonical SDFs are resolved lazily only for ligand pairs whose `pocket_qcov` is
positive. With `PLINDER_OFFLINE=true` or `PLINDER_OFFLINE_MODE=true`, the same calls
use the normal PLINDER cache and report a missing shard instead of accessing the
network. The prefetch function covers alignment shards only; run reconstruction once
on the online node as well if the annotation files and any positive-pocket ligand
archives have not already been cached.

```python
from pathlib import Path

import pandas as pd

from plinder.data.utils.annotations.save_utils import (
    SystemReconstructionOptions,
    SystemReconstructionOutputs,
    save_reconstructed_system,
)

annotation = pd.read_parquet("annotation_table.parquet")
row = annotation.query("system_id == '2y4i__1__1.B__1.E_1.F'").iloc[0]
save_reconstructed_system(
    "pdb_00002y4i_xyz-enrich.cif.gz",
    row,
    outputs=SystemReconstructionOutputs(
        system_cif=Path("system.cif"),
        receptor_cif=Path("receptor.cif"),
        # Leave an output as None when it should not be written.
        sequences_fasta=None,
    ),
    options=SystemReconstructionOptions(
        system_waters="interacting",  # "none", "interacting", or "all"
        receptor_waters="none",
        system_include_other_protein_chains=False,
        system_include_other_ligand_chains=False,
    ),
)
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

- `cluster`: the cluster algorithm used
  - `communities`: clusters derived from community detection algorithm
  - `components`: clusters derived from disconnected component of similarity graph
- `directed`: type of graph used for cluster input
  - `False`: undirected
  - `True`: directed
- `metric`: the similarity metrics used for generating the clusters
  - `pli_qcov`: Protein-ligand interaction similarity between aligned ligand-binding region (pocket) residues of two systems.
  - `pli_unique_qcov`: Protein-ligand interaction similarity between aligned pocket residues of two systems, taking only unique interaction type into consideration.
  - `pocket_fident`: Pocket region sequence identity of the ligand-binding (pocket) region of a system to a (possibly non-pocket) region of another system.
  - `pocket_fident_qcov`: Sequence identity between ligand binding region (pocket) of two systems.
  - `pocket_lddt`: Structural similarity between ligand-binding region (pocket) of a system to any region (possibly non-pocket) of another system.
  - `pocket_lddt_qcov`: Structural similarity between ligand-binding region (pocket) two systems.
  - `pocket_qcov`: Query coverage between ligand-binding region of two systems.
  - `protein_fident_max`: Local sequence identity between components of two systems, aggregated by max score across all pairs of protein chains or ligand chains.
  - `protein_fident_qcov_max`: Global protein sequence identity between components of two systems multiplied by query system coverage, aggregated by max score across all pairs of protein chains or ligand chains.
  - `protein_fident_qcov_weighted_max`: Global protein sequence identity between components of two systems, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
  - `protein_fident_qcov_weighted_sum`: Global protein sequence identity between components of two systems, aggregated by length-weighted max score across all pairs of protein chains or ligand chains.
  - `protein_fident_weighted_max`: Local sequence identity between components of two systems, aggregated by length-weighted max score across all pairs of protein chains or ligand chains.
  - `protein_fident_weighted_sum`: Local sequence identity between components of two systems, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
  - `protein_lddt_max`: Local structural similarity between chains of two systems, aggregated by max score across all pairs of protein chains or ligand chains.
  - `protein_lddt_qcov_max`: Global protein structural similarity multiplied by query system coverage, aggregated by max score across all pairs of protein chains or ligand chains.
  - `protein_lddt_qcov_weighted_max`: Global protein structural similarity multiplied by query system coverage, aggregated by length-weighted max score across all pairs of protein chains or ligand chains.
  - `protein_lddt_qcov_weighted_sum`: Global protein structural similarity multiplied by query system coverage, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
  - `protein_lddt_weighted_max`: Local structural similarity between chains of two systems, aggregated by length-weighted max score across all pairs of protein chains or ligand chains.
  - `protein_lddt_weighted_sum`: Local structural similarity between chains of two systems, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
  - `protein_qcov_weighted_sum`: Global protein query coverage, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
  - `protein_seqsim_max`: Global protein sequence similarity between components of two systems, aggregated by max score across all pairs of protein chains or ligand chains.
  - `protein_seqsim_qcov_max`: Global protein sequence similarity between components of two systems multiplied by query system coverage, aggregated by max score across all pairs of protein chains or ligand chains.
  - `protein_seqsim_qcov_weighted_max`: Global protein sequence similarity between components of two systems multiplied by query system coverage, aggregated by length-weighted max score across all pairs of protein chains or ligand chains.
  - `protein_seqsim_qcov_weighted_sum`: Global protein sequence similarity between components of two systems multiplied by query system coverage, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
  - `protein_seqsim_weighted_max`: Global protein sequence similarity between components of two systems, aggregated by length-weighted max score across all pairs of protein chains or ligand chains.
  - `protein_seqsim_weighted_sum`: Global protein sequence similarity between components of two systems, aggregated by length-weighted sum of scores across mapped protein or ligand chains.
- `threshold`: similarity threshold in percent.
  - ...
  - `50`
  - `70`
  - `95`
  - `100`

### Splits (`splits/`)

This directory contains split files and the configs used to generate them.

- `split.parquet`: listing the split category for each system
- `split.yaml`: the config used to generate the split

:::{list-table} `split.parquet`
:widths: 10 5 30
:header-rows: 1

- - Name
  - Type
  - Description
- - system_id
  - str
  - The PLINDER system ID
- - split
  - str
  - Split category: either `train` (training set), `test` (test set),`val` (training set) or `removed` (removed for de-leaking purposes)
- - cluster
  - str
  - Cluster label used in sampling test set
- - cluster_for_val_split
  - str
  - Cluster label used in sampling validation set.
- - uniqueness
  - str
  - system label used to remove redundant systems from the split
- - system_pass_validation_criteria
  - bool
  - does as system pass the crystal quality for test?
- - system_pass_statistics_criteria
  - bool
  - does a system fit the statistics criteria for test?
- - system_proper_num_ligand_chains
  - int
  - number of ligand entries in a system that are not classified as ion or artifact (i.e "proper" ligands)
- - system_proper_pocket_num_residues
  - int
  - total number of pocket residues that are within 6 Å distance to a "proper" ligand(s) in a system
- - system_proper_num_interactions
  - int
  - total number of PLI interactions to a "proper" ligand(s) in a system
- - system_proper_ligand_max_molecular_weight
  - float
  - maximum molecular weight of the "proper" ligand(s) in a system
- - system_has_binding_affinity
  - bool
  - does the system have a ligand with an annotated binding affinity?
- - system_has_apo_or_pred
  - bool
  - does the system have either `apo` or `pred` structure linked?
:::

The content of `split.yaml` is described below:

```bash
split:
  graph_configs: # Similarity graph configuration
  - metric: pli_unique_qcov # Metric used to generate the base graph from which all partitioning is done.
    threshold: 30 # Threshold used to generate the base graph from which all partitioning is done.
    depth: 1 # Depth at which the neighbors are defined.
  - metric: protein_seqsim_weighted_sum # Same as above
    threshold: 30 # Same as above
    depth: 1 # Same as above
  mms_unique_quality_count: 3 # How many unique congeneric IDs passing quality to consider as MMS
  ligand_cluster_metric: Tanimoto_similarity_max # which metric to use for ligand clusters (these are added to test from removed if they are different from train/val and corresponding leaked systems are removed from train/val)
  ligand_cluster_threshold: 50 # Which threshold to use for ligand clusters

  ligand_cluster_cluster: components # Which cluster to use for ligand clusters
  test_cluster_cluster: communities # What kind of cluster to use for sampling test
  test_cluster_metric: pli_unique_qcov # Metric to use for sampling representatives from each test cluster
  test_cluster_threshold: 50  # Threshold to use for sampling representatives from each test cluster
  test_cluster_directed: false # Directed to use for sampling representatives from each test cluster
  num_test_representatives: 2 # Max number of representatives from each test cluster
  num_per_entry_pdb_id_and_unique_ccd_codes: 1 # Max number of systems to choose per entry pdb id and unique ccd codes
  min_test_cluster_size: 5 # Test should not be singletons
  min_test_leakage_count: 30  # Test should not be too unique
  max_test_leakage_count: 1000 # Test should not be in too big communities or cause too many train cases to be removed
  max_removed_fraction: 0.2 # Maximum fraction of systems that can be removed due to test set selection
  num_test: 1000 # test set size
  val_cluster_cluster: components # What kind of cluster to use for sampling val
  val_cluster_metric: pocket_qcov # Metric to use for splitting train and val
  val_cluster_threshold: 50  # Threshold to use for splitting train and val
  val_cluster_directed: false # Directed to use for splitting train and val
  num_val_representatives: 3 # Max number of representatives from each val cluster
  min_val_cluster_size: 30  # Val should not be singletons
  num_val: 1000  # Val set size
  min_max_pli: # Test/val should not have too few or too many interactions
  - 3
  - 50
  min_max_pocket: # Test/val should not have too few or too many pocket residues
  - 5
  - 100
  min_max_ligand: # Test/val should not have too small or too large ligands
  - 200
  - 800
  test_additional_criteria: # Priority columns to use for scoring systems with a weight attached to each column
  - - system_pass_validation_criteria # Indicator of whether a system is passing validation criteria
    - ==
    - 'True'
  - - system_pass_statistics_criteria # Indicator of whether a system is passing statistic criteria
    - ==
    - 'True'
  - - biounit_num_ligands # Number of ligands in the biounit.
    - <=
    - 20
  priority_columns:
    system_ligand_has_cofactor: -40.0
    leakage_count: -1.0
```

### Linked structures (`linked_structures/`)

This directory contains the linked apo and predicted structures for PLINDER systems. These structures are intended to be used for augmenting the PLINDER dataset, eg. for flexible docking or pocket prediction purposes.
The files are grouped into zipped subdirectories by using `two_char_code` of the system.
Each unzipped subdirectory contains `pred` and `apo` subfolders that in turn contain folders named by `system_id`.
Inside each `apo/{system_id}` and `pred/{system_id}` folder is another directory containing a superposed system: `{source_id}_{chain_id}/superposed.cif`, where `{source_id}` and `{chain_id}` for apo systems is `pdb_id` with a source chain identifier, and for predicted structures, `{source_id}` is `uniprot_id` used in AF2DB with a chain identifier set to `A`.

### Linked systems (`links/`)

This directory contains parquet files linking PLINDER systems to their apo and predicted structures in `linked_structures/`.

:::{list-table} `{apo|pred}_links.parquet`
:widths: 10 5 30
:header-rows: 1

- - Name
  - Type
  - Description
- - reference_system_id
  - str
  - The PLINDER system ID
- - id
  - str
  - The PDB or AF2DB (for `apo` and `pred`, respectively) `{source_id}_{chain_id}` tag
- - pocket_fident
  - float
  - sequence identity for pocket residues
- - pocket_lddt
  - float
  - Local Distance Difference Test (lDDT) score for the pocket residue alpha carbons as returned by Foldseek.
- - protein_fident_qcov_weighted_sum
  - float
  - Sum of fident \* qcov for all templates, weighted by the number of residues in the template
- - protein_fident_weighted_sum
  - float
  - Sum of fident for all templates, weighted by the number of residues in the template
- - protein_lddt_weighted_sum
  - float
  - Sum of lDDT for all residues, weighted by the number of residues in the template
- - target_id
  - str
  - apo or pred stucture `{source_id}` tag
- - sort_score
  - float
  - Score used to sort linked structures. This is resolution for apos and plddt for preds.
- - receptor_file
  - str
  - intermediate aligned linked receptor file path
- - ligand_files
  - str
  - intermediate file path for ligands used in calculations
- - num_reference_ligands
  - int
  - number of ligands in reference structure
- - num_model_ligands
  - int
  - number of ligands in model structure
- - num_reference_proteins
  - int
  - number of protein chains in reference structure
- - num_model_proteins
  - int
  - number of protein chains in model structure
- - fraction_reference_ligands_mapped
  - float
  - Fraction of reference ligands that were successfully mapped to model ligands
- - fraction_model_ligands_mapped
  - float
  - Fraction of model ligands that were successfully mapped to reference ligands
- - lddt_pli_ave
  - float
  - Average lDDT score for protein-ligand interactions
- - lddt_pli_wave
  - float
  - Weighted average lDDT score for protein-ligand interactions
- - bisy_rmsd_ave
  - float
  - Average binding-site superposed symmetry-corrected RMSD between reference and model ligands
- - bisy_rmsd_wave
  - float
  - Weighted average binding-site superposed symmetry-corrected RMSD between reference and model ligands
- - lddt_lp_ave
  - float
  - Average lDDT score for ligand poses
- - lddt_lp_wave
  - float
  - Weighted average lDDT score for ligand poses
- - fraction_reference_proteins_mapped
  - float
  - Fraction of reference protein chains with corresponding model chains
- - fraction_model_proteins_mapped
  - float
  - Fraction of model protein chains mapped to corresponding reference chains
- - lddt
  - float
  - Global lDDT score calculated over all atoms in the structure
- - bb_lddt
  - float
  - Global lDDT score calculated over backbone atoms (N, CA, C, O) in the structure
- - per_chain_lddt_ave
  - float
  - Average per-chain lDDT score calculated over all atoms
- - per_chain_bb_lddt_ave
  - float
  - Average per-chain lDDT score calculated over backbone atoms (N, CA, C, O)
:::

(miscellaneous-target)=

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

#### Raw annotation parts (`raw_entries/`)

During ingest this directory contains one `{pdb_id}.parquet` annotation part and one per-entry directory containing `entry_chains.parquet`, `entry_source.parquet`, and canonical ligand SDFs, grouped by `two_char_code`. The join step consolidates the ligand-level parts into `index/annotation_table.parquet`, the normalized chain rows into `index/entry_chains.parquet`, and one pinned source revision per PDB into `index/entry_sources.parquet`; entry JSON archives are not produced.

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

#### Mapped protein alignments (`alignments/`, V3)

V3 release artifacts contain mapped Foldseek and MMseqs hits rather than the much
larger materialized pairwise score table:

```bash
|-- search_db=holo
|   |-- alignment_type=foldseek
|   |   |-- shard={two_char_code}.parquet
|   |-- alignment_type=mmseqs
|       |-- shard={two_char_code}.parquet
|-- search_db=apo
|-- search_db=pred
```

Rows are sorted by query entry, target entry, and mapped chain identifiers before
being written with Zstandard compression and bounded row groups. The public
`reconstruct_similarity_scores()` API uses Parquet filters to calculate a requested
system cross-product without materializing or downloading the global score table.

#### Protein similarity dataset (`scores/`, V2 and V3 ingest intermediate)

These tables contain the protein or pocket similarity scores used for clustering.
They remain the V2 release format and a local V3 ingest intermediate, but are not a
V3 release artifact after cluster columns have been merged into the annotation
parquet.

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

- - Name
  - Type
  - Description
- - query_system
  - str
  - The PLINDER system ID of query system
- - target_system
  - str
  - The PLINDER system ID of target system
- - protein_mapping
  - str
  - Chain mapping between query system and target system
- - protein_mapper
  - str
  - Alignment method used for mapping.
- - similarity
  - int
  - Similarity metric of interest
- - source
  - str
  - Source of similarity metric. It could either be `foldseek`, `mmseqs` or `both`
- - metric
  - str
  - Similarity metric of interest
- - mapping
  - str
  - Local region mapping between query system and target system
- - search_db
  - str
  - Search database type. Could be `apo`, `holo` or `pred`
:::
