---
sd_hide_title: true
---

# Dataset

PLINDER, the **Protein & Ligand INteraction Dataset and Evaluation Resource**,
extracts protein-ligand systems and protein-protein interfaces from the PDB.

## Release layout

PLINDER releases are identified by the month of their PDB snapshot and a
release number within that month. The file layout is:

:::{note}
Historical releases remain available at `gs://plinder/2024-04/v0`,
`gs://plinder/2024-04/v1`, and `gs://plinder/2024-06/v2`. Their layout and APIs
differ from the current release.
:::

```text
<release-month>/
└── <release-number>/
    ├── index/
    │   ├── annotation_table.parquet
    │   ├── system_validation.parquet
    │   ├── entry_chains.parquet
    │   ├── entry_biounit_chains.parquet
    │   ├── entry_metadata.parquet
    │   ├── entry_sources.parquet
    │   ├── interface_annotation_table.parquet
    │   ├── alignment_chain_lookup.parquet
    │   ├── linked_apo_structures.parquet
    │   ├── ligand_pocket_membership.parquet
    │   ├── ligand_pocket_representatives.parquet
    │   ├── ligand_clusters.parquet
    │   ├── ligand_mmp_pairs.parquet
    │   ├── interface_half_representatives.parquet
    │   ├── interface_membership.parquet
    │   ├── interface_representatives.parquet
    │   └── interface_clusters.parquet
    ├── ligand_archives/
    │   ├── {two_char_code}.parquet
    │   └── manifest.json
    ├── alignments/
    │   └── search_db=holo/
    │       └── alignment_type={foldseek,mmseqs}/
    │           └── shard={two_char_code}.parquet
    ├── ligand_scores/
    │   └── {fragment}.parquet
    ├── interface_scores/
    │   ├── shard={two_char_code}.parquet
    │   └── shard={two_char_code}.json
    ├── ligand_sampling/
    │   ├── set_cover/metric={metric}/threshold={threshold}.parquet
    │   └── directed_set_cover/metric={metric}/threshold={threshold}.parquet
    ├── interface_sampling/
    │   └── directed_set_cover/metric=interface_qcov/threshold={threshold}.parquet
    ├── search_databases/
    │   ├── manifest.json
    │   ├── holo_foldseek/
    │   └── holo_mmseqs/
    └── exports/
        ├── ligand_similarity_scores.parquet
        └── interface_similarity_scores.parquet
```

`plinder_download` downloads the index and representative-cover tables by
default. It asks before downloading ligand archives, alignments, score files,
complete similarity tables, and search databases; `--yes` downloads every
group. APIs fetch missing optional files when they need them unless offline
mode is enabled. Source PDB mmCIFs are separate: they are fetched per PDB entry
during reconstruction and are not part of the bulk download.

(annotation-table-target)=
(annotation-tables-index)=

## Release tables

The following are all registered release tables available through
`query_table()`. The main tables of interest are `annotation`, which has one
row per ligand, and `interface_annotations`, which has one row per
protein-protein interface. The remaining tables separate entry,
chain, validation, cluster, and representative data so clients only load it
when needed. `query_table()` can request columns across related tables in a
single query and selects the required relationships from those columns.

- `annotation`: ligand annotations and reconstructable system IDs;
- `system_validation`: ligand and pocket validation summaries for each
  system;
- `entry_chains`: one polymer chain in a PDB entry;
- `entry_biounit_chains`: one chain instance in a biological assembly;
- `entry_metadata`: experimental and entry-validation metadata;
- `entry_sources`: the source mmCIF revision used for the release;
- `interface_annotations`: one protein-chain interface;
- `alignment_chain_lookup`: protein-chain search identifiers and
  residue mappings;
- `linked_apo_structures`: ranked apo chains linked to holo ligands;
- `ligand_pocket_membership`: ligand-to-pocket-representative
  assignments;
- `ligand_pocket_representatives`: the receptor chains, pocket
  residues, and protein-ligand interactions for each selected representative;
- `ligand_clusters`: ligand cover labels, representative
  indicators, and coverage statistics;
- `ligand_mmp_pairs`: matched molecular pairs over unique ligand SMILES;
- `interface_membership`: interface-to-representative assignments;
- `interface_representatives`: representative full interfaces;
- `interface_half_representatives`: representative interface sides;
- `interface_clusters`: protein-interface cover labels.

The checked-in column reference below is generated from these release tables.

:::{include} table.html
:::

### Table relationships

The tables can be read directly with any Parquet reader. The identifiers used
to combine them are:

- `entry_pdb_id` for entry metadata and source revisions;
- the pair (`entry_pdb_id`, `chain_asym_id`) between `entry_biounit_chains`
  and `entry_chains`;
- `system_id` for ligand systems, system validation, protein interfaces, and
  interface clusters;
- `ligand_id` for ligand annotations, pocket-representative membership, and
  ligand clusters;
- `representative_ligand_id` and `representative_system_id` for the compact
  representative tables.

All registered relationships point from a table to rows that are unique on the
join columns. Adding entry metadata, system validation, or representative
membership therefore does not multiply the rows of the starting table. See
{doc}`/examples/2_query_filter_index` for examples that select the required
tables automatically.

(protein-interface-reference)=

## Protein interfaces

`interface_annotations` has one row for an unordered pair of
protein-chain instances in a biological assembly. Its `system_id` has the form
`<pdb>__<assembly>__<chain-instance-1>--<chain-instance-2>`.
Assembly-chain IDs such as `1.A` and `2.A` distinguish repeated copies of the
same asymmetric-unit chain.

The table records both chain-instance IDs, resolved residue numbers and indices
for each side, the number of contacting residue pairs, and PRODIGY-cryst
features. Full deposited sequences and source taxonomies are in the chain and
entry tables.

`interface_representatives` stores the full-interface representatives,
while `interface_half_representatives` stores individual interface
sides. The latter is useful when one wants diverse protein surfaces without
requiring both sides of the same complex. `interface_clusters` provides
cover columns for `query_table()`, while detailed assignments are under
`interface_sampling/`.

(structure-asset-reference)=

## Structures and reconstruction

PLINDER stores one asymmetric-unit SDF for each ligand chain in
`ligand_archives/{two_char_code}.parquet`. Assembly copies are not stored because
they have the same conformation under a rigid-body transform. The archive stores
the SDF together with the PDB entry ID and asymmetric-unit ligand ID.

`PlinderSystem.reconstruct()` rebuilds selected system structures from the
deposited PDB mmCIF, the ligand rows, and `entry_biounit_chains.parquet`.
Accessing `PlinderSystem.system_cif` or `PlinderSystem.receptor_cif` performs
reconstruction when the complete system or ligand-free receptor mmCIF is not
already cached. Reconstruction options can include additional chains from the
same PDB entry.

`PlinderInterface.reconstruct()` similarly rebuilds a protein-protein interface,
and `PlinderInterface.interface_cif` returns its standard mmCIF path. The source
revision recorded in `entry_sources.parquet` is fetched from the
[wwPDB versioned archive](https://www.wwpdb.org/ftp/pdb-versioned-ftp-site) and
cached separately from the bulk release. Source PDB mmCIFs are therefore not
included by `plinder_download`; `download_pdb_mmcifs()` can prefetch selected
source revisions for offline reconstruction.

The reconstructed mmCIFs contain the atom, sequence, component, assembly, and
bond information needed to read them as self-contained PDBx/mmCIF files.
Asymmetric ligand SDFs preserve the curated bond orders. Assembly-coordinate
SDFs use the same chemistry with the recorded rigid-body transforms applied.

### Linked apo chains

`linked_apo_structures` associates a holo ligand system with ranked
deposited apo protein chains. A candidate must pass the configured pocket and
whole-chain similarity requirements for every proper ligand pocket in the holo
system. Candidates from the same PDB entry are excluded.

Ranking prefers chains with no nearby ligand-like components, then ion-only,
artifact, and other-ligand contacts. Resolution and similarity break later
ties. The release stores the assembly chain instance that was scored;
`PlinderSystem.reconstruct_linked_apo()` rebuilds its coordinates from the
recorded source entry and assembly membership.

(similarity-score-reference)=

## Similarity scores

### Protein alignments

The release includes mapped Foldseek and MMseqs chain-level alignments grouped
by query PDB code. These alignments provide the chain and residue mappings used
to calculate protein sequence and structure scores, ligand-pocket and
protein-ligand interaction scores, and protein-interface scores.

### Ligand similarities

`exports/ligand_similarity_scores.parquet` stores the complete directed
ligand-pair table derived from the corresponding protein alignments and ligand
comparisons. Each row contains `pocket_qcov`, `pocket_fident_qcov`, `pli_qcov`,
and `sucos_shape`. Swapping query and target can change the values. Symmetric
Tanimoto similarities calculated from 1,024-bit Morgan fingerprints with
radius 2 (ECFP4) over unique canonical SMILES remain in `ligand_scores/`, which
is used to build the reciprocal chemical covers.

The Tanimoto rows contain `query_ligand_id`, `target_ligand_id`, and
`tanimoto_similarity_ecfp4_1024`; both ID columns refer to `ligand_smiles_id`
values in `annotation`.

### Protein-interface similarities

`exports/interface_similarity_scores.parquet` contains the complete directed
interface table with `iface1_qcov`, `iface2_qcov`, and their combined
`similarity` value. Swapping query and target can change the coverage.

(representative-cover-reference)=

## Representative covers

- `ligand_sampling/set_cover/` contains an undirected set cover for reciprocal
  Tanimoto similarity;
- `ligand_sampling/directed_set_cover/` contains directed covers for pocket,
  interaction, and pocket-weighted ligand 3D metrics;
- `interface_sampling/directed_set_cover/` contains directed covers for protein
  interfaces.

Each metric has files named `metric={metric}/threshold={threshold}.parquet`.
The release also includes one row per ligand in
`index/ligand_clusters.parquet` and one row per interface to
`index/interface_clusters.parquet`. `query_table()` joins these columns to the
corresponding annotation table when requested. Tanimoto columns end in
`__ligand__set_cover`; directional ligand columns end in
`__ligand__directed_set_cover`.

The Tanimoto cover uses only reciprocal edges meeting the requested threshold.
For a directional cover, representative selection first maximizes residual gain
at the requested threshold. If uncovered nodes remain above a threshold of 50,
the algorithm uses 50-percent edges before assigning remaining singletons. The
detailed file records `representative_selection_threshold` and
`assignment_threshold`, so downstream selection can distinguish strict and
fallback assignments. It also records selection order, marginal gain, and each
node's potential coverage count and fraction.

## Custom-scoring databases

`search_databases/` contains Foldseek structure and MMseqs sequence databases
for PLINDER receptor and interface chains, together with the chain mappings
used by custom scoring.

## Matched molecular pairs

`ligand_mmp_pairs` contains compact mmpdb transformations over
unique canonical ligand SMILES. Rows identify the two `ligand_smiles_id` values,
both SMILES, the transformation and shared core, cut count, heavy-atom counts,
and the fraction of each ligand contained in the shared core.
