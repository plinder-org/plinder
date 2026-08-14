---
sd_hide_title: true
---

# Dataset

PLINDER, the **Protein & Ligand INteraction Dataset and Evaluation Resource**,
publishes ligand systems and protein-protein interfaces from the same PDB ingest
while keeping their distinct row grains and coordinate APIs.

## Release layout

PLINDER releases are addressed by an ingest month and a release number within
that month. The public layout is:

```text
<ingest-month>/
└── <release-number>/
    ├── index/
    │   ├── annotation_table.parquet
    │   ├── entry_chains.parquet
    │   ├── entry_biounit_chains.parquet
    │   ├── entry_metadata.parquet
    │   ├── entry_sources.parquet
    │   ├── interface_annotation_table.parquet
    │   ├── alignment_chain_lookup.parquet
    │   ├── linked_apo_structures.parquet
    │   ├── ligand_pocket_membership.parquet
    │   ├── ligand_pocket_representatives.parquet
    │   ├── ligand_mmp_pairs.parquet
    │   ├── interface_half_representatives.parquet
    │   ├── interface_membership.parquet
    │   └── interface_representatives.parquet
    ├── ligand_archives/
    │   ├── {two_char_code}.parquet
    │   └── manifest.json
    ├── alignments/
    │   └── search_db=holo/
    │       └── alignment_type={foldseek,mmseqs}/
    │           └── shard={two_char_code}.parquet
    ├── ligand_scores/
    ├── interface_scores/
    ├── ligand_clusters/
    ├── ligand_sampling/
    ├── interface_clusters/
    ├── interface_sampling/
    ├── search_databases/
    └── exports/
```

`plinder_download` downloads the index tables, clustering diagnostics, and
sampling tables by default. It asks before downloading the larger ligand,
alignment, score, export, and search-database groups; `--yes` downloads every
group. APIs fetch missing optional artifacts when they need them unless offline
mode is enabled. Source PDB mmCIFs are separate: they are fetched per PDB entry
during reconstruction and are not part of the bulk download.

(annotation-table-target)=

## Release tables

The annotation table has one row per ligand, not one row per system. Values
whose natural grain is an entry, chain, interface, or representative are stored
once in narrower tables:

- `annotation_table.parquet`: ligand annotations and reconstructable system IDs;
- `entry_chains.parquet`: one polymer chain in a PDB entry;
- `entry_biounit_chains.parquet`: one chain instance in a biological assembly;
- `entry_metadata.parquet`: experimental and entry-validation metadata;
- `entry_sources.parquet`: the source mmCIF revision used during ingest;
- `interface_annotation_table.parquet`: one protein-chain interface;
- `alignment_chain_lookup.parquet`: protein-chain search identifiers and residue mappings;
- `linked_apo_structures.parquet`: ranked deposited apo chains linked to holo systems;
- `ligand_pocket_membership.parquet`: ligand-to-pocket-representative assignments;
- `ligand_pocket_representatives.parquet`: receptor, pocket, and interaction payloads for ligand-pocket representatives;
- `ligand_mmp_pairs.parquet`: matched molecular pairs over unique ligand SMILES;
- `interface_membership.parquet`: interface-to-representative assignments;
- `interface_representatives.parquet`: representative full interfaces;
- `interface_half_representatives.parquet`: representative interface sides.

The checked-in column reference below is generated from these release tables.

:::{include} table.html
:::

### Querying tables

`query_table()` reads only the requested columns and rows. Optional joins are
explicit and are limited to relationships that preserve the base table's row
grain. This makes entry metadata or pocket membership available to a ligand
query without duplicating or dropping ligand rows.

```python
from plinder.core import query_table

ligands = query_table(
    "annotation",
    columns=[
        "system_id",
        "ligand_id",
        "ligand_ccd_code",
        "entry_resolution",
        "representative_ligand_id",
    ],
    joins=["entry_metadata", "ligand_pocket_membership"],
    filters=[
        ("entry_resolution", "<=", 2.5),
        ("system_pass_validation_criteria", "==", True),
    ],
)
```

Registered table names, row grains, and keys are available through
`plinder.core.RELEASE_TABLES`. `PlinderRelease` resolves the corresponding local
paths and can point at either the configured release or an explicit local copy.

```python
from pathlib import Path

from plinder.core import PlinderRelease, query_table

release = PlinderRelease(data_dir=Path("/data/plinder-release"))
entries = query_table(
    "entry_metadata",
    columns=["entry_pdb_id", "entry_release_date", "entry_resolution"],
    release=release,
)
```

## Protein interfaces

PLINDER publishes protein-protein interfaces in a dedicated table because their
natural row is a chain pair, not a ligand. `interface_annotation_table.parquet`
has one row for an unordered pair of protein-chain instances in a biological
assembly. Its `system_id` has the form
`<pdb>__<assembly>__<chain-instance-1>--<chain-instance-2>`.

The annotation records both chains, the resolved residues on each side, the
number of contacting residue pairs, and PRODIGY-cryst features. Join entry
metadata and representative membership when selecting a working set:

```python
from plinder.core import query_table

interfaces = query_table(
    "interface_annotations",
    columns=[
        "system_id",
        "entry_pdb_id",
        "interface_chain_1",
        "interface_chain_2",
        "interface_num_contact_residue_pairs",
        "prodigy_label",
        "prodigy_probability_bio",
        "entry_release_date",
        "entry_resolution",
        "entry_source_taxonomy_ids",
        "representative_system_id",
    ],
    joins=["entry_metadata", "interface_membership"],
    filters=[
        ("interface_num_contact_residue_pairs", ">=", 10),
        ("entry_resolution", "<=", 3.0),
        ("prodigy_label", "==", "BIO"),
    ],
)
```

`PlinderInterface` is the coordinate-level interface API. It expands the
deposited assembly and retains exactly the annotated two-chain pair. In-memory
chain keys use the assembly-instance IDs from the table, so repeated copies of
one asymmetric-unit chain remain distinct.

```python
from plinder.core import PlinderInterface

interface = PlinderInterface(system_id=interfaces.iloc[0]["system_id"])

annotation = interface.annotation
sequences = interface.sequences
complex_atoms = interface.atom_array
chain_atoms = interface.chain_structures
interface_atoms = interface.interface_structure
interface_masks = interface.interface_residue_masks
interface_cif = interface.interface_cif
```

`sequences` contains the full deposited polymer sequence for each side;
`atom_array` contains resolved coordinates. Each value in
`interface_residue_masks` is an atom mask over `atom_array`, derived from the
stored resolved-residue indices. The written mmCIF is self-contained and
contains only the two protein chains. The API does not assign a receptor and a
ligand orientation because the annotated interface is unordered.

An explicit `release=PlinderRelease(data_dir=...)` keeps reconstructed files
under that release root. A caller may also supply `source_mmcif`, but its
resolved residue numbering must match the release annotation; otherwise the
interface-residue properties raise a descriptive error instead of selecting
different residues.

`interface_representatives.parquet` stores the full-interface representatives,
while `interface_half_representatives.parquet` stores individual interface
sides. The latter is useful when one wants diverse protein surfaces without
requiring both sides of the same complex. Detailed set-cover assignments live
under `interface_sampling/`.

## Structure assets and reconstruction

PLINDER stores one canonical asymmetric-unit SDF for each ligand chain in
`ligand_archives/{two_char_code}.parquet`. Assembly copies are not stored because
they have the same conformation under a rigid-body transform. `PlinderSystem`
extracts only the SDFs required for the requested system.

System and receptor mmCIFs are rebuilt from the deposited PDB mmCIF, the ligand
rows, and `entry_biounit_chains.parquet`. The source revision recorded in
`entry_sources.parquet` is fetched from the
[wwPDB versioned archive](https://www.wwpdb.org/ftp/pdb-versioned-ftp-site) and
cached under the configured PLINDER directory. An explicit `source_mmcif`
supplied to `PlinderSystem` takes precedence.

```python
from pathlib import Path

from plinder.core import PlinderSystem

system = PlinderSystem(system_id="2y4i__1__1.B__1.E_1.F")
system_cif = Path(system.system_cif)
receptor_cif = Path(system.receptor_cif)
canonical_sdfs = system.canonical_ligand_sdfs
assembly_sdfs = system.ligand_sdfs
```

The reconstructed mmCIFs contain the atom, sequence, component, assembly, and
bond information needed to read them as self-contained PDBx/mmCIF files.
Canonical ligand SDFs preserve the curated bond orders; `ligand_sdfs` writes the
corresponding assembly coordinates.

For an offline workflow, fetch the required source files on an online node
before enabling `PLINDER_OFFLINE=true`:

```python
from plinder.core import download_pdb_mmcifs

download_pdb_mmcifs(["2y4i", "1a3b"])
```

Only the requested PDB entries are fetched. An offline request for an absent
source file or release artifact raises an error with its expected cache path.

### Linked apo chains

`linked_apo_structures.parquet` associates a holo system with ranked deposited
apo protein chains. A candidate must pass the configured pocket and whole-chain
similarity requirements for every proper ligand pocket in the holo system.
Candidates from the same PDB entry are excluded.

Ranking prefers chains with no nearby ligand-like components, then ion-only,
artifact, and other-ligand contacts. Resolution and similarity break later
ties. The release stores the exact assembly chain instance that was scored; it
does not store a copied or pre-fitted coordinate file.

```python
from plinder.core import PlinderSystem

system = PlinderSystem(system_id="2y4i__1__1.B__1.E_1.F")
links = system.linked_apo_structures

# The highest-ranked link is used when no ID is supplied.
apo_cif = system.reconstruct_linked_apo()

# Fit the selected apo chain to a holo receptor chain when desired.
fitted_apo_cif = system.superpose_linked_apo(reference_chain="1.B")
```

For a multichain receptor, `reference_chain` is required for fitting so that a
chain is never selected arbitrarily.

## Similarity artifacts

### Ligand similarities

`ligand_scores/` stores sharded BulkTanimoto edges over unique canonical-SMILES
nodes. Every edge at or above the configured minimum is retained.

`exports/all_sucos_shape_pocket_qcov.parquet` stores the complete published
ligand-level SuCOS/pocket-coverage export, including values below the clustering
threshold.

### Protein alignments

The release publishes mapped Foldseek and MMseqs hits rather than the complete
pairwise protein-score table. Alignment rows are ordered by query and target and
stored in deterministic PDB two-character shards.

`reconstruct_similarity_scores()` filters those shards to a requested system
cross-product and calculates directed ligand-level scores. Positive-pocket pairs
load their canonical ligand SDFs for the gated 3D metrics.

```python
from plinder.core.scores import (
    prefetch_similarity_alignments,
    reconstruct_similarity_scores,
)

queries = ["2y4i__1__1.B__1.E_1.F"]
targets = ["6cex__1__1.D__1.M"]
prefetch_similarity_alignments(queries)
scores = reconstruct_similarity_scores(queries, targets)
```

Prefetching downloads only the Foldseek/MMseqs shards required by the query PDB
entries. Run reconstruction once while online as well if the annotation rows or
positive-pocket ligand archives are not already cached.

Protein-interface coverage can be rebuilt for a bounded interface cross-product
from the same mapped alignments:

```python
from plinder.core.scores import reconstruct_interface_similarity_scores

interface_scores = reconstruct_interface_similarity_scores(
    query_interface_ids=["7cm8__1__1.A--2.A"],
    target_interface_ids=["7cma__1__1.A--1.B"],
)
```

The result is directional: swapping query and target can change interface
coverage.

## Representative covers

Connectivity components are build-time helpers and are not public cluster
labels. The release publishes greedy representative covers:

- `ligand_sampling/set_cover/` contains an undirected set cover for reciprocal
  Tanimoto similarity;
- `ligand_sampling/directed_set_cover/` contains directed covers for pocket,
  interaction, and pocket-weighted ligand 3D metrics;
- `interface_sampling/directed_set_cover/` contains directed covers for protein
  interfaces.

Each metric has files named `metric={metric}/threshold={threshold}.parquet`.
Ligand cover labels are also merged into the annotation table. Tanimoto columns
end in `__ligand__set_cover`; directional columns end in
`__ligand__directed_set_cover`.

The Tanimoto cover uses only reciprocal edges meeting the requested threshold.
For a directional cover, representative selection first maximizes residual gain
at the requested threshold. If uncovered nodes remain above a threshold of 50,
the algorithm uses 50-percent edges before assigning remaining singletons. The
detailed file records `representative_selection_threshold` and
`assignment_threshold`, so downstream selection can distinguish strict and
fallback assignments. It also records selection order, marginal gain, and each
node's potential coverage count and fraction.

## Matched molecular pairs

`index/ligand_mmp_pairs.parquet` contains compact mmpdb transformations over
unique canonical ligand SMILES. Rows identify the two `ligand_smiles_id` values,
both SMILES, the transformation and shared core, cut count, heavy-atom counts,
and the fraction of each ligand contained in the shared core.

```python
from plinder.core import query_table

pairs = query_table(
    "ligand_mmp_pairs",
    columns=[
        "ligand_smiles_id_1",
        "ligand_smiles_id_2",
        "transformation",
        "shared_core_smiles",
        "ligand_1_shared_core_fraction",
        "ligand_2_shared_core_fraction",
    ],
    filters=[("ligand_1_shared_core_fraction", ">=", 0.5)],
)
```
