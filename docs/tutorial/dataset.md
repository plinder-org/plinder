# Dataset tutorial

## Downloading a release

PLINDER data is published in `gs://plinder`. A release is identified by:

- `PLINDER_RELEASE`: the ingest month in `YYYY-MM` form;
- `PLINDER_RELEASE_NUMBER`: the numbered release within that month.

Install the package and download the compact index, clustering diagnostics, and
representative tables:

```bash
plinder_download --release 2026-07 --release-number 1
```

The command asks before downloading larger groups such as ligand archives,
alignments, scores, exports, and custom-scoring search databases. Pass `--yes`
to download all groups without prompts. If a group is skipped, an online API
call fetches the specific artifact it needs later.

Files can also be copied directly:

```console
$ export PLINDER_RELEASE=2026-07
$ export PLINDER_RELEASE_NUMBER=1
$ mkdir -p ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_RELEASE_NUMBER}/
$ gsutil -m cp -r \
    "gs://plinder/${PLINDER_RELEASE}/${PLINDER_RELEASE_NUMBER}/index" \
    ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_RELEASE_NUMBER}/
```

:::{note}
The preprint releases remain available under `gs://plinder/2024-04/v1` and
`gs://plinder/2024-04/v0`. Their layout and APIs differ from the current
release described here.
:::

## Querying ligand and entry data

`annotation_table.parquet` has one row per ligand. Start with `query_table()` so
only the selected columns and rows are read:

```python
from plinder.core import query_table

ligands = query_table(
    "annotation",
    columns=[
        "entry_pdb_id",
        "system_id",
        "ligand_id",
        "ligand_ccd_code",
        "ligand_smiles",
        "entry_resolution",
    ],
    joins=["entry_metadata"],
    filters=[("entry_resolution", "<=", 2.5)],
)

print(ligands.head())
```

Entry metadata is stored once in `entry_metadata.parquet`; the explicit join
adds it without changing the ligand row count. Other useful starting tables are
`entry_chains`, `entry_biounit_chains`, `interface_annotations`,
`linked_apo_structures`, and `ligand_mmp_pairs`.

For a fully local release, pass an explicit `PlinderRelease`:

```python
from pathlib import Path

from plinder.core import PlinderRelease, query_table

release = PlinderRelease(data_dir=Path("/data/plinder-release"))
chains = query_table(
    "entry_chains",
    columns=["entry_pdb_id", "chain_asym_id", "chain_receptor_type"],
    filters=[("chain_receptor_type", "==", "protein")],
    release=release,
)
```

## Selecting and loading protein interfaces

Protein interfaces have their own one-row-per-chain-pair table. Query it
directly instead of starting from ligand annotations:

```python
from plinder.core import query_table

interfaces = query_table(
    "interface_annotations",
    columns=[
        "system_id",
        "interface_chain_1",
        "interface_chain_2",
        "interface_chain_1_residue_indices",
        "interface_chain_2_residue_indices",
        "interface_num_contact_residue_pairs",
        "prodigy_label",
        "entry_resolution",
        "entry_source_taxonomy_ids",
        "representative_system_id",
    ],
    joins=["entry_metadata", "interface_membership"],
    filters=[
        ("interface_num_contact_residue_pairs", ">=", 10),
        ("entry_resolution", "<=", 3.0),
    ],
)

print(interfaces.head())
```

The chain pair is unordered. `interface_chain_1` and `interface_chain_2` are
stable assembly-instance IDs, not receptor and ligand roles. Entry-level source
organisms come from `entry_metadata`; full sequences come from `entry_chains`.

Load one row as a reconstructed two-chain complex:

```python
from plinder.core import PlinderInterface

interface = PlinderInterface(system_id=interfaces.iloc[0]["system_id"])

print(interface.chains)
print(interface.sequences)

complex_atoms = interface.atom_array
side_1 = interface.chain_structures[interface.chains[0]]
side_2 = interface.chain_structures[interface.chains[1]]
interface_atoms = interface.interface_structure
interface_cif = interface.interface_cif
```

These objects are Biotite `AtomArray` instances. The two masks in
`interface.interface_residue_masks` select the annotated contact surface from
`complex_atoms`, which is useful for residue-level featurization:

```python
side_1_mask = interface.interface_residue_masks[interface.chains[0]]
side_2_mask = interface.interface_residue_masks[interface.chains[1]]

side_1_interface_atoms = complex_atoms[side_1_mask]
side_2_interface_atoms = complex_atoms[side_2_mask]
```

To use a fixed list in a training pipeline, store the selected `system_id`
values and create `PlinderInterface` objects in the dataset's `__getitem__`.
This leaves atom selection, tensor conversion, and batching under the model's
control rather than imposing one protein-interface representation.

## Reconstructing a system

Choose a `system_id` from the ligand query and create a `PlinderSystem`:

```python
from pathlib import Path

from plinder.core import PlinderSystem

system = PlinderSystem(system_id="2y4i__1__1.B__1.E_1.F")

system_cif = Path(system.system_cif)
receptor_cif = Path(system.receptor_cif)
canonical_sdfs = system.canonical_ligand_sdfs
assembly_sdfs = system.ligand_sdfs
```

The first coordinate request fetches the source PDB mmCIF revision recorded by
the release. It then writes a self-contained system or receptor mmCIF locally.
Canonical SDFs are extracted from the ligand archive; assembly SDFs use the
coordinates of the selected biological-assembly instances.

`plinder_download` does not download the complete source PDB archive. To prepare
specific systems for an offline machine, fetch their unique PDB entries first:

```python
from plinder.core import download_pdb_mmcifs

download_pdb_mmcifs(["2y4i", "6cex"])
```

After the required release tables, source mmCIFs, and ligand archives are
cached, set `PLINDER_OFFLINE=true`.

## Reconstructing a linked apo chain

Some holo systems have ranked deposited apo-chain links:

```python
links = system.linked_apo_structures

if not links.empty:
    apo_cif = system.reconstruct_linked_apo()
    print(links[["linked_structure_id", "rank"]])
```

The output stays in the deposited apo coordinates. Fit it to a holo receptor
chain only when your application needs that frame:

```python
fitted_apo_cif = system.superpose_linked_apo(reference_chain="1.B")
```

If the holo receptor has exactly one protein chain, the reference chain can be
inferred. Multichain receptors require an explicit `reference_chain`.

(cluster-target)=

## Inspecting representative covers

Representative files live below `ligand_sampling/` and
`interface_sampling/`. For example, inspect the directional pocket cover at a
70-percent requested threshold:

```python
import pandas as pd

from plinder.core import PlinderRelease

release = PlinderRelease()
sampling_dir = release.fetch("ligand_sampling")
cover = pd.read_parquet(
    sampling_dir
    / "directed_set_cover"
    / "metric=pocket_qcov"
    / "threshold=70.parquet"
)

print(
    cover[
        [
            "ligand_id",
            "centroid_ligand_id",
            "similarity_to_centroid",
            "assignment_threshold",
            "label",
        ]
    ].head()
)
```

The detailed assignment table tells you whether a ligand was assigned at the
requested threshold or during the 50-percent fallback pass. The annotation
table also contains the compact label and centroid indicator columns for direct
filtering.

Tanimoto uses the undirected path instead:

```python
tanimoto_cover = pd.read_parquet(
    sampling_dir
    / "set_cover"
    / "metric=tanimoto_similarity_ecfp4_1024"
    / "threshold=70.parquet"
)
```

## Querying matched molecular pairs

The matched-molecular-pair table is another release index table:

```python
from plinder.core import query_table

pairs = query_table(
    "ligand_mmp_pairs",
    columns=[
        "ligand_smiles_id_1",
        "ligand_smiles_id_2",
        "transformation",
        "shared_core_smiles",
    ],
    filters=[("shared_core_num_heavy_atoms", ">=", 10)],
)

print(pairs.head())
```
