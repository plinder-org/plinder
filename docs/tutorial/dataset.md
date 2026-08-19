# Getting started

The `plinder` Python package is the main entry point for downloading release
files, querying annotations, and reconstructing structures.

## Install PLINDER

```bash
pip install plinder
```

A release is identified by the month of the PDB release and a number within that month:

```bash
export PLINDER_RELEASE=2026-07
export PLINDER_RELEASE_NUMBER=1
```

Download the index and representative-cover tables:

```bash
plinder_download --release 2026-07 --release-number 1
```

The downloader asks separately about large assets such as ligand archives,
similarity scores, alignments, and custom-scoring search databases. When
`PLINDER_OFFLINE` is unset, API calls can fetch an omitted artifact when it is
first needed.

:::{note}
The preprint releases remain available under v0-v2 of `gs://plinder/2024-04/`.
They are no longer actively supported and their layout/APIs differ from the release described here.
:::

## Query ligands

The two main tables are `annotation_table.parquet`, which has one row per
ligand, and `interface_annotation_table.parquet`, which has one row per
protein-protein interface. Use `query_table()` to read only the columns and rows
needed for an analysis.

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
    filters=[("entry_resolution", "<=", 2.5)],
)

print(ligands.head())
```

## Reconstruct a ligand system

Create a `PlinderSystem` from a complete `system_id`:

```python
from plinder.core import PlinderSystem

system = PlinderSystem(system_id="2y4i__1__1.B__1.E_1.F")

system_annotations = system.system
receptor_atoms = system.receptor_structure
ligand_atoms = system.ligand_structures
system_cif = system.system_cif
ligand_sdfs = system.ligand_sdfs
```

The first coordinate request fetches the deposited PDB mmCIF revision recorded
by the release and writes a self-contained system locally. Asymmetric ligand
SDFs come from the release ligand archive; `ligand_sdfs` contains copies in the
selected assembly coordinates.

## Query and reconstruct a protein interface

Protein interfaces are queried independently of ligand annotations:

```python
from plinder.core import PlinderInterface, query_table

interfaces = query_table(
    "interface_annotations",
    columns=[
        "system_id",
        "interface_chain_1",
        "interface_chain_2",
        "interface_num_contact_residue_pairs",
        "entry_resolution",
    ],
    filters=[
        ("interface_num_contact_residue_pairs", ">=", 10),
        ("entry_resolution", "<=", 3.0),
    ],
)

interface = PlinderInterface(system_id=interfaces.iloc[0]["system_id"])
interface_atoms = interface.atom_array
interface_cif = interface.interface_cif
```

`interface_atoms` is a Biotite `AtomArray`. `interface_cif` is the path to the
written two-chain mmCIF.

## Continue with a guide

The {doc}`/examples/index` page groups complete workflows for data access,
structure reconstruction, similarity analysis, and custom scoring.

See {doc}`/dataset` for the release layout and table reference, or
{doc}`/api/index` for the generated Python API reference.
