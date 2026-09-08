# Getting started

The `plinder` Python package is the main entry point for downloading release
files, querying annotations, and reconstructing structures.

## Install PLINDER

```bash
pip install plinder
```

A release is identified by the month of its PDB snapshot and a number within
that month:

```bash
export PLINDER_RELEASE=2026-07
export PLINDER_RELEASE_NUMBER=1
```

## Continue with a runnable guide

Start with {doc}`/examples/1_download`, then use
{doc}`/examples/2_query_filter_index` to select ligand or protein-interface
records. The remaining {doc}`/examples/index` guides cover coordinate access,
similarity analysis, and custom scoring.

Use the {doc}`/dataset` page as the release and table reference, and
{doc}`/api/index` for individual Python signatures.
