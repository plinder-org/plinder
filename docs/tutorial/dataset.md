# Getting started

The `plinder` Python package is the main entry point for downloading release
files, querying annotations, and reconstructing structures.

PLINDER downloads the files needed by each API call automatically and reuses
them on subsequent calls. You can start exploring straight away.

## Install PLINDER

```bash
pip install plinder
```

A release is identified by the month of its PDB snapshot. The current release
(`2026-09`) is the default; select another with:

```bash
export PLINDER_RELEASE=2026-09
```

## Continue with a runnable guide

Start with {doc}`/examples/2_query_filter_index` to select ligand or protein-interface
records. The remaining {doc}`/examples/index` guides cover coordinate access,
similarity analysis, and custom scoring.

For a local copy to use offline, see {doc}`/examples/1_download`.

Use the {doc}`/dataset` page as the release and table reference, and
{doc}`/api/index` for individual Python signatures.
