# Getting started

The `plinder` Python package is the main entry point for downloading release
files, querying annotations, and reconstructing structures.

## Install PLINDER

```bash
pip install plinder
```

Releases download from Cloudflare R2
(`https://plinderdata.org/<release-month>/`) with checksum verification and automatic retries;
no Google Cloud SDK is required. `PLINDER_MIRROR_URL` selects another HTTP(S)
mirror serving the same release layout and manifest.

A release is identified by the month of its PDB snapshot. The current release
(`2026-09`) is the default; select another with:

```bash
export PLINDER_RELEASE=2026-09
```

## Continue with a runnable guide

Start with {doc}`/examples/1_download`, then use
{doc}`/examples/2_query_filter_index` to select ligand or protein-interface
records. The remaining {doc}`/examples/index` guides cover coordinate access,
similarity analysis, and custom scoring.

Use the {doc}`/dataset` page as the release and table reference, and
{doc}`/api/index` for individual Python signatures.
