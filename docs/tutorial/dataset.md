# Dataset tutorial

## Getting the data

:::{todo}
- link the data location
- explain what a bucket is
- reference installation guide of `gsutil`
- give CLI commands for getting the the parquet file
:::

## Understanding the directory structure

The directory downloaded from the bucket has the following structure:

```bash
2024-04/                     # the "`plinder` release" (`PLINDER_RELEASE`)
|-- v1                       # the "`plinder` iteration" (`PLINDER_ITERATION`)
|   |-- clusters             # pre-calculated cluster labels derived from the protein similarity dataset
|   |-- dbs                  # TSVs containing the raw files and IDs in the foldseek and mmseqs sub-databases
|   |-- entries              # raw annotations prior to consolidation (split by `two_char_code` and zipped)
|   |-- fingerprints         # index mapping files for the ligand similarity dataset
|   |-- index                # consolidated tabular annotations
|   |-- leakage              # leakage results
|   |-- ligand_scores        # ligand similarity parquet dataset
|   |-- ligands              # ligand data expanded from entries for computing similarity
|   |-- mmp                  # MMP and MMS data
|   |-- scores               # protein similarity parquet dataset
|   |-- splits               # split files and the configs used to generate them (if available)
|   |-- systems              # structure files for all systems (split by `two_char_code` and zipped)
```

:::{todo}
- explain most important parts of directory structure
- explain split into subdirectories of chunks
- mention location of central `.parquet` file
:::

## Interpreting the table

:::{todo}
- explain what a parquet file is and link to parquet reference
- explain how to open it (pandas and maybe some CLI tool?)
- show an excerpt of it (most important columns, e.g. the system name column, must be included)
- explain how to to interpret the columns in general
:::

## Inspecting the entries

:::{todo}
- transition from system name in table to structure entries
- mention where the structure files are
- mention available formats
- show prominent excerpt of structure file?
:::