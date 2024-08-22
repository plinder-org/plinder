![plinder](https://github.com/user-attachments/assets/05088c51-36c8-48c6-a7b2-8a69bd40fb44)

<div align="center">
    <h1>The Protein Ligand INteractions Dataset and Evaluation Resource</h1>
</div>

---

[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/plinder-org/plinder/blob/master/LICENSE.txt)
[![publish](https://github.com/plinder-org/plinder/actions/workflows/main.yaml/badge.svg)](https://github.com/plinder-org/plinder/pkgs/container/plinder)
[![website](https://img.shields.io/badge/website-plinder-blue.svg)](https://www.plinder.sh/)
[![bioarXiv](https://img.shields.io/badge/bioarXiv-2024.07.17.603955v1-blue.svg)](https://www.biorxiv.org/content/10.1101/2024.07.17.603955v1)
[![docs](https://github.com/plinder-org/plinder/actions/workflows/docs.yaml/badge.svg)](https://plinder-org.github.io/plinder/)
[![coverage](https://github.com/plinder-org/plinder/raw/python-coverage-comment-action-data/badge.svg)](https://github.com/plinder-org/plinder/tree/python-coverage-comment-action-data)

![overview](https://github.com/user-attachments/assets/39d251b1-8114-4242-b9fc-e0cce900d22f)

# 📚 About

![plinder_banner](https://github.com/user-attachments/assets/43d129f2-3bb6-4903-81fa-182c351c64b6)

**plinder**, short for **p**rotein **l**igand **in**teractions **d**ataset and **e**valuation **r**esource, is a dataset and resource for training and evaluation of protein-ligand docking algorithms.

It is a comprehensive, annotated, high quality dataset, including:

- \> 400k PLI systems across > 11k SCOP domains and > 50k unique small molecules
- 500+ annotations for each system, including protein and ligand properties, quality, matched molecular series and more
- Automated curation pipeline to keep up with the PDB
- 14 PLI metrics and over 20 billion similarity scores
- Unbound \(_apo_\) and _predicted_ Alphafold2 structures linked to _holo_ systems
- `train-val-test` splits and ability to tune splitting based on the learning task
- Robust evaluation harness to simplify and standard performance comparison between models

The `plinder` project is a community effort, launched by the University of Basel, SIB Swiss Institute of Bioinformatics, VantAI, NVIDIA, MIT CSAIL, and will be regularly updated.
We highly welcome contributions!
If you find `plinder` useful, please see the citation file for details on how to cite.

To accelerate community adoption, PLINDER will be used as the field’s new Protein-Ligand interaction dataset standard as part of an exciting competition at the upcoming 2024 [Machine Learning in Structural Biology (MLSB)](https://mlsb.io/) Workshop at NeurIPS, one of the fields’ premiere academic gatherings, which will be announced shortly.

# 👨💻 Getting Started

Please use a virtual environment for the `plinder` project.
We recommend the [miniforge](https://github.com/conda-forge/miniforge) environment manager.

**NOTE**: We currently only support a Linux environment. `plinder`
uses `openstructure` for some of its functionality and is available
from the `aivant` conda channel using `conda install aivant::openstructure`,
but it is only built targeting Linux architectures. We also depend on
`networkit>=11.0` which as of the time of writing does not build from source
cleanly on MacOS. For MacOS users, please see the relevant
[docker](#package-publishing) resources below.

## Install plinder

The `plinder` package can be obtained from GitHub:

    git clone https://github.com/plinder-org/plinder.git
    cd plinder
    mamba env create -f environment.yml
    mamba activate plinder
    pip install .

Or with a development installation:

    cd plinder
    pip install -e '.[dev]'

# 📚 Documentation

The documentation can be found [here](https://plinder-org.github.io/plinder/). We also provide a quick "Getting started" below.


# ⬇️ Getting the dataset

We provide 2 ways of interacting with `plinder`:

1. A python API: using the `plinder.core` API, you can transparently and lazily
download and interact with most of the components of the dataset. **(WIP)**

2. Via raw files on bucket: if you prefer to use the dataset directly, you can fetch it using [`gsutil`](https://cloud.google.com/storage/docs/gsutil) from google cloud storage. **(recommended)**


If you go with route 2, see below sections.

## Downloading the dataset

    export PLINDER_RELEASE=2024-04 # Current release
    export PLINDER_ITERATION=v1 # Current iteration
    gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/* ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/

**NOTE!**: the version used for the preprint is `gs://plinder/2024-04/v1`, while the current version we plan to use for the MLSB challenge is `gs://plinder/2024-06/v2`.

This yields the following structure, with the `systems`, `splits`, and `index/annotation_table.parquet` being most important for direct usage and the rest containing useful annotations and medadata.


```bash
2024-06/                     # The "`plinder` release" (`PLINDER_RELEASE`)
|-- v2                       # The "`plinder` iteration" (`PLINDER_ITERATION`)
|   |-- systems              # Actual structure files for all systems (split by `two_char_code` and zipped)
|   |-- splits               # List of system ids in a .parquet and each split  the configs used to generate them (if available)
|   |-- clusters             # Pre-calculated cluster labels derived from the protein similarity dataset
|   |-- entries              # Raw annotations prior to consolidation (split by `two_char_code` and zipped)
|   |-- fingerprints         # Index mapping files for the ligand similarity dataset
|   |-- index                # Consolidated tabular annotations
|   |-- leakage              # Leakage results
|   |-- ligand_scores        # Ligand similarity parquet dataset
|   |-- ligands              # Ligand data expanded from entries for computing similarity
|   |-- linked_structures    # Linked structures
|   |-- links                # Table of linked (apo and pred) structures and their respective similarity scores.
|   |-- mmp                  # Matched Molecular Series/Pair data
|   |-- scores               # Extended protein similarity parquet dataset
```
**Note:** We added a new sub-directory to `links` to v2 which contains table of linked (apo and pred) structures and their respective similarity scores.


## Unpacking the structure files:

Identical to the PDB NextGen Archive, we split the structures into subdirectories of chunks to make loading and querying speed palatable.

The structure files can be found in the subfolder `plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems`. To unpack the structures do the following

```bash
cd plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/systems; for i in ls *zip; do unzip $i; done
```

This will yield directories such as `7zzh__1__1.A_4.A__4.C`, which is what we call a `plinder system id` of the form `<pdbid>__<biounit>__<chain_ids_of_receptor>__<chain_ids_of_ligand>`.
The directory contains `mmcif`, `pdb` and `sdf` file formats as well as some additional metadata (like chain mapping and fastas).

The naming of the directory is by `system` - since a given PDB entry can contain multiple protein-ligand complexes, we assign a plinder system id to every system. These `system id`s are also used for retrieving the splits, as shown in the next step

### Accessing the splits files:

`plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/splits` contains an index for splits contained in a single parquet file. The most current split is `gs://plinder/2024-04/v1/splits/plinder-pl50.parquet` containing the pl50 split from the preprint. (Again: note this will be shortly updated to v2 and the v2 split will be used for the MLSB leaderboard)

```python
>>> import pandas as pd
>>> df = pd.read_parquet("gs://plinder/2024-04/v1/splits/plinder-pl50.parquet")
>>> df.head()
               system_id  split  cluster cluster_for_val_split
0      3grt__1__1.A__1.B  train  c100718                    c0
1  3grt__1__1.A_2.A__1.C  train  c196491                    c0
2      3grt__1__2.A__2.B  train  c100727                    c0
3  3grt__1__1.A_2.A__2.C  train  c234445                    c0
4      1grx__1__1.A__1.B  train  c186691                  c154
```

With the `system_id` contained in these split files, you can load the respective train, val & test splits **after unzipping** the `systems` directory. E.g. as shown in the Dataframe above, `~/.local/share/plinder/2024-04/v1/systems/1grx__1__1.A__1.B/system.cif` will contain the full mmcif of the system. We also provide cif files of seperated receptors (`*/receptor.cif`) and ligands (`*/ligand_files/*.sdf`) as well as pdb files (`*/receptor.pdb`) but **strongly encourage cif**, pdb is considered a [legacy file format](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif).


Note: a non-redundant and single-ligand smaller subset of this version was used to train diffdock in the paper and is available at 2024-04/v0.

The folder also contains a `.yaml` which is the config used to generate the split and can be ignored unless you want to reproduce the splits.


## 🔢 Plinder versions

We version the `plinder` dataset with two controls:

- `PLINDER_RELEASE`: the month stamp of the last RCSB sync
- `PLINDER_ITERATION`: value that enables iterative development within a release

We version the `plinder` application using an automated semantic
versioning scheme based on the `git` commit history.
The `plinder.data` package is responsible for generating a dataset
release and the `plinder.core` package makes it easy to interact
with the dataset.

Changelog:
- 2024-06/v2:
    - Improved SDF saving to handle some bond order issues
    - Updated system definition to be more stable and independent of PLIP
    - Added binding affinities from BindingDB
    - Annotated all crystal contacts
    - Improved covalency detection
    - Added linked apo/pred structures to v2/links and v2/linked_structures
    - Added statistics requirement in the split, and lots of deduplication and enrichment for unique ligand.

- 2024-04/v1 (Current): Version with redundancy removal by protein pocket and ligand similarity.
- 2024-04/v0: Version used to re-train DiffDock in the paper, with redundancy removal based on \<pdbid\>\_\<ligand ccd codes\>

## 🏅 Gold standard benchmark sets

As part of `plinder` resource we also provide train, validation and test splits that are curated to minimize the information leakage based on protein-ligand interaction similarity. In addition, we have prioritized the systems that has a linked experimental `apo` structure or matched molecular series to support realistic inference scenarios for hit discovery and optimization.
Finally, a particular care is taken for test set that is further prioritized to contain high quality structures to provide unambiguous ground-truths for performance benchmarking.

![test_stratification](https://github.com/user-attachments/assets/5bb96534-f939-42b5-bf85-5ac3a71aa324)

Moreover, as we enticipate this resource to be used for benchmarking a wide range of methods, including those simultaneously predicting protein structure (aka. co-folding) or those generating novel ligand structures, we further stratified test (by novel ligand, pocket, protein or all) to cover a wide range of tasks.

## ⚖️ Evaluation harness

See the [`plinder.eval`](docs/eval.md) docs for more details.

## 📦 Dataloader

The `plinder.data.loader` package contains a `PyTorch` dataloader for the dataset using the `atom3d` format. It is an example of using the `plinder.core` API
to implement a dataloader, but is not the only way to use the dataset.

**Note**: The dataloader requires both `torch` and `atom3d` to be installed. You use the `[loader]` dependency block when installing `plinder`:

    pip install .[loader]

## 📡 Future work

We are currently working on the following:

- Implementing the Dataloader
- Establishing a leaderboard
- Improving the documentation and examples

# 👨💻 Code organization

This code is split into 4 sub-packages

- `plinder.core`: core data structures for interacting with and loading the dataset.
- `plinder.data`: core code for generating the dataset
- `plinder.eval`: evaluation harness for the dataset that takes as an input predicted and ground truth structures in a pre-determined folder structure and returns a leaderboard-ready set of entries
- `plinder.methods`: implementations of the methods in the leaderboard that leverage plinder-primitives for training & running

# 💽 Dataset Generation

![workflow](https://github.com/user-attachments/assets/cde72643-5fdf-4998-8719-216d0cef2706)

See the [End-to-end pipeline](docs/data/process.md) description for technical details about the dataset generation.

# 📝 Examples & documentation

Package documentation, including API documentation, [example notebooks](examples/), and supplementary guides, are made available.

# ⚙️ Dev guide

To develop and test changes to the source code, please use a development installation:

## Dev mode install

    git clone https://github.com/plinder-org/plinder.git
    # OR
    git clone git@github.com:plinder-org/plinder.git
    cd plinder
    mamba env create -f environment.yml
    mamba activate plinder
    pip install -e '.[dev]'

Please install pre-commit hooks (the same checks will run in CI):

    pre-commit install

## Test suite

Test linting checks:

    tox -e lint

Run typing checks:

    tox -e type

Run the test suite:

    tox -e test

We lint with `ruff`. See `tox.ini` and `.pre-commit-config.yaml` for details.

## Debugging

In order to change log levels in `plinder`, please set:

    export PLINDER_LOG_LEVEL=10

## Contributing

This is a community effort and as such we highly encourage contributions.

## Package publishing

We publish the `plinder` project as a docker container, to ensure the highest
level of compatibility with non-Linux platforms. See the relevant docker resources
here for more details:

- `docker-compose.yml`: defines a `base` image, the `plinder` "app" and a `test` container
- `dockerfiles/base/`: contains the files for the `base` image
- `dockerfiles/main/`: contains the files for the `plinder` "app" image

The CI workflow will automatically semver bump the `plinder` version and publish
the `plinder` image to the GitHub Container Registry on merges to `main`. Control
over the version bumping semantics is handled by inspecting the commit history
since the previous release:

- If `bumpversion skip` is present in the commit message, the version will not be bumped
- If `bumpversion major` is present in the commit message, the major version will be bumped
- If `bumpversion minor` is present in the commit message, the minor version will be bumped
- If `bumpversion patch` is present in the commit message (or nothing is found), the patch version will be bumped

**NOTE**: The CI workflow will use the **most recent** match in the commit history to make its decision.

# 📃 Publications

Durairaj, Janani, Yusuf Adeshina, Zhonglin Cao, Xuejin Zhang, Vladas Oleinikovas, Thomas Duignan, Zachary McClure, Xavier Robin, Emanuele Rossi, Guoqing Zhou, Srimukh Prasad Veccham, Clemens Isert, Yuxing Peng, Prabindh Sundareson, Mehmet Akdel, Gabriele Corso, Hannes Stärk, Zachary Wayne Carpenter, Michael M. Bronstein, Emine Kucukbenli, Torsten Schwede, Luca Naef. 2024. “PLINDER: The Protein-Ligand Interactions Dataset and Evaluation Resource.”
[bioRxiv](https://doi.org/10.1101/2024.07.17.603955)
[ICML'24 ML4LMS](https://openreview.net/forum?id=7UvbaTrNbP)
