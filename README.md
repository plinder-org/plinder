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

The `plinder` project is a community effort, launched by the University of Basel, SIB Swiss Institute of Bioinformatics, VantAI, NVIDIA, MIT CSAIL, and will be regularly updated
We highly welcome contributions!
If you find `plinder` useful, please see the citation file for details on how to cite.

To accelerate community adoption, PLINDER will be used as the field’s new Protein-Ligand interaction dataset standard as part of an exciting competition at the upcoming 2024 [Machine Learning in Structural Biology (MLSB)](https://mlsb.io/) Workshop at NeurIPS, one of the fields’ premiere academic gatherings, which will be announced shortly.

# 👨💻 Getting Started

Please use a virtual environment for the `plinder` project.
We recommend the [miniforge](https://github.com/conda-forge/miniforge) environment manager.

**NOTE**: We currently only support a Linux environment. `plinder`
uses `openstructure` for some of its functionality and is available
from the `aivant` conda channel using `conda install aivant::openstructure`, but it is only built targeting Linux architectures.
For MacOS users, please see the relevant [docker](#package-publishing) resources below.

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

# ⬇️ Getting the dataset

Using the `plinder.core` API, you can transparently and lazily
download and interact with most of the components of the dataset.
However, if you prefer to use the dataset directly, you can
fetch it using [`gsutil`](https://cloud.google.com/storage/docs/gsutil).
For portions of the code that still assume
an already available pre-existing dataset, you will need to use `gsutil`.

To download the manifest of available versions:

    gsutil -m cp -r gs://plinder/manifest.md .

Then you can choose to download a particular README version or download the entire dataset with:

    export PLINDER_RELEASE=2024-06
    export PLINDER_ITERATION=v2
    gsutil -m cp -r gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/* ~/.local/share/plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/

View the hierarchy of a particular version of the dataset with:

    gsutil ls gs://plinder/2024-06/v2/

Note that `gsutil` also supports glob patterns like:

    gsutil ls 'gs://plinder/2024-04/v1/fingerprints/*'

This will list the files at the given path in cloud storage.

**NOTE**: We keep the default `PLINDER_RELEASE` and `PLINDER_ITERATION` in the source code up-to-date
with the latest version of the dataset, so if you plan to use a different version, be sure to set
the environment variables accordingly.

## 🔢 Plinder versions

We version the `plinder` dataset with two controls:

- `PLINDER_RELEASE`: the month stamp of the last RCSB sync
- `PLINDER_ITERATION`: value that enables iterative development within a release

We version the `plinder` application using an automated semantic
versioning scheme based on the `git` commit history.
The `plinder.data` package is responsible for generating a dataset
release and the `plinder.core` package makes it easy to interact
with the dataset.

## 🏅 Gold standard benchmark sets

As part of `plinder` resource we also provide train, validation and test splits that are curated to minimize the information leakage based on protein-ligand interaction similarity. In addition, we have prioritized the systems that has a linked experimental `apo` structure or matched molecular series to support realistic inference scenarios for hit discovery and optimization.
Finally, a particular care is taken for test set that is further prioritized to contain high quality structures to provide unambiguous ground-truths for performance benchmarking.

![test_stratification](https://github.com/user-attachments/assets/5bb96534-f939-42b5-bf85-5ac3a71aa324)

Moreover, as we enticipate this resource to be used for benchmarking a wide range of methods, including those simultaneously predicting protein structure (aka. co-folding) or those generating novel ligand structures, we further stratified test (by novel ligand, pocket, protein or all) to cover a wide range of tasks.

Our latest test split contains:

| Novel   |   # of systems | # of high quality |  stratification criteria |
|:--------|---------------:|------------------:|:---------------:|
| pocket  | 5206 | 5203 | PLI shared < 50 _&_  Pocket shared lDDT < 0.5 |
| ligand  | 2395 | 2395 | ECFP4 fingerprint similarity < 0.3 |
| protein |  983 |  983 | Protein Seq. Sim. < 0.3 _&_ Protein lDDT > 0.7 |
| all     |  268 |  268 | all of the above |
| none    |    0 |    0 | none of the above |

## ⚖️ Evaluation harness

See the [`plinder.eval`](docs/eval/README.md) docs for more details.

## 📦 Dataloader

The `plinder.data.loader` package contains a `PyTorch` dataloader for the dataset using the `atom3d` format. It is an example of using the `plinder.core` API
to implement a dataloader, but is not the only way to use the dataset.

**Note**: The dataloader requires both `torch` and `atom3d` to be installed. You use the `[loader]` dependency block when installing `plinder`:

    pip install .[loader]

## ℹ️ Filters & Annotations

See the [`plinder.data`](docs/data/README.md) docs for more details.

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

See the [End-to-end pipeline](docs/data/README.md) description for technical details about the dataset generation.

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

**NOTE**: The CI workflow will use the __most recent__ match in the commit history to make its decision.

# 📃 Publications
Durairaj, Janani, Yusuf Adeshina, Zhonglin Cao, Xuejin Zhang, Vladas Oleinikovas, Thomas Duignan, Zachary McClure, Xavier Robin, Emanuele Rossi, Guoqing Zhou, Srimukh Prasad Veccham, Clemens Isert, Yuxing Peng, Prabindh Sundareson, Mehmet Akdel, Gabriele Corso, Hannes Stärk, Zachary Wayne Carpenter, Michael M. Bronstein, Emine Kucukbenli, Torsten Schwede, Luca Naef. 2024. “PLINDER: The Protein-Ligand Interactions Dataset and Evaluation Resource.”
[bioRxiv](https://doi.org/10.1101/2024.07.17.603955)
[ICML'24 ML4LMS](https://openreview.net/forum?id=7UvbaTrNbP)
