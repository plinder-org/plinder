# Python API tutorial

**Work in progress**

## Installation

### Getting Plinder

Due to dependencies that are not installable via `pip`, `plinder` is currently not
available at PyPI.
You can download the official
[_GitHub_ repository](https://github.com/plinder-org/plinder/)
instead, for example via `git`.

```console
$ git clone https://github.com/plinder-org/plinder.git
```

### Creating the Conda environment

The most convenient way to install the aforementioned extra dependencies is a _Conda_
environment.
If you have not _Conda_ installed yet, we recommend its installation via
[miniforge](https://github.com/conda-forge/miniforge).
Afterwards the environment can be created from the `environment.yml` in the local
repository clone.

:::{note}
We currently only support a Linux environment.
`plinder` uses `openstructure` for some of its functionality and is available from the
`aivant` conda channel using `conda install aivant::openstructure`, but it is only built
targeting Linux architectures.
For Windows and MacOS users, please see the relevant
[_Docker_](#docker-target) resources.
:::

```console
$ mamba env create -f environment.yml
$ mamba activate plinder
```

### Installing plinder

Now `plinder` can be installed into the created environment:

```console
$ pip install .
```

(docker-target)=
### Alternative: Using a Docker container

We also publish the `plinder` project as a docker container as alternative to the
_Conda_-based installation, to ensure the highest level of compatibility with
non-Linux platforms.
See the relevant docker resources here for more details:

- `docker-compose.yml`: defines a `base` image, the `plinder` "app" and a `test`
  container
- `dockerfiles/base/`: contains the files for the `base` image
- `dockerfiles/main/`: contains the files for the `plinder` "app" image

## Overview

The user-faced subpackage of `plinder` is `plinder.core`.

## Example Usage

### Configure dataset environment variable
> We need to set environment variables to point to the release and iteration of choice. For the sake of demonstartion, this will be set to point to a smaller toy example dataset, which are `PLINDER_RELEASE=2024-06` and `PLINDER_ITERATION=toy`.

>NOTE!: the version used for the preprint is `PLINDER_RELEASE=2024-04` and `PLINDER_ITERATION=v1`, while the current version with updated annotations to be used for the MLSB challenge is`PLINDER_RELEASE=2024-06` and `PLINDER_ITERATION=v2`. You could do this directly from a shell terminal with `export PLINDER_RELEASE=2024-04 && export PLINDER_ITERATION=toy` or do it with python with `os.environ`.

```python
from __future__ import annotations
import os
from pathlib import Path

release = "2024-04"
iteration = "tutorial"
os.environ["PLINDER_RELEASE"] = release
os.environ["PLINDER_ITERATION"] = iteration
os.environ["PLINDER_REPO"] =  str(Path.home()/"plinder-org/plinder")
os.environ["PLINDER_LOCAL_DIR"] =  str(Path.home()/".local/share/plinder")
version = f"{release}/{iteration}"
```

### Get config

```python
import plinder.core.utils.config

cfg = plinder.core.utils.config.get_config()
print(f"""
local cache directory: {cfg.data.plinder_dir}
remote data directory: {cfg.data.plinder_remote}
""")
data_dir = Path(cfg.data.plinder_dir)
```

### Inspect annotation file
```python
from plinder.core.index.utils import get_plindex
get_plindex(cfg=cfg)
```

### Query protein similarity
```python
from plinder.core.scores.protein import query_protein_similarity
query_protein_similarity(
    search_db="apo",
    filters=[("similarity", ">", "50")]
)
```

### Load plinder system data onject from system_id
```python
from plinder.core.system.utils import load_systems
load_systems(
    system_ids=["5dax__1__1.A__1.B_1.C_1.D", "3to9__2__2.A__2.B"]
)
```

### Inspect manifest table
```python
from plinder.core.index.utils import get_manifest
get_manifest(cfg=cfg)
```

### Load entries
```python
from plinder.core.index.utils import load_entries
entries = load_entries(
    two_char_codes="lp"

)
entries["2lp7"].keys()
```
