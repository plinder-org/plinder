# Python API tutorial

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