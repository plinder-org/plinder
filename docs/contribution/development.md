# Development

## Installation

For package development, the installation procedure is more complex than the one for the
end user

### Getting `plinder`

For development you need a clone of the official
[_GitHub_ repository](https://github.com/plinder-org/plinder/).

```console
$ git clone https://github.com/plinder-org/plinder.git
```

### Creating the Conda environment

The data generation pipeline (`plinder.data`) requires a few tools that are only
available via _Conda_ (mmseqs2, foldseek, reduce).
If you have not _Conda_ installed yet, we recommend its installation via
[miniforge](https://github.com/conda-forge/miniforge).

```console
$ mamba env create -f environment.yml
$ mamba activate plinder
```

### Installing `plinder`

The base install covers data generation and the core library (numpy 2 compatible):

```console
$ pip install -e ".[dev]"
```

### Evaluation scoring (optional)

`plinder.eval` runs the [OpenStructure](https://openstructure.org/) command-line
actions for ligand and protein-interface evaluation. OpenStructure 2.12.0 or
newer is installed from Bioconda by the repository's `environment.yml`; it is
not a PyPI dependency. Install the optional evaluation dependencies with:

```console
$ pip install -e ".[eval]"
```

:::{note}
The `eval` extra installs PoseBusters. OpenStructure is Conda-only and
is installed by `mamba env create -f environment.yml` above.
Data generation (`plinder.data`) does **not** require OpenStructure and
works with numpy 2.

For the full data pipeline, additional dependencies are needed:

```console
$ pip install -r requirements_data.txt
```

This includes Linux pytorch (for the loader) and pipeline-specific tools.
For Windows and MacOS users, please see the relevant
[_Docker_](#docker-target) resources.
:::

### Enabling Pre-commit hooks

Please install pre-commit hooks, that will run the same code quality checks as the CI:

```console
$ pre-commit install
```

(docker-target)=
### Alternative: Using a Docker container

We also publish the `plinder` project as a
[docker container](https://github.com/plinder-org/plinder/pkgs/container/plinder)
as alternative to the _Conda_-based installation, to ensure the highest level of
compatibility with non-Linux platforms.
See the relevant docker resources here for more details:

- `docker-compose.yml`: defines a `base` image, the `plinder` "app" and a `test`
  container
- `dockerfiles/base/`: contains the files for the `base` image
- `dockerfiles/main/`: contains the files for the `plinder` "app" image

## Testing and linting

`plinder` uses [`tox`](https://tox.wiki) for running tests, type checks and code linting
(with [`ruff`](https://docs.astral.sh/ruff/)).

```console
$ tox -e test
$ tox -e type
$ tox -e lint
```

See `tox.ini` and `.pre-commit-config.yaml` for details.

## Debugging

In order to change log levels in `plinder`, please set:

```console
export PLINDER_LOG_LEVEL=10
```
