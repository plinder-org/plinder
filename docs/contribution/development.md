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

The data generation pipeline (`plinder.data`) and evaluation require a few tools
that are only available via _Conda_ (mmseqs2, foldseek, OpenStructure).
If you have not _Conda_ installed yet, we recommend its installation via
[miniforge](https://github.com/conda-forge/miniforge).

```console
$ mamba env create -f environment.yml
$ mamba activate plinder
```

### Installing `plinder`

All Python dependencies are declared in `pyproject.toml` and locked in `uv.lock`.
With [uv](https://docs.astral.sh/uv/), install `plinder` in editable mode with the
complete `dev` dependency group (all extras used in CI, including CPU-only pytorch
on Linux and the git-only pipeline packages) into the active Conda environment:

```console
$ UV_PROJECT_ENVIRONMENT="$CONDA_PREFIX" uv sync --inexact
```

CI runs the same command with `--locked`. Without the Conda environment, a plain `uv sync`
installs into `.venv`; tests that need the Conda-only tools will then fail.

With pip, the base install covers data generation and the core library:

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

The full data pipeline also needs the git-only packages in the `dev`
dependency group, which only `uv sync` installs.
:::

### Enabling Pre-commit hooks

Please install pre-commit hooks, that will run the same code quality checks as the CI:

```console
$ pre-commit install
```

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
