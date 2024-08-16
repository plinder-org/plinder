# Contributor guide

The PLINDER project is a community effort, launched by the University of Basel,
SIB Swiss Institute of Bioinformatics, VantAI, NVIDIA, MIT CSAIL, and will be regularly
updated.
We highly welcome contributions!

## Development installation

```console
$ git clone git@github.com:plinder-org/plinder.git
$ cd plinder
$ mamba env create -f environment.yml
$ mamba activate plinder
$ pip install -e '.[dev]'
```

Please install pre-commit hooks (the same checks will run in CI):

```console
$ pre-commit install
```

## Test suite

Test linting checks:

```console
$ tox -e lint
```

Run typing checks:

```console
$ tox -e type
```

Run the test suite:

```console
$ tox -e test
```

We lint with `ruff`.
See `tox.ini` and `.pre-commit-config.yaml` for details.

## Debugging

In order to change log levels in `plinder`, please set:

```console
export PLINDER_LOG_LEVEL=10
```

:::{todo}
- more details for each command
- add detailed information on other `plinder` subpackages aprt from `core`
- split into multiple documents
:::
