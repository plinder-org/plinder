# Development

## Installation

For the scope of the entire PLINDER workflow non-Python dependencies are required

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
