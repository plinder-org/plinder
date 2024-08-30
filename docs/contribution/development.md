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

The `plinder` subpackages beside `plinder.core` require dependencies that are not
installable via `pip`.
The most convenient way to install the aforementioned extra dependencies is a _Conda_
environment.
If you have not _Conda_ installed yet, we recommend its installation via
[miniforge](https://github.com/conda-forge/miniforge).
Afterwards the environment can be created from the `environment.yml` in the local
repository clone.

:::{note}
We currently only support a Linux environment.
`plinder.data` uses a number of dependencies which are not simply pip-installable.
`openstructure` is for some of its functionality and is available from the
`aivant` conda channel using `conda install aivant::openstructure`, but it is only built
targeting Linux architectures. Additionally, `networkit>=11.0`, which at the time of writing,
does not install cleanly on MacOS, along with a number of dependencies which are referenced
by a GitHub link directly, make a pip-installable package problematic. These additional
dependencies can be installed by running:

```console
$ pip install -r requirements_data.txt
```

`plinder.eval` also relies on `openstructure` for metrics
calculations. For Windows and MacOS users, please see the relevant
[_Docker_](#docker-target) resources.
:::

```console
$ mamba env create -f environment.yml
$ mamba activate plinder
```

### Installing `plinder`

Now `plinder` can be installed into the created environment:

```console
$ pip install -e .
```

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
