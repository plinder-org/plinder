# Documentation

The *PLINDER* documentation resides in the `docs` directory.
The documents are written in *Markdown* files rendered with [*Sphinx*](https://www.sphinx-doc.org) and [*MyST*](https://myst-parser.readthedocs.io).
The [Python API Tutorial](/tutorial/api) uses a [*Jupyter*](https://jupyter.org/)
notebook rendered with [*MyST-NB*](https://myst-nb.readthedocs.io)

## Building

Building the documentation requires some extra requirements that are specified
in the `pyproject.toml`.

```console
$ pip install -e ".[docs]"
```

To build the documentation run

```console
$ sphinx-build docs <output directory>
```
