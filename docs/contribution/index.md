# Contributor guide

The PLINDER project is a community effort, launched by the University of Basel,
SIB Swiss Institute of Bioinformatics, Proxima (formerly VantAI), NVIDIA, MIT CSAIL,
and will be regularly updated.
We highly welcome contributions!

This guide introduces how to maintain and improve `plinder` as a developer.

## Code organization

The code is split into four subpackages:

- `plinder.core`: Provides core data structures for interacting with
  and loading the dataset.
  Parts of it are exposed as the [public API](/api/index).
- `plinder.data`: Contains core code for generating the PLINDER dataset.
- `plinder.eval`: Provides an evaluation harness that compares predicted and
  ground-truth structures and produces leaderboard-ready results.
- `plinder.methods`: Implements methods in the leaderboard that leverage
  PLINDER primitives for training and running.

:::{toctree}
:maxdepth: 1
:hidden:

pipeline
development
documentation
release
:::
