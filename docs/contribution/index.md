# Contributor guide

The PLINDER project is a community effort, launched by the University of Basel,
SIB Swiss Institute of Bioinformatics, VantAI, NVIDIA, MIT CSAIL, and will be regularly
updated.
We highly welcome contributions!

This guide gives an introduction about how to maintain and improve `plinder` as
developer

# Code organization

This code is split into 4 sub-packages

- `plinder.core`: Provides core data structures for interacting with
  and loading the dataset.
  Parts of it are exposed as the [public API](/api/index).
- `plinder.data`: Contains core code for generating the PLINDER dataset.
- `plinder.eval`: Offers evaluation harness for the dataset that takes as an input
  predicted and ground truth structures in a pre-determined folder structure an
  returns a leaderboard-ready set of entries.
  Parts of it are user-faced via [CLI scripts](/evaluation).
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
