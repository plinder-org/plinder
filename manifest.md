plinder.data versions
---------------------

The latest version of this `manifest.md` is available [here](https://storage.googleapis.com/plinder/manifest.md).

The `plinder.data` dataset itself comes in releases based on when the seed data was synced from the RCSB.
It is additionally versioned by the iteration of the pipeline that produced the dataset.
The numbers reported in the manuscripts for ML4LMS and NeurIPS can be derived from files contained
within the `2024-04/v1` version of the dataset at `gs://plinder/2024-04/v1`.

Here is a high level summary of the pertinent auxiliary assets you can peruse prior to committing
to downloading an entire dataset (note that a full dataset is large, >~1TB):

- `2024-04/v1`
  * `readme` -- [`gs://plinder/2024-04/v1/README.md`](https://storage.googleapis.com/plinder/2024-04/v1/README.md)
  * `source` -- [`gs://plinder/2024-04/v1/plinder.tar.gz`](https://storage.googleapis.com/plinder/2024-04/v1/plinder.tar.gz)


Each `readme` contains details of the corresponding contents of the dataset and each `source` artifact
contains the full source code used to generate and interact with the dataset.
