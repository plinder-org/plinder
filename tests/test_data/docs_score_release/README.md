# Documentation score examples

This is a bounded slice of the September 2026 PLINDER release for executing
`similarity_and_representatives.ipynb` and `representatives_and_leakage.ipynb`
in CI and documentation deployment. It contains real published rows, not
synthetic scores. The manifest records the selection, row counts, and SHA-256
of each generated file.

Cluster IDs and representative assignments are copied from the full release.
They are not recalculated on this subset. Pair-score files retain both endpoints
within the example entries (and the ligand cover representative's entry).
The matched-molecular-pair table retains rows with either endpoint in the
example ligand SMILES IDs.

Only these two notebooks use `PLINDER_SCORE_EXAMPLE_DATA_DIR`. Without it,
they query the configured full release. Other guides continue to use the
configured release and their usual structure/prediction downloads. This slice
is not a complete release and must not be used for benchmark construction.

To regenerate from a complete, local September release, run on a compute node:

```bash
python scripts/prepare_docs_score_release.py /path/to/PLINDER-2026-09 /tmp/docs-score-release
```

Review the manifest and copy the generated files into this directory. Then run
`tests/core/test_docs_score_release.py` and execute the two notebooks with
`PLINDER_SCORE_EXAMPLE_DATA_DIR` set to this directory's absolute path.
