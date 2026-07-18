# Batched PDB ingest on Slurm

The entry jobs read the managed NextGen mmCIF and validation archives directly.
They do not stage inputs or create symlinks.

Set site-specific locations in the submission environment:

```bash
export PLINDER_ENV_ROOT=/path/to/plinder/environment
export PLINDER_PDB_NEXTGEN_ROOT=/path/to/pdb_nextgen/data/entries/divided
export PLINDER_VALIDATION_ROOT=/path/to/validation_reports
export OUTPUT_ROOT=/path/to/ingest-output
```

Prepare shared reference data once on a node with internet access:

```bash
python scripts/prepare_ingest_reference_data.py "${OUTPUT_ROOT}" --threads 4
```

If the BindingDB TSV was downloaded but its derived `affinity.json` still needs
to be built, it can be transformed on an offline Slurm node:

```bash
sbatch \
  --output="${OUTPUT_ROOT}/logs/transform-bindingdb-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT \
  scripts/slurm/transform_bindingdb.sbatch "${OUTPUT_ROOT}"
```

For a full release, first create a size-balanced manifest. This scans only local
metadata and does not download structure or validation files:

```bash
python -m plinder.data.pipeline.ingest manifest \
  "${PLINDER_PDB_NEXTGEN_ROOT}" \
  "${PLINDER_VALIDATION_ROOT}" \
  "${OUTPUT_ROOT}/manifests/pdb-nextgen-all.txt" \
  --batch-size 100 --threads 8
```

The neighboring `.txt.json` report contains the exact entry and batch counts,
compressed-CIF size distribution, and largest inputs; `.txt.entries.tsv` retains
the full size inventory for sampling and runtime analysis. Add `--check-validation`
only when an up-front XML availability audit is worth the extra metadata I/O.
Ordinary ingest records exact XML availability per entry. Each contiguous 100-line
slice is balanced by compressed CIF bytes. Entries within each slice run smallest
first, so a 30-minute timeout leaves the expensive tail for the retry rather than
discarding time before many small entries are completed.

Submit the complete one-CPU array, replacing `LAST_BATCH_INDEX` with
`batch_count - 1` from the report. No array throttle is imposed here; Slurm can
use all concurrency allowed by available capacity and fair-share:

```bash
sbatch \
  --array=0-LAST_BATCH_INDEX \
  --output="${OUTPUT_ROOT}/logs/ingest-batch-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_PDB_NEXTGEN_ROOT,PLINDER_VALIDATION_ROOT \
  scripts/slurm/ingest_pdb_batch.sbatch \
  "${OUTPUT_ROOT}/manifests/pdb-nextgen-all.txt" "${OUTPUT_ROOT}" 100
```

The batch runner records progress atomically after every PDB ID. Resubmitting the
same array skips entries whose per-entry metrics and outputs are complete. This
makes it safe to retry the full array under the longer QoS; normally only timed-out
or failed entries do work:

```bash
sbatch \
  --qos=6hours \
  --mem=64G \
  --array=0-LAST_BATCH_INDEX \
  --output="${OUTPUT_ROOT}/logs/ingest-batch-retry-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_PDB_NEXTGEN_ROOT,PLINDER_VALIDATION_ROOT \
  scripts/slurm/ingest_pdb_batch.sbatch \
  "${OUTPUT_ROOT}/manifests/pdb-nextgen-all.txt" "${OUTPUT_ROOT}" 100
```

Both array passes force both Plinder offline-mode variables. A missing source CIF
or required shared reference fails the affected entry without a network request.
Per-entry timings are written to sharded `metrics/ingest-one-<pdb_id>.json` files;
whole-process resource use is written once per array task to
`metrics/ingest-batch-<job_id>-<batch_index>.time-v.txt`.
