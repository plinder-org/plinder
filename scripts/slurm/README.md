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

## Collate the V3 annotation index

After entry ingest is complete, freeze the two-character code list, inventory
the exact per-entry outputs in parallel, and atomically merge the inventories.
The plan fails if any materialized entry is missing an annotation, chain,
biological-assembly-chain, source, or ligand Parquet:

```bash
sbatch \
  --output="${OUTPUT_ROOT}/logs/collate-plan-start-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/collate_v3_shards.sbatch plan-start "${OUTPUT_ROOT}"

# Read code_count from index/.staging/v3_collation/plan-build.json. With a
# batch size of four, LAST_PLAN_BATCH_INDEX is ceil(code_count / 4) - 1.
sbatch \
  --array=0-LAST_PLAN_BATCH_INDEX --cpus-per-task=4 --mem=8G \
  --output="${OUTPUT_ROOT}/logs/collate-plan-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/collate_v3_shards.sbatch plan-shard "${OUTPUT_ROOT}" 4

sbatch \
  --dependency=afterok:PLAN_ARRAY_JOB_ID \
  --output="${OUTPUT_ROOT}/logs/collate-plan-finish-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/collate_v3_shards.sbatch plan-finish "${OUTPUT_ROOT}"
```

When augmenting an existing ligand release with interface-only ingest, set
`PLINDER_COLLATE_INTERFACES_ONLY=true` on `plan-start`. This preserves the
installed `index/annotation_table.parquet` byte-for-byte and collates only the
shared chain, entry, source, and interface tables. The setting is frozen in the
plan, so later stages do not need the environment variable.

Read `index/.staging/v3_collation/plan.json` after that job succeeds. With a
batch size of four, set `LAST_CODE_BATCH_INDEX` to
`ceil(code_count / 4) - 1`, then submit the unthrottled shard array:

```bash
sbatch \
  --array=0-LAST_CODE_BATCH_INDEX \
  --output="${OUTPUT_ROOT}/logs/collate-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/collate_v3_shards.sbatch shard "${OUTPUT_ROOT}" 4
```

Each shard is atomic and resumable. Once the array succeeds, merge and validate
the planned shards with a larger single task:

```bash
sbatch \
  --qos=6hours \
  --cpus-per-task=4 --mem=48G \
  --output="${OUTPUT_ROOT}/logs/collate-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT,PLINDER_COLLATE_MEMORY_LIMIT=40GB \
  scripts/slurm/collate_v3_shards.sbatch finalize "${OUTPUT_ROOT}"
```

Finalization writes the local `index/*.parquet` files only after validating row
counts, keys, cross-table references, and ligand scoreability. Each shard
verifies its raw inputs before producing a frozen output; finalization validates
and merges those outputs without rescanning every per-entry source file. It
performs no upload or external release operation.

## Build publishable protein-search shards

Protein scoring uses the union of protein chains that are holo ligand receptors
or members of `index/interface_annotation_table.parquet`; unrelated chains
retained only for biological-unit context are not searched. Freeze that query
universe before starting any arrays:

```bash
python -m plinder.data.pipeline.score plan "${OUTPUT_ROOT}" --max-seqs 10000
```

The resulting `manifests/protein_scoring_plan.json` gives `query_count` and
`shard_count`. Both Foldseek and MMseqs searches are filtered by E-value and
coverage, with the minimum sequence identity explicitly set to zero.

Cache `pdb_seqres.txt.gz` on an online node and set `PLINDER_SEQRES_PATH` to it.
Run `createdb` on a host with enough memory and fast access to the managed PDB
tree; this reads the original NextGen divided mmCIF archive directly. The
command derives a TSV of only the protein-containing PDB entries in the frozen
plan, avoiding a recursive walk over the full archive:

```bash
python -m plinder.data.pipeline.score create-dbs "${OUTPUT_ROOT}" \
  --cif-root "${PLINDER_PDB_NEXTGEN_ROOT}" \
  --seqres-path "${PLINDER_SEQRES_PATH}" \
  --threads 4
```

If the login host is memory-limited, run the same resumable operation on Slurm.
Each backend is built and validated under node-local `/scratch`, then staged and
swapped into shared storage. A retry reuses either backend whose installed
`createdb` output already completed:

```bash
sbatch \
  --qos=6hours --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/protein-create-dbs-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT,PLINDER_PDB_NEXTGEN_ROOT,PLINDER_SEQRES_PATH \
  scripts/slurm/score_v3.sbatch create-dbs "${OUTPUT_ROOT}"
```

On Slurm, select the V3 protein chains, cluster the Foldseek and MMseqs targets
at exactly 100% identity and 100% bidirectional coverage, and index the
representative targets. This job is offline and uses node-local `/scratch` for
all clustering and index working files:

```bash
sbatch \
  --qos=6hours --cpus-per-task=16 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/protein-make-sub-dbs-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch make-sub-dbs "${OUTPUT_ROOT}"
```

The 10,000-hit search cap therefore applies to unique representatives rather
than being consumed by identical chains. Foldseek's cluster search and the
MMseqs expansion workflow both realign expanded members before conversion, so
the resulting Parquets still contain chain-level scores and identifiers.

Foldseek and MMseqs use separate arrays so their query-batch sizes can be tuned
independently. Pilot before submitting the complete arrays. For a Foldseek
batch size of 5,000, `LAST_FOLDSEEK_INDEX` is
`ceil(query_count / 5000) - 1`:

```bash
sbatch \
  --qos=6hours --array=0-LAST_FOLDSEEK_INDEX \
  --cpus-per-task=128 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/protein-foldseek-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch search-foldseek "${OUTPUT_ROOT}" 5000

sbatch \
  --qos=30min --array=0-LAST_MMSEQS_INDEX \
  --cpus-per-task=128 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/protein-mmseqs-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch search-mmseqs "${OUTPUT_ROOT}" 5000
```

Foldseek uses the longer QoS up front because its expanded hit sets and raw
checkpoints are larger. For MMseqs, resubmit only incomplete batches with
`--qos=6hours`; per-query raw checkpoints make those retries skip outputs
already installed by the first pass.

The Foldseek release output stores one scalar LDDT value per hit, not
per-residue `lddtfull`. Full aligned strings remain generation intermediates;
release shards retain only scalar search fields and compact pocket-relevant
residue-number/identity mappings needed to reconstruct system and ligand scores.

Mapping is separate from derived scoring so completed searches remain usable
when downstream score generation needs retries. Each task owns one PDB
two-character query shard, writes temporary per-query mappings only on
node-local scratch, and atomically publishes one Parquet per available backend.
With one shard per task, `LAST_SHARD_INDEX` is `shard_count - 1`:

```bash
sbatch \
  --qos=30min --array=0-LAST_SHARD_INDEX --cpus-per-task=2 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/protein-map-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch map "${OUTPUT_ROOT}" 1

# Pack each canonical ASU ligand once for bulk scoring and distribution.
sbatch \
  --qos=30min --array=0-LAST_LIGAND_PACK_INDEX --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/ligand-pack-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch pack-ligands "${OUTPUT_ROOT}" 50

sbatch \
  --output="${OUTPUT_ROOT}/logs/ligand-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch finalize-ligands "${OUTPUT_ROOT}"

sbatch \
  --output="${OUTPUT_ROOT}/logs/protein-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch finalize-alignments "${OUTPUT_ROOT}"
```

The same mapped protein alignments also drive protein-interface scoring. Plan
one two-character query shard per array task, run the shards independently,
and then publish the compact all-vs-all table:

```bash
sbatch \
  --qos=30min --cpus-per-task=1 --mem=8G \
  --output="${OUTPUT_ROOT}/logs/interface-score-plan-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch plan-interface-scores "${OUTPUT_ROOT}" 1

# Read batch_count from manifests/interface_scoring_plan.json, then:
sbatch \
  --qos=30min --array=0-LAST_INTERFACE_INDEX --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/interface-score-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch score-interface-shards "${OUTPUT_ROOT}" 1

sbatch \
  --qos=6hours --cpus-per-task=8 --mem=64G \
  --output="${OUTPUT_ROOT}/logs/interface-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT,PLINDER_DUCKDB_MEMORY_LIMIT=56GB \
  scripts/slurm/score_v3.sbatch finalize-interface-scores "${OUTPUT_ROOT}"
```

Each winning direct or swapped chain assignment retains `iface1_qcov` and
`iface2_qcov` before multiplication. They are the coverages of the query
interface's canonical first and second chain, respectively, under that winning
assignment. The compact release file
`exports/all_interface_qcov.parquet` contains the query and target interface
IDs, both side coverages, and the final directional 0--100 similarity. Positive
scores below the lowest clustering threshold are retained.

Derived per-ligand system scores are generation intermediates used for graph
clustering, not release artifacts. Before scattering them, estimate work from
the compact alignment hits and greedily spread expensive queries across fixed
50-query batches within each two-character alignment shard:

```bash
sbatch \
  --qos=6hours --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/protein-score-plan-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch plan-score-batches "${OUTPUT_ROOT}" 50

# Read score_batch_count from manifests/protein_scoring_plan.json, then:
sbatch \
  --qos=30min --array=0-LAST_SCORE_INDEX --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/protein-score-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch score "${OUTPUT_ROOT}" 50
```

Query-system safety caps default to 30 protein receptor chains and 30 proper
ligand chains. Larger holo systems remain available as targets. Override the
query caps with `PLINDER_MAX_QUERY_PROTEIN_CHAINS` and
`PLINDER_MAX_QUERY_PROPER_LIGAND_CHAINS` when planning.

This first pass calculates only protein and pocket scores. For every proper,
3D-scoreable ligand pair with positive pocket query coverage it also records a
full-precision candidate row; it does not load an SDF or run shape alignment.
Retry only any batches that exceed the 30-minute QoS with `--qos=6hours`.
If a batch exceeds memory, freeze its still-missing PDB IDs into a Parquet with
a unique `pdb_id` column and retry them independently:

```bash
export PLINDER_SCORE_PDB_MANIFEST="${OUTPUT_ROOT}/manifests/protein_scoring_retry.parquet"
sbatch \
  --qos=6hours --array=0-LAST_RETRY_INDEX --cpus-per-task=4 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/protein-score-pdb-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT,PLINDER_SCORE_PDB_MANIFEST \
  scripts/slurm/score_v3.sbatch score-pdbs "${OUTPUT_ROOT}" 1
```

Each array task then owns one PDB and still reuses completed score checkpoints.

After every protein-score batch has completed, first consolidate the per-PDB
candidate checkpoints into two-character query shards. With a collation batch
size of four, use `ceil(shard_count / 4) - 1` as
`LAST_CANDIDATE_SHARD_INDEX`. This prevents the global planner from opening one
small Parquet for every PDB entry:

```bash
sbatch \
  --qos=30min --array=0-LAST_CANDIDATE_SHARD_INDEX --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-candidates-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch collate-ligand-3d-candidates "${OUTPUT_ROOT}" 4
```

Then deduplicate the sharded candidate rows into canonical ASU ligand pairs,
retain every pair with positive pocket coverage, and stream deterministic
query-local batches ordered by ligand size:

```bash
sbatch \
  --qos=6hours --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-plan-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch plan-ligand-3d "${OUTPUT_ROOT}" 30000

# Read ligand_3d_batch_count from protein_scoring_plan.json, then:
sbatch \
  --qos=30min --array=0-LAST_LIGAND_3D_INDEX --cpus-per-task=1 --mem=16G \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-score-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch score-ligand-3d "${OUTPUT_ROOT}" 30000
```

The pair scorer bulk-loads only the required canonical SDF records from the
packed Parquets and calculates shape, color, and SuCOS once per unique directed
ASU ligand pair. If any coarse batches fail, split only those batches into
smaller, resumable retry units and require every retry before continuing:

```bash
sbatch \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-retry-plan-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch plan-ligand-3d-retries "${OUTPUT_ROOT}" 500

# Read ligand_3d_retry_batch_count from protein_scoring_plan.json, then:
sbatch \
  --qos=6hours --array=0-LAST_LIGAND_3D_RETRY_INDEX \
  --cpus-per-task=1 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-retry-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch score-ligand-3d-retry "${OUTPUT_ROOT}" 500

sbatch \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-retry-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch finalize-ligand-3d-retries "${OUTPUT_ROOT}"
```

Repartition the complete balanced outputs by query PDB shard before
expanding them back to system-level rows. With a collation batch size of four,
`LAST_LIGAND_3D_SHARD_INDEX` is
`ceil(ligand_3d_query_shard_count / 4) - 1`:

```bash
sbatch \
  --qos=30min --array=0-LAST_LIGAND_3D_SHARD_INDEX --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-collate-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch collate-ligand-3d "${OUTPUT_ROOT}" 4

sbatch \
  --qos=30min --array=0-LAST_LIGAND_MERGE_INDEX --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/ligand-3d-merge-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch merge-ligand-3d "${OUTPUT_ROOT}" 50

sbatch \
  --qos=6hours --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/protein-score-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch finalize-scores "${OUTPUT_ROOT}"
```

The merge computes `sucos_shape_pocket_qcov` from full-precision pocket
coverage before converting retained similarities to the 0–100 integer schema.

The complete `sucos_shape_pocket_qcov` edge table is a release artifact. Export
it directly from the full-precision candidate and canonical-pair shards so that
scores below the 30% clustering threshold are retained. Use the same query-shard
array bound as ligand-3D collation, then validate and concatenate the shards:

```bash
sbatch \
  --qos=30min --array=0-LAST_LIGAND_3D_SHARD_INDEX \
  --cpus-per-task=4 --mem=32G \
  --output="${OUTPUT_ROOT}/logs/sucos-export-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch export-sucos-shards "${OUTPUT_ROOT}" 4

sbatch \
  --qos=6hours --cpus-per-task=8 --mem=64G \
  --output="${OUTPUT_ROOT}/logs/sucos-export-finalize-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch finalize-sucos-export "${OUTPUT_ROOT}"
```

The default final path is
`exports/all_sucos_shape_pocket_qcov.parquet`; set
`PLINDER_SUCOS_EXPORT_DIR` and `PLINDER_SUCOS_EXPORT_OUTPUT` to override the
temporary shard directory or final release location.

Targeted collation repairs deliberately leave `index/collation.json` in
`requires_downstream_repair` state and remove the stale nonredundant index.
After repairing fingerprints and affected score shards, rerun the clustering
and final-index sequence below. Finalization now rejects cluster artifacts that
do not exactly cover the current eligible ligand universe and marks the repair
complete only after recreating `annotation_table_nonredundant.parquet`.

Build ligand-level reciprocal-minimum components and greedy centroid
communities only after score finalization.
The planning job reports exact array counts for the default metrics and the
30, 50, 70, 90, and 100 thresholds:

```bash
sbatch \
  --qos=30min --cpus-per-task=1 --mem=8G \
  --output="${OUTPUT_ROOT}/logs/cluster-plan-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch plan-clusters "${OUTPUT_ROOT}" 20
```

First read `symmetric_fragment_batch_count` from that log. Each task copies 20
regular score shards to node-local scratch and hash-partitions canonical ligand
pairs while retaining both directional maxima:

```bash
sbatch \
  --qos=6hours --array=0-LAST_SYMMETRIC_FRAGMENT_INDEX \
  --cpus-per-task=4 --mem=64G \
  --output="${OUTPUT_ROOT}/logs/symmetric-fragment-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch symmetric-edge-fragments "${OUTPUT_ROOT}" 1

sbatch \
  --qos=6hours --array=0-LAST_SYMMETRIC_EDGE_INDEX \
  --cpus-per-task=4 --mem=64G \
  --output="${OUTPUT_ROOT}/logs/symmetric-edge-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch symmetric-edge-shards "${OUTPUT_ROOT}" 1
```

The second array takes the minimum of the two directional maxima. Pairs with
only one stored direction retain that direction but have a null reciprocal
score: public reciprocal-minimum clustering ignores them, while the optional
directed cover can still use them. Then read
`component_reduction_batch_count` from the plan and reduce one compact edge
shard per task:

```bash
sbatch \
  --qos=6hours --array=0-LAST_COMPONENT_REDUCTION_INDEX \
  --cpus-per-task=4 --mem=64G \
  --output="${OUTPUT_ROOT}/logs/component-reduction-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch component-reductions "${OUTPUT_ROOT}" 1

sbatch \
  --qos=6hours --cpus-per-task=8 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/component-merge-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch merge-components "${OUTPUT_ROOT}"
```

The merge publishes exact connected components of the reciprocal-minimum graph
without retaining a full in-memory 30%-threshold graph. It also creates
any-direction connectivity partitions used by the optional directed sampling
cover. Finally, read `community_batch_count` from the plan log and scatter one
metric/threshold per task:

```bash
COMMUNITY_JOB=$(sbatch --parsable \
  --qos=6hours --array=0-LAST_COMMUNITY_INDEX \
  --cpus-per-task=8 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/community-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch communities "${OUTPUT_ROOT}" 1)
```

Communities use deterministic greedy centroid cover followed by reassignment
to the highest-scoring selected centroid. Every member meets the threshold in
both directions to its centroid; two non-centroid members need not meet it
directly. For Tanimoto, the symmetric score is used directly.

Run `directed-covers` with the
`directed_cover_batch_count` from the plan. A centroid covers query ligand `Q`
when `score(Q -> centroid)` meets the threshold. These labels are required by
final index enrichment and the detailed assignments also support training-set
sampling:

```bash
DIRECTED_COVER_JOB=$(sbatch --parsable \
  --qos=6hours --array=0-LAST_DIRECTED_COVER_INDEX \
  --cpus-per-task=8 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/directed-cover-%A-%a.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch directed-covers "${OUTPUT_ROOT}" 1)
```
Set comma-delimited `PLINDER_CLUSTER_METRICS` or
`PLINDER_CLUSTER_THRESHOLDS` consistently on every clustering command to run a
subset. Because Slurm itself treats commas inside `--export` as separators,
define comma-delimited values in the submitting shell and export their names:

```bash
export PLINDER_CLUSTER_METRICS='pocket_qcov,pli_qcov'
export PLINDER_CLUSTER_THRESHOLDS='100,90,70,50,30'
sbatch \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT,PLINDER_CLUSTER_METRICS,PLINDER_CLUSTER_THRESHOLDS \
  scripts/slurm/score_v3.sbatch plan-clusters "${OUTPUT_ROOT}" 1
```

Protein-interface clustering uses the same commands and thresholds in a
separate namespace. Set these variables for every command in the sequence:

```bash
export PLINDER_CLUSTER_ENTITY_TYPE=interface
export PLINDER_CLUSTER_METRICS=interface_qcov
export PLINDER_CLUSTER_THRESHOLDS='100,90,70,50,30'
```

The interface plan reads `interface_scores/shard=*.parquet`. Reciprocal
components and communities are written below `interface_clusters/`; the
directional centroid cover is written below `interface_sampling/`. Final index
enrichment adds `interface_qcov__THRESHOLD__component`, `__community`, and
`__directed_set_cover` columns to `index/interface_annotation_table.parquet`.
The ligand and interface plans and artifacts never share cache paths.
Run both entity sequences before `finalize-index`; finalization rejects a
non-empty interface annotation table when its interface clusters are absent.

After component, community, and directed-cover branches complete, validate every
published artifact and write `ligand_clusters/stats.parquet` (or
`interface_clusters/stats.parquet`) plus `stats.json`.
This gate checks artifact coverage, duplicate and null labels, consistent ligand
counts, and monotonic component counts across thresholds. Community and cover
counts are reported but are not required to be monotonic:

```bash
STATS_JOB=$(sbatch --parsable \
  --dependency="afterok:${COMMUNITY_JOB}:${DIRECTED_COVER_JOB}" \
  --qos=6hours --cpus-per-task=8 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/cluster-stats-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch cluster-stats "${OUTPUT_ROOT}")
```

Only after that audit succeeds, merge the ligand-level annotations and cluster
IDs into the index. Cluster files are loaded incrementally and progress with an
ETA is logged, avoiding a giant concatenated long-form table:

```bash
sbatch \
  --dependency="afterok:${STATS_JOB}" \
  --qos=6hours --cpus-per-task=8 --mem=128G \
  --output="${OUTPUT_ROOT}/logs/finalize-index-%j.out" \
  --export=ALL,PLINDER_ENV_ROOT,PLINDER_REPO_ROOT \
  scripts/slurm/score_v3.sbatch finalize-index "${OUTPUT_ROOT}"
```

The release scoring artifacts are
`alignments/search_db=holo/alignment_type=*/shard=*.parquet`, their validated
manifest, and `exports/all_sucos_shape_pocket_qcov.parquet`. Per-PDB raw search
files and the other derived score datasets are generation intermediates; mapped
per-PDB files exist only transiently on node-local scratch.
