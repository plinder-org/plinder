#!/usr/bin/env bash
#SBATCH --job-name=plinder-batch
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail
umask 002

MANIFEST=${1:?"usage: sbatch ingest_pdb_batch.sbatch MANIFEST OUTPUT_ROOT BATCH_SIZE"}
OUTPUT_ROOT=${2:?"usage: sbatch ingest_pdb_batch.sbatch MANIFEST OUTPUT_ROOT BATCH_SIZE"}
BATCH_SIZE=${3:?"usage: sbatch ingest_pdb_batch.sbatch MANIFEST OUTPUT_ROOT BATCH_SIZE"}
BATCH_INDEX=${SLURM_ARRAY_TASK_ID:-0}
JOB_ID=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}
ENV_ROOT=${PLINDER_ENV_ROOT:?"set PLINDER_ENV_ROOT to the ingest Conda environment"}
: "${PLINDER_PDB_NEXTGEN_ROOT:?set PLINDER_PDB_NEXTGEN_ROOT}"
: "${PLINDER_VALIDATION_ROOT:?set PLINDER_VALIDATION_ROOT}"
REPO_ROOT=${PLINDER_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}

export PATH="${ENV_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${ENV_ROOT}/lib:${LD_LIBRARY_PATH:-}"
if [[ -n "${REPO_ROOT}" && -d "${REPO_ROOT}/src/plinder" ]]; then
    export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"
fi
export PLINDER_OFFLINE=true
export PLINDER_OFFLINE_MODE=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p "${OUTPUT_ROOT}/metrics"
echo "host=$(hostname) job_id=${JOB_ID} batch_index=${BATCH_INDEX} batch_size=${BATCH_SIZE}"
echo "PLINDER offline mode is enabled; source CIFs and reference data must exist"

ARGS=(
    "${MANIFEST}"
    "${OUTPUT_ROOT}"
    --batch-size "${BATCH_SIZE}"
    --batch-index "${BATCH_INDEX}"
)
if [[ "${FORCE:-0}" == "1" ]]; then
    ARGS+=(--force)
fi

/usr/bin/time -v \
    -o "${OUTPUT_ROOT}/metrics/ingest-batch-${JOB_ID}-${BATCH_INDEX}.time-v.txt" \
    python -m plinder.data.pipeline.ingest batch "${ARGS[@]}"
