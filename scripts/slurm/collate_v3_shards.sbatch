#!/usr/bin/env bash
#SBATCH --job-name=plinder-collate
#SBATCH --partition=scicore,bigmem
#SBATCH --qos=30min
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail
umask 002

MODE=${1:-shard}
case "${MODE}" in
    plan|finalize)
        DATA_DIR=${2:?'usage: sbatch collate_v3_shards.sbatch MODE DATA_DIR'}
        BATCH_SIZE=0
        ;;
    shard)
        DATA_DIR=${2:?'usage: sbatch collate_v3_shards.sbatch shard DATA_DIR BATCH_SIZE'}
        BATCH_SIZE=${3:?'usage: sbatch collate_v3_shards.sbatch shard DATA_DIR BATCH_SIZE'}
        ;;
    *)
        # Preserve the original DATA_DIR BATCH_SIZE interface.
        DATA_DIR=${1:?'usage: sbatch collate_v3_shards.sbatch DATA_DIR BATCH_SIZE'}
        BATCH_SIZE=${2:?'usage: sbatch collate_v3_shards.sbatch DATA_DIR BATCH_SIZE'}
        MODE=shard
        ;;
esac
BATCH_INDEX=${SLURM_ARRAY_TASK_ID:-0}
JOB_ID=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}
ENV_ROOT=${PLINDER_ENV_ROOT:?'set PLINDER_ENV_ROOT to the ingest Conda environment'}
REPO_ROOT=${PLINDER_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
SCRATCH_ROOT=${PLINDER_COLLATE_SCRATCH:-${SLURM_TMPDIR:-/scratch/${USER}/plinder-collate-${SLURM_JOB_ID:-local}}}

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

mkdir -p "${DATA_DIR}/metrics" "${SCRATCH_ROOT}/${JOB_ID}-${BATCH_INDEX}"
echo "host=$(hostname) job_id=${JOB_ID} mode=${MODE} batch_index=${BATCH_INDEX} batch_size=${BATCH_SIZE}"

TIME_OUTPUT="${DATA_DIR}/metrics/collate-${MODE}-${JOB_ID}-${BATCH_INDEX}.time-v.txt"
MEMORY_LIMIT=${PLINDER_COLLATE_MEMORY_LIMIT:-7GB}
if [[ "${MODE}" == "plan" ]]; then
    /usr/bin/time -v -o "${TIME_OUTPUT}" \
        python -m plinder.data.pipeline.collate plan "${DATA_DIR}" \
            --threads "${SLURM_CPUS_PER_TASK:-1}"
elif [[ "${MODE}" == "finalize" ]]; then
    /usr/bin/time -v -o "${TIME_OUTPUT}" \
        python -m plinder.data.pipeline.collate finalize "${DATA_DIR}" \
            --threads "${SLURM_CPUS_PER_TASK:-1}" \
            --memory-limit "${MEMORY_LIMIT}" \
            --scratch-dir "${SCRATCH_ROOT}/${JOB_ID}-${BATCH_INDEX}"
else
    /usr/bin/time -v -o "${TIME_OUTPUT}" \
        python -m plinder.data.pipeline.collate shard "${DATA_DIR}" \
            --batch-index "${BATCH_INDEX}" \
            --batch-size "${BATCH_SIZE}" \
            --threads "${SLURM_CPUS_PER_TASK:-1}" \
            --memory-limit "${MEMORY_LIMIT}" \
            --scratch-dir "${SCRATCH_ROOT}/${JOB_ID}-${BATCH_INDEX}"
fi
