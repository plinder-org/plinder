#!/usr/bin/env bash
#SBATCH --job-name=plinder-bindingdb
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

set -euo pipefail
umask 002

OUTPUT_ROOT=${1:?"usage: sbatch transform_bindingdb.sbatch OUTPUT_ROOT"}
ENV_ROOT=${PLINDER_ENV_ROOT:?"set PLINDER_ENV_ROOT to the ingest Conda environment"}
REPO_ROOT=${PLINDER_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
RAW_BINDINGDB="${OUTPUT_ROOT}/dbs/affinity/BindingDB_All.tsv"

if [[ ! -s "${RAW_BINDINGDB}" ]]; then
    echo "missing pre-downloaded BindingDB TSV: ${RAW_BINDINGDB}" >&2
    exit 1
fi

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
echo "host=$(hostname) job_id=${SLURM_JOB_ID:-none} output_root=${OUTPUT_ROOT}"
echo "PLINDER offline mode is enabled; transforming the existing BindingDB TSV"

/usr/bin/time -v \
    -o "${OUTPUT_ROOT}/metrics/transform-bindingdb.time-v.txt" \
    python -c \
    'from pathlib import Path; from plinder.data.pipeline.io import download_affinity_data; download_affinity_data(data_dir=Path(__import__("sys").argv[1]), force_update=False)' \
    "${OUTPUT_ROOT}"

ls -lh "${OUTPUT_ROOT}/dbs/affinity/affinity.json"
