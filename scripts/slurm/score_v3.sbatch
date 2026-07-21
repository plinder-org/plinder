#!/usr/bin/env bash
#SBATCH --job-name=plinder-score
#SBATCH --partition=scicore,bigmem
#SBATCH --qos=30min
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail
umask 002

MODE=${1:?"usage: sbatch score_v3.sbatch MODE DATA_DIR [BATCH_SIZE]"}
DATA_DIR=${2:?"usage: sbatch score_v3.sbatch MODE DATA_DIR [BATCH_SIZE]"}
BATCH_SIZE=${3:-1}
BATCH_INDEX=${SLURM_ARRAY_TASK_ID:-0}
JOB_ID=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}
ENV_ROOT=${PLINDER_ENV_ROOT:?"set PLINDER_ENV_ROOT to the ingest Conda environment"}
REPO_ROOT=${PLINDER_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
LOCAL_SCRATCH=${PLINDER_SCORE_SCRATCH:-${SLURM_TMPDIR:-/scratch/${USER}/plinder-score-${SLURM_JOB_ID:-local}}}
TASK_SCRATCH=${LOCAL_SCRATCH}/${MODE}-${BATCH_INDEX}

export PATH="${ENV_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${ENV_ROOT}/lib:${LD_LIBRARY_PATH:-}"
if [[ -n "${REPO_ROOT}" && -d "${REPO_ROOT}/src/plinder" ]]; then
    export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"
fi
export PLINDER_OFFLINE=true
export PLINDER_OFFLINE_MODE=true
export TMPDIR="${TASK_SCRATCH}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p "${DATA_DIR}/metrics" "${TASK_SCRATCH}"
echo "host=$(hostname) job_id=${JOB_ID} mode=${MODE} batch_index=${BATCH_INDEX} batch_size=${BATCH_SIZE} scratch=${TASK_SCRATCH}"
echo "PLINDER offline mode is enabled; all sources and reference data must exist"

CLI_MODE=${MODE}
if [[ "${MODE}" == "search-foldseek" || "${MODE}" == "search-mmseqs" ]]; then
    CLI_MODE=search
fi
ARGS=("${CLI_MODE}" "${DATA_DIR}")
case "${MODE}" in
    plan)
        ARGS+=(--max-seqs "${PLINDER_MAX_SEQS:-10000}")
        ;;
    create-dbs)
        : "${PLINDER_PDB_NEXTGEN_ROOT:?set PLINDER_PDB_NEXTGEN_ROOT}"
        : "${PLINDER_SEQRES_PATH:?set PLINDER_SEQRES_PATH}"
        ARGS+=(
            --cif-root "${PLINDER_PDB_NEXTGEN_ROOT}"
            --seqres-path "${PLINDER_SEQRES_PATH}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}/createdb"
        )
        if [[ "${FORCE:-0}" == "1" ]]; then
            ARGS+=(--force)
        fi
        ;;
    make-sub-dbs)
        ARGS+=(
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
        )
        ;;
    plan-score-batches)
        ARGS+=(
            --batch-size "${BATCH_SIZE}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
        )
        ;;
    plan-ligand-3d)
        ARGS+=(
            --batch-size "${BATCH_SIZE}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-8GB}"
        )
        ;;
    plan-ligand-3d-retries)
        ARGS+=(--batch-size "${BATCH_SIZE}")
        ;;
    plan-clusters)
        ARGS+=(
            --source-batch-size "${BATCH_SIZE}"
            --community-batch-size "${PLINDER_COMMUNITY_BATCH_SIZE:-1}"
        )
        ;;
    drop-score-queries)
        : "${PLINDER_DROPPED_QUERY_MANIFEST:?set PLINDER_DROPPED_QUERY_MANIFEST}"
        ARGS+=(--pdb-manifest "${PLINDER_DROPPED_QUERY_MANIFEST}")
        ;;
    finalize-sucos-export)
        ARGS+=(
            --source-dir "${PLINDER_SUCOS_EXPORT_DIR:-${DATA_DIR}/exports/all_sucos_shape_pocket_qcov}"
            --output "${PLINDER_SUCOS_EXPORT_OUTPUT:-${DATA_DIR}/exports/all_sucos_shape_pocket_qcov.parquet}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}"
        )
        ;;
    finalize-alignments|finalize-ligands|finalize-scores|finalize-ligand-3d-retries|finalize-index|merge-components|cluster-stats)
        ;;
    search-foldseek|search-mmseqs|map|pack-ligands|score|score-pdbs|collate-ligand-3d-candidates|score-ligand-3d|score-ligand-3d-retry|collate-ligand-3d|merge-ligand-3d|export-sucos-shards|collate-alignments|collate-score-partitions|component-reductions|communities)
        ARGS+=(
            --batch-index "${BATCH_INDEX}"
            --batch-size "${BATCH_SIZE}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
        )
        if [[ "${FORCE:-0}" == "1" ]]; then
            ARGS+=(--force)
        fi
        if [[ "${MODE}" == "search-foldseek" || "${MODE}" == "search-mmseqs" ]]; then
            ARGS+=(--alignment-type "${MODE#search-}")
        fi
        if [[ "${MODE}" == "score-pdbs" ]]; then
            : "${PLINDER_SCORE_PDB_MANIFEST:?set PLINDER_SCORE_PDB_MANIFEST}"
            ARGS+=(--pdb-manifest "${PLINDER_SCORE_PDB_MANIFEST}")
        fi
        if [[ "${MODE}" == "export-sucos-shards" ]]; then
            ARGS+=(
                --output-dir "${PLINDER_SUCOS_EXPORT_DIR:-${DATA_DIR}/exports/all_sucos_shape_pocket_qcov}"
                --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-16GB}"
            )
        fi
        ;;
    *)
        echo "unsupported mode: ${MODE}" >&2
        exit 2
        ;;
esac

if [[ "${MODE}" == "plan-clusters" || "${MODE}" == "component-reductions" || "${MODE}" == "merge-components" || "${MODE}" == "communities" || "${MODE}" == "cluster-stats" ]]; then
    if [[ -n "${PLINDER_CLUSTER_METRICS:-}" ]]; then
        IFS=',' read -r -a CLUSTER_METRICS <<< "${PLINDER_CLUSTER_METRICS}"
        for METRIC in "${CLUSTER_METRICS[@]}"; do
            ARGS+=(--metric "${METRIC}")
        done
    fi
    if [[ -n "${PLINDER_CLUSTER_THRESHOLDS:-}" ]]; then
        IFS=',' read -r -a CLUSTER_THRESHOLDS <<< "${PLINDER_CLUSTER_THRESHOLDS}"
        for THRESHOLD in "${CLUSTER_THRESHOLDS[@]}"; do
            ARGS+=(--threshold "${THRESHOLD}")
        done
    fi
fi

/usr/bin/time -v \
    -o "${DATA_DIR}/metrics/protein-${MODE}-${JOB_ID}-${BATCH_INDEX}.time-v.txt" \
    python -m plinder.data.pipeline.score "${ARGS[@]}"
