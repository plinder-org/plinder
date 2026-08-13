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
SEARCH_DB=${PLINDER_SEARCH_DB:-holo}
ENV_ROOT=${PLINDER_ENV_ROOT:?"set PLINDER_ENV_ROOT to the ingest Conda environment"}
REPO_ROOT=${PLINDER_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
LOCAL_SCRATCH=${PLINDER_SCORE_SCRATCH:-${SLURM_TMPDIR:-${TMPDIR:-/scratch/${USER}/plinder-score-${SLURM_JOB_ID:-local}}}}
TASK_SCRATCH=${LOCAL_SCRATCH}/${MODE}-${BATCH_INDEX}
SCRATCH_MONITOR_PID=

cleanup_task_scratch() {
    if [[ -n "${SCRATCH_MONITOR_PID}" ]]; then
        kill "${SCRATCH_MONITOR_PID}" 2>/dev/null || true
        wait "${SCRATCH_MONITOR_PID}" 2>/dev/null || true
    fi
    if [[ -n "${TASK_SCRATCH}" && "${TASK_SCRATCH}" == "${LOCAL_SCRATCH}"/* ]]; then
        rm -rf -- "${TASK_SCRATCH}"
    fi
}
trap cleanup_task_scratch EXIT

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
echo "host=$(hostname) job_id=${JOB_ID} mode=${MODE} search_db=${SEARCH_DB} batch_index=${BATCH_INDEX} batch_size=${BATCH_SIZE} scratch=${TASK_SCRATCH}"
echo "PLINDER offline mode is enabled; all sources and reference data must exist"

monitor_scratch_usage() {
    local output=${DATA_DIR}/metrics/protein-${MODE}-${JOB_ID}-${BATCH_INDEX}.scratch.tsv
    local interval=${PLINDER_SCRATCH_MONITOR_INTERVAL_SECONDS:-5}
    local scratch_kib
    printf 'unix_seconds\tscratch_kib\n' > "${output}"
    while [[ -d "${TASK_SCRATCH}" ]]; do
        scratch_kib=$(du -sk "${TASK_SCRATCH}" 2>/dev/null | cut -f1) || scratch_kib=0
        printf '%(%s)T\t%s\n' -1 "${scratch_kib:-0}" >> "${output}"
        sleep "${interval}"
    done
}

if [[ "${PLINDER_MONITOR_SCRATCH:-0}" == "1" ]]; then
    monitor_scratch_usage &
    SCRATCH_MONITOR_PID=$!
fi

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
            --search-db "${SEARCH_DB}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
        )
        ;;
    refresh-alignment-lookup)
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
            --max-query-protein-chains "${PLINDER_MAX_QUERY_PROTEIN_CHAINS:-30}"
            --max-query-proper-ligand-chains "${PLINDER_MAX_QUERY_PROPER_LIGAND_CHAINS:-30}"
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
    plan-interface-scores)
        ARGS+=(--batch-size "${BATCH_SIZE}")
        ;;
    plan-interface-score-repair)
        ARGS+=(--batch-size "${BATCH_SIZE}")
        if [[ -n "${PLINDER_INTERFACE_REPAIR_QUERY_MANIFEST:-}" ]]; then
            ARGS+=(
                --query-manifest
                "${PLINDER_INTERFACE_REPAIR_QUERY_MANIFEST}"
            )
        fi
        ;;
    plan-ligand-pocket-scores)
        ARGS+=(
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}"
            --max-query-protein-chains "${PLINDER_MAX_QUERY_PROTEIN_CHAINS:-30}"
            --max-query-proper-ligand-chains "${PLINDER_MAX_QUERY_PROPER_LIGAND_CHAINS:-30}"
        )
        ;;
    plan-score-repair)
        : "${PLINDER_AFFECTED_PDB_MANIFEST:?set PLINDER_AFFECTED_PDB_MANIFEST}"
        ARGS+=(
            --affected-manifest "${PLINDER_AFFECTED_PDB_MANIFEST}"
            --batch-size "${BATCH_SIZE}"
            --target-batch-size "${PLINDER_TARGET_REPAIR_BATCH_SIZE:-100}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}"
        )
        if [[ -n "${PLINDER_ADDITIONAL_FULL_QUERY_MANIFEST:-}" ]]; then
            ARGS+=(
                --additional-full-query-manifest
                "${PLINDER_ADDITIONAL_FULL_QUERY_MANIFEST}"
            )
        fi
        if [[ -n "${PLINDER_SCORE_REPAIR_MANIFEST:-}" ]]; then
            ARGS+=(--output "${PLINDER_SCORE_REPAIR_MANIFEST}")
        fi
        ;;
    plan-bounded-score-repair)
        : "${PLINDER_BOUNDED_QUERY_MANIFEST:?set PLINDER_BOUNDED_QUERY_MANIFEST}"
        ARGS+=(
            --pdb-manifest "${PLINDER_BOUNDED_QUERY_MANIFEST}"
            --batch-size "${BATCH_SIZE}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}"
        )
        ;;
    plan-score-repair-ligand-3d)
        : "${PLINDER_SCORE_REPAIR_MANIFEST:?set PLINDER_SCORE_REPAIR_MANIFEST}"
        ARGS+=(
            --repair-manifest "${PLINDER_SCORE_REPAIR_MANIFEST}"
            --batch-size "${BATCH_SIZE}"
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-8GB}"
        )
        ;;
    finalize-score-repair-queries)
        : "${PLINDER_SCORE_REPAIR_MANIFEST:?set PLINDER_SCORE_REPAIR_MANIFEST}"
        ARGS+=(
            --repair-manifest "${PLINDER_SCORE_REPAIR_MANIFEST}"
            --max-new-drops "${PLINDER_MAX_NEW_REPAIR_DROPS:-2000}"
        )
        ;;
    validate-score-repair)
        : "${PLINDER_SCORE_REPAIR_MANIFEST:?set PLINDER_SCORE_REPAIR_MANIFEST}"
        ARGS+=(--repair-manifest "${PLINDER_SCORE_REPAIR_MANIFEST}")
        ;;
    plan-clusters)
        ARGS+=(
            --source-batch-size "${BATCH_SIZE}"
            --cover-batch-size "${PLINDER_COVER_BATCH_SIZE:-1}"
            --symmetric-bucket-count "${PLINDER_SYMMETRIC_EDGE_BUCKET_COUNT:-64}"
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
    finalize-interface-scores)
        ARGS+=(
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}"
        )
        if [[ -n "${PLINDER_INTERFACE_SCORE_OUTPUT:-}" ]]; then
            ARGS+=(--output "${PLINDER_INTERFACE_SCORE_OUTPUT}")
        fi
        ;;
    finalize-interface-score-repair)
        ARGS+=(
            --threads "${SLURM_CPUS_PER_TASK:-1}"
            --scratch-dir "${TASK_SCRATCH}"
            --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}"
        )
        ;;
    finalize-alignments|finalize-ligands|finalize-scores|finalize-ligand-3d-retries|finalize-index|merge-components|cluster-stats)
        ;;
    search-foldseek|search-mmseqs|map|pack-ligands|score|score-pdbs|repair-score-pdbs|repair-score-target-remainders|repair-pack-ligands|repair-candidates|repair-ligand-3d|repair-merge-ligand-3d|repair-score-shards|repair-sucos-shards|collate-ligand-3d-candidates|score-ligand-3d|score-ligand-3d-retry|collate-ligand-3d|merge-ligand-3d|export-sucos-shards|score-ligand-pocket-shards|materialize-ligand-3d-candidates|score-interface-shards|score-interface-repair-batches|collate-alignments|collate-score-partitions|symmetric-edge-fragments|symmetric-edge-shards|component-reductions|set-covers|directed-covers)
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
        if [[ "${MODE}" == "search-foldseek" || "${MODE}" == "search-mmseqs" || "${MODE}" == "map" || "${MODE}" == "score" || "${MODE}" == "score-pdbs" ]]; then
            ARGS+=(--search-db "${SEARCH_DB}")
        fi
        if [[ "${MODE}" == "score-pdbs" ]]; then
            : "${PLINDER_SCORE_PDB_MANIFEST:?set PLINDER_SCORE_PDB_MANIFEST}"
            ARGS+=(--pdb-manifest "${PLINDER_SCORE_PDB_MANIFEST}")
        fi
        if [[ "${MODE}" == repair-* ]]; then
            : "${PLINDER_SCORE_REPAIR_MANIFEST:?set PLINDER_SCORE_REPAIR_MANIFEST}"
            ARGS+=(--repair-manifest "${PLINDER_SCORE_REPAIR_MANIFEST}")
        fi
        if [[ "${MODE}" == "export-sucos-shards" || "${MODE}" == "repair-sucos-shards" ]]; then
            ARGS+=(
                --output-dir "${PLINDER_SUCOS_EXPORT_DIR:-${DATA_DIR}/exports/all_sucos_shape_pocket_qcov}"
                --memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-16GB}"
            )
        fi
        if [[ "${MODE}" == "score-interface-shards" || "${MODE}" == "score-interface-repair-batches" || "${MODE}" == "score-ligand-pocket-shards" || "${MODE}" == "materialize-ligand-3d-candidates" ]]; then
            ARGS+=(--memory-limit "${PLINDER_DUCKDB_MEMORY_LIMIT:-32GB}")
        fi
        if [[ "${MODE}" == "score-interface-shards" || "${MODE}" == "score-interface-repair-batches" ]]; then
            ARGS+=(--side-coverage-buckets "${PLINDER_INTERFACE_SIDE_COVERAGE_BUCKETS:-32}")
        fi
        if [[ "${MODE}" == "score-interface-repair-batches" ]]; then
            ARGS+=(
                --repair-batches-per-task
                "${PLINDER_INTERFACE_REPAIR_BATCHES_PER_TASK:-1}"
            )
            if [[ "${FORCE:-0}" == "1" ]]; then
                ARGS+=(--force)
            fi
        fi
        if [[ "${MODE}" == "score-ligand-pocket-shards" && -n "${PLINDER_SCORE_SHARDS:-}" ]]; then
            IFS=',' read -r -a SCORE_SHARDS <<< "${PLINDER_SCORE_SHARDS}"
            for SHARD in "${SCORE_SHARDS[@]}"; do
                ARGS+=(--shard "${SHARD}")
            done
        fi
        ;;
    *)
        echo "unsupported mode: ${MODE}" >&2
        exit 2
        ;;
esac

if [[ "${MODE}" == "plan-clusters" || "${MODE}" == "symmetric-edge-fragments" || "${MODE}" == "symmetric-edge-shards" || "${MODE}" == "component-reductions" || "${MODE}" == "merge-components" || "${MODE}" == "set-covers" || "${MODE}" == "directed-covers" || "${MODE}" == "cluster-stats" ]]; then
    ARGS+=(--entity-type "${PLINDER_CLUSTER_ENTITY_TYPE:-ligand}")
fi

if [[ "${MODE}" == "plan-clusters" || "${MODE}" == "component-reductions" || "${MODE}" == "merge-components" || "${MODE}" == "set-covers" || "${MODE}" == "directed-covers" || "${MODE}" == "cluster-stats" ]]; then
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

set +e
/usr/bin/time -v \
    -o "${DATA_DIR}/metrics/protein-${MODE}-${JOB_ID}-${BATCH_INDEX}.time-v.txt" \
    python -m plinder.data.pipeline.score "${ARGS[@]}"
COMMAND_STATUS=$?
set -e

if [[ -n "${SCRATCH_MONITOR_PID}" ]]; then
    kill "${SCRATCH_MONITOR_PID}" 2>/dev/null || true
    wait "${SCRATCH_MONITOR_PID}" 2>/dev/null || true
    SCRATCH_MONITOR_PID=
fi
exit "${COMMAND_STATUS}"
