#!/bin/bash
# AnDOM2 — cleanup old temporary search files and expired batch jobs
# Run via cron: 0 3 * * * /var/www/AnDOM2/cleanup.sh >> /var/log/andom2_cleanup.log 2>&1

set -e

OUTPUT_DIR="${ANDOM_OUTPUT_DIR:-/app/output}"
TMP_DIR="$OUTPUT_DIR/tmp"
RESULTS_DIR="$OUTPUT_DIR/results"

# ── configurable retention ────────────────────────────────────────────────────
TMP_DAYS="${ANDOM_TMP_DAYS:-1}"        # MMseqs2/Foldseek tmp dirs: 1 day
JOB_DAYS="${ANDOM_JOB_DAYS:-7}"        # Batch job results: 7 days
BENCHMARK_DAYS="${ANDOM_BM_DAYS:-30}"  # Benchmark JSON files: 30 days

echo "=== AnDOM2 cleanup $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "Retention: tmp=${TMP_DAYS}d  jobs=${JOB_DAYS}d  benchmarks=${BENCHMARK_DAYS}d"

# ── 1. MMseqs2 / Foldseek tmp dirs ───────────────────────────────────────────
if [ -d "$TMP_DIR" ]; then
    before=$(find "$TMP_DIR" -mindepth 1 -maxdepth 1 | wc -l)
    find "$TMP_DIR" -mindepth 1 -maxdepth 1 \
        ! -name "latest" \
        -mtime +${TMP_DAYS} \
        -exec rm -rf {} + 2>/dev/null || true
    after=$(find "$TMP_DIR" -mindepth 1 -maxdepth 1 | wc -l)
    echo "tmp/: removed $((before - after)) dirs (${before} → ${after})"
else
    echo "tmp/: directory not found, skipping"
fi

# ── 2. Batch job results ──────────────────────────────────────────────────────
if [ -d "$RESULTS_DIR" ]; then
    removed=0
    for job_dir in "$RESULTS_DIR"/*/; do
        [ -d "$job_dir" ] || continue
        job_json="$job_dir/job.json"
        # only remove completed/failed/cancelled jobs older than JOB_DAYS
        if [ -f "$job_json" ]; then
            status=$(python3 -c "import json,sys; d=json.load(open('$job_json')); print(d.get('status',''))" 2>/dev/null || echo "")
            if [[ "$status" =~ ^(done|failed|cancelled)$ ]]; then
                if find "$job_dir" -maxdepth 0 -mtime +${JOB_DAYS} | grep -q .; then
                    rm -rf "$job_dir"
                    removed=$((removed + 1))
                fi
            fi
        fi
    done
    echo "results/: removed $removed expired job dirs"
fi

# ── 3. Benchmark JSON files ───────────────────────────────────────────────────
if [ -d "$RESULTS_DIR" ]; then
    bm_removed=$(find "$RESULTS_DIR" -maxdepth 1 -name "benchmark_*.json" \
        -mtime +${BENCHMARK_DAYS} -delete -print | wc -l)
    echo "results/: removed $bm_removed old benchmark files"
fi

# ── 4. Orphaned queryDB files in project root ────────────────────────────────
# These are created by MMseqs2 in the working dir — move them to tmp
PROJECT_ROOT="${ANDOM_PROJECT_ROOT:-/app}"
orphans=0
for f in queryDB queryDB.dbtype queryDB_h queryDB_h.dbtype queryDB_h.index \
          queryDB.index queryDB.lookup queryDB.source \
          resultDB resultDB.dbtype resultDB.index \
          query_input.fasta query_struct.pdb \
          seq_results.tsv struct_results.tsv; do
    if [ -f "$PROJECT_ROOT/$f" ]; then
        rm -f "$PROJECT_ROOT/$f"
        orphans=$((orphans + 1))
    fi
done
[ $orphans -gt 0 ] && echo "root/: removed $orphans orphaned MMseqs2 files"

echo "=== cleanup done ==="
