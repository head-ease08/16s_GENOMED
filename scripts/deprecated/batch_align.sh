#!/usr/bin/env bash
# Batch-align every R1/R2 pair found flat in one directory (no per-sample
# subfolders) against a reference, via standalone_align.sh. Pairs files by
# stripping the R1/R2 token from the filename.
#
# Usage: batch_align.sh <reads_dir> <reference.fasta> <out_bam_dir> [jobs] [threads_per_job]
#   jobs default: 4 parallel alignments
#   threads_per_job default: nproc / jobs
set -euo pipefail

READS_DIR=${1:?Usage: batch_align.sh <reads_dir> <reference.fasta> <out_bam_dir> [jobs] [threads_per_job]}
REFERENCE=${2:?missing reference.fasta}
OUT_DIR=${3:?missing out_bam_dir}
JOBS=${4:-4}
THREADS=${5:-$(( $(nproc) / JOBS > 0 ? $(nproc) / JOBS : 1 ))}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT_DIR"

shopt -s nullglob nocaseglob
r1_files=("$READS_DIR"/*R1*.f*q*.gz)
shopt -u nullglob nocaseglob
if [ ${#r1_files[@]} -eq 0 ]; then
    echo "No *R1*.f*q*.gz files found in $READS_DIR" >&2
    exit 1
fi

align_one() {
    r1="$1"; reference="$2"; out_dir="$3"; threads="$4"
    base=$(basename "$r1")
    r2_base=$(echo "$base" | sed -E 's/R1/R2/; s/r1/r2/')
    r2="$(dirname "$r1")/$r2_base"
    sample=$(echo "$base" | sed -E 's/[._-]?[Rr]1.*$//')

    if [ ! -f "$r2" ]; then
        echo "WARNING: no R2 match for $base (looked for $r2_base) -- skipping" >&2
        return 0
    fi
    "$SCRIPT_DIR/standalone_align.sh" "$reference" "$r1" "$r2" \
        "$out_dir/${sample}.bam" "$sample" "$threads"
}
export -f align_one
export SCRIPT_DIR

echo "Found ${#r1_files[@]} R1 files in $READS_DIR, aligning with $JOBS parallel jobs x $THREADS threads each" >&2

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${r1_files[@]}" | parallel -j "$JOBS" align_one {} "$REFERENCE" "$OUT_DIR" "$THREADS"
else
    for r1 in "${r1_files[@]}"; do
        align_one "$r1" "$REFERENCE" "$OUT_DIR" "$THREADS" &
        while [ "$(jobs -r -p | wc -l)" -ge "$JOBS" ]; do wait -n; done
    done
    wait
fi

echo "Done. BAMs in $OUT_DIR" >&2
