#!/usr/bin/env bash
# Reference-free insert size via fastp's R1/R2 overlap analysis. Only valid
# for fragments shorter than R1_len + R2_len (fragment must overlap) -- use
# this for merge-able regions (V1_V2, V9), not the long non-overlapping ones.
#
# Usage: insert_size_fastp.sh <reads_dir> <out_dir> [jobs] [threads_per_job]
#   reads_dir: flat dir of R1/R2 pairs (e.g. one region's demux output)
set -euo pipefail

READS_DIR=${1:?Usage: insert_size_fastp.sh <reads_dir> <out_dir> [jobs] [threads_per_job]}
OUT_DIR=${2:?missing out_dir}
JOBS=${3:-4}
THREADS=${4:-2}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT_DIR"

shopt -s nullglob nocaseglob
r1_files=("$READS_DIR"/*R1*.f*q*.gz)
shopt -u nullglob nocaseglob
if [ ${#r1_files[@]} -eq 0 ]; then
    echo "No *R1*.f*q*.gz files found in $READS_DIR" >&2
    exit 1
fi

extract_one() {
    r1="$1"; out_dir="$2"; threads="$3"
    base=$(basename "$r1")
    r2_base=$(echo "$base" | sed -E 's/R1/R2/; s/r1/r2/')
    r2="$(dirname "$r1")/$r2_base"
    sample=$(echo "$base" | sed -E 's/[._-]?[Rr]1.*$//')
    [ -f "$r2" ] || { echo "WARNING: no R2 match for $base -- skipping" >&2; return 0; }

    json="$out_dir/${sample}.fastp.json"
    fastp -i "$r1" -I "$r2" -w "$threads" \
        --json "$json" --html /dev/null \
        > "$out_dir/${sample}.fastp.log" 2>&1

    python3 "$SCRIPT_DIR/fastp_json_to_insert_sizes.py" "$json" "$out_dir/${sample}.insert_sizes.txt"
    n=$(wc -l < "$out_dir/${sample}.insert_sizes.txt")
    [ "$n" -eq 0 ] && echo "WARNING: $sample - 0 insert sizes from fastp overlap" >&2
}
export -f extract_one
export SCRIPT_DIR

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${r1_files[@]}" | parallel -j "$JOBS" extract_one {} "$OUT_DIR" "$THREADS"
else
    for r1 in "${r1_files[@]}"; do
        extract_one "$r1" "$OUT_DIR" "$THREADS" &
        while [ "$(jobs -r -p | wc -l)" -ge "$JOBS" ]; do wait -n; done
    done
    wait
fi

Rscript "$SCRIPT_DIR/plot_insert_size_hist.R" "$OUT_DIR"
echo "Done. Histograms in $OUT_DIR" >&2
