#!/usr/bin/env bash
# Inverse of flatten_raw_samples.sh: turns flat {sample}_R1.fq.gz / _R2.fq.gz
# files into the "one subfolder per sample" layout (RAW_DIR/{sample}/*_R1*.fq.gz)
# that the Snakefile's discover_samples() requires. Symlinks, doesn't copy.
#
# Usage: nest_raw_samples.sh <flat_dir> <nested_out_dir>
set -euo pipefail

FLAT_DIR=${1:?Usage: nest_raw_samples.sh <flat_dir> <nested_out_dir>}
OUT_DIR=${2:?missing nested_out_dir}

mkdir -p "$OUT_DIR"

shopt -s nullglob nocaseglob
r1_files=("$FLAT_DIR"/*[Rr]1*.f*q*.gz)
shopt -u nullglob nocaseglob

n=0
for r1 in "${r1_files[@]}"; do
    base=$(basename "$r1")
    r2_base=$(echo "$base" | sed -E 's/[Rr]1/R2/')
    r2="$FLAT_DIR/$r2_base"
    sample=$(echo "$base" | sed -E 's/[._-]?[Rr]1.*$//')

    if [ ! -f "$r2" ]; then
        echo "WARNING: $sample - no matching R2 for $base, skipping" >&2
        continue
    fi

    sample_dir="$OUT_DIR/$sample"
    mkdir -p "$sample_dir"
    ln -sf "$(realpath "$r1")" "$sample_dir/${base}"
    ln -sf "$(realpath "$r2")" "$sample_dir/${r2_base}"
    n=$((n+1))
done

echo "Nested $n samples into $OUT_DIR" >&2
