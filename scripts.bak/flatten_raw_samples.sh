#!/usr/bin/env bash
# Adapter for the "one subfolder per sample" raw layout (RAW_DIR/{sample}/*_R1*.fq.gz,
# same as Snakefile's discover_samples()) -> flat {sample}_R1.fq.gz / _R2.fq.gz dir
# that demux_regions.sh / batch_align.sh / standalone_align.sh expect. Symlinks,
# doesn't copy -- raw fastqs can be huge.
#
# Usage: flatten_raw_samples.sh <raw_dir> <flat_out_dir>
set -euo pipefail

RAW_DIR=${1:?Usage: flatten_raw_samples.sh <raw_dir> <flat_out_dir>}
OUT_DIR=${2:?missing flat_out_dir}

mkdir -p "$OUT_DIR"

n=0
for d in "$RAW_DIR"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")

    shopt -s nullglob nocaseglob
    r1_candidates=("$d"*[Rr]1*.f*q*.gz)
    r2_candidates=("$d"*[Rr]2*.f*q*.gz)
    shopt -u nullglob nocaseglob

    if [ ${#r1_candidates[@]} -eq 0 ] || [ ${#r2_candidates[@]} -eq 0 ]; then
        echo "WARNING: $sample - missing R1 or R2, skipping" >&2
        continue
    fi

    ln -sf "$(realpath "${r1_candidates[0]}")" "$OUT_DIR/${sample}_R1.fq.gz"
    ln -sf "$(realpath "${r2_candidates[0]}")" "$OUT_DIR/${sample}_R2.fq.gz"
    n=$((n+1))
done

echo "Linked $n samples into $OUT_DIR" >&2
