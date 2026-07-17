#!/usr/bin/env bash
# Extract per-sample TLEN (insert size) values from BAMs and render histograms.
# Usage: insert_size_histograms.sh [bam_dir] [out_dir]
#   bam_dir default: results/qc_align/bam
#   out_dir default: results/qc_align/insert_size_hist
set -euo pipefail

BAM_DIR=${1:-results/qc_align/bam}
OUT_DIR=${2:-results/qc_align/insert_size_hist}
JOBS=${JOBS:-$(nproc)}

mkdir -p "$OUT_DIR"

shopt -s nullglob
bams=("$BAM_DIR"/*.bam)
shopt -u nullglob
if [ ${#bams[@]} -eq 0 ]; then
    echo "No BAMs found in $BAM_DIR" >&2
    exit 1
fi

extract_one() {
    bam="$1"
    out_dir="$2"
    sample=$(basename "$bam" .bam)
    # first-in-pair, both mates mapped, not secondary/supplementary.
    # Deliberately NOT filtering on the "proper pair" flag (-f 2): short
    # 16S amplicons routinely produce dovetailed/overlapping mates that
    # bwa-mem2 won't mark proper-paired even though both reads map fine,
    # which silently zeroed out this file for amplicon data.
    out="$out_dir/${sample}.insert_sizes.txt"
    samtools view -f 64 -F 2828 "$bam" \
        | awk '{v=$9; if (v<0) v=-v; if (v>0) print v}' > "$out"
    n=$(wc -l < "$out")
    [ "$n" -eq 0 ] && echo "WARNING: $sample - 0 insert sizes extracted (check samtools flagstat)" >&2
}
export -f extract_one

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${bams[@]}" | parallel -j "$JOBS" extract_one {} "$OUT_DIR"
else
    for bam in "${bams[@]}"; do
        extract_one "$bam" "$OUT_DIR" &
        while [ "$(jobs -r -p | wc -l)" -ge "$JOBS" ]; do wait -n; done
    done
    wait
fi

Rscript "$(dirname "$0")/plot_insert_size_hist.R" "$OUT_DIR"

echo "Histograms written to $OUT_DIR"
