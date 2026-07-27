#!/usr/bin/env bash
# Pivot per-region cutadapt logs from demux_regions.sh (or the Snakemake
# demux_region rule) into a PRIMER x SAMPLE table of read pairs kept.
#
# Expects <logs_dir>/<region>/<sample>.log, e.g.:
#   results/region/V1_V2/sample1.log        (demux_regions.sh output dir)
#   logs/demux_region/V1_V2/sample1.log      (Snakemake log dir)
#
# Usage: primer_sample_table.sh <logs_dir> <output_tsv>
set -euo pipefail

LOGS_DIR=${1:?Usage: primer_sample_table.sh <logs_dir> <output_tsv>}
OUTPUT=${2:?missing output_tsv}

regions=()
for d in "$LOGS_DIR"/*/; do
    [ -d "$d" ] || continue
    regions+=("$(basename "$d")")
done
[ ${#regions[@]} -gt 0 ] || { echo "No region subdirectories found in $LOGS_DIR" >&2; exit 1; }

samples=()
for log in "$LOGS_DIR/${regions[0]}"/*.log; do
    [ -f "$log" ] || continue
    samples+=("$(basename "$log" .log)")
done
[ ${#samples[@]} -gt 0 ] || { echo "No .log files found in $LOGS_DIR/${regions[0]}" >&2; exit 1; }

pairs_written() {
    local log="$1"
    [ -f "$log" ] || { echo "NA"; return; }
    grep -E '^Pairs written' "$log" | awk -F: '{sub(/\(.*/,"",$2); gsub(/[, ]/,"",$2); print $2}'
}

{
    printf 'Sample'
    for region in "${regions[@]}"; do printf '\t%s' "$region"; done
    printf '\n'

    for sample in "${samples[@]}"; do
        printf '%s' "$sample"
        for region in "${regions[@]}"; do
            n=$(pairs_written "$LOGS_DIR/$region/${sample}.log")
            printf '\t%s' "${n:-NA}"
        done
        printf '\n'
    done
} > "$OUTPUT"

echo "Done. ${#samples[@]} samples x ${#regions[@]} regions -> $OUTPUT" >&2
