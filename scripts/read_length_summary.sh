#!/usr/bin/env bash
# Average read length after DADA2 filtering, R1+R2 combined.
# Usage: read_length_summary.sh <qc_dir> <output_tsv> <sample1> [sample2 ...]
set -euo pipefail

QC_DIR=${1:?Usage: read_length_summary.sh <qc_dir> <output_tsv> <sample...>}
OUTPUT=${2:?Usage: read_length_summary.sh <qc_dir> <output_tsv> <sample...>}
shift 2

echo -e "Sample\tAverage_read_length" > "$OUTPUT"
for sample in "$@"; do
    avg=$(zcat "${QC_DIR}/${sample}_R1.fq.gz" "${QC_DIR}/${sample}_R2.fq.gz" 2>/dev/null \
        | awk 'NR%4==2 {sum+=length($0); n++} END {if(n>0) printf "%.1f", sum/n; else print "NA"}')
    [ -z "$avg" ] && avg="NA"
    echo -e "${sample}\t${avg}" >> "$OUTPUT"
done
