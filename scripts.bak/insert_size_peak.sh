#!/usr/bin/env bash
# Insert size peak (mode of the TLEN distribution) per sample, no Picard.
# Usage: insert_size_peak.sh <bam_dir> <output_tsv>
set -euo pipefail

BAM_DIR=${1:?Usage: insert_size_peak.sh <bam_dir> <output_tsv>}
OUTPUT=${2:?Usage: insert_size_peak.sh <bam_dir> <output_tsv>}

echo -e "Sample\tInsert_size_peak" > "$OUTPUT"
for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    peak=$(samtools view -f 2 -F 256 "$bam" \
        | awk '$9>0{print $9}' \
        | sort -n | uniq -c | sort -k1,1nr | head -1 | awk '{print $2}')
    [ -z "$peak" ] && peak="NA"
    echo -e "${sample}\t${peak}" >> "$OUTPUT"
done
