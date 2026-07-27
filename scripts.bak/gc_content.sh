#!/usr/bin/env bash
# Usage: gc_content.sh <bam_dir> <output_tsv>
BAM_DIR=${1:?Usage: gc_content.sh <bam_dir> <output_tsv>}
OUTPUT=${2:?Usage: gc_content.sh <bam_dir> <output_tsv>}

process_one() {
    local bam=$1
    local sample=$(basename "$bam" .bam)

    # GC% для R1 (флаг 64 = first in pair)
    gc_r1=$(samtools view -F 4 -f 64 -@ 2 "$bam" | head -1000000 | awk '
        {
            seq = $10
            gc = gsub(/[GCgc]/, "", seq)
            total_chars = length($10)
            if (total_chars > 0) {
                sum += gc / total_chars * 100
                n++
            }
        }
        END { if (n > 0) printf "%.2f", sum / n; else print "NA" }')

    # GC% для R2 (флаг 128 = second in pair)
    gc_r2=$(samtools view -F 4 -f 128 -@ 2 "$bam" | head -1000000 | awk '
        {
            seq = $10
            gc = gsub(/[GCgc]/, "", seq)
            total_chars = length($10)
            if (total_chars > 0) {
                sum += gc / total_chars * 100
                n++
            }
        }
        END { if (n > 0) printf "%.2f", sum / n; else print "NA" }')

    gc_mean=$(awk -v r1="$gc_r1" -v r2="$gc_r2" 'BEGIN {
        if (r1 == "NA" || r2 == "NA") print "NA"
        else printf "%.2f", (r1 + r2) / 2
    }')

    echo -e "${sample}\t${gc_r1}\t${gc_r2}\t${gc_mean}"
}

export -f process_one

echo -e "Sample\tGC_R1_%\tGC_R2_%\tGC_Mean_%" > "$OUTPUT"
ls "${BAM_DIR}"/*.bam | parallel -j 8 process_one {} >> "$OUTPUT"
