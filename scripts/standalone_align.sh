#!/usr/bin/env bash
# Standalone bwa-mem2 alignment, no Snakemake / no fixed data/raw layout.
# Give it any R1/R2 fastq pair + any reference FASTA, get a sorted,
# duplicate-marked, indexed BAM back. Builds the bwa-mem2 index itself
# if it's missing (one-time cost, cached next to the reference).
#
# Usage:
#   standalone_align.sh <reference.fasta> <r1.fq.gz> <r2.fq.gz> <output.bam> [sample] [threads]
#
# Example:
#   scripts/standalone_align.sh references/silva_for_alignment.fasta \
#       ~/data/S1_R1.fq.gz ~/data/S1_R2.fq.gz out/S1.bam S1 8
set -euo pipefail

REFERENCE=${1:?Usage: standalone_align.sh <reference.fasta> <r1.fq.gz> <r2.fq.gz> <output.bam> [sample] [threads]}
R1=${2:?missing r1.fq.gz}
R2=${3:?missing r2.fq.gz}
OUTPUT_BAM=${4:?missing output.bam}
SAMPLE=${5:-$(basename "$R1" | sed -E 's/[._]R?1.*$//')}
THREADS=${6:-$(nproc)}

for f in "$REFERENCE" "$R1" "$R2"; do
    [ -f "$f" ] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

command -v bwa-mem2 >/dev/null || { echo "ERROR: bwa-mem2 not on PATH" >&2; exit 1; }
command -v samtools  >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }

mkdir -p "$(dirname "$OUTPUT_BAM")"

# build the bwa-mem2 index next to the reference if it's not there yet
if [ ! -f "${REFERENCE}.bwt.2bit.64" ]; then
    echo "[standalone_align] no index found for $REFERENCE, building one (one-off)..." >&2
    bwa-mem2 index "$REFERENCE"
fi

TMP_PREFIX="$(dirname "$OUTPUT_BAM")/${SAMPLE}.tmpsort"
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

echo "[standalone_align] aligning $SAMPLE ($R1 + $R2) against $REFERENCE -> $OUTPUT_BAM" >&2

bwa-mem2 mem -M -t "$THREADS" -R "$RG" "$REFERENCE" "$R1" "$R2" \
    | samtools fixmate -m -u -@ 2 - - \
    | samtools sort -u -@ 4 -m 2G -T "$TMP_PREFIX" - \
    | samtools markdup -@ 4 --write-index - "${OUTPUT_BAM}##idx##${OUTPUT_BAM}.bai"

echo "[standalone_align] done. flagstat:" >&2
samtools flagstat "$OUTPUT_BAM" >&2
