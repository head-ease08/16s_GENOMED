#!/usr/bin/env bash
# bwa-mem2 alignment of DADA2-filtered reads against the SILVA reference,
# for alignment-based QC metrics (mapped %, insert size, coverage, GC, dup %).
# Fused pipe, no intermediate BAM: bwa-mem2 mem -> fixmate -> sort -> markdup.
#
# Usage: align_silva.sh <reference> <r1.fq.gz> <r2.fq.gz> <output.bam> <sample> <threads>
set -euo pipefail

REFERENCE="$1"
R1="$2"
R2="$3"
OUTPUT_BAM="$4"
SAMPLE="$5"
THREADS="$6"

TMP_PREFIX="$(dirname "$OUTPUT_BAM")/${SAMPLE}.tmpsort"
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

bwa-mem2 mem -M -t "$THREADS" -R "$RG" "$REFERENCE" "$R1" "$R2" \
    | samtools fixmate -m -u -@ 2 - - \
    | samtools sort -u -@ 4 -m 2G -T "$TMP_PREFIX" - \
    | samtools markdup -@ 4 --write-index - "${OUTPUT_BAM}##idx##${OUTPUT_BAM}.bai"
