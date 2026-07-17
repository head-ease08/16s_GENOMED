#!/usr/bin/env bash
# Collapse a multi-species 16S reference (RDP/SILVA-style, thousands of
# near-identical strain records) into cluster-representative sequences via
# vsearch, then build a bwa-mem2 index on the result. Fixes two problems
# alignment-based insert-size QC hits on raw multi-strain databases:
#   - mates scattering onto different near-duplicate contigs (TLEN=0)
#   - forced softclip when the one hand-picked reference diverges from the
#     sample's actual organism in the amplicon's variable region
# Clustering keeps enough distinct references for real diversity while
# merging near-duplicate strains into one contig, so mates consistently
# land on the same contig with low divergence. Default 97% is the classic
# 16S species-level OTU threshold -- tested against 99%/97%/95%/90% on a
# synthetic near-duplicate reference: 99% was too strict to actually merge
# same-species near-duplicates (mate-scatter persisted), 97% and below
# collapsed them fully and fixed it.
#
# Usage: cluster_reference.sh <reference.fasta> <out_clustered.fasta> [identity]
#   identity default: 0.97
set -euo pipefail

REFERENCE=${1:?Usage: cluster_reference.sh <reference.fasta> <out_clustered.fasta> [identity]}
OUT=${2:?missing out_clustered.fasta}
IDENTITY=${3:-0.97}

command -v vsearch >/dev/null || { echo "ERROR: vsearch not on PATH (conda install -c bioconda vsearch)" >&2; exit 1; }
command -v bwa-mem2 >/dev/null || { echo "ERROR: bwa-mem2 not on PATH" >&2; exit 1; }

[ -f "$REFERENCE" ] || { echo "ERROR: file not found: $REFERENCE" >&2; exit 1; }

echo "[cluster_reference] clustering $REFERENCE at ${IDENTITY} identity..." >&2
vsearch --cluster_fast "$REFERENCE" \
    --id "$IDENTITY" \
    --centroids "$OUT" \
    --threads "$(nproc)" \
    --fasta_width 0

echo "[cluster_reference] $(grep -c '^>' "$REFERENCE") input seqs -> $(grep -c '^>' "$OUT") cluster representatives" >&2

echo "[cluster_reference] building bwa-mem2 index..." >&2
bwa-mem2 index "$OUT"

echo "[cluster_reference] done: $OUT (+ index)" >&2
