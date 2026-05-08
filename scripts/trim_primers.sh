#!/usr/bin/env bash
set -euo pipefail

FWD="$1"
REV="$2"
FWD_RC="$3"
REV_RC="$4"
R1_IN="$5"
R2_IN="$6"
R1_OUT="$7"
R2_OUT="$8"
ERROR_RATE="$9"
MIN_LEN="${10}"
THREADS="${11}"

cutadapt \
    -g "file:$FWD" \
    -a "file:$REV_RC" \
    -G "file:$REV" \
    -A "file:$FWD_RC" \
    -n 2 \
    -e "$ERROR_RATE" \
    --discard-untrimmed \
    --minimum-length "$MIN_LEN" \
    -j "$THREADS" \
    -o "$R1_OUT" \
    -p "$R2_OUT" \
    "$R1_IN" "$R2_IN"
