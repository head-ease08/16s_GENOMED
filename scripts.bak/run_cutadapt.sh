#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="$1"
OUTPUT_DIR="$2"
FWD_FASTA="$3"
REV_FASTA="$4"

THREADS="${THREADS:-4}"
MIN_LEN="${MIN_LEN:-50}"
ERROR_RATE="${ERROR_RATE:-0.1}"

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

revcomp_fasta() {
    awk 'BEGIN{RS=">";ORS=""} NR>1 {
            n=index($0,"\n"); header=substr($0,1,n-1);
            seq=substr($0,n+1); gsub("\n","",seq);
            print ">"header"_RC\n"seq"\n";
         }' "$1" \
    | awk 'BEGIN{
            map["A"]="T"; map["T"]="A"; map["G"]="C"; map["C"]="G";
            map["N"]="N"; map["R"]="Y"; map["Y"]="R"; map["S"]="S";
            map["W"]="W"; map["K"]="M"; map["M"]="K"; map["B"]="V";
            map["V"]="B"; map["D"]="H"; map["H"]="D";
         }
         /^>/ { print; next }
         {
            rc=""; n=length($0);
            for (i=n; i>=1; i--) {
                base=toupper(substr($0,i,1));
                rc = rc (base in map ? map[base] : base);
            }
            print rc;
         }' > "$2"
}

FWD_RC_FASTA="$TMP_DIR/fwd_rc.fasta"
REV_RC_FASTA="$TMP_DIR/rev_rc.fasta"
revcomp_fasta "$FWD_FASTA" "$FWD_RC_FASTA"
revcomp_fasta "$REV_FASTA" "$REV_RC_FASTA"

shopt -s nullglob
R1_FILES=( "$INPUT_DIR"/*_R1*.fq.gz "$INPUT_DIR"/*_R1*.fq )
shopt -u nullglob

SUMMARY="$OUTPUT_DIR/cutadapt_summary.tsv"
printf 'sample\tinput_pairs\toutput_pairs\tkept_pct\n' > "$SUMMARY"

for r1 in "${R1_FILES[@]}"; do
    base="$(basename "$r1")"
    r2="${r1/_R1/_R2}"
    sample="${base%%_R1*}"

    out_r1="$OUTPUT_DIR/${sample}_R1.fq.gz"
    out_r2="$OUTPUT_DIR/${sample}_R2.fq.gz"
    log="$OUTPUT_DIR/logs/${sample}.log"

    cutadapt \
        -g "file:$FWD_FASTA" \
        -a "file:$REV_RC_FASTA" \
        -G "file:$REV_FASTA" \
        -A "file:$FWD_RC_FASTA" \
        -n 2 \
        -e "$ERROR_RATE" \
        --discard-untrimmed \
        --minimum-length "$MIN_LEN" \
        -j "$THREADS" \
        -o "$out_r1" \
        -p "$out_r2" \
        "$r1" "$r2" \
        > "$log" 2>&1

    total_pairs=$(grep -E '^Total read pairs processed' "$log" \
                  | awk -F: '{gsub(/[, ]/,"",$2); print $2}')
    written_pairs=$(grep -E '^Pairs written' "$log" \
                    | awk -F: '{sub(/\(.*/,"",$2); gsub(/[, ]/,"",$2); print $2}')

    kept_pct="NA"
    if [[ -n "${total_pairs:-}" && "${total_pairs:-0}" -gt 0 ]]; then
        kept_pct=$(awk -v a="$written_pairs" -v b="$total_pairs" \
                   'BEGIN{printf "%.2f", 100*a/b}')
    fi

    printf '%s\t%s\t%s\t%s\n' \
        "$sample" "${total_pairs:-NA}" "${written_pairs:-NA}" "$kept_pct" \
        >> "$SUMMARY"
done
