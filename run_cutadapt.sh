#!/usr/bin/env bash
#
# run_cutadapt.sh — remove 16S primers from paired-end Illumina reads.
#
# Usage:
#   ./run_cutadapt.sh <input_dir> <output_dir> <fwd_primers.fasta> <rev_primers.fasta>
#
# Arguments:
#   input_dir         directory with raw paired-end FASTQ files
#                     (expects *_R1*.fastq(.gz) and *_R2*.fastq(.gz))
#   output_dir        directory to write trimmed reads, logs, and summary
#   fwd_primers.fasta FASTA with forward primer(s); supports a pool for
#                     multiplex protocols
#   rev_primers.fasta FASTA with reverse primer(s)
#
# Outputs:
#   <output_dir>/*.fastq.gz                 primer-trimmed reads
#   <output_dir>/logs/<sample>.log          per-sample cutadapt log
#   <output_dir>/cutadapt_summary.tsv       sample, in_pairs, out_pairs, kept_pct
#
# Requirements: cutadapt (>=3.0).

set -euo pipefail

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <fwd_primers.fasta> <rev_primers.fasta>" >&2
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
FWD_FASTA="$3"
REV_FASTA="$4"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input_dir does not exist: $INPUT_DIR" >&2
    exit 1
fi

if [[ ! -f "$FWD_FASTA" ]]; then
    echo "ERROR: forward primers FASTA not found: $FWD_FASTA" >&2
    exit 1
fi

if [[ ! -f "$REV_FASTA" ]]; then
    echo "ERROR: reverse primers FASTA not found: $REV_FASTA" >&2
    exit 1
fi

# Tunables
THREADS="${THREADS:-4}"          # threads per cutadapt invocation
MIN_LEN="${MIN_LEN:-50}"         # minimum read length after trimming
ERROR_RATE="${ERROR_RATE:-0.1}"  # cutadapt error tolerance (fraction)

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------------------------------------------------------------------------
# Build reverse-complement FASTAs (for 3' read-through trimming on short
# amplicons where R1 reads into the reverse primer and vice versa).
# Supports IUPAC ambiguity codes.
# ---------------------------------------------------------------------------

revcomp_fasta() {
    # $1 = input FASTA, $2 = output FASTA
    awk 'BEGIN{RS=">";ORS=""} NR>1 {
            n=index($0,"\n");
            header=substr($0,1,n-1);
            seq=substr($0,n+1);
            gsub("\n","",seq);
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

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

command -v cutadapt >/dev/null 2>&1 \
    || { echo "ERROR: cutadapt not found in PATH" >&2; exit 1; }

echo "cutadapt version: $(cutadapt --version)"
echo "Forward primers: $(grep -c '^>' "$FWD_FASTA")"
echo "Reverse primers: $(grep -c '^>' "$REV_FASTA")"

# ---------------------------------------------------------------------------
# Discover paired samples
# ---------------------------------------------------------------------------

shopt -s nullglob
R1_FILES=( "$INPUT_DIR"/*_R1*.fastq.gz "$INPUT_DIR"/*_R1*.fastq )
shopt -u nullglob

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no *_R1*.fastq(.gz) files found in $INPUT_DIR" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Main loop: one sample at a time
# ---------------------------------------------------------------------------

SUMMARY="$OUTPUT_DIR/cutadapt_summary.tsv"
printf 'sample\tinput_pairs\toutput_pairs\tkept_pct\n' > "$SUMMARY"

for r1 in "${R1_FILES[@]}"; do
    base="$(basename "$r1")"
    r2="${r1/_R1/_R2}"
    sample="${base%%_R1*}"

    if [[ ! -f "$r2" ]]; then
        echo "WARN: missing R2 for $r1, skipping" >&2
        continue
    fi

    out_r1="$OUTPUT_DIR/${sample}_R1.fastq.gz"
    out_r2="$OUTPUT_DIR/${sample}_R2.fastq.gz"
    log="$OUTPUT_DIR/logs/${sample}.log"

    echo "Processing: $sample"

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

    # Pull total/written read counts from cutadapt log
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

echo
echo "Done. Summary: $SUMMARY"
column -t -s $'\t' "$SUMMARY"
