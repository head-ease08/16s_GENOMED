#!/usr/bin/env bash
# Standalone V-region demultiplexing, no Snakemake. Splits every R1/R2 pair
# found flat in a directory into per-region trimmed fastq (primers stripped),
# same cutadapt logic and primer table as `rule demux_region` in the Snakefile
# (see REGION_PRIMERS there). Reads not matching a region's primer pair are
# discarded for that region (--discard-untrimmed) -- normal, a given read
# only belongs to one region.
#
# Usage: demux_regions.sh <reads_dir> <out_dir> [jobs]
#   out_dir layout: <out_dir>/<region>/<sample>_R1.fq.gz + _R2.fq.gz
set -euo pipefail

READS_DIR=${1:?Usage: demux_regions.sh <reads_dir> <out_dir> [jobs]}
OUT_DIR=${2:?missing out_dir}
JOBS=${3:-4}

# region -> "fwd_primer rev_primer" (must match Snakefile REGION_PRIMERS)
REGIONS="V1_V2 V9 V3_V4 V4_V5 V6_V8"
fwd_V1_V2="AGAGTTTGATCMTGGCTCAG";  rev_V1_V2="GGACCGTGTCTCAGTTCCAG"
fwd_V9="TGCCACGGTGAATACGTTCC";     rev_V9="CCTTGTTACGACTTCACCCCA"
fwd_V3_V4="CCTACGGGNGGCWGCAG";     rev_V3_V4="GGACTACHVGGGTATCTAATCC"
fwd_V4_V5="GGAGGGTGCAAGCGTTAATC"; rev_V4_V5="TTAACCTTGCGGCCGTACTC"
fwd_V6_V8="CGGTGGAGCATGTGGTTTAA"; rev_V6_V8="AGTTGCAGACTCCAATCCGG"

shopt -s nullglob nocaseglob
r1_files=("$READS_DIR"/*R1*.f*q*.gz)
shopt -u nullglob nocaseglob
if [ ${#r1_files[@]} -eq 0 ]; then
    echo "No *R1*.f*q*.gz files found in $READS_DIR" >&2
    exit 1
fi

demux_one() {
    r1="$1"; region="$2"; fwd="$3"; rev="$4"; out_dir="$5"
    base=$(basename "$r1")
    r2_base=$(echo "$base" | sed -E 's/R1/R2/; s/r1/r2/')
    r2="$(dirname "$r1")/$r2_base"
    sample=$(echo "$base" | sed -E 's/[._-]?[Rr]1.*$//')
    [ -f "$r2" ] || { echo "WARNING: no R2 match for $base -- skipping" >&2; return 0; }

    region_dir="$out_dir/$region"
    mkdir -p "$region_dir"
    cutadapt -g "$fwd" -G "$rev" --discard-untrimmed -e 0.1 \
        --minimum-length 50 -j 1 \
        -o "$region_dir/${sample}_R1.fq.gz" -p "$region_dir/${sample}_R2.fq.gz" \
        "$r1" "$r2" > "$region_dir/${sample}.log" 2>&1
}
export -f demux_one

jobs=()
for region in $REGIONS; do
    fwd_var="fwd_${region}"; rev_var="rev_${region}"
    fwd="${!fwd_var}"; rev="${!rev_var}"
    for r1 in "${r1_files[@]}"; do
        jobs+=("$r1|$region|$fwd|$rev")
    done
done

echo "Demuxing ${#r1_files[@]} samples x 5 regions = ${#jobs[@]} cutadapt runs, $JOBS parallel" >&2

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${jobs[@]}" | parallel -j "$JOBS" --colsep '\|' demux_one {1} {2} {3} {4} "$OUT_DIR"
else
    for j in "${jobs[@]}"; do
        IFS='|' read -r r1 region fwd rev <<< "$j"
        demux_one "$r1" "$region" "$fwd" "$rev" "$OUT_DIR" &
        while [ "$(jobs -r -p | wc -l)" -ge "$JOBS" ]; do wait -n; done
    done
    wait
fi

echo "Done. Per-region reads in $OUT_DIR/<region>/" >&2
for region in $REGIONS; do
    n=$(find "$OUT_DIR/$region" -iname '*_R1.fq.gz' 2>/dev/null | wc -l)
    echo "  $region: $n samples" >&2
done
