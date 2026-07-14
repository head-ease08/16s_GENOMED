#!/usr/bin/env bash
# DADA2 SILVA training-set FASTA headers are bare taxonomy strings
# (">Kingdom;Phylum;...;Species"), no accession — many entries share an
# identical lineage string. bwa-mem2 needs unique contig names, so rewrite
# headers to ref<N> and keep the original lineage in a side TSV for lookup
# after alignment (taxon_mapping_composition.py joins on this file).
#
# Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>
set -euo pipefail

IN=${1:?Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>}
OUT=${2:?Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>}
MAP=${3:?Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>}

: > "$MAP"

zcat "$IN" | awk -v map="$MAP" '
/^>/ {
    n++
    lineage = substr($0, 2)
    id = "ref" n
    print id "\t" lineage >> map
    print ">" id
    next
}
{ print }
' > "$OUT"
