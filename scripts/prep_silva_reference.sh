#!/usr/bin/env bash
# DADA2 SILVA training-set FASTA headers are bare taxonomy strings
# (">Kingdom;Phylum;...;Species"), no accession — hundreds of near-identical
# entries share the same lineage. Aligning against all of them scatters a
# sample's reads thinly across many near-duplicate contigs: coverage per
# contig rounds to ~0x and mates often land on two different copies of "the
# same" organism, so proper-pair/insert-size comes out empty. Keep only the
# FIRST sequence seen per unique lineage string (one representative contig
# per organism) so reads from the same organism actually pile up on one
# contig. Also rewrites headers to ref<N> (bwa-mem2 needs unique names) and
# keeps id->lineage in a side TSV for lookup after alignment
# (taxon_mapping_composition.py joins on this file).
#
# Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>
set -euo pipefail

IN=${1:?Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>}
OUT=${2:?Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>}
MAP=${3:?Usage: prep_silva_reference.sh <silva_training_fasta.gz> <out_fasta> <out_map_tsv>}

: > "$MAP"

zcat "$IN" | awk -v map="$MAP" '
/^>/ {
    lineage = substr($0, 2)
    if (lineage in seen) {
        keep = 0
        next
    }
    seen[lineage] = 1
    n++
    id = "ref" n
    print id "\t" lineage >> map
    print ">" id
    keep = 1
    next
}
keep { print }
' > "$OUT"

echo "=== $(grep -c '^>' "$OUT") unique-lineage contigs written to $OUT ===" >&2
