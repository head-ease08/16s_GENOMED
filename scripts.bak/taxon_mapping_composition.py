#!/usr/bin/env python3
"""
Per-sample mapped-read composition: samtools idxstats gives read count per
reference contig, contig id is joined back to its SILVA lineage string via
the map TSV written by prep_silva_reference.sh, counts are summed per lineage
and turned into %. Cross-check against the DADA2-classifier-based
abundance_table_silva.csv (independent method: alignment vs naive-Bayes).

Usage: taxon_mapping_composition.py <sample.bam> <map.tsv> <output_tsv>
"""
import sys
import subprocess
from collections import defaultdict


def main():
    if len(sys.argv) != 4:
        print("Usage: taxon_mapping_composition.py <sample.bam> <map.tsv> <output_tsv>")
        sys.exit(1)

    bam, map_tsv, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
    sample = bam.split("/")[-1].split(".")[0]

    id_to_lineage = {}
    with open(map_tsv) as f:
        for line in f:
            ref_id, lineage = line.rstrip("\n").split("\t", 1)
            id_to_lineage[ref_id] = lineage

    idxstats = subprocess.run(
        ["samtools", "idxstats", bam], check=True, capture_output=True, text=True
    ).stdout

    # Group by Genus+Species (last two ranks), not the full lineage string.
    # SILVA has inconsistent higher-rank naming across entries for the same
    # organism (e.g. "Proteobacteria" vs "Pseudomonadota" -- a historical
    # phylum rename that both spellings still appear under in the DB), so
    # grouping on the full string leaves the same species split across
    # multiple rows with different counts instead of one summed total.
    counts = defaultdict(int)
    total_mapped = 0
    for line in idxstats.splitlines():
        if not line:
            continue
        ref_id, _length, mapped, _unmapped = line.split("\t")
        mapped = int(mapped)
        if mapped == 0:
            continue
        lineage = id_to_lineage.get(ref_id, "unknown")
        ranks = [r for r in lineage.rstrip(";").split(";") if r]
        taxon = " ".join(ranks[-2:]) if len(ranks) >= 2 else (ranks[0] if ranks else "unknown")
        counts[taxon] += mapped
        total_mapped += mapped

    with open(out_tsv, "w") as out:
        out.write("Sample\tTaxon\tReads\tPct\n")
        if total_mapped == 0:
            out.write(f"{sample}\tNA\t0\tNA\n")
            return
        for taxon, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            pct = n / total_mapped * 100
            out.write(f"{sample}\t{taxon}\t{n}\t{pct:.2f}\n")


if __name__ == "__main__":
    main()
