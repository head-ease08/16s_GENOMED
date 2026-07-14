#!/usr/bin/env python3
"""
Total_Bases / Bases_0x / Mean_Coverage / Median_Coverage for one sample,
restricted to reference contigs that received >=1 mapped read.

The SILVA reference has tens of thousands of entries; a sample's reads only
ever hit a handful of them. Running mosdepth over the whole reference would
drown every stat in zero-coverage bases from contigs that were never in this
sample, making Mean_Coverage ~0 and Bases_0x ~100% regardless of how well the
detected organisms are actually covered.

Earlier version tried to subset the BAM (samtools view -b bam <contigs>) and
then trim the header to match. That's wrong: BAM alignment records store
their reference as a numeric index into the @SQ list, not by name. Dropping
@SQ lines renumbers everything after them, so every read's refID silently
points at the wrong contig (or an out-of-range one) after reheadering --
that's what made `samtools index` fail. No BAM surgery needed: mosdepth's
`-b <bed>` restricts BOTH per-base and summary output to the listed regions
directly on the original, untouched, already-indexed BAM.

Usage: coverage_metrics.py <sample.bam> <output_tsv>
"""
import sys
import os
import subprocess
import gzip
import numpy as np


def main():
    if len(sys.argv) != 3:
        print("Usage: coverage_metrics.py <sample.bam> <output_tsv>")
        sys.exit(1)

    bam, out_tsv = sys.argv[1], sys.argv[2]
    sample = os.path.basename(bam).split(".")[0]
    workdir = os.path.dirname(out_tsv) or "."
    os.makedirs(workdir, exist_ok=True)
    bed_path = os.path.join(workdir, f"{sample}.covered.bed")
    prefix = os.path.join(workdir, f"{sample}.mosdepth")

    idxstats = subprocess.run(
        ["samtools", "idxstats", bam], check=True, capture_output=True, text=True
    ).stdout
    covered = []
    for line in idxstats.splitlines():
        if not line:
            continue
        ref_name, length, mapped, _unmapped = line.split("\t")
        if int(mapped) > 0:
            covered.append((ref_name, int(length)))

    if not covered:
        with open(out_tsv, "w") as f:
            f.write("Sample\tTotal_Bases\tBases_0x\tMean_Coverage\tMedian_Coverage\n")
            f.write(f"{sample}\tNA\tNA\tNA\tNA\n")
        return

    with open(bed_path, "w") as f:
        for ref_name, length in covered:
            f.write(f"{ref_name}\t0\t{length}\n")

    subprocess.run(["mosdepth", "-x", "-t", "4", "-b", bed_path, prefix, bam], check=True)

    per_base = prefix + ".per-base.bed.gz"
    depths, weights = [], []
    with gzip.open(per_base, "rt") as f:
        for line in f:
            _, start, end, cov = line.rstrip("\n").split("\t")
            depths.append(int(cov))
            weights.append(int(end) - int(start))

    depths = np.array(depths)
    weights = np.array(weights)
    total_bases = int(weights.sum())
    bases_0x = int(weights[depths == 0].sum())
    mean_cov = float(np.average(depths, weights=weights))

    order = np.argsort(depths)
    d_sorted, w_sorted = depths[order], weights[order]
    cum = np.cumsum(w_sorted)
    median_cov = float(d_sorted[np.searchsorted(cum, total_bases / 2.0)])

    with open(out_tsv, "w") as f:
        f.write("Sample\tTotal_Bases\tBases_0x\tMean_Coverage\tMedian_Coverage\n")
        f.write(f"{sample}\t{total_bases}\t{bases_0x}\t{mean_cov:.2f}\t{median_cov:.2f}\n")

    for ext in (
        ".mosdepth.global.dist.txt",
        ".mosdepth.region.dist.txt",
        ".mosdepth.summary.txt",
        ".per-base.bed.gz",
        ".per-base.bed.gz.csi",
        ".regions.bed.gz",
        ".regions.bed.gz.csi",
    ):
        fp = prefix + ext
        if os.path.exists(fp):
            os.remove(fp)
    os.remove(bed_path)


if __name__ == "__main__":
    main()
