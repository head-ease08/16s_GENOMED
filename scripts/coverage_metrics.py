#!/usr/bin/env python3
"""
Total_Bases / Bases_0x / Mean_Coverage / Median_Coverage for one sample,
restricted to reference contigs that received >=1 mapped read.

The SILVA reference has tens of thousands of entries; a sample's reads only
ever hit a handful of them. Running mosdepth over the whole reference would
drown every stat in zero-coverage bases from contigs that were never in this
sample, making Mean_Coverage ~0 and Bases_0x ~100% regardless of how well the
detected organisms are actually covered. So: subset the BAM to covered
contigs first, then run mosdepth on that subset.

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
    subset_bam = os.path.join(workdir, f"{sample}.covered.bam")
    prefix = os.path.join(workdir, f"{sample}.mosdepth")

    idxstats = subprocess.run(
        ["samtools", "idxstats", bam], check=True, capture_output=True, text=True
    ).stdout
    covered = [
        line.split("\t")[0]
        for line in idxstats.splitlines()
        if line and int(line.split("\t")[2]) > 0
    ]

    if not covered:
        with open(out_tsv, "w") as f:
            f.write("Sample\tTotal_Bases\tBases_0x\tMean_Coverage\tMedian_Coverage\n")
            f.write(f"{sample}\tNA\tNA\tNA\tNA\n")
        return

    with open(subset_bam, "wb") as f:
        subprocess.run(["samtools", "view", "-b", bam] + covered, check=True, stdout=f)
    subprocess.run(["samtools", "index", subset_bam], check=True)
    subprocess.run(["mosdepth", "-x", "-t", "4", prefix, subset_bam], check=True)

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
        ".mosdepth.summary.txt",
        ".per-base.bed.gz",
        ".per-base.bed.gz.csi",
    ):
        fp = prefix + ext
        if os.path.exists(fp):
            os.remove(fp)
    os.remove(subset_bam)
    if os.path.exists(subset_bam + ".bai"):
        os.remove(subset_bam + ".bai")


if __name__ == "__main__":
    main()
