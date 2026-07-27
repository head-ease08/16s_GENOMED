#!/usr/bin/env python3
"""
Total_Bases / Bases_0x / Mean_Coverage / Median_Coverage for one sample,
restricted to reference contigs that received >=1 mapped read.

The SILVA reference has thousands of entries; a sample's reads only ever hit
a handful of them. Running mosdepth over the whole reference would drown
every stat in zero-coverage bases from contigs that were never touched.

Two earlier approaches both turned out wrong:
  1. `samtools view -b bam <contigs>` then hand-trim the header's @SQ lines:
     BAM records store their reference as a numeric index into the @SQ list,
     not by name. Dropping @SQ lines renumbers everything after them, so
     every read's refID silently points at the wrong (or out-of-range)
     contig -- `samtools index` correctly rejects the result.
  2. `mosdepth -b <bed-of-covered-contigs>`: this only restricts the
     *.regions.bed.gz summary. *.per-base.bed.gz is STILL computed over
     every contig in the BAM header regardless -- mosdepth just collapses
     each untouched contig into one big zero-depth row, so the row count
     looks small but Total_Bases still sums the entire multi-thousand-
     contig reference.

The only correct way to shrink the header AND keep every read's reference
binding consistent is to rebuild the BAM from SAM text: samtools resolves
RNAME by name (not index) when parsing text, so trimming the @SQ lines in
the header and feeding only the matching body lines back through
`samtools view -b` makes samtools recompute correct refIDs itself as part
of that conversion. No manual reheadering.

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
    covered_names = []
    for line in idxstats.splitlines():
        if not line:
            continue
        ref_name, _length, mapped, _unmapped = line.split("\t")
        if int(mapped) > 0:
            covered_names.append(ref_name)

    if not covered_names:
        with open(out_tsv, "w") as f:
            f.write("Sample\tTotal_Bases\tBases_0x\tMean_Coverage\tMedian_Coverage\n")
            f.write(f"{sample}\tNA\tNA\tNA\tNA\n")
        return

    covered_set = set(covered_names)
    header = subprocess.run(
        ["samtools", "view", "-H", bam], check=True, capture_output=True, text=True
    ).stdout
    trimmed_lines = []
    for line in header.splitlines():
        if line.startswith("@SQ"):
            sn = next((f[3:] for f in line.split("\t") if f.startswith("SN:")), None)
            if sn not in covered_set:
                continue
        trimmed_lines.append(line)
    trimmed_header = "\n".join(trimmed_lines) + "\n"

    body = subprocess.run(
        ["samtools", "view", bam] + covered_names, check=True, capture_output=True, text=True
    ).stdout

    subprocess.run(
        ["samtools", "view", "-b", "-o", subset_bam, "-"],
        input=trimmed_header + body,
        text=True,
        check=True,
    )
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
