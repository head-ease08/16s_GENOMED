#!/usr/bin/env python3
"""
Assemble the final per-sample QC table from track.csv (DADA2 filterAndTrim
counts) + the alignment-based metrics TSVs (flagstat, GC, insert size,
adapter dimers, coverage).

Usage: build_qc_table.py --track results/track.csv --read_length ... \
    --flagstat ... --gc ... --insert_size ... --dimers ... --coverage ... \
    --out results/qc_align/final_qc_metrics.tsv
"""
import argparse
import csv
from collections import OrderedDict

FINAL_COLS = [
    "Sample",
    "Raw total sequences",
    "Reads passed filters",
    "Reads mapped",
    "Alignment reads_%",
    "Average read length (after filtering)",
    "Insert size peak",
    "Reads duplicated_%",
    "Adapter dimers_%",
    "GC_R1_%",
    "GC_R2_%",
    "GC_Mean_%",
    "Total_Bases",
    "Bases_0x",
    "Mean_Coverage",
    "Median_Coverage",
]


def load_tsv(path, key="Sample"):
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[row[key]] = row
    return out


def load_track(path):
    """track.csv is comma-separated, written by scripts/create_summary.R."""
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            out[row["sample"]] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", required=True)
    ap.add_argument("--read_length", required=True)
    ap.add_argument("--flagstat", required=True)
    ap.add_argument("--gc", required=True)
    ap.add_argument("--insert_size", required=True)
    ap.add_argument("--dimers", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    track = load_track(args.track)
    read_length = load_tsv(args.read_length)
    flagstat = load_tsv(args.flagstat)
    gc = load_tsv(args.gc)
    insert_size = load_tsv(args.insert_size)
    dimers = load_tsv(args.dimers)
    coverage = load_tsv(args.coverage)

    all_samples = sorted(
        set(track) | set(read_length) | set(flagstat) | set(gc)
        | set(insert_size) | set(dimers) | set(coverage)
    )

    with open(args.out, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(FINAL_COLS)

        for s in all_samples:
            row = OrderedDict.fromkeys(FINAL_COLS, "NA")
            row["Sample"] = s

            t = track.get(s, {})
            row["Raw total sequences"] = t.get("input", "NA")
            row["Reads passed filters"] = t.get("filtered", "NA")

            rl = read_length.get(s, {})
            row["Average read length (after filtering)"] = rl.get("Average_read_length", "NA")

            fs = flagstat.get(s, {})
            row["Reads mapped"] = fs.get("Reads_mapped", "NA")
            row["Alignment reads_%"] = fs.get("Alignment_reads_%", "NA")
            row["Reads duplicated_%"] = fs.get("Reads_duplicated_%", "NA")

            ins = insert_size.get(s, {})
            row["Insert size peak"] = ins.get("Insert_size_peak", "NA")

            dm = dimers.get(s, {})
            row["Adapter dimers_%"] = dm.get("Adapter_dimers_%", "NA")

            g = gc.get(s, {})
            row["GC_R1_%"] = g.get("GC_R1_%", "NA")
            row["GC_R2_%"] = g.get("GC_R2_%", "NA")
            row["GC_Mean_%"] = g.get("GC_Mean_%", "NA")

            cov = coverage.get(s, {})
            row["Total_Bases"] = cov.get("Total_Bases", "NA")
            row["Bases_0x"] = cov.get("Bases_0x", "NA")
            row["Mean_Coverage"] = cov.get("Mean_Coverage", "NA")
            row["Median_Coverage"] = cov.get("Median_Coverage", "NA")

            writer.writerow([row[c] for c in FINAL_COLS])

    print(f"=== Done. Output: {args.out} ({len(all_samples)} samples) ===")


if __name__ == "__main__":
    main()
