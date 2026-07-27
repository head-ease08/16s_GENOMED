#!/usr/bin/env python3
"""
Per-sample Taxon x Regions x Reads x Pct table, built from the ASV-level
abundance tables that `region_dada2_all` writes per V-region
(results/region/<region>/abundance_table_{silva,rdp}.csv --
scripts/make_abundance_table.R). ASVs are collapsed to Genus+Species per
region, then summed across regions for each sample -- Reads is the total
across all regions, Regions lists which ones actually contributed, Pct is
of that sample's total reads (all taxa, all regions).

Usage:
    region_taxon_summary.py --region_dir results/region --db silva \
        --out taxon_by_sample.tsv
"""
import argparse
import csv
import glob
import os
from collections import defaultdict


def taxon_name(row):
    genus = (row.get("Genus") or "").strip()
    species = (row.get("Species") or "").strip()
    if not genus or genus.upper() == "NA":
        return "na"
    return f"{genus} {species}" if species and species.upper() != "NA" else genus


def load_region_table(path):
    """-> {sample: {taxon: reads}}, reading only raw count columns (not _pct)."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        non_sample = {"asv_id", "sequence", "Kingdom", "Phylum", "Class",
                       "Order", "Family", "Genus", "Species"}
        sample_cols = [c for c in fieldnames if c not in non_sample and not c.endswith("_pct")]

        per_sample = defaultdict(lambda: defaultdict(int))
        for row in reader:
            taxon = taxon_name(row)
            for s in sample_cols:
                reads = int(float(row[s] or 0))
                if reads:
                    per_sample[s][taxon] += reads
        return per_sample


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region_dir", default="results/region")
    ap.add_argument("--db", default="silva", choices=["silva", "rdp"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pattern = os.path.join(args.region_dir, "*", f"abundance_table_{args.db}.csv")
    tables = sorted(glob.glob(pattern))
    if not tables:
        raise SystemExit(f"No abundance_table_{args.db}.csv found under {args.region_dir}/*/")

    # sample -> taxon -> {"reads": int, "regions": set}
    combined = defaultdict(lambda: defaultdict(lambda: {"reads": 0, "regions": set()}))
    for path in tables:
        region = os.path.basename(os.path.dirname(path))
        per_sample = load_region_table(path)
        for sample, taxa in per_sample.items():
            for taxon, reads in taxa.items():
                entry = combined[sample][taxon]
                entry["reads"] += reads
                entry["regions"].add(region)

    with open(args.out, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["Sample", "Taxon", "Reads", "Regions", "Pct"])
        for sample in sorted(combined):
            taxa = combined[sample]
            total = sum(t["reads"] for t in taxa.values())
            for taxon, entry in sorted(taxa.items(), key=lambda kv: -kv[1]["reads"]):
                pct = entry["reads"] / total * 100 if total else 0.0
                regions = ",".join(sorted(entry["regions"]))
                writer.writerow([sample, taxon, entry["reads"], regions, f"{pct:.2f}"])

    print(f"=== Done. {len(tables)} region tables -> {args.out} ===")


if __name__ == "__main__":
    main()
