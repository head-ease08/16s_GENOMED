#!/usr/bin/env python3
"""
Long-format sample × species × pct table from a per-region abundance CSV.

Aggregates ASV-level rows to species (Genus + Species), sums counts per
(sample, species), recomputes pct against per-sample total (across ALL ASVs
in the input, not just resolved-to-species — so pcts are comparable to the
original _pct columns).

Usage:
    python3 species_by_sample.py <abundance_table.csv> <out.tsv> [--include-unresolved]

--include-unresolved keeps rows where Species is NA (labeled by lowest
non-NA rank, e.g. "Neisseria sp." or "Neisseriaceae unclassified").
Without the flag: only species-level hits.
"""

import argparse
import re
import sys
import pandas as pd

SAMPLE_RE = re.compile(r"RnDD_260713_(\d+)_RnDD-16S_n0$")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--include-unresolved", action="store_true")
    ap.add_argument("--min-pct", type=float, default=0.0,
                    help="drop rows below this pct (default 0 = keep all non-zero)")
    args = ap.parse_args()

    df = pd.read_csv(args.infile)

    # identify sample count columns (integer id, no _pct suffix)
    sample_cols = {c: int(SAMPLE_RE.match(c).group(1))
                   for c in df.columns
                   if not c.endswith("_pct") and SAMPLE_RE.match(c)}
    if not sample_cols:
        sys.exit("no sample columns matched")

    # build a species label
    def label(row):
        g, s = row.get("Genus"), row.get("Species")
        if pd.notna(g) and pd.notna(s):
            return f"{g} {s}"
        if not args.include_unresolved:
            return None
        # fallback: lowest non-NA rank
        for rank in ["Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]:
            v = row.get(rank)
            if pd.notna(v):
                return f"{v} (unresolved)"
        return "unclassified"

    df["_species"] = df.apply(label, axis=1)
    df = df[df["_species"].notna()].copy()

    # per-sample totals BEFORE aggregation (from full original table, to keep
    # pcts comparable to make_abundance_table.R output)
    original = pd.read_csv(args.infile)
    totals = {sid: original[col].sum() for col, sid in sample_cols.items()}

    # aggregate counts by species
    agg = df.groupby("_species")[list(sample_cols.keys())].sum()
    agg.columns = [sample_cols[c] for c in agg.columns]  # rename to sample ids

    # melt to long
    long = agg.reset_index().melt(id_vars="_species", var_name="sample",
                                   value_name="reads")
    long = long.rename(columns={"_species": "species"})
    long["sample"] = long["sample"].astype(int)
    long["total_reads_in_sample"] = long["sample"].map(totals)
    long["pct"] = long["reads"] / long["total_reads_in_sample"] * 100

    long = long[long["reads"] > 0]
    long = long[long["pct"] >= args.min_pct]
    long = long.sort_values(["sample", "pct"], ascending=[True, False])

    long[["sample", "species", "reads", "pct"]].to_csv(
        args.outfile, sep="\t", index=False, float_format="%.3f"
    )
    print(f"wrote {len(long)} rows to {args.outfile}")
    print(f"samples: {sorted(long['sample'].unique())}")
    print(f"species: {long['species'].nunique()}")


if __name__ == "__main__":
    main()
