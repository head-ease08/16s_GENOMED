#!/usr/bin/env python3
# Expand fastp's insert_size.histogram (index=bp, value=count) into a flat
# one-value-per-line file, same format plot_insert_size_hist.R expects.
# Usage: fastp_json_to_insert_sizes.py <report.json> <out.txt>
import json
import sys

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: fastp_json_to_insert_sizes.py <report.json> <out.txt>")
    json_path, out_path = sys.argv[1], sys.argv[2]

    with open(json_path) as f:
        report = json.load(f)

    hist = report.get("insert_size", {}).get("histogram", [])
    with open(out_path, "w") as out:
        for size, count in enumerate(hist):
            if count:
                out.write((str(size) + "\n") * count)

if __name__ == "__main__":
    main()
