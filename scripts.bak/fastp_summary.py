#!/usr/bin/env python3
"""
Parse fastp *_fastp.json files into one TSV — same fields/paths
QC_pipe_frag/scripts/build_final_table.py reads from fastp JSON:
  Raw total sequences            <- summary.before_filtering.total_reads
  Reads passed filters           <- filtering_result.passed_filter_reads
  Average read length (after..)  <- summary.after_filtering.read1_mean_length
  Reads duplicated_%              <- duplication.rate * 100

Usage: fastp_summary.py <fastp_json_dir> <output_tsv>
Expects files named {sample}_fastp.json in the given directory.
"""
import sys
import os
import json
import glob


def fastp_value(j, *path, default="NA"):
    cur = j
    for k in path:
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return default
    return cur


def main():
    if len(sys.argv) != 3:
        print("Usage: fastp_summary.py <fastp_json_dir> <output_tsv>")
        sys.exit(1)

    json_dir, out_tsv = sys.argv[1], sys.argv[2]

    with open(out_tsv, "w") as out:
        out.write(
            "Sample\tRaw total sequences\tReads passed filters\t"
            "Average read length (after filtering)\tReads duplicated_%\n"
        )
        for f in sorted(glob.glob(os.path.join(json_dir, "*_fastp.json"))):
            sample = os.path.basename(f).replace("_fastp.json", "")
            with open(f) as fh:
                j = json.load(fh)

            raw_total = fastp_value(j, "summary", "before_filtering", "total_reads")
            passed = fastp_value(j, "filtering_result", "passed_filter_reads")
            avg_len = fastp_value(j, "summary", "after_filtering", "read1_mean_length")
            dup = fastp_value(j, "duplication", "rate", default=None)
            dup_pct = "NA"
            if dup is not None:
                try:
                    dup_pct = f"{float(dup) * 100:.2f}"
                except (TypeError, ValueError):
                    pass

            out.write(f"{sample}\t{raw_total}\t{passed}\t{avg_len}\t{dup_pct}\n")

    print(f"=== Done: {out_tsv} ===")


if __name__ == "__main__":
    main()
