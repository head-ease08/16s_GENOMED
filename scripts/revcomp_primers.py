IUPAC_RC = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

def revcomp_fasta(src, dst):
    with open(src) as f, open(dst, "w") as out:
        for rec in f.read().lstrip(">").split(">"):
            if not rec.strip():
                continue
            header, *lines = rec.splitlines()
            seq = "".join(lines)
            out.write(f">{header}_RC\n{seq[::-1].translate(IUPAC_RC)}\n")

revcomp_fasta(snakemake.input.fwd, snakemake.output.fwd_rc)
revcomp_fasta(snakemake.input.rev, snakemake.output.rev_rc)
