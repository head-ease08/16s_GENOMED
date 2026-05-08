from pathlib import Path

RAW_DIR = "data/raw"
QC_DIR = "data/qc"
TRIMMED_DIR = "data/trimmed"
PRIMERS_DIR = "primers"

SAMPLES, = glob_wildcards(f"{RAW_DIR}/{{sample}}_R1.fastq.gz")

import os
for s in SAMPLES:
    r2 = f"{RAW_DIR}/{s}_R2.fastq.gz"
    assert os.path.exists(r2), f"Нет пары для {s}: {r2}"

IUPAC_RC = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

rule revcomp_primers:
    input:
        fwd = PRIMERS_DIR + "/fwd_primers.fasta",
        rev = PRIMERS_DIR + "/rev_primers.fasta",
    output:
        fwd_rc = temp(PRIMERS_DIR + "/fwd_rc.fasta"),
        rev_rc = temp(PRIMERS_DIR + "/rev_rc.fasta"),
    run:
        def revcomp_fasta(src, dst):
            with open(src) as f, open(dst, "w") as out:
                for rec in f.read().lstrip(">").split(">"):
                    if not rec.strip():
                        continue
                    header, *lines = rec.splitlines()
                    seq = "".join(lines)
                    out.write(f">{header}_RC\n{seq[::-1].translate(IUPAC_RC)}\n")
        revcomp_fasta(input.fwd, output.fwd_rc)
        revcomp_fasta(input.rev, output.rev_rc)


rule trim_primers:
    input:
        r1      = RAW_DIR + "/{sample}_R1.fastq.gz",
        r2      = RAW_DIR + "/{sample}_R2.fastq.gz",
        fwd     = PRIMERS_DIR + "/fwd_primers.fasta",
        rev     = PRIMERS_DIR + "/rev_primers.fasta",
        fwd_rc  = PRIMERS_DIR + "/fwd_rc.fasta",
        rev_rc  = PRIMERS_DIR + "/rev_rc.fasta",
    output:
        r1 = TRIMMED_DIR + "/{sample}_R1.fastq.gz",
        r2 = TRIMMED_DIR + "/{sample}_R2.fastq.gz",
    log:
        "logs/cutadapt/{sample}.log",
    params:
        threads    = 4,
        min_len    = 50,
        error_rate = 0.1,
    shell:
        """
        cutadapt \
            -g file:{input.fwd} \
            -a file:{input.rev_rc} \
            -G file:{input.rev} \
            -A file:{input.fwd_rc} \
            -n 2 \
            -e {params.error_rate} \
            --discard-untrimmed \
            --minimum-length {params.min_len} \
            -j {params.threads} \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule quality_per_sample:
    input:
        r1 = expand(f"{RAW_DIR}/{{sample}}_R1.fastq.gz", sample=SAMPLES),
        r2 = expand(f"{RAW_DIR}/{{sample}}_R2.fastq.gz", sample=SAMPLES),
    output:
        r1_pdf = f"{QC_DIR}/per_sample_forward.pdf",
        r2_pdf = f"{QC_DIR}/per_sample_reverse.pdf",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/quality_per_sample.R"

rule quality_aggregated:
    input:
        r1 = expand(f"{RAW_DIR}/{{sample}}_R1.fastq.gz", sample=SAMPLES),
        r2 = expand(f"{RAW_DIR}/{{sample}}_R2.fastq.gz", sample=SAMPLES),
    output:
        r1_pdf = f"{QC_DIR}/aggregated_forward.pdf",
        r2_pdf = f"{QC_DIR}/aggregated_reverse.pdf",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/quality_aggregated.R"

rule filter_reads:
    input:
        r1 = TRIMMED_DIR + "/{sample}_R1.fastq.gz",
        r2 = TRIMMED_DIR + "/{sample}_R2.fastq.gz",
    output:
        r1 = QC_DIR + "/{sample}_R1.fastq.gz",
        r2 = QC_DIR + "/{sample}_R2.fastq.gz",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/filter_reads.R"
