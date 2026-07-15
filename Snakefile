from pathlib import Path
import os
import glob

RAW_DIR     = config.get("run_dir", "data/raw")
QC_DIR      = "data/qc"
TRIMMED_DIR = "data/trimmed"
PRIMERS_DIR = "primers"
REF_DIR     = "references"

# Alignment-based QC (rules at the bottom of this file, target: qc_align_all).
# Reuses the SILVA training-set FASTA already required for assign_taxonomy_silva.
SILVA_TRAINSET = REF_DIR + "/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
QC_ALIGN_DIR   = "results/qc_align"
BAM_DIR        = QC_ALIGN_DIR + "/bam"
FLAGSTAT_DIR   = QC_ALIGN_DIR + "/flagstat"
DIMERS_DIR     = QC_ALIGN_DIR + "/adapter_dimers"
COV_DIR        = QC_ALIGN_DIR + "/coverage"
COMP_DIR       = QC_ALIGN_DIR + "/composition"

# Per-region DADA2 (rules at the bottom, target: region_dada2_all).
# Multiplexed multi-V-region primer pool (see primers/*.fasta) can't go
# through one pooled DADA2 run: learnErrors assumes one error model for one
# amplicon length, and 2x151bp reads physically can't overlap-merge every
# region's amplicon. V3_V4 (~407bp insert), V4_V5 (~330bp), V6_V8 (~341bp)
# are too long for 2x151bp to ever meet in the middle regardless of truncLen
# -- those run DADA2 forward-read-only (no merge step). V9 (~132bp insert)
# clearly merges. V1_V2 (~280bp insert) is borderline once primer removal is
# accounted for (see truncLen note below) -- kept in the merge bucket for
# now, but check its actual merge rate in results/region/V1_V2/merged/ once
# this runs; if it's near-zero, move it to the SE bucket like the others.
#
# truncLen note: demux_region strips the primer with cutadapt BEFORE
# filterAndTrim ever sees the read, so the truncLen budget is against the
# post-primer-trim length (measured ~130bp for a 151bp raw read, not 151),
# and filterAndTrim discards any read shorter than truncLen outright rather
# than soft-trimming it. 125/120 leaves a small quality-trim margin under
# that 130bp ceiling.
REGION_PRIMERS = {
    "V1_V2": {"fwd": "AGAGTTTGATCMTGGCTCAG",  "rev": "GGACCGTGTCTCAGTTCCAG",    "truncLen": (125, 120), "merge": True},
    "V9":    {"fwd": "TGCCACGGTGAATACGTTCC",  "rev": "CCTTGTTACGACTTCACCCCA",  "truncLen": (125, 120), "merge": True},
    "V3_V4": {"fwd": "CCTACGGGNGGCWGCAG",     "rev": "GGACTACHVGGGTATCTAATCC", "truncLen": 125,        "merge": False},
    "V4_V5": {"fwd": "GGAGGGTGCAAGCGTTAATC",  "rev": "TTAACCTTGCGGCCGTACTC",   "truncLen": 125,        "merge": False},
    "V6_V8": {"fwd": "CGGTGGAGCATGTGGTTTAA",  "rev": "AGTTGCAGACTCCAATCCGG",   "truncLen": 125,        "merge": False},
}
ALL_REGIONS    = list(REGION_PRIMERS.keys())
MERGE_REGIONS  = [r for r, c in REGION_PRIMERS.items() if c["merge"]]
SE_REGIONS     = [r for r, c in REGION_PRIMERS.items() if not c["merge"]]
REGION_DIR     = "results/region"
SILVA_SPECIES  = REF_DIR + "/silva_v138.2_assignSpecies.fa.gz"
RDP_TRAINSET   = REF_DIR + "/rdp_19_toSpecies_trainset.fa.gz"

FILTER_INPUT = RAW_DIR if config.get("skip_trimming") else TRIMMED_DIR

# Sample discovery: RAW_DIR/{sample_dir}/*_R1*.fq.gz + *_R2*.fq.gz.
# Same logic as QC_pipe_frag/Snakefile discover_samples() — one subfolder
# per sample, filenames inside need not follow {sample}_R1.fq.gz.
SAMPLE_DIR_GLOB = config.get("sample_dir_glob", "*")
R1_GLOB = config.get("r1_glob", "*[Rr]1*.f*q*.gz")
R2_GLOB = config.get("r2_glob", "*[Rr]2*.f*q*.gz")


def discover_samples():
    samples = {}
    for d in sorted(glob.glob(os.path.join(RAW_DIR, SAMPLE_DIR_GLOB))):
        if not os.path.isdir(d):
            continue
        name = os.path.basename(d)
        r1 = glob.glob(os.path.join(d, R1_GLOB))
        r2 = glob.glob(os.path.join(d, R2_GLOB))
        if not r1:
            print(f"WARNING: {name} - no R1, skipping")
            continue
        if not r2:
            print(f"WARNING: {name} - no R2, skipping")
            continue
        samples[name] = {"R1": r1[0], "R2": r2[0]}
    return samples


SAMPLES_DICT = discover_samples()
SAMPLES = list(SAMPLES_DICT.keys())

if not SAMPLES:
    raise ValueError(f"No samples found in {RAW_DIR} matching {SAMPLE_DIR_GLOB}")


rule all:
    input:
        "results/track.csv",
        "results/seqtab_nochim.rds",
        "results/taxa/silva/taxa_species.rds",
        "results/taxa/rdp/taxa_species.rds",
        "results/abundance_table_silva.csv",
        "results/abundance_table_rdp.csv",
        "results/tree/tree.pdf",
        "results/plots/composition_silva.pdf",
        "results/plots/composition_rdp.pdf",


rule revcomp_primers:
    input:
        fwd = PRIMERS_DIR + "/fwd_primers.fasta",
        rev = PRIMERS_DIR + "/rev_primers.fasta",
    output:
        fwd_rc = temp(PRIMERS_DIR + "/fwd_rc.fasta"),
        rev_rc = temp(PRIMERS_DIR + "/rev_rc.fasta"),
    log:
        "logs/revcomp_primers.log",
    script:
        "scripts/revcomp_primers.py"


rule trim_primers:
    input:
        r1     = lambda wc: SAMPLES_DICT[wc.sample]["R1"],
        r2     = lambda wc: SAMPLES_DICT[wc.sample]["R2"],
        fwd    = PRIMERS_DIR + "/fwd_primers.fasta",
        rev    = PRIMERS_DIR + "/rev_primers.fasta",
        fwd_rc = PRIMERS_DIR + "/fwd_rc.fasta",
        rev_rc = PRIMERS_DIR + "/rev_rc.fasta",
    output:
        r1 = TRIMMED_DIR + "/{sample}_R1.fq.gz",
        r2 = TRIMMED_DIR + "/{sample}_R2.fq.gz",
    log:
        "logs/cutadapt/{sample}.log",
    params:
        threads    = config.get("threads", 4),
        min_len    = config.get("min_len", 50),
        error_rate = config.get("error_rate", 0.1),
    shell:
        "bash scripts/trim_primers.sh "
        "{input.fwd} {input.rev} {input.fwd_rc} {input.rev_rc} "
        "{input.r1} {input.r2} {output.r1} {output.r2} "
        "{params.error_rate} {params.min_len} {params.threads} "
        "> {log} 2>&1"


rule quality_per_sample:
    input:
        r1 = [SAMPLES_DICT[s]["R1"] for s in SAMPLES],
        r2 = [SAMPLES_DICT[s]["R2"] for s in SAMPLES],
    output:
        r1_pdf = f"{QC_DIR}/per_sample_forward.pdf",
        r2_pdf = f"{QC_DIR}/per_sample_reverse.pdf",
    log:
        "logs/quality_per_sample.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/quality_per_sample.R"


rule quality_aggregated:
    input:
        r1 = [SAMPLES_DICT[s]["R1"] for s in SAMPLES],
        r2 = [SAMPLES_DICT[s]["R2"] for s in SAMPLES],
    output:
        r1_pdf = f"{QC_DIR}/aggregated_forward.pdf",
        r2_pdf = f"{QC_DIR}/aggregated_reverse.pdf",
    log:
        "logs/quality_aggregated.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/quality_aggregated.R"


rule filter_reads:
    input:
        r1 = FILTER_INPUT + "/{sample}_R1.fq.gz",
        r2 = FILTER_INPUT + "/{sample}_R2.fq.gz",
    output:
        r1    = QC_DIR + "/{sample}_R1.fq.gz",
        r2    = QC_DIR + "/{sample}_R2.fq.gz",
        stats = "results/filter_stats/{sample}.rds",
    log:
        "logs/filter_reads/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/filter_reads.R"


rule error_correction:
    input:
        r1 = expand(QC_DIR + "/{sample}_R1.fq.gz", sample=SAMPLES),
        r2 = expand(QC_DIR + "/{sample}_R2.fq.gz", sample=SAMPLES),
    output:
        r1_rds = "results/err_forward.rds",
        r2_rds = "results/err_reverse.rds",
        r1_pdf = "results/plots/error_model_forward.pdf",
        r2_pdf = "results/plots/error_model_reverse.pdf",
    log:
        "logs/error_correction.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/error_correction_model.R"


rule dereplication:
    input:
        r1 = QC_DIR + "/{sample}_R1.fq.gz",
        r2 = QC_DIR + "/{sample}_R2.fq.gz",
    output:
        r1_rds = "results/derep/{sample}_R1.rds",
        r2_rds = "results/derep/{sample}_R2.rds",
    log:
        "logs/dereplication/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dereplication.R"


rule dada2_inference:
    input:
        r1_rds     = "results/derep/{sample}_R1.rds",
        r2_rds     = "results/derep/{sample}_R2.rds",
        r1_err_rds = "results/err_forward.rds",
        r2_err_rds = "results/err_reverse.rds",
    output:
        r1_rds = "results/dada/{sample}_R1.rds",
        r2_rds = "results/dada/{sample}_R2.rds",
    log:
        "logs/dada2_inference/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dada2_inference.R"


rule merge_reads:
    input:
        r1_rds_derep = "results/derep/{sample}_R1.rds",
        r2_rds_derep = "results/derep/{sample}_R2.rds",
        r1_rds_dada  = "results/dada/{sample}_R1.rds",
        r2_rds_dada  = "results/dada/{sample}_R2.rds",
    output:
        merged_reads = "results/merged/{sample}.rds",
    log:
        "logs/merge_reads/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/merge_reads.R"


rule make_seqtable:
    input:
        merged_reads = expand("results/merged/{sample}.rds", sample=SAMPLES),
    output:
        sequence_table = "results/seqtab.rds",
    log:
        "logs/make_seqtable.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_seqtable.R"


rule remove_chimera:
    input:
        sequence_table = "results/seqtab.rds",
    output:
        seq_tab_nochim = "results/seqtab_nochim.rds",
    log:
        "logs/remove_chimera.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/remove_chimera.R"


rule create_summary:
    input:
        filter_stats  = expand("results/filter_stats/{sample}.rds", sample=SAMPLES),
        dada_fwd      = expand("results/dada/{sample}_R1.rds", sample=SAMPLES),
        mergers       = expand("results/merged/{sample}.rds", sample=SAMPLES),
        seqtab_nochim = "results/seqtab_nochim.rds",
    output:
        summary = "results/track.csv",
    log:
        "logs/create_summary.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/create_summary.R"


rule assign_taxonomy_silva:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        silva         = REF_DIR + "/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
    output:
        taxa = "results/taxa/silva/taxa.rds",
    log:
        "logs/assign_taxonomy_silva.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/assign_taxa.R"


rule add_species_silva:
    input:
        taxa  = "results/taxa/silva/taxa.rds",
        silva = REF_DIR + "/silva_v138.2_assignSpecies.fa.gz",
    output:
        taxa = "results/taxa/silva/taxa_species.rds",
    log:
        "logs/add_species_silva.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/add_species.R"


rule assign_taxonomy_rdp:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        silva         = REF_DIR + "/rdp_19_toSpecies_trainset.fa.gz",
    output:
        taxa = "results/taxa/rdp/taxa_species.rds",
    log:
        "logs/assign_taxonomy_rdp.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/assign_taxa.R"


rule plot_composition_silva:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        taxa_species  = "results/taxa/silva/taxa_species.rds",
    output:
        plot_pdf = "results/plots/composition_silva.pdf",
    log:
        "logs/plot_composition_silva.log",
    params:
        samples  = SAMPLES,
        top_n    = config.get("top_n", 20),
        width    = config.get("plot_width", 14),
        height   = config.get("plot_height", 7),
        metadata = config.get("metadata", "sra_metadata.csv"),
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/plot_composition.R"


rule plot_composition_rdp:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        taxa_species  = "results/taxa/rdp/taxa_species.rds",
    output:
        plot_pdf = "results/plots/composition_rdp.pdf",
    log:
        "logs/plot_composition_rdp.log",
    params:
        samples  = SAMPLES,
        top_n    = config.get("top_n", 20),
        width    = config.get("plot_width", 14),
        height   = config.get("plot_height", 7),
        metadata = config.get("metadata", "sra_metadata.csv"),
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/plot_composition.R"


rule make_abundance_table_silva:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        taxa_species  = "results/taxa/silva/taxa_species.rds",
    output:
        abundance_table = "results/abundance_table_silva.csv",
    log:
        "logs/make_abundance_table_silva.log",
    params:
        samples = SAMPLES,
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_abundance_table.R"


rule make_abundance_table_rdp:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        taxa_species  = "results/taxa/rdp/taxa_species.rds",
    output:
        abundance_table = "results/abundance_table_rdp.csv",
    log:
        "logs/make_abundance_table_rdp.log",
    params:
        samples = SAMPLES,
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_abundance_table.R"


rule build_tree:
    input:
        seqtab_nochim = "results/seqtab_nochim.rds",
        taxa_species  = "results/taxa/silva/taxa_species.rds",
    output:
        tree_rds = "results/tree/tree.rds",
        tree_nwk = "results/tree/tree.nwk",
        tree_pdf = "results/tree/tree.pdf",
    log:
        "logs/build_tree.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/build_tree.R"


# =====================================================================
# Alignment-based QC (opt-in — not in `rule all`, run explicitly:
#   snakemake --cores N --use-conda qc_align_all
# Maps DADA2-filtered reads (data/qc/{sample}_R[12].fq.gz) against the
# SILVA training-set FASTA with bwa-mem2, then derives generic
# alignment/coverage QC metrics the same way QC_pipe_frag does for WGS.
# Reference headers there are bare taxonomy strings with no accession, and
# many entries share an identical lineage — prep_silva_reference.sh gives
# every record a unique contig id and keeps id->lineage in a side TSV so
# bwa-mem2 gets unique reference names and taxon_mapping_composition.py can
# still recover taxonomy after the fact.
# =====================================================================

rule prep_silva_reference:
    input:
        SILVA_TRAINSET,
    output:
        fasta = REF_DIR + "/silva_for_alignment.fasta",
        map   = REF_DIR + "/silva_for_alignment.map.tsv",
    log:
        "logs/prep_silva_reference.log",
    conda:
        "envs/qc_align.yaml"
    shell:
        "bash scripts/prep_silva_reference.sh {input} {output.fasta} {output.map} "
        "> {log} 2>&1"


rule build_silva_index:
    input:
        REF_DIR + "/silva_for_alignment.fasta",
    output:
        REF_DIR + "/silva_for_alignment.fasta.bwt.2bit.64",
    log:
        "logs/build_silva_index.log",
    conda:
        "envs/qc_align.yaml"
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"


rule align_silva:
    input:
        r1  = QC_DIR + "/{sample}_R1.fq.gz",
        r2  = QC_DIR + "/{sample}_R2.fq.gz",
        ref = REF_DIR + "/silva_for_alignment.fasta",
        idx = REF_DIR + "/silva_for_alignment.fasta.bwt.2bit.64",
    output:
        bam = BAM_DIR + "/{sample}.bam",
        bai = BAM_DIR + "/{sample}.bam.bai",
    log:
        "logs/align_silva/{sample}.log",
    threads: config.get("align_threads", 8)
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        mkdir -p {BAM_DIR}
        bash scripts/align_silva.sh {input.ref} {input.r1} {input.r2} \
            {output.bam} {wildcards.sample} {threads} > {log} 2>&1
        """


rule flagstat:
    input:
        bam = BAM_DIR + "/{sample}.bam",
    output:
        txt = FLAGSTAT_DIR + "/{sample}.flagstat.txt",
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        mkdir -p {FLAGSTAT_DIR}
        samtools flagstat {input.bam} > {output.txt}
        """


rule flagstat_summary:
    input:
        expand(FLAGSTAT_DIR + "/{sample}.flagstat.txt", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/alignment_summary.tsv",
    run:
        with open(output.tsv, "w") as out:
            out.write("Sample\tReads_mapped\tAlignment_reads_%\tReads_duplicated_%\n")
            for f in input:
                sample = os.path.basename(f).replace(".flagstat.txt", "")
                total = mapped = pct = dup = "NA"
                with open(f) as fh:
                    lines = fh.read().splitlines()
                for line in lines:
                    if " in total " in line:
                        total = line.split()[0]
                    if " duplicates" in line:
                        dup = line.split()[0]
                    if " mapped (" in line and "primary" not in line and "mate" not in line:
                        parts = line.split()
                        mapped = parts[0]
                        for p in parts:
                            if "%" in p:
                                pct = p.replace("(", "").replace("%", "")
                                break
                dup_pct = "NA"
                try:
                    dup_pct = f"{float(dup) / float(total) * 100:.2f}"
                except (ValueError, ZeroDivisionError):
                    pass
                out.write(f"{sample}\t{mapped}\t{pct}\t{dup_pct}\n")


rule gc_content:
    input:
        expand(BAM_DIR + "/{sample}.bam", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/gc_content.tsv",
    conda:
        "envs/qc_align.yaml"
    shell:
        "bash scripts/gc_content.sh {BAM_DIR} {output.tsv}"


rule insert_size_summary:
    input:
        expand(BAM_DIR + "/{sample}.bam", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/insert_size_summary.tsv",
    conda:
        "envs/qc_align.yaml"
    shell:
        "bash scripts/insert_size_peak.sh {BAM_DIR} {output.tsv}"


rule read_length_summary:
    input:
        expand(QC_DIR + "/{sample}_R1.fq.gz", sample=SAMPLES),
        expand(QC_DIR + "/{sample}_R2.fq.gz", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/read_length_summary.tsv",
    params:
        samples = SAMPLES,
    shell:
        "bash scripts/read_length_summary.sh {QC_DIR} {output.tsv} {params.samples}"


# SG-protocol linked adapter (Illumina universal, see QC_pipe_frag/config.yaml
# protocols.SG.cutadapt_linked) — override with --config cutadapt_linked_adapter=...
# if a different library-prep adapter was used.
rule adapter_dimers:
    input:
        r1 = lambda wc: SAMPLES_DICT[wc.sample]["R1"],
    output:
        log  = DIMERS_DIR + "/{sample}.log",
        json = DIMERS_DIR + "/{sample}.json",
    params:
        linked = config.get(
            "cutadapt_linked_adapter",
            "ACACTCTTCCCTACACGACGCTCTCCGATCTTT...GATCGGAAGAGCACACGTCTGAACTCCAGTC",
        ),
    threads: config.get("threads", 4)
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        mkdir -p {DIMERS_DIR}
        cutadapt \
            -g "{params.linked}" \
            --discard-untrimmed \
            --minimum-length 20 \
            -o /dev/null \
            --json {output.json} \
            -j {threads} \
            {input.r1} \
            > {output.log} 2>&1
        """


rule adapter_dimers_summary:
    input:
        expand(DIMERS_DIR + "/{sample}.log", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/adapter_dimers_summary.tsv",
    shell:
        """
        echo -e "Sample\\tTotal_reads\\tReads_with_both_adapters\\tAdapter_dimers\\tAdapter_dimers_%" > {output.tsv}
        for log in {input}; do
            sample=$(basename "$log" .log)
            total=$(grep "Total reads processed:" "$log" | awk '{{gsub(",","",$NF); print $NF}}')
            with_adapters=$(grep "Reads with adapters:" "$log" | awk '{{gsub(",","",$4); print $4}}')
            dimers=$(grep "Reads that were too short:" "$log" | awk '{{gsub(",","",$NF); print $NF}}')
            pct=$(awk -v d="$dimers" -v t="$total" 'BEGIN {{if(t>0) printf "%.4f", d/t*100; else print "NA"}}')
            echo -e "${{sample}}\\t${{total}}\\t${{with_adapters}}\\t${{dimers}}\\t${{pct}}" >> {output.tsv}
        done
        """


rule coverage_metrics:
    input:
        bam = BAM_DIR + "/{sample}.bam",
        bai = BAM_DIR + "/{sample}.bam.bai",
    output:
        tsv = COV_DIR + "/{sample}_coverage.tsv",
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        mkdir -p {COV_DIR}
        python scripts/coverage_metrics.py {input.bam} {output.tsv}
        """


rule coverage_summary:
    input:
        expand(COV_DIR + "/{sample}_coverage.tsv", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/coverage_summary.tsv",
    shell:
        """
        head -1 {input[0]} > {output.tsv}
        for f in {input}; do
            tail -n +2 "$f" >> {output.tsv}
        done
        """


rule taxon_mapping_composition:
    input:
        bam = BAM_DIR + "/{sample}.bam",
        bai = BAM_DIR + "/{sample}.bam.bai",
        map = REF_DIR + "/silva_for_alignment.map.tsv",
    output:
        tsv = COMP_DIR + "/{sample}_mapped_composition.tsv",
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        mkdir -p {COMP_DIR}
        python scripts/taxon_mapping_composition.py {input.bam} {input.map} {output.tsv}
        """


rule mapped_composition_summary:
    input:
        expand(COMP_DIR + "/{sample}_mapped_composition.tsv", sample=SAMPLES),
    output:
        tsv = QC_ALIGN_DIR + "/mapped_composition_summary.tsv",
    shell:
        """
        head -1 {input[0]} > {output.tsv}
        for f in {input}; do
            tail -n +2 "$f" >> {output.tsv}
        done
        """


rule build_qc_table:
    input:
        track       = "results/track.csv",
        read_length = QC_ALIGN_DIR + "/read_length_summary.tsv",
        flagstat    = QC_ALIGN_DIR + "/alignment_summary.tsv",
        gc          = QC_ALIGN_DIR + "/gc_content.tsv",
        insert_size = QC_ALIGN_DIR + "/insert_size_summary.tsv",
        dimers      = QC_ALIGN_DIR + "/adapter_dimers_summary.tsv",
        coverage    = QC_ALIGN_DIR + "/coverage_summary.tsv",
    output:
        tsv = QC_ALIGN_DIR + "/final_qc_metrics.tsv",
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        python scripts/build_qc_table.py \
            --track {input.track} \
            --read_length {input.read_length} \
            --flagstat {input.flagstat} \
            --gc {input.gc} \
            --insert_size {input.insert_size} \
            --dimers {input.dimers} \
            --coverage {input.coverage} \
            --out {output.tsv}
        """


rule qc_align_all:
    input:
        QC_ALIGN_DIR + "/final_qc_metrics.tsv",
        QC_ALIGN_DIR + "/mapped_composition_summary.tsv",


# =====================================================================
# Per-region DADA2 (opt-in — not in `rule all`, run explicitly:
#   snakemake --cores N --use-conda region_dada2_all
# Demultiplexes raw reads by V-region primer pair, then runs a separate
# DADA2 pipeline per region with a region-appropriate truncLen. V1_V2/V9
# merge normally; V3_V4/V4_V5/V6_V8 can't physically overlap at 2x151bp
# (see REGION_PRIMERS comment above) so those run forward-read-only.
# =====================================================================

rule demux_region:
    input:
        r1 = lambda wc: SAMPLES_DICT[wc.sample]["R1"],
        r2 = lambda wc: SAMPLES_DICT[wc.sample]["R2"],
    output:
        r1 = REGION_DIR + "/{region}/trimmed/{sample}_R1.fq.gz",
        r2 = REGION_DIR + "/{region}/trimmed/{sample}_R2.fq.gz",
    params:
        fwd = lambda wc: REGION_PRIMERS[wc.region]["fwd"],
        rev = lambda wc: REGION_PRIMERS[wc.region]["rev"],
    log:
        "logs/demux_region/{region}/{sample}.log",
    conda:
        "envs/qc_align.yaml"
    shell:
        """
        mkdir -p $(dirname {output.r1})
        cutadapt -g {params.fwd} -G {params.rev} --discard-untrimmed -e 0.1 \
            --minimum-length 50 -j 4 \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} > {log} 2>&1
        """


# --- filter_reads: PE (merge regions) vs SE (forward-only regions) ---

rule region_filter_reads:
    input:
        r1 = REGION_DIR + "/{region}/trimmed/{sample}_R1.fq.gz",
        r2 = REGION_DIR + "/{region}/trimmed/{sample}_R2.fq.gz",
    output:
        r1    = REGION_DIR + "/{region}/qc/{sample}_R1.fq.gz",
        r2    = REGION_DIR + "/{region}/qc/{sample}_R2.fq.gz",
        stats = REGION_DIR + "/{region}/filter_stats/{sample}.rds",
    params:
        truncLen = lambda wc: REGION_PRIMERS[wc.region]["truncLen"],
    wildcard_constraints:
        region = "|".join(MERGE_REGIONS),
    log:
        "logs/region_filter_reads/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/region_filter_reads.R"


rule region_filter_reads_se:
    input:
        r1 = REGION_DIR + "/{region}/trimmed/{sample}_R1.fq.gz",
    output:
        r1    = REGION_DIR + "/{region}/qc/{sample}_R1.fq.gz",
        stats = REGION_DIR + "/{region}/filter_stats/{sample}.rds",
    params:
        truncLen = lambda wc: REGION_PRIMERS[wc.region]["truncLen"],
    wildcard_constraints:
        region = "|".join(SE_REGIONS),
    log:
        "logs/region_filter_reads/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/region_filter_reads_se.R"


# --- error model: PE vs SE ---

rule region_error_correction:
    input:
        r1 = expand(REGION_DIR + "/{{region}}/qc/{sample}_R1.fq.gz", sample=SAMPLES),
        r2 = expand(REGION_DIR + "/{{region}}/qc/{sample}_R2.fq.gz", sample=SAMPLES),
    output:
        r1_rds = REGION_DIR + "/{region}/err_forward.rds",
        r2_rds = REGION_DIR + "/{region}/err_reverse.rds",
        r1_pdf = REGION_DIR + "/{region}/plots/error_model_forward.pdf",
        r2_pdf = REGION_DIR + "/{region}/plots/error_model_reverse.pdf",
    wildcard_constraints:
        region = "|".join(MERGE_REGIONS),
    log:
        "logs/region_error_correction/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/error_correction_model.R"


rule region_error_correction_se:
    input:
        r1 = expand(REGION_DIR + "/{{region}}/qc/{sample}_R1.fq.gz", sample=SAMPLES),
    output:
        r1_rds = REGION_DIR + "/{region}/err_forward.rds",
        r1_pdf = REGION_DIR + "/{region}/plots/error_model_forward.pdf",
    wildcard_constraints:
        region = "|".join(SE_REGIONS),
    log:
        "logs/region_error_correction/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/error_correction_model_se.R"


# --- dereplication: PE vs SE ---

rule region_dereplication:
    input:
        r1 = REGION_DIR + "/{region}/qc/{sample}_R1.fq.gz",
        r2 = REGION_DIR + "/{region}/qc/{sample}_R2.fq.gz",
    output:
        r1_rds = REGION_DIR + "/{region}/derep/{sample}_R1.rds",
        r2_rds = REGION_DIR + "/{region}/derep/{sample}_R2.rds",
    wildcard_constraints:
        region = "|".join(MERGE_REGIONS),
    log:
        "logs/region_dereplication/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dereplication.R"


rule region_dereplication_se:
    input:
        r1 = REGION_DIR + "/{region}/qc/{sample}_R1.fq.gz",
    output:
        r1_rds = REGION_DIR + "/{region}/derep/{sample}_R1.rds",
    wildcard_constraints:
        region = "|".join(SE_REGIONS),
    log:
        "logs/region_dereplication/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dereplication_se.R"


# --- dada inference: PE vs SE ---

rule region_dada2_inference:
    input:
        r1_rds     = REGION_DIR + "/{region}/derep/{sample}_R1.rds",
        r2_rds     = REGION_DIR + "/{region}/derep/{sample}_R2.rds",
        r1_err_rds = REGION_DIR + "/{region}/err_forward.rds",
        r2_err_rds = REGION_DIR + "/{region}/err_reverse.rds",
    output:
        r1_rds = REGION_DIR + "/{region}/dada/{sample}_R1.rds",
        r2_rds = REGION_DIR + "/{region}/dada/{sample}_R2.rds",
    wildcard_constraints:
        region = "|".join(MERGE_REGIONS),
    log:
        "logs/region_dada2_inference/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dada2_inference.R"


rule region_dada2_inference_se:
    input:
        r1_rds     = REGION_DIR + "/{region}/derep/{sample}_R1.rds",
        r1_err_rds = REGION_DIR + "/{region}/err_forward.rds",
    output:
        r1_rds = REGION_DIR + "/{region}/dada/{sample}_R1.rds",
    wildcard_constraints:
        region = "|".join(SE_REGIONS),
    log:
        "logs/region_dada2_inference/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/dada2_inference_se.R"


# --- seqtable: merge (PE) vs forward-only (SE) ---

rule region_merge_reads:
    input:
        r1_rds_derep = REGION_DIR + "/{region}/derep/{sample}_R1.rds",
        r2_rds_derep = REGION_DIR + "/{region}/derep/{sample}_R2.rds",
        r1_rds_dada  = REGION_DIR + "/{region}/dada/{sample}_R1.rds",
        r2_rds_dada  = REGION_DIR + "/{region}/dada/{sample}_R2.rds",
    output:
        merged_reads = REGION_DIR + "/{region}/merged/{sample}.rds",
    wildcard_constraints:
        region = "|".join(MERGE_REGIONS),
    log:
        "logs/region_merge_reads/{region}/{sample}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/merge_reads.R"


rule region_make_seqtable:
    input:
        merged_reads = expand(REGION_DIR + "/{{region}}/merged/{sample}.rds", sample=SAMPLES),
    output:
        sequence_table = REGION_DIR + "/{region}/seqtab.rds",
    wildcard_constraints:
        region = "|".join(MERGE_REGIONS),
    log:
        "logs/region_make_seqtable/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_seqtable.R"


rule region_make_seqtable_se:
    input:
        dada_fwd = expand(REGION_DIR + "/{{region}}/dada/{sample}_R1.rds", sample=SAMPLES),
    output:
        sequence_table = REGION_DIR + "/{region}/seqtab.rds",
    wildcard_constraints:
        region = "|".join(SE_REGIONS),
    log:
        "logs/region_make_seqtable/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_seqtable_se.R"


# --- everything past here is region-agnostic (same format seqtab either way) ---

rule region_remove_chimera:
    input:
        sequence_table = REGION_DIR + "/{region}/seqtab.rds",
    output:
        seq_tab_nochim = REGION_DIR + "/{region}/seqtab_nochim.rds",
    log:
        "logs/region_remove_chimera/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/remove_chimera.R"


rule region_assign_taxonomy_silva:
    input:
        seqtab_nochim = REGION_DIR + "/{region}/seqtab_nochim.rds",
        silva         = SILVA_TRAINSET,
    output:
        taxa = REGION_DIR + "/{region}/taxa/silva/taxa.rds",
    log:
        "logs/region_assign_taxonomy_silva/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/assign_taxa.R"


rule region_add_species_silva:
    input:
        taxa  = REGION_DIR + "/{region}/taxa/silva/taxa.rds",
        silva = SILVA_SPECIES,
    output:
        taxa = REGION_DIR + "/{region}/taxa/silva/taxa_species.rds",
    log:
        "logs/region_add_species_silva/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/add_species.R"


rule region_assign_taxonomy_rdp:
    input:
        seqtab_nochim = REGION_DIR + "/{region}/seqtab_nochim.rds",
        silva         = RDP_TRAINSET,
    output:
        taxa = REGION_DIR + "/{region}/taxa/rdp/taxa_species.rds",
    log:
        "logs/region_assign_taxonomy_rdp/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/assign_taxa.R"


rule region_abundance_table_silva:
    input:
        seqtab_nochim = REGION_DIR + "/{region}/seqtab_nochim.rds",
        taxa_species  = REGION_DIR + "/{region}/taxa/silva/taxa_species.rds",
    output:
        abundance_table = REGION_DIR + "/{region}/abundance_table_silva.csv",
    params:
        samples = SAMPLES,
    log:
        "logs/region_abundance_table_silva/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_abundance_table.R"


rule region_abundance_table_rdp:
    input:
        seqtab_nochim = REGION_DIR + "/{region}/seqtab_nochim.rds",
        taxa_species  = REGION_DIR + "/{region}/taxa/rdp/taxa_species.rds",
    output:
        abundance_table = REGION_DIR + "/{region}/abundance_table_rdp.csv",
    params:
        samples = SAMPLES,
    log:
        "logs/region_abundance_table_rdp/{region}.log",
    conda:
        "envs/dada2.yaml"
    script:
        "scripts/make_abundance_table.R"


rule region_dada2_all:
    input:
        expand(REGION_DIR + "/{region}/abundance_table_silva.csv", region=ALL_REGIONS),
        expand(REGION_DIR + "/{region}/abundance_table_rdp.csv", region=ALL_REGIONS),
