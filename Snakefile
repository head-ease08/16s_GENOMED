from pathlib import Path
import os
import glob

RAW_DIR     = config.get("run_dir", "data/raw")
QC_DIR      = "data/qc"
TRIMMED_DIR = "data/trimmed"
PRIMERS_DIR = "primers"
REF_DIR     = "references"

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
