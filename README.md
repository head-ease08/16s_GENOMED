# 16S rRNA Amplicon Sequencing Pipeline

Snakemake + DADA2 pipeline: paired-end Illumina 16S reads → chimera-free ASV table with dual taxonomy (SILVA + RDP), abundance tables, composition plots, and phylogenetic tree.

---

## Requirements

| Tool | Min version | Notes |
|---|---|---|
| Snakemake | >= 7.0 | pipeline manager |
| conda / mamba | any | environment management |
| cutadapt | >= 3.0 | `-j` threads and `file:` primer pools require >= 3.0 |
| Python | >= 3.8 | Snakemake dependency |
| R packages | — | installed automatically via `envs/dada2.yaml` |

```bash
# install Snakemake
conda install -c bioconda -c conda-forge snakemake mamba

# install cutadapt
conda install -c bioconda "cutadapt>=3.0"
```

---

## Project structure

```
16s_GM_local/
├── Snakefile
├── README.md
├── .gitignore
│
├── primers/
│   ├── fwd_primers.fasta        # forward primer(s) in FASTA format
│   └── rev_primers.fasta        # reverse primer(s)
│
├── references/                  # not tracked by git — place DB files here
│   ├── silva_nr99_v138.2_toSpecies_trainset.fa.gz
│   ├── silva_v138.2_assignSpecies.fa.gz
│   └── rdp_19_toSpecies_trainset.fa.gz
│
├── data/
│   ├── raw/                     # input FASTQ (read-only)
│   ├── trimmed/                 # cutadapt output
│   └── qc/                      # filterAndTrim output
│
├── results/
│   ├── filter_stats/            # per-sample filter stats (.rds)
│   ├── derep/                   # dereplication objects (.rds)
│   ├── dada/                    # DADA2 denoised objects (.rds)
│   ├── merged/                  # merged paired reads (.rds)
│   ├── plots/                   # error model + composition PDFs
│   ├── taxa/
│   │   ├── silva/               # taxa.rds, taxa_species.rds
│   │   └── rdp/                 # taxa_species.rds
│   ├── tree/                    # tree.rds, tree.nwk, tree.pdf
│   ├── seqtab.rds
│   ├── seqtab_nochim.rds
│   ├── abundance_table_silva.csv
│   ├── abundance_table_rdp.csv
│   └── track.csv
│
├── logs/                        # per-rule log files
│   ├── cutadapt/
│   ├── filter_reads/
│   ├── dereplication/
│   ├── dada2_inference/
│   ├── merge_reads/
│   └── ...
│
├── scripts/                     # R and shell scripts (Snakemake-driven)
├── util/                        # standalone utility scripts
└── envs/
    └── dada2.yaml               # conda environment spec
```

---

## Input format

Raw reads go in `data/raw/`. Required naming pattern:

```
{sample}_R1.fq.gz
{sample}_R2.fq.gz
```

`{sample}` must not contain underscores before `_R1` / `_R2`. Example:
```
SRR123456_R1.fq.gz   SRR123456_R2.fq.gz
patient01_R1.fq.gz   patient01_R2.fq.gz
```

The pipeline asserts at startup that every `_R1.fq.gz` has a matching `_R2.fq.gz`.

---

## Reference databases

Place files in `references/`. Download from Zenodo (SILVA v138.2 / RDP 19):

```bash
mkdir -p references

# SILVA — genus-level training set (assignTaxonomy)
wget -P references/ https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz

# SILVA — species assignment (addSpecies)
wget -P references/ https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz

# RDP 19 — genus-to-species training set (assignTaxonomy)
wget -P references/ https://zenodo.org/records/14169026/files/rdp_19_toSpecies_trainset.fa.gz
```

---

## Primer files

Each file is a standard FASTA. Multiple primers per file are supported (cutadapt `file:` pools).

```
# primers/fwd_primers.fasta  (example: V3-V4 341F)
>341F
CCTACGGGNGGCWGCAG

# primers/rev_primers.fasta  (example: V3-V4 805R)
>805R
GACTACHVGGGTATCTAATCC
```

IUPAC ambiguity codes are handled. The pipeline generates reverse complements automatically (`revcomp_primers` rule).

---

## Running the pipeline

### Dry run (check DAG, no execution)
```bash
snakemake --cores 8 --use-conda -n
```

### Full run
```bash
snakemake --cores 8 --use-conda
```

### Skip primer trimming (reads already trimmed)
```bash
snakemake --cores 8 --use-conda --config skip_trimming=true
```
When `skip_trimming=true`, `filter_reads` reads directly from `data/raw/` instead of `data/trimmed/`.

### Run only specific targets
```bash
# quality plots only
snakemake --cores 4 --use-conda data/qc/per_sample_forward.pdf

# up to chimera removal
snakemake --cores 8 --use-conda results/seqtab_nochim.rds

# SILVA taxonomy only (skip RDP + tree + plots)
snakemake --cores 8 --use-conda results/taxa/silva/taxa_species.rds

# single sample through dereplication
snakemake --cores 4 --use-conda results/derep/sample01_R1.rds
```

### Resume after failure
```bash
snakemake --cores 8 --use-conda --rerun-incomplete
```

### Cluster / HPC submission
```bash
snakemake --cores 64 --use-conda \
    --cluster "sbatch -c {threads} --mem={resources.mem_mb}M" \
    --jobs 32
```

---

## Config flags

Passed via `--config key=value` or a `--configfile config.yaml`.

| Key | Default | Description |
|---|---|---|
| `skip_trimming` | `false` | Skip cutadapt; read from `data/raw/` directly |
| `threads` | `4` | Threads for cutadapt per sample |
| `min_len` | `50` | Minimum read length after primer trimming |
| `error_rate` | `0.1` | Cutadapt max error rate for primer matching |
| `top_n` | `20` | Top N genera shown in composition plots |
| `plot_width` | `14` | Composition plot width (inches) |
| `plot_height` | `7` | Composition plot height (inches) |
| `metadata` | `sra_metadata.csv` | CSV with `Run` and `Sample Name` columns for plot grouping |

Example `config.yaml`:
```yaml
threads: 8
min_len: 100
error_rate: 0.15
top_n: 30
metadata: my_metadata.csv
```

```bash
snakemake --cores 8 --use-conda --configfile config.yaml
```

---

## Pipeline stages

| # | Snakemake rule | Script | Output |
|---|---|---|---|
| 1 | `revcomp_primers` | `revcomp_primers.py` | `primers/fwd_rc.fasta`, `primers/rev_rc.fasta` (temp) |
| 2 | `trim_primers` | `trim_primers.sh` (cutadapt) | `data/trimmed/{sample}_R[12].fq.gz` |
| 3 | `quality_per_sample` | `quality_per_sample.R` | `data/qc/per_sample_forward.pdf`, `per_sample_reverse.pdf` |
| 4 | `quality_aggregated` | `quality_aggregated.R` | `data/qc/aggregated_forward.pdf`, `aggregated_reverse.pdf` |
| 5 | `filter_reads` | `filter_reads.R` | `data/qc/{sample}_R[12].fq.gz`, `results/filter_stats/{sample}.rds` |
| 6 | `error_correction` | `error_correction_model.R` | `results/err_forward.rds`, `err_reverse.rds`, error model PDFs |
| 7 | `dereplication` | `dereplication.R` | `results/derep/{sample}_R[12].rds` |
| 8 | `dada2_inference` | `dada2_inference.R` | `results/dada/{sample}_R[12].rds` |
| 9 | `merge_reads` | `merge_reads.R` | `results/merged/{sample}.rds` |
| 10 | `make_seqtable` | `make_seqtable.R` | `results/seqtab.rds` |
| 11 | `remove_chimera` | `remove_chimera.R` | `results/seqtab_nochim.rds` |
| 12 | `create_summary` | `create_summary.R` | `results/track.csv` |
| 13 | `assign_taxonomy_silva` | `assign_taxa.R` | `results/taxa/silva/taxa.rds` |
| 14 | `add_species_silva` | `add_species.R` | `results/taxa/silva/taxa_species.rds` |
| 15 | `assign_taxonomy_rdp` | `assign_taxa.R` | `results/taxa/rdp/taxa_species.rds` |
| 16 | `plot_composition_silva` | `plot_composition.R` | `results/plots/composition_silva.pdf` |
| 17 | `plot_composition_rdp` | `plot_composition.R` | `results/plots/composition_rdp.pdf` |
| 18 | `make_abundance_table_silva` | `make_abundance_table.R` | `results/abundance_table_silva.csv` |
| 19 | `make_abundance_table_rdp` | `make_abundance_table.R` | `results/abundance_table_rdp.csv` |
| 20 | `build_tree` | `build_tree.R` | `results/tree/tree.rds`, `tree.nwk`, `tree.pdf` |

---

## Key tunable parameters (inside scripts)

These are hardcoded in the R scripts — edit directly or parameterize via Snakemake `params:` if needed.

| Parameter | Script | Default | Notes |
|---|---|---|---|
| `truncLen` | `filter_reads.R` | `c(240, 200)` | Set based on quality plots; must retain overlap region |
| `maxEE` | `filter_reads.R` | `c(2, 2)` | Lower = stricter, fewer reads pass |
| `truncQ` | `filter_reads.R` | `2` | Truncate at first base below this Q score |
| `rm.phix` | `filter_reads.R` | `TRUE` | Remove PhiX spike-in reads |
| `MAX_CONSIST` | `error_correction_model.R` | `20` | Max EM iterations for error model |
| `pool` | `dada2_inference.R` | `FALSE` | `TRUE` or `"pseudo"` increases sensitivity, much slower |
| `minOverlap` | `merge_reads.R` | `12` | Minimum overlap for R1/R2 merge |
| `maxMismatch` | `merge_reads.R` | `0` | Mismatches allowed in merge overlap |

---

## Outputs

### `results/track.csv`
Read counts at each pipeline stage per sample. Use to diagnose where reads are lost.

| Column | Description |
|---|---|
| `sample` | Sample name |
| `input` | Reads entering `filterAndTrim` |
| `filtered` | Reads passing quality filter |
| `denoised` | Unique sequences after DADA2 denoising |
| `merged` | Reads successfully merged (R1 + R2) |
| `nonchim` | Reads remaining after chimera removal |

### `results/abundance_table_silva.csv` / `abundance_table_rdp.csv`
One row per ASV. Columns: `asv_id`, `sequence`, taxonomy columns (`Kingdom` → `Species`), raw count per sample, relative abundance (%) per sample.

### `results/taxa/silva/taxa_species.rds` / `results/taxa/rdp/taxa_species.rds`
R matrix: rows = ASV sequences, columns = taxonomic ranks. Load in R with `readRDS()`.

### `results/seqtab_nochim.rds`
Chimera-free sequence table. Rows = samples, columns = ASV sequences.

### `results/tree/`
- `tree.nwk` — Newick format phylogenetic tree (NJ/GTR via DECIPHER)
- `tree.rds` — `phylo` object for use with phyloseq
- `tree.pdf` — annotated tree visualization

### `results/plots/`
- `error_model_forward.pdf` / `error_model_reverse.pdf` — learned error model diagnostics
- `composition_silva.pdf` / `composition_rdp.pdf` — stacked bar charts, top N genera by V-region group

---

## Logs

Every rule writes to `logs/<rule>/<sample>.log` or `logs/<rule>.log`. Check these first on failure:

```bash
# example: why did filter_reads fail for sample01?
cat logs/filter_reads/sample01.log

# error correction log
cat logs/error_correction.log
```

---

## Utility scripts (`util/`)

These are standalone scripts, not part of the Snakemake DAG. Run manually from the project root.

### `util/uniq.R`
Detects genera present in mock community samples (Zymo / ZIEL) using `sra_metadata.csv` for SRR→sample mapping.

```bash
Rscript util/uniq.R
```

Reads `results/abundance_table_silva.csv`.

### `util/uniq2.R`
Same as `uniq.R` but reads SRR lists from `/tmp/{mock}_srr.txt` files instead of metadata CSV.

```bash
# create SRR list files first
echo "SRR123456" > /tmp/zymo_srr.txt
Rscript util/uniq2.R
```

### `util/u_tsv.R`
Compares detected genera across Zymo/ZIEL mock samples for both SILVA and RDP databases. Outputs `results/mock_detected_taxa.tsv`.

```bash
Rscript util/u_tsv.R
```

Output columns: `mock`, `db`, `Genus`, `Species`.

---

## Interactive notebook

`scripts/dada2_paired_with_tree.Rmd` — RMarkdown notebook that runs the full pipeline interactively (no Snakemake). Useful for exploratory analysis and debugging. Expects reads in `./` matching `_R1*.fq` pattern. Run in RStudio or:

```bash
Rscript -e "rmarkdown::render('scripts/dada2_paired_with_tree.Rmd')"
```

---

## Environment

The conda environment `envs/dada2.yaml` provides all R dependencies. Snakemake creates it automatically on first run with `--use-conda`.

To create manually:
```bash
conda env create -f envs/dada2.yaml
conda activate dada2
```

Packages included: `dada2 >= 1.28`, `ggplot2`, `tidyr`, `DECIPHER`, `ape`, `ggtree`.

---

## Troubleshooting

**`No paired R2 for {sample}`** — R2 file missing or naming mismatch. Check `data/raw/` contains both `_R1.fq.gz` and `_R2.fq.gz` for every sample.

**All reads lost at `filter_reads`** — `truncLen` too aggressive. Inspect quality plots and lower values. R1 + R2 truncated lengths must sum to > amplicon length + required overlap.

**Low merge rate** — overlap too short (reduce `truncLen`) or `minOverlap` too high. Check amplicon expected length.

**`learnErrors` convergence warning** — increase `MAX_CONSIST` in `error_correction_model.R`, or ensure enough reads (>= 1M bases recommended per read direction).

**SILVA/RDP rule fails with path error** — reference files must be in `references/` with exact filenames shown in the [Reference databases](#reference-databases) section.

**Composition plot fails with `sra_metadata.csv` not found** — provide path via `--config metadata=/path/to/file.csv`, or skip this rule and use abundance tables directly.
