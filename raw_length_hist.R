#!/usr/bin/env Rscript
# Read-length histograms for raw FASTQs. One PDF, one facet per sample.
#
# Usage:
#   Rscript raw_length_hist.R data/raw all_lengths.tsv raw_lengths.pdf R1
#     args: <raw_dir> <input_tsv_cache> <output_pdf> <R1|R2>
#
# If <input_tsv_cache> exists, it's read directly (skip re-scanning FASTQs).
# Otherwise, walks data/raw/<sample>/*R1*.fq.gz (or R2) and builds it.

suppressMessages({
    library(ggplot2)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[1] %||% "data/raw"
cache   <- args[2] %||% "all_lengths_R1.tsv"
out_pdf <- args[3] %||% "raw_lengths.pdf"
which_r <- args[4] %||% "R1"

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

if (!file.exists(cache)) {
    cat("Building length table from FASTQs...\n")
    pattern <- paste0("*", which_r, "*.fq.gz")
    fqs <- Sys.glob(file.path(raw_dir, "*", pattern))
    if (length(fqs) == 0) stop("No FASTQs found under ", raw_dir)
    rows <- list()
    for (f in fqs) {
        sample <- basename(dirname(f))
        cat("  ", sample, "\n")
        # awk: длина каждой 2-й строки в каждой 4-строке FASTQ
        cmd <- sprintf(
            "zcat %s | awk 'NR%%4==2 {print length($0)}' | sort -n | uniq -c",
            shQuote(f))
        con <- pipe(cmd, "r")
        lines <- readLines(con); close(con)
        for (ln in lines) {
            p <- strsplit(trimws(ln), "\\s+")[[1]]
            rows[[length(rows)+1]] <- data.frame(
                sample = sample, length = as.integer(p[2]), count = as.integer(p[1])
            )
        }
    }
    df <- do.call(rbind, rows)
    write.table(df, cache, sep="\t", row.names=FALSE, quote=FALSE)
    cat("Cached to ", cache, "\n")
} else {
    cat("Reading cache: ", cache, "\n")
    df <- read.table(cache, header=TRUE, sep="\t")
}

df$sample <- factor(df$sample, levels = sort(unique(df$sample)))

# Общий overview: все samples на одном plot'е (mean-normalized frequencies)
p_overview <- df %>%
    group_by(sample) %>%
    mutate(pct = 100 * count / sum(count)) %>%
    ggplot(aes(x = length, y = pct, group = sample, color = sample)) +
    geom_line(alpha = 0.5) +
    labs(title = paste("Raw read length distribution -", which_r),
         x = "read length (bp)", y = "% of reads in sample") +
    theme_minimal() +
    theme(legend.position = "none")

# Faceted: одна panel на sample
p_facets <- ggplot(df, aes(x = length, y = count)) +
    geom_col(width = 1) +
    facet_wrap(~ sample, scales = "free_y", ncol = 4) +
    labs(title = paste("Per-sample length histogram -", which_r),
         x = "read length (bp)", y = "reads") +
    theme_minimal(base_size = 8)

pdf(out_pdf, width = 14, height = 10)
print(p_overview)
print(p_facets)
dev.off()

cat("Wrote ", out_pdf, "\n")

# Печатаем summary в консоль
cat("\n=== per-sample length summary ===\n")
df %>%
    group_by(sample) %>%
    summarize(
        n_reads    = sum(count),
        min_len    = min(length),
        median_len = round(sum(rep(length, count) > 0) / 2),  # approx
        mode_len   = length[which.max(count)],
        max_len    = max(length),
        pct_at_max = round(100 * max(count) / sum(count), 1)
    ) %>%
    print(n = 100)
