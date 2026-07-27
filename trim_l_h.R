#!/usr/bin/env Rscript
# Read-length histograms for FASTQs in a FLAT directory (e.g. data/trimmed/).
# Filenames like: SAMPLE_R1.fq.gz  or  SAMPLE_R2.fq.gz
#
# Usage:
#   Rscript trimmed_length_hist.R data/trimmed all_lengths_trimmed_R1.tsv trimmed_lengths_R1.pdf R1
#   Rscript trimmed_length_hist.R data/trimmed all_lengths_trimmed_R2.tsv trimmed_lengths_R2.pdf R2

suppressMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "data/trimmed"
cache   <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "all_lengths_trimmed_R1.tsv"
out_pdf <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "trimmed_lengths.pdf"
which_r <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "R1"

if (!file.exists(cache)) {
    cat("Building length table from FASTQs in flat dir...\n")
    fqs <- Sys.glob(file.path(in_dir, paste0("*_", which_r, ".fq.gz")))
    if (length(fqs) == 0) stop("No files matched: ", file.path(in_dir, paste0("*_", which_r, ".fq.gz")))

    all_rows <- list()
    for (f in fqs) {
        # filename → sample name (strip _R1.fq.gz / _R2.fq.gz)
        sample <- sub(paste0("_", which_r, "\\.fq\\.gz$"), "", basename(f))
        cat("  ", sample, "\n")
        cmd <- sprintf(
            "zcat %s | awk 'NR%%4==2 {print length($0)}' | sort -n | uniq -c",
            shQuote(f))
        con <- pipe(cmd, "r")
        lines <- readLines(con); close(con)
        for (ln in lines) {
            p <- strsplit(trimws(ln), "\\s+")[[1]]
            all_rows[[length(all_rows)+1]] <- data.frame(
                sample = sample, length = as.integer(p[2]), count = as.integer(p[1]),
                stringsAsFactors = FALSE
            )
        }
    }
    df <- do.call(rbind, all_rows)
    write.table(df, cache, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Cached to ", cache, "\n")
} else {
    cat("Reading cache: ", cache, "\n")
    first <- readLines(cache, n = 1)
    has_header <- startsWith(first, "sample\t")
    if (has_header) {
        df <- read.table(cache, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        df <- read.table(cache, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE,
                         col.names = c("sample", "length", "count"))
    }
}

df$sample <- factor(df$sample, levels = sort(unique(df$sample)))

# per-sample % normalization
totals <- tapply(df$count, df$sample, sum)
df$pct <- 100 * df$count / totals[as.character(df$sample)]

p_overview <- ggplot(df, aes(x = length, y = pct, group = sample, color = sample)) +
    geom_line(alpha = 0.5) +
    geom_vline(xintercept = 140, linetype = "dashed", color = "red") +
    annotate("text", x = 122, y = Inf, label = "truncLen=140", hjust = 0, vjust = 1.5,
             color = "red", size = 3) +
    labs(title = paste("Post-cutadapt read length distribution -", which_r),
         subtitle = "reads to the LEFT of the red line are discarded by filterAndTrim",
         x = "read length (bp)", y = "% of reads in sample") +
    theme_minimal() +
    theme(legend.position = "none")

p_facets <- ggplot(df, aes(x = length, y = count)) +
    geom_col(width = 1) +
    geom_vline(xintercept = 140, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_wrap(~ sample, scales = "free_y", ncol = 4) +
    labs(title = paste("Per-sample post-cutadapt length histogram -", which_r),
         subtitle = "red dashed line = truncLen=140 cutoff",
         x = "read length (bp)", y = "reads") +
    theme_minimal(base_size = 8)

pdf(out_pdf, width = 14, height = 10)
print(p_overview)
print(p_facets)
dev.off()

cat("Wrote ", out_pdf, "\n\n")

cat("=== per-sample summary (post-cutadapt) ===\n")
summary_df <- do.call(rbind, lapply(split(df, df$sample), function(x) {
    tot <- sum(x$count)
    below_140 <- sum(x$count[x$length < 140])
    data.frame(
        sample     = as.character(x$sample[1]),
        n_reads    = tot,
        mode_len   = x$length[which.max(x$count)],
        min_len    = min(x$length),
        max_len    = max(x$length),
        below_140  = below_140,
        pct_below_140 = round(100 * below_140 / tot, 1),
        stringsAsFactors = FALSE
    )
}))
print(summary_df, row.names = FALSE)
