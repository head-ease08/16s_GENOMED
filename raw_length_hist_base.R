#!/usr/bin/env Rscript
# Read-length histograms for raw FASTQs, base R only (ggplot2 for plotting).
# Usage: Rscript raw_length_hist_base.R <raw_dir> <cache.tsv> <out.pdf> <R1|R2>

suppressMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "data/raw"
cache   <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "all_lengths_R1.tsv"
out_pdf <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "raw_lengths.pdf"
which_r <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "R1"

if (!file.exists(cache)) {
    cat("Building length table from FASTQs...\n")
    pattern <- paste0("*", which_r, "*.fq.gz")
    fqs <- Sys.glob(file.path(raw_dir, "*", pattern))
    if (length(fqs) == 0) stop("No FASTQs found under ", raw_dir)

    all_rows <- list()
    for (f in fqs) {
        sample <- basename(dirname(f))
        cat("  ", sample, "\n")
        cmd <- sprintf(
            "zcat %s | awk 'NR%%4==2 {print length($0)}' | sort -n | uniq -c",
            shQuote(f))
        con <- pipe(cmd, "r")
        lines <- readLines(con); close(con)
        for (ln in lines) {
            p <- strsplit(trimws(ln), "\\s+")[[1]]
            all_rows[[length(all_rows)+1]] <- data.frame(
                sample = sample,
                length = as.integer(p[2]),
                count  = as.integer(p[1]),
                stringsAsFactors = FALSE
            )
        }
    }
    df <- do.call(rbind, all_rows)
    # exclude non-sample directories that happen to have R1 files (e.g. fastp_output)
    df <- df[!df$sample %in% c("fastp_output"), , drop = FALSE]
    write.table(df, cache, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Cached to ", cache, "\n")
} else {
    cat("Reading cache: ", cache, "\n")
    # autodetect: if first field of line 1 isn't "sample", treat as no-header
    first <- readLines(cache, n = 1)
    has_header <- startsWith(first, "sample\t")
    df <- read.table(cache, header = has_header, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE,
                     col.names = if (!has_header) c("sample","length","count") else NULL)
    df <- df[!df$sample %in% c("fastp_output"), , drop = FALSE]
    cat("Columns found: ", paste(colnames(df), collapse=", "), "\n")
    cat("Rows: ", nrow(df), "\n")
    print(head(df, 3))
    if (!all(c("sample", "length", "count") %in% colnames(df))) {
        stop("Cache file malformed. Delete it and re-run: rm ", cache)
    }
    df$length <- as.integer(df$length)
    df$count  <- as.integer(df$count)
}

df$sample <- factor(df$sample, levels = sort(unique(df$sample)))

# per-sample % normalization (base R)
totals <- tapply(df$count, df$sample, sum)
df$pct <- 100 * df$count / totals[as.character(df$sample)]

p_overview <- ggplot(df, aes(x = length, y = pct, group = sample, color = sample)) +
    geom_line(alpha = 0.5) +
    labs(title = paste("Raw read length distribution -", which_r),
         x = "read length (bp)", y = "% of reads in sample") +
    theme_minimal() +
    theme(legend.position = "none")

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

cat("Wrote ", out_pdf, "\n\n")

# Summary в консоль, base R
cat("=== per-sample summary ===\n")
summary_df <- do.call(rbind, lapply(split(df, df$sample), function(x) {
    data.frame(
        sample     = as.character(x$sample[1]),
        n_reads    = sum(x$count),
        mode_len   = x$length[which.max(x$count)],
        pct_mode   = round(100 * max(x$count) / sum(x$count), 1),
        min_len    = min(x$length),
        max_len    = max(x$length),
        stringsAsFactors = FALSE
    )
}))
print(summary_df, row.names = FALSE)
