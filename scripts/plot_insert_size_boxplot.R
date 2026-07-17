#!/usr/bin/env Rscript
# Cross-sample insert size comparison: amplicon libraries have tight,
# narrow-range per-sample distributions (fixed primer positions), so a grid
# of near-identical single-peak histograms is a worse comparison view than
# one boxplot with all samples side by side, sorted by median. Also writes a
# numeric summary TSV (n / median / IQR / mean / sd per sample) since a plot
# alone won't tell you if e.g. one sample's peak is exactly 101 vs 103.
#
# Usage: plot_insert_size_boxplot.R <dir_with_insert_sizes_txt>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_insert_size_boxplot.R <dir>")
out_dir <- args[1]

files <- list.files(out_dir, pattern = "\\.insert_sizes\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No *.insert_sizes.txt files found in ", out_dir)

samples <- sub("\\.insert_sizes\\.txt$", "", basename(files))
data <- setNames(lapply(files, function(f) {
    x <- scan(f, what = numeric(), quiet = TRUE)
    x[is.finite(x) & x > 0]
}), samples)

medians <- sapply(data, function(x) if (length(x)) median(x) else NA)
ord <- order(medians, na.last = TRUE)
samples_ord <- samples[ord]
data_ord <- data[samples_ord]

# summary TSV
summary_df <- do.call(rbind, lapply(samples_ord, function(s) {
    x <- data[[s]]
    if (length(x) == 0) {
        data.frame(sample = s, n = 0, min = NA, q1 = NA, median = NA,
                   q3 = NA, max = NA, mean = NA, sd = NA)
    } else {
        q <- quantile(x, c(0.25, 0.75), names = FALSE)
        data.frame(sample = s, n = length(x), min = min(x), q1 = q[1],
                   median = median(x), q3 = q[2], max = max(x),
                   mean = round(mean(x), 1), sd = round(sd(x), 2))
    }
}))
write.csv(summary_df, file.path(out_dir, "insert_size_summary.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)

# boxplot, width scales with sample count so labels don't collide
n <- length(samples_ord)
png(file.path(out_dir, "insert_size_boxplot.png"),
    width = max(900, 60 * n), height = 600, res = 120)
par(mar = c(max(6, 0.5 * max(nchar(samples_ord))), 4, 3, 1))
boxplot(data_ord, names = samples_ord, las = 2, col = "steelblue",
        outline = TRUE, pch = 20, cex = 0.4,
        main = "Insert size by sample", ylab = "Insert size (bp)")
dev.off()

cat("Wrote insert_size_boxplot.png and insert_size_summary.csv to", out_dir, "\n")
