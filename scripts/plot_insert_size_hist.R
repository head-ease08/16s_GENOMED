#!/usr/bin/env Rscript
# Render insert-size (TLEN) histograms from *.insert_sizes.txt files.
# Usage: plot_insert_size_hist.R <dir_with_insert_sizes_txt>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_insert_size_hist.R <dir>")
out_dir <- args[1]

files <- list.files(out_dir, pattern = "\\.insert_sizes\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No *.insert_sizes.txt files found in ", out_dir)

samples <- sub("\\.insert_sizes\\.txt$", "", basename(files))
data <- setNames(lapply(files, function(f) {
    x <- scan(f, what = numeric(), quiet = TRUE)
    x[is.finite(x) & x > 0]
}), samples)

# clip at 99th percentile (pooled) so a few outliers don't blow out the x-axis
pooled <- unlist(data, use.names = FALSE)
xmax <- if (length(pooled)) quantile(pooled, 0.99, names = FALSE) else 1000

# per-sample PNGs
for (s in samples) {
    x <- data[[s]]
    png(file.path(out_dir, paste0(s, ".insert_size_hist.png")), width = 900, height = 600, res = 120)
    if (length(x) == 0) {
        plot.new(); title(paste(s, "- no data"))
    } else {
        hist(x[x <= xmax], breaks = 100, col = "steelblue", border = NA,
             main = paste("Insert size -", s), xlab = "Insert size (bp)")
        abline(v = median(x), col = "red", lty = 2)
        legend("topright", legend = sprintf("median = %d", as.integer(median(x))),
               col = "red", lty = 2, bty = "n")
    }
    dev.off()
}

# combined grid, one panel per sample, shared x-axis for easy comparison
n <- length(samples)
ncol <- ceiling(sqrt(n))
nrow <- ceiling(n / ncol)
png(file.path(out_dir, "insert_size_hist_all.png"),
    width = 400 * ncol, height = 300 * nrow, res = 120)
par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))
for (s in samples) {
    x <- data[[s]]
    if (length(x) == 0) {
        plot.new(); title(paste(s, "- no data")); next
    }
    hist(x[x <= xmax], breaks = 100, col = "steelblue", border = NA,
         main = s, xlab = "Insert size (bp)", xlim = c(0, xmax))
    abline(v = median(x), col = "red", lty = 2)
}
dev.off()

cat("Wrote", n, "per-sample PNGs and insert_size_hist_all.png to", out_dir, "\n")
