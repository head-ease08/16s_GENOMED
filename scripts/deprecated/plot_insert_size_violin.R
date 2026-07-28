#!/usr/bin/env Rscript
# Horizontal violin plot, one row per sample, on a shared dynamic x-axis.
# Unlike the boxplot, this shows the actual distribution SHAPE per sample --
# a sample with a bimodal mix (e.g. some reads at a short junk length, some
# at the real amplicon length) shows up as two humps instead of one huge
# IQR box that hides the bimodality. Base R only (density() + polygon()),
# no extra package dependency.
#
# Usage: plot_insert_size_violin.R <dir_with_insert_sizes_txt>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_insert_size_violin.R <dir>")
out_dir <- args[1]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), "insert_size_sample_order.R"))

files <- list.files(out_dir, pattern = "\\.insert_sizes\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No *.insert_sizes.txt files found in ", out_dir)

samples_raw <- sub("\\.insert_sizes\\.txt$", "", basename(files))
data <- setNames(lapply(files, function(f) {
    x <- scan(f, what = numeric(), quiet = TRUE)
    x[is.finite(x) & x > 0]
}), samples_raw)

ord_info <- sample_order_and_labels(samples_raw)
samples <- samples_raw[ord_info$order]
labels  <- ord_info$labels
data <- data[samples]

# top = sample 1 (reverse so row 1 in plot space, at the top, is sample 1)
data   <- rev(data)
labels <- rev(labels)
n <- length(samples)

medians <- sapply(data, function(x) if (length(x)) median(x) else NA)

pooled <- unlist(data, use.names = FALSE)
xlim <- if (length(pooled)) quantile(pooled, c(0.005, 0.995), names = FALSE) else c(0, 1)
if (diff(xlim) == 0) xlim <- xlim + c(-5, 5)

curves <- lapply(data, function(x) {
    if (length(x) == 0) return(NULL)
    if (length(unique(x)) < 2) {
        density(c(x, x + 1e-6, x - 1e-6), bw = 0.3, n = 512)
    } else {
        density(x, n = 512)
    }
})

png(file.path(out_dir, "insert_size_violin.png"),
    width = 1400, height = max(900, 45 * n), res = 130)
par(mar = c(4, 5, 7, 2))
plot(NA, xlim = xlim, ylim = c(0.3, n + 0.7), yaxt = "n",
     xlab = "Insert size (bp)", ylab = "Sample",
     main = "Insert size by sample (violin)")
axis(2, at = seq_len(n), labels = labels, las = 1, cex.axis = 0.8)

half_width <- 0.42
for (i in seq_len(n)) {
    d <- curves[[i]]
    if (is.null(d)) {
        text(mean(xlim), i, "no data", cex = 0.7, col = "grey40")
        next
    }
    scaled <- d$y / max(d$y) * half_width
    polygon(c(d$x, rev(d$x)), i + c(scaled, -rev(scaled)),
            col = "steelblue", border = "grey30", lwd = 0.5)
    if (!is.na(medians[i])) {
        segments(medians[i], i - half_width, medians[i], i + half_width,
                 col = "red", lwd = 1.5)
    }
}

# median value ticked on its own axis (top), same as the boxplot version
axis(3, at = medians, labels = round(medians, 1), las = 2, cex.axis = 0.6,
     col = "red", col.axis = "red", tick = TRUE)
dev.off()

cat("Wrote insert_size_violin.png to", out_dir, "\n")
