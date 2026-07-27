#!/usr/bin/env Rscript
# Render insert-size (TLEN) histograms from *.insert_sizes.txt files.
# Usage: plot_insert_size_hist.R <dir_with_insert_sizes_txt>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_insert_size_hist.R <dir>")
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
samples <- samples_raw[ord_info$order]      # sample names, in display order
labels  <- ord_info$labels                  # e.g. "1","2",...,"28" (or full name if pattern didn't match)
label_of <- setNames(labels, samples)

for (s in samples) {
    cat(sprintf("%-40s n=%d %s\n", s, length(data[[s]]),
                if (length(data[[s]])) paste("range", min(data[[s]]), "-", max(data[[s]])) else ""))
}

# Per-sample clip at that SAME sample's own 99th percentile -- not a pooled
# one. A pooled/global xmax silently emptied out any sample whose scale
# differed from the rest (e.g. one sample mostly short junk fragments while
# others are normal-length): x[x <= global_xmax] could drop ~all of that
# sample's own points, rendering a blank-looking plot despite the txt file
# having data.
plot_hist <- function(x, main, xlim = NULL) {
    if (length(x) == 0) {
        plot.new(); title(paste(main, "- no data")); return(invisible())
    }
    own_xmax <- quantile(x, 0.99, names = FALSE)
    xv <- x[x <= own_xmax]
    if (length(xv) == 0) xv <- x  # degenerate: fall back to unclipped
    if (length(unique(xv)) < 2) {
        # single repeated value -- hist() can choke on this, draw a bar manually.
        # Keep the axis tight around the value (not a fixed 0..N range) so the
        # bar is wide enough on-screen to actually be visible, not a 1px hairline.
        center <- xv[1]
        span <- max(5, center * 0.05)
        plot_xlim <- if (!is.null(xlim)) xlim else c(center - span, center + span)
        bar_w <- diff(plot_xlim) * 0.02
        plot(1, type = "n", xlim = plot_xlim,
             ylim = c(0, length(xv)), xlab = "Insert size (bp)", ylab = "Frequency", main = main)
        rect(center - bar_w, 0, center + bar_w, length(xv), col = "steelblue", border = NA)
    } else if (is.null(xlim)) {
        hist(xv, breaks = min(100, length(unique(xv))), col = "steelblue", border = NA,
             main = main, xlab = "Insert size (bp)")
    } else {
        hist(xv, breaks = min(100, length(unique(xv))), col = "steelblue", border = NA,
             main = main, xlab = "Insert size (bp)", xlim = xlim)
    }
    abline(v = median(x), col = "red", lty = 2)
}

# per-sample PNGs -- filename keeps the real sample name (traceability),
# on-image title shows just the short numeric label
for (s in samples) {
    x <- data[[s]]
    png(file.path(out_dir, paste0(s, ".insert_size_hist.png")), width = 900, height = 600, res = 120)
    plot_hist(x, paste("Insert size - Sample", label_of[[s]]))
    if (length(x) > 0) {
        legend("topright", legend = sprintf("median = %d, n = %d", as.integer(median(x)), length(x)),
               col = "red", lty = 2, bty = "n")
    }
    dev.off()
}

# combined grid, one panel per sample, in sample-number order (1,2,3...
# not alphabetical) -- each panel auto-scales to its own data
n <- length(samples)
ncol <- ceiling(sqrt(n))
nrow <- ceiling(n / ncol)
png(file.path(out_dir, "insert_size_hist_all.png"),
    width = 400 * ncol, height = 300 * nrow, res = 120)
par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))
for (s in samples) {
    plot_hist(data[[s]], label_of[[s]])
}
dev.off()

# second grid: same panels, but on a SHARED x-axis so peak position is
# directly comparable across samples (like the boxplot, but with the actual
# distribution shape visible, not just the IQR box). Range is computed
# dynamically from the real pooled data (0.5th-99.5th pctile, not a fixed
# guess) so it's wide enough to not blank out a sample that sits at a
# different scale than the rest -- it'll just look narrow in that panel,
# same tradeoff a shared axis always has, but never empty.
pooled <- unlist(data, use.names = FALSE)
if (length(pooled) > 0) {
    shared_xlim <- quantile(pooled, c(0.005, 0.995), names = FALSE)
    if (diff(shared_xlim) == 0) shared_xlim <- shared_xlim + c(-5, 5)
    png(file.path(out_dir, "insert_size_hist_all_shared_axis.png"),
        width = 400 * ncol, height = 300 * nrow, res = 120)
    par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))
    for (s in samples) {
        plot_hist(data[[s]], label_of[[s]], xlim = shared_xlim)
    }
    dev.off()
}

cat("Wrote", n, "per-sample PNGs, insert_size_hist_all.png (own-scale) and",
    "insert_size_hist_all_shared_axis.png (shared scale) to", out_dir, "\n")
