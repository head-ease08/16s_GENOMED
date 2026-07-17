#!/usr/bin/env Rscript
# One plot per region: all samples' insert-size distributions overlaid as
# density curves on a shared axis, colored by sample number (gradient, not
# a 28-entry legend). Lets you see at a glance whether peaks line up across
# samples or a subset drifted -- more informative than a grid of near-
# identical panels, and shows actual shape unlike the boxplot.
#
# Usage: plot_insert_size_density.R <dir_with_insert_sizes_txt>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_insert_size_density.R <dir>")
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
n <- length(samples)

data <- data[samples]
data <- data[sapply(data, length) > 0]  # skip empty samples, nothing to draw
if (length(data) == 0) stop("All samples empty -- nothing to plot")
labels <- labels[samples %in% names(data)]

# one density (or, for near-constant data, a narrow synthetic bump so it's
# still visible) per sample
curves <- lapply(data, function(x) {
    if (length(unique(x)) < 2) {
        # density() errors on zero-variance input -- fake a narrow spike
        d <- density(c(x, x + 1e-6, x - 1e-6), bw = 0.3)
    } else {
        d <- density(x)
    }
    d
})

pooled <- unlist(data, use.names = FALSE)
xlim <- quantile(pooled, c(0.005, 0.995), names = FALSE)
if (diff(xlim) == 0) xlim <- xlim + c(-5, 5)
ylim <- c(0, max(sapply(curves, function(d) max(d$y))))

pal <- hcl.colors(length(curves), "Viridis")

png(file.path(out_dir, "insert_size_density_overlay.png"), width = 1100, height = 700, res = 120)
par(mar = c(4, 4, 3, 6), xpd = FALSE)
plot(NA, xlim = xlim, ylim = ylim, xlab = "Insert size (bp)", ylab = "Density",
     main = paste("Insert size distribution -", basename(out_dir), "-", length(curves), "samples"))
for (i in seq_along(curves)) {
    lines(curves[[i]], col = pal[i], lwd = 1.5)
}

# color legend: gradient bar labelled with sample numbers, not a 28-row key
usr <- par("usr")
legend_x <- usr[2] + 0.02 * diff(usr[1:2])
legend_y <- seq(usr[3], usr[4], length.out = length(curves) + 1)
par(xpd = TRUE)
for (i in seq_along(curves)) {
    rect(legend_x, legend_y[i], legend_x + 0.03 * diff(usr[1:2]), legend_y[i + 1],
         col = pal[i], border = NA)
}
label_idx <- unique(round(seq(1, length(curves), length.out = min(10, length(curves)))))
text(legend_x + 0.05 * diff(usr[1:2]), (legend_y[label_idx] + legend_y[label_idx + 1]) / 2,
     labels = labels[label_idx], cex = 0.7, adj = 0)
dev.off()

cat("Wrote insert_size_density_overlay.png (", length(curves), "samples) to", out_dir, "\n")
