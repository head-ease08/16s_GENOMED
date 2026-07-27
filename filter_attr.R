#!/usr/bin/env Rscript
# Attribute filterAndTrim losses to individual filters.
#
# Method: run baseline config, then run with each filter turned OFF individually.
# The delta in kept reads = number of reads that ONLY that filter was removing.
#
# Usage:
#   Rscript filter_attribution.R <sample_R1.fq.gz> <sample_R2.fq.gz>
#
# Or point it at data/trimmed/{sample}_R{1,2}.fq.gz directly.

suppressMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: filter_attribution.R R1.fq.gz R2.fq.gz")
r1 <- args[1]; r2 <- args[2]

tmp_r1 <- tempfile(fileext = ".fq.gz")
tmp_r2 <- tempfile(fileext = ".fq.gz")

run_cfg <- function(name, ...) {
    args <- list(...)
    # baseline
    baseline <- list(truncLen = c(120, 120), maxN = 0,
                     maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE)
    cfg <- modifyList(baseline, args)
    out <- do.call(filterAndTrim, c(
        list(r1, tmp_r1, r2, tmp_r2,
             compress = TRUE, multithread = TRUE, verbose = FALSE),
        cfg))
    cat(sprintf("%-25s in=%7d out=%7d keep=%5.1f%%\n",
                name, out[1,1], out[1,2], 100 * out[1,2] / out[1,1]))
    invisible(out[1,2])
}

cat("=== leave-one-out on filterAndTrim ===\n")
cat(sprintf("R1: %s\nR2: %s\n\n", basename(r1), basename(r2)))

base <- run_cfg("baseline")
cat("\n--- turn off each filter individually ---\n")
no_ee   <- run_cfg("no maxEE",    maxEE   = c(Inf, Inf))
no_eeR2 <- run_cfg("no maxEE on R2 only", maxEE = c(2, Inf))
no_trQ  <- run_cfg("no truncQ",   truncQ  = 0)
no_len  <- run_cfg("no truncLen (0 keeps all len)", truncLen = c(0, 0))
no_n    <- run_cfg("no maxN",     maxN    = Inf)
no_phix <- run_cfg("no rm.phix",  rm.phix = FALSE)

cat("\n--- interpretation ---\n")
cat("For each row: (kept when OFF) - (baseline kept) = reads THIS filter alone was killing.\n")
cat("Reads killed by:\n")
cat(sprintf("  maxEE (both):   %+7d\n", no_ee - base))
cat(sprintf("  maxEE (R2):     %+7d  <- expected biggest for 2x151 Illumina\n", no_eeR2 - base))
cat(sprintf("  truncQ:         %+7d\n", no_trQ - base))
cat(sprintf("  truncLen<120:   %+7d\n", no_len - base))
cat(sprintf("  maxN:           %+7d\n", no_n - base))
cat(sprintf("  rm.phix:        %+7d\n", no_phix - base))

cat("\n--- some candidate configs ---\n")
run_cfg("softR2",           maxEE = c(2, 5))
run_cfg("softR2 + short",   maxEE = c(2, 5), truncLen = c(115, 110))
run_cfg("R1 only",          maxEE = c(2, Inf), truncLen = c(115, 0))
