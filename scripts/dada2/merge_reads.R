log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

dada_r1  <- readRDS(snakemake@input$r1_rds_dada)
derep_r1 <- readRDS(snakemake@input$r1_rds_derep)
dada_r2  <- readRDS(snakemake@input$r2_rds_dada)
derep_r2 <- readRDS(snakemake@input$r2_rds_derep)

# upstream NULL means this sample/region had 0 reads survive filtering --
# propagate NULL through rather than calling mergePairs() on it.
merged <- if (is.null(dada_r1) || is.null(dada_r2)) {
    NULL
} else {
    mergePairs(
        dada_r1, derep_r1, dada_r2, derep_r2,
        minOverlap  = 12,
        maxMismatch = 0,
        verbose     = TRUE
    )
}

saveRDS(merged, snakemake@output$merged_reads)
