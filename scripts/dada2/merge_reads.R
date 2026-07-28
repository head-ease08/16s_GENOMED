log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    mergePairs(
        readRDS(snakemake@input$r1_rds_dada),
        readRDS(snakemake@input$r1_rds_derep),
        readRDS(snakemake@input$r2_rds_dada),
        readRDS(snakemake@input$r2_rds_derep),
        minOverlap  = 12,
        maxMismatch = 0,
        verbose     = TRUE
    ),
    snakemake@output$merged_reads
)
