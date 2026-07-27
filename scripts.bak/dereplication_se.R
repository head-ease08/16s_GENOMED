log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    derepFastq(snakemake@input$r1, n = 1e7, verbose = TRUE),
    snakemake@output$r1_rds
)
