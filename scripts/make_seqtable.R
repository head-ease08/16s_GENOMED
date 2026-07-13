log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    makeSequenceTable(lapply(snakemake@input$merged_reads, readRDS)),
    snakemake@output$sequence_table
)
