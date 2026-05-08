library(dada2)

saveRDS(
    makeSequenceTable(lapply(snakemake@input$merged_reads, readRDS)),
    snakemake@output$sequence_table
)