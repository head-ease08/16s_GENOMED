log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

# No merge step for regions where R1+R2 can't physically overlap (amplicon
# longer than R1+R2 combined) -- ASV sequences here are forward-read-only.
saveRDS(
    makeSequenceTable(lapply(snakemake@input$dada_fwd, readRDS)),
    snakemake@output$sequence_table
)
