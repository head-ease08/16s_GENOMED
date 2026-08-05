log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

# No merge step for regions where R1+R2 can't physically overlap (amplicon
# longer than R1+R2 combined) -- ASV sequences here are forward-read-only.
# drop samples that had 0 reads survive filtering upstream (NULL placeholder)
dada_fwd <- lapply(snakemake@input$dada_fwd, readRDS)
dada_fwd <- dada_fwd[!vapply(dada_fwd, is.null, logical(1))]

saveRDS(
    makeSequenceTable(dada_fwd),
    snakemake@output$sequence_table
)
