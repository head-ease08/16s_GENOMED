log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

# No merge step for regions where R1+R2 can't physically overlap (amplicon
# longer than R1+R2 combined) -- ASV sequences here are forward-read-only.
files <- snakemake@input$dada_fwd
names(files) <- sub("_R1\\.rds$", "", basename(files))

# drop samples that had 0 reads survive filtering upstream (NULL placeholder)
dada_fwd <- lapply(files, readRDS)
dada_fwd <- dada_fwd[!vapply(dada_fwd, is.null, logical(1))]

saveRDS(
    makeSequenceTable(dada_fwd),
    snakemake@output$sequence_table
)
