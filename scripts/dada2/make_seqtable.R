log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

files <- snakemake@input$merged_reads
names(files) <- sub("\\.rds$", "", basename(files))

# drop samples that had 0 reads survive filtering upstream (NULL placeholder)
merged <- lapply(files, readRDS)
merged <- merged[!vapply(merged, is.null, logical(1))]

saveRDS(
    makeSequenceTable(merged),
    snakemake@output$sequence_table
)
