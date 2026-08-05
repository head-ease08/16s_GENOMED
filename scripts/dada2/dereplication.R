log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)
library(ShortRead)

# derepFastq errors on a fastq with 0 reads (empty gzip, upstream filter
# dropped everything for this sample/region) -- write NULL through instead
# so downstream rules can detect and skip the sample.
derep_or_null <- function(f) {
    if (countFastq(f)$records == 0) return(NULL)
    derepFastq(f, n = 1e7, verbose = TRUE)
}

saveRDS(derep_or_null(snakemake@input$r1), snakemake@output$r1_rds)
saveRDS(derep_or_null(snakemake@input$r2), snakemake@output$r2_rds)
