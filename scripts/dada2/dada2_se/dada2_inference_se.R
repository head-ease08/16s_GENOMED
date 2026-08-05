log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

derep <- readRDS(snakemake@input$r1_rds)

# derep step writes NULL when a sample/region had 0 reads survive filtering --
# propagate NULL through rather than calling dada() on it.
result <- if (is.null(derep)) {
    NULL
} else {
    dada(derep, err = readRDS(snakemake@input$r1_err_rds), pool = FALSE, multithread = TRUE)
}

saveRDS(result, snakemake@output$r1_rds)
