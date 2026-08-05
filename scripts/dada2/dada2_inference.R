log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

# derep step writes NULL when a sample/region had 0 reads survive filtering --
# propagate NULL through rather than calling dada() on it.
dada_or_null <- function(derep_rds, err_rds) {
    derep <- readRDS(derep_rds)
    if (is.null(derep)) return(NULL)
    dada(derep, err = readRDS(err_rds), pool = FALSE, multithread = TRUE)
}

saveRDS(dada_or_null(snakemake@input$r1_rds, snakemake@input$r1_err_rds), snakemake@output$r1_rds)
saveRDS(dada_or_null(snakemake@input$r2_rds, snakemake@input$r2_err_rds), snakemake@output$r2_rds)
