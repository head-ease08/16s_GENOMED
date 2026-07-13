log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    assignTaxonomy(
        readRDS(snakemake@input$seqtab_nochim),
        snakemake@input$silva,
        multithread = TRUE,
        tryRC       = TRUE
    ),
    snakemake@output$taxa
)
