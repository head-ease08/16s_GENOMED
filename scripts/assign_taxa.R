library(dada2)

saveRDS(
    assignTaxonomy(
        readRDS(snakemake@input$seqtab_nochim),
        readRDS(snakemake@input$silva),
        multithread = TRUE,
        tryRC       = TRUE
    ),
    snakemake@output$taxa
)