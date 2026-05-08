library(dada2)
saveRDS(
    assignTaxonomy(
        readRDS(snakemake@input$taxa),
        readRDS(snakemake@input$silva),
        tryRC       = TRUE
    ),
    snakemake@output$taxa
)