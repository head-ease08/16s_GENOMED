log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    addSpecies(
        readRDS(snakemake@input$taxa),
        snakemake@input$silva,
        tryRC = TRUE
    ),
    snakemake@output$taxa
)
