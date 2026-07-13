log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    removeBimeraDenovo(
        readRDS(snakemake@input$sequence_table),
        method      = "consensus",
        multithread = TRUE,
        verbose     = TRUE
    ),
    snakemake@output$seq_tab_nochim
)
