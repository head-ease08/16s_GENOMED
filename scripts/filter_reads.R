log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

saveRDS(
    filterAndTrim(
        snakemake@input$r1,  snakemake@output$r1,
        snakemake@input$r2,  snakemake@output$r2,
        truncLen    = c(140, 141),
        maxN        = 0,
        maxEE       = c(2, 2),
        truncQ      = 2,
        rm.phix     = TRUE,
        compress    = TRUE,
        multithread = TRUE
    ),
    snakemake@output$stats
)
