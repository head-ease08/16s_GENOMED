library(dada2)

saveRDS(
    dada(readRDS(snakemake@input$r1_rds), err = readRDS(snakemake@input$r1_err_rds),
         pool = FALSE, multithread = TRUE),
    snakemake@output$r1_rds
)
saveRDS(
    dada(readRDS(snakemake@input$r2_rds), err = readRDS(snakemake@input$r2_err_rds),
         pool = FALSE, multithread = TRUE),
    snakemake@output$r2_rds
)
