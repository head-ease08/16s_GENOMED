library(dada2)

saveRDS(
    derepFastq(snakemake@input$r1, n = 1e7, verbose = TRUE),
    snakemake@output$r1_rds
)
saveRDS(
    derepFastq(snakemake@input$r2, n = 1e7, verbose = TRUE),
    snakemake@output$r2_rds
)
forward_derep <- lapply(snakemake@input$r1_rds, readRDS)
names(forward_derep) <- sub("_R1\\.rds$", "", basename(snakemake@input$r1_rds))

reverse_derep <- lapply(snakemake@input$r2_rds, readRDS)
names(reverse_derep) <- sub("_R2\\.rds$", "", basename(snakemake@input$r2_rds)
