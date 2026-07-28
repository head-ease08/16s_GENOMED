log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)
library(ggplot2)

ggsave(snakemake@output$r1_pdf,
       plotQualityProfile(snakemake@input$r1, aggregate = TRUE, n = 100000),
       width = 10, height = 6)

ggsave(snakemake@output$r2_pdf,
       plotQualityProfile(snakemake@input$r2, aggregate = TRUE, n = 100000),
       width = 10, height = 6)
