library(dada2)

ggsave(snakemake@output$r1_pdf,
       plotQualityProfile(snakemake@input$r1, aggregate = TRUE, n = 100000),
       width = 10, height = 6)

ggsave(snakemake@output$r2_pdf,
       plotQualityProfile(snakemake@input$r2, aggregate = TRUE, n = 100000),
       width = 10, height = 6)
