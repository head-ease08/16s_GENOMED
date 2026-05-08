library(dada2)

ggsave(snakemake@output$r1_pdf,
       plotQualityProfile(snakemake@input$r1, n = 50000),
       width = 10, height = 6)

ggsave(snakemake@output$r2_pdf,
       plotQualityProfile(snakemake@input$r2, n = 50000),
       width = 10, height = 6)
