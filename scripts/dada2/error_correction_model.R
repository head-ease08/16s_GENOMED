log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)
library(ggplot2)

forward_err <- learnErrors(snakemake@input$r1, multithread = TRUE, MAX_CONSIST = 20)
reverse_err <- learnErrors(snakemake@input$r2, multithread = TRUE, MAX_CONSIST = 20)

saveRDS(forward_err, snakemake@output$r1_rds)
saveRDS(reverse_err, snakemake@output$r2_rds)

ggsave(snakemake@output$r1_pdf,
       plotErrors(forward_err, nominalQ = TRUE),
       width = 10, height = 6)

ggsave(snakemake@output$r2_pdf,
       plotErrors(reverse_err, nominalQ = TRUE),
       width = 10, height = 6)
