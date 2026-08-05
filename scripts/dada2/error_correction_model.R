log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)
library(ggplot2)
library(ShortRead)

# some samples/regions had 0 reads survive upstream filtering, leaving an
# empty placeholder gzip -- learnErrors chokes on those, so drop them.
nonempty <- function(files) Filter(function(f) countFastq(f)$records > 0, files)

forward_err <- learnErrors(nonempty(snakemake@input$r1), multithread = TRUE, MAX_CONSIST = 20)
reverse_err <- learnErrors(nonempty(snakemake@input$r2), multithread = TRUE, MAX_CONSIST = 20)

saveRDS(forward_err, snakemake@output$r1_rds)
saveRDS(reverse_err, snakemake@output$r2_rds)

ggsave(snakemake@output$r1_pdf,
       plotErrors(forward_err, nominalQ = TRUE),
       width = 10, height = 6)

ggsave(snakemake@output$r2_pdf,
       plotErrors(reverse_err, nominalQ = TRUE),
       width = 10, height = 6)
