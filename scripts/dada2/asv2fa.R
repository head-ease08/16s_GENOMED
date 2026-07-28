log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seqtab <- readRDS(snakemake@input$seqtab_nochim)

seqs <- colnames(seqtab)
size <- colSums(seqtab)
ids  <- paste0("ASV", seq_along(seqs), ";size=", size)

fa <- file(snakemake@output$fasta, open = "wt")
for (i in seq_along(seqs)) {
    writeLines(paste0(">", ids[i]), fa)
    writeLines(seqs[i], fa)
}
close(fa)
