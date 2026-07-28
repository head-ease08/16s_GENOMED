log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seqtab <- readRDS(snakemake@input$seqtab_nochim)
seqtab <- seqtab[, colSums(seqtab) > 0, drop = FALSE]

seqs    <- colnames(seqtab)
size    <- colSums(seqtab)
asv_id  <- paste0("ASV", seq_along(seqs))
headers <- paste0(asv_id, ";size=", size)

fa <- file(snakemake@output$fasta, open = "wt")
for (i in seq_along(seqs)) {
    writeLines(paste0(">", headers[i]), fa)
    writeLines(seqs[i], fa)
}
close(fa)

# asv_id <-> sequence map, so downstream steps can trace a pooled-abundance
# vsearch cluster label (ASV<N>) back to the real per-sample seqtab column.
write.table(data.frame(asv_id = asv_id, sequence = seqs),
            snakemake@output$asv_map, sep = "\t", quote = FALSE, row.names = FALSE)
