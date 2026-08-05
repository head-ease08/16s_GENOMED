log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seqtab <- readRDS(snakemake@input$seqtab_nochim)
seqtab <- seqtab[, colSums(seqtab) > 0, drop = FALSE]

# colnames() on a 0-column matrix with no dimnames returns NULL, not
# character(0) -- passing NULL into data.frame() below silently drops the
# whole column instead of producing an empty one, so coerce explicitly.
seqs <- as.character(colnames(seqtab))
size <- colSums(seqtab)

# paste0("ASV", integer(0)) returns "ASV" (length 1), not character(0) --
# guard the 0-ASV case explicitly so asv_id/headers stay empty too.
asv_id  <- if (length(seqs)) paste0("ASV", seq_along(seqs)) else character(0)
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
