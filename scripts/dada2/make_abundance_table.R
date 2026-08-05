log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

seqtab  <- readRDS(snakemake@input$seqtab_nochim)
taxa    <- readRDS(snakemake@input$taxa_species)

counts <- t(seqtab)

shared <- intersect(rownames(counts), rownames(taxa))
counts <- counts[shared, , drop = FALSE]
taxa   <- taxa[shared, , drop = FALSE]

taxa <- taxa[, !duplicated(colnames(taxa), fromLast = TRUE), drop = FALSE]

relabund <- sweep(counts, 2, colSums(counts), "/") * 100
colnames(relabund) <- paste0(colnames(counts), "_pct")

# paste0("ASV", integer(0)) returns "ASV" (length 1), not character(0) --
# guard the 0-row case explicitly.
asv_id <- if (nrow(counts)) paste0("ASV", seq_len(nrow(counts))) else character(0)

abundance <- cbind(
    asv_id   = asv_id,
    sequence = shared,
    as.data.frame(taxa),
    as.data.frame(counts),
    as.data.frame(relabund)
)

write.csv(abundance, snakemake@output$abundance_table, row.names = FALSE)
