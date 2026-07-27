log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

seqtab  <- readRDS(snakemake@input$seqtab_nochim)
taxa    <- readRDS(snakemake@input$taxa_species)
samples <- snakemake@params$samples

rownames(seqtab) <- samples

counts <- t(seqtab)

shared <- intersect(rownames(counts), rownames(taxa))
counts <- counts[shared, , drop = FALSE]
taxa   <- taxa[shared, , drop = FALSE]

taxa <- taxa[, !duplicated(colnames(taxa), fromLast = TRUE), drop = FALSE]

relabund <- sweep(counts, 2, colSums(counts), "/") * 100
colnames(relabund) <- paste0(colnames(counts), "_pct")

abundance <- cbind(
    asv_id   = paste0("ASV", seq_len(nrow(counts))),
    sequence = shared,
    as.data.frame(taxa),
    as.data.frame(counts),
    as.data.frame(relabund)
)

write.csv(abundance, snakemake@output$abundance_table, row.names = FALSE)
