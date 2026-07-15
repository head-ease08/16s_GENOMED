log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seqtab <- readRDS(snakemake@input$seqtab_nochim)
taxa   <- readRDS(snakemake@input$taxa_species)

# taxa rows and seqtab columns are both keyed by ASV sequence
taxa <- taxa[colnames(seqtab), , drop = FALSE]

is_bacteria <- taxa[, "Kingdom"] %in% c("Bacteria", "Archaea") &
    !(taxa[, "Family"] %in% c("Mitochondria")) &
    !(taxa[, "Order"] %in% c("Chloroplast"))
is_bacteria[is.na(is_bacteria)] <- FALSE

seqtab_kingdom <- seqtab[, is_bacteria, drop = FALSE]
cat("Kingdom/organelle filter: kept", ncol(seqtab_kingdom), "of", ncol(seqtab), "ASVs\n")

# drop ASVs seen in fewer than 2 samples -- random noise/index-hopping rarely
# repeats across samples, real biology usually does
present_in_samples <- colSums(seqtab_kingdom > 0)
seqtab_clean <- seqtab_kingdom[, present_in_samples >= 2, drop = FALSE]
cat("Prevalence>=2 filter: kept", ncol(seqtab_clean), "of", ncol(seqtab_kingdom), "ASVs\n")

saveRDS(seqtab_clean, snakemake@output$seqtab_clean)
