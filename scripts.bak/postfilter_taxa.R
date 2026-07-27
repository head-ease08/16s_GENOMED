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

# drop ASVs with fewer than 10 reads total across the whole run -- catches
# singleton/near-singleton noise the earlier filters don't touch
MIN_ABUNDANCE <- 10
total_reads <- colSums(seqtab_kingdom)
seqtab_abund <- seqtab_kingdom[, total_reads >= MIN_ABUNDANCE, drop = FALSE]
cat("Abundance>=", MIN_ABUNDANCE, " filter: kept ", ncol(seqtab_abund), " of ",
    ncol(seqtab_kingdom), " ASVs\n", sep = "")

# drop ASVs seen in fewer than 2 samples -- random noise/index-hopping rarely
# repeats across samples, real biology usually does
present_in_samples <- colSums(seqtab_abund > 0)
seqtab_clean <- seqtab_abund[, present_in_samples >= 2, drop = FALSE]
cat("Prevalence>=2 filter: kept", ncol(seqtab_clean), "of", ncol(seqtab_abund), "ASVs\n")

saveRDS(seqtab_clean, snakemake@output$seqtab_clean)
