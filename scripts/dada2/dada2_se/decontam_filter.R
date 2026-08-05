log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(decontam)

seqtab <- readRDS(snakemake@input$seqtab_nochim)
neg_samples <- snakemake@params$neg_samples

is_neg <- rownames(seqtab) %in% neg_samples
cat("Negative controls found in seqtab:", sum(is_neg), "of", length(neg_samples), "expected\n")

contam <- isContaminant(seqtab, method = "prevalence", neg = is_neg)
cat("Contaminant ASVs:", sum(contam$contaminant, na.rm = TRUE), "of", nrow(contam), "\n")

# Drop contaminant ASVs (columns) and the negative-control samples themselves
# (rows) -- blanks did their job identifying contaminants, they aren't real
# biological samples and don't belong in the final abundance table.
seqtab_clean <- seqtab[!is_neg, !contam$contaminant, drop = FALSE]

saveRDS(seqtab_clean, snakemake@output$seqtab_decontam)
saveRDS(contam, snakemake@output$contam_stats)
