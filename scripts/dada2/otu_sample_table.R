log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seqtab <- readRDS(snakemake@input$seqtab_nochim)
seqtab <- seqtab[, colSums(seqtab) > 0, drop = FALSE]

# 0 ASVs survived upstream -- vsearch writes a 0-byte otu_table.tsv in that
# case, which read.table can't parse at all ("no lines available in input").
if (file.size(snakemake@input$otu_table) == 0) {
    write.table(seqtab, snakemake@output$otu_sample_table,
                sep = "\t", quote = FALSE, col.names = NA)
    quit(save = "no", status = 0)
}

asv_map <- read.table(snakemake@input$asv_map, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
seq_of_asv <- setNames(asv_map$sequence, asv_map$asv_id)

# vsearch --cluster_size --otutabout: rows = OTU (centroid asv_id), columns
# = input asv_id, each column has exactly one nonzero row -- that's the
# ASV's cluster assignment. This pooled table only encodes membership here,
# never used for actual counts (those came from summing all samples
# together in asv2fa.R, which is not what LULU expects).
membership <- read.table(snakemake@input$otu_table, header = TRUE, row.names = 1,
                          sep = "\t", check.names = FALSE, comment.char = "")
otu_of_asv <- apply(membership, 2, function(col) rownames(membership)[which(col > 0)[1]])

real_counts <- t(seqtab)
rownames(real_counts) <- names(seq_of_asv)[match(rownames(real_counts), seq_of_asv)]
real_counts <- real_counts[names(otu_of_asv), , drop = FALSE]

otu_table_real <- rowsum(real_counts, group = otu_of_asv[rownames(real_counts)])

write.table(otu_table_real, snakemake@output$otu_sample_table,
            sep = "\t", quote = FALSE, col.names = NA)
