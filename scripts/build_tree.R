log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)
library(DECIPHER)
library(ape)
library(ggtree)

seqs <- getSequences(readRDS(snakemake@input$seqtab_nochim))
names(seqs) <- seqs

alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)

tree_dend <- Treeline(alignment,
                      method  = "NJ",
                      model   = "GTR",
                      verbose = FALSE)

tree <- as.phylo(as.hclust(tree_dend))

saveRDS(tree, snakemake@output$tree_rds)
write.tree(tree, snakemake@output$tree_nwk)

taxa <- as.data.frame(readRDS(snakemake@input$taxa_species))
taxa <- taxa[, !duplicated(colnames(taxa))]
taxa <- cbind(label = rownames(taxa), taxa)

taxa$tip_label <- ifelse(
    !is.na(taxa$Species),
    paste(taxa$Genus, taxa$Species),
    ifelse(!is.na(taxa$Genus), taxa$Genus,
    ifelse(!is.na(taxa$Family), taxa$Family, taxa$Order))
)

p <- ggtree(tree) %<+% taxa +
    geom_tiplab(aes(label = tip_label, color = Phylum), size = 2.5, align = TRUE) +
    theme_tree2() +
    theme(legend.position = "bottom")

ggsave(snakemake@output$tree_pdf, p, width = 12, height = 14)
