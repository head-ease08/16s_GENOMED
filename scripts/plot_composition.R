log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)
library(ggplot2)
library(tidyr)

seqtab   <- readRDS(snakemake@input$seqtab_nochim)
taxa     <- readRDS(snakemake@input$taxa_species)
samples  <- snakemake@params$samples
top_n    <- snakemake@params$top_n
meta_csv <- snakemake@params$metadata

rownames(seqtab) <- samples

meta      <- read.csv(meta_csv, check.names = FALSE, stringsAsFactors = FALSE)
srr_col   <- "Run"
name_col  <- "Sample Name"

extract_vregion <- function(x) {
    m <- regmatches(x, regexpr("V\\d+(-V\\d+)*$", x))
    ifelse(length(m) == 0, x, m)
}

meta$vregion <- sapply(meta[[name_col]], extract_vregion)

srr_to_vregion <- setNames(meta$vregion, meta[[srr_col]])

common <- intersect(samples, names(srr_to_vregion))
seqtab <- seqtab[common, , drop = FALSE]
group  <- srr_to_vregion[common]

taxa_df <- as.data.frame(taxa)
taxa_df <- taxa_df[, !duplicated(colnames(taxa_df), fromLast = TRUE), drop = FALSE]

genus <- ifelse(
    !is.na(taxa_df$Genus),    taxa_df$Genus,
    ifelse(!is.na(taxa_df$Family), paste0("unknown_", taxa_df$Family),
    ifelse(!is.na(taxa_df$Order),  paste0("unknown_", taxa_df$Order),
                                   "unassigned"))
)
names(genus) <- rownames(taxa_df)

counts <- t(seqtab)
shared <- intersect(rownames(counts), names(genus))
counts <- counts[shared, , drop = FALSE]
genus  <- genus[shared]

vregions     <- sort(unique(group))
agg_counts   <- sapply(vregions, function(vr) {
    cols <- names(group)[group == vr]
    cols <- intersect(cols, colnames(counts))
    if (length(cols) == 1) counts[, cols]
    else rowSums(counts[, cols, drop = FALSE])
})
colnames(agg_counts) <- vregions

genus_counts <- rowsum(agg_counts, group = genus)
relabund     <- sweep(genus_counts, 2, colSums(genus_counts), "/") * 100

sorted   <- sort(rowSums(relabund), decreasing = TRUE)
top_taxa <- names(sorted)[seq_len(min(top_n, length(sorted)))]
top_taxa <- top_taxa[top_taxa %in% rownames(relabund)]

rest  <- rownames(relabund)[!rownames(relabund) %in% top_taxa]
other <- if (length(rest) > 0) colSums(relabund[rest, , drop = FALSE]) else rep(0, ncol(relabund))

plot_mat        <- rbind(relabund[top_taxa, , drop = FALSE], other = other)
df              <- as.data.frame(plot_mat)
df$taxon        <- factor(rownames(df), levels = rev(c(top_taxa, "other")))
df_long         <- pivot_longer(df, cols = -taxon, names_to = "vregion", values_to = "pct")
df_long$vregion <- factor(df_long$vregion, levels = vregions)

p <- ggplot(df_long, aes(x = vregion, y = pct, fill = taxon)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(x = NULL, y = "Relative abundance (%)", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(
        axis.text.x        = element_text(angle = 45, hjust = 1),
        legend.position    = "right",
        panel.grid.major.x = element_blank()
    )

ggsave(snakemake@output$plot_pdf, p,
       width  = snakemake@params$width,
       height = snakemake@params$height)
