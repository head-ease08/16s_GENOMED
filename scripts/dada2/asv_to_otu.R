log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Per-ASV composition of each OTU's vsearch cluster: which ASVs fell into
# which OTU, and each ASV's own (pre-clustering) taxonomy call -- a
# diagnostic/audit view alongside the consensus taxonomy in otu_taxonomy.R.
strip_size <- function(x) sub(";size=.*$", "", x)

empty_out <- data.frame(OTU_ID = character(0), ASV_ID = character(0),
                         Taxonomy = character(0), Sequence = character(0),
                         stringsAsFactors = FALSE)

# 0 ASVs survived upstream -- vsearch writes a 0-byte cluster_map.uc in that
# case, which read.table can't parse at all ("no lines available in input").
if (file.size(snakemake@input$uc) == 0) {
    write.table(empty_out, snakemake@output$composition, sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 0)
}

asv_map <- read.table(snakemake@input$asv_map, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
seq_of_asv <- setNames(toupper(asv_map$sequence), asv_map$asv_id)

taxa <- readRDS(snakemake@input$taxa_species)
taxa <- taxa[, !duplicated(colnames(taxa), fromLast = TRUE), drop = FALSE]
tax_string <- apply(taxa, 1, function(row) paste(row[!is.na(row)], collapse = ";"))

uc <- read.table(snakemake@input$uc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
uc <- uc[uc$V1 %in% c("S", "H"), ]

asv_id <- strip_size(uc$V9)
otu_id <- ifelse(uc$V1 == "S", asv_id, strip_size(uc$V10))
seqs   <- seq_of_asv[asv_id]

out <- data.frame(
    OTU_ID   = otu_id,
    ASV_ID   = asv_id,
    Taxonomy = unname(tax_string[seqs]),
    Sequence = unname(seqs),
    stringsAsFactors = FALSE
)
out <- out[order(out$OTU_ID, out$ASV_ID), ]

write.table(out, snakemake@output$composition, sep = "\t", row.names = FALSE, quote = FALSE)
