log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

read_fasta <- function(path) {
    lines <- readLines(path)
    id_idx <- grep("^>", lines)
    ids    <- sub("^>", "", lines[id_idx])
    ids    <- sub(";size=.*$", "", ids)
    ends   <- c(id_idx[-1] - 1, length(lines))
    seqs   <- vapply(seq_along(id_idx), function(i) {
        paste(lines[(id_idx[i] + 1):ends[i]], collapse = "")
    }, character(1))
    setNames(toupper(seqs), ids)
}

curated <- read.table(snakemake@input$curated_table, header = TRUE, row.names = 1,
                       sep = "\t", check.names = FALSE, comment.char = "")
centroid_seqs <- read_fasta(snakemake@input$centroids)
taxa          <- readRDS(snakemake@input$taxa_species)

otu_seqs <- centroid_seqs[rownames(curated)]

taxa <- taxa[, !duplicated(colnames(taxa), fromLast = TRUE), drop = FALSE]
taxa_matched <- taxa[otu_seqs, , drop = FALSE]
rownames(taxa_matched) <- rownames(curated)

out <- cbind(
    otu_id   = rownames(curated),
    sequence = otu_seqs,
    as.data.frame(taxa_matched),
    curated
)

write.table(out, snakemake@output$otu_taxonomy, sep = "\t", quote = FALSE, row.names = FALSE)
