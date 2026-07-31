log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Per-rank consensus taxonomy across all members of an OTU's vsearch cluster,
# weighted by each member's total (pooled) abundance -- not just the
# centroid ASV that vsearch greedily picked as cluster seed. At each rank
# (Kingdom -> Species), the label carrying >= MAJORITY_THRESHOLD of the
# WEIGHT AMONG MEMBERS CLASSIFIED AT THAT RANK wins (members with no call
# at this rank don't count toward the denominator, so an unresolved
# centroid can't itself block consensus among the members that do have a
# call); once no label clears the threshold at a rank, that rank and all
# deeper ranks are left unclassified (NA).
MAJORITY_THRESHOLD <- 0.51

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

strip_size <- function(x) sub(";size=.*$", "", x)
parse_size <- function(x) as.numeric(sub("^.*;size=", "", x))

curated <- read.table(snakemake@input$curated_table, header = TRUE, row.names = 1,
                       sep = "\t", check.names = FALSE, comment.char = "")
centroid_seqs <- read_fasta(snakemake@input$centroids)

taxa <- readRDS(snakemake@input$taxa_species)
taxa <- taxa[, !duplicated(colnames(taxa), fromLast = TRUE), drop = FALSE]
ranks <- colnames(taxa)

asv_map <- read.table(snakemake@input$asv_map, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
seq_of_asv <- setNames(toupper(asv_map$sequence), asv_map$asv_id)

uc <- read.table(snakemake@input$uc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
uc <- uc[uc$V1 %in% c("S", "H"), ]

member_id     <- strip_size(uc$V9)
member_weight <- parse_size(uc$V9)
# S rows are the centroid of their own (new) cluster; H rows point at their
# centroid via column 10.
centroid_id <- ifelse(uc$V1 == "S", member_id, strip_size(uc$V10))

consensus_taxonomy <- function(otu_id) {
    in_cluster <- centroid_id == otu_id
    members <- member_id[in_cluster]
    weights <- member_weight[in_cluster]
    tax_mat <- taxa[seq_of_asv[members], , drop = FALSE]

    out <- setNames(character(length(ranks)), ranks)
    broken <- FALSE
    for (r in ranks) {
        if (broken) {
            out[r] <- NA
            next
        }
        tab <- tapply(weights, tax_mat[, r], sum)
        tab <- tab[!is.na(names(tab))]
        classified_weight <- sum(tab)
        if (length(tab) == 0 || max(tab) / classified_weight < MAJORITY_THRESHOLD) {
            out[r] <- NA
            broken <- TRUE
        } else {
            out[r] <- names(tab)[which.max(tab)]
        }
    }
    out
}

otu_ids <- rownames(curated)
consensus <- t(vapply(otu_ids, consensus_taxonomy, character(length(ranks))))
rownames(consensus) <- otu_ids
colnames(consensus) <- ranks

otu_seqs <- centroid_seqs[otu_ids]

out <- cbind(
    otu_id   = otu_ids,
    sequence = otu_seqs,
    as.data.frame(consensus, stringsAsFactors = FALSE),
    curated
)

write.table(out, snakemake@output$otu_taxonomy, sep = "\t", quote = FALSE, row.names = FALSE)
