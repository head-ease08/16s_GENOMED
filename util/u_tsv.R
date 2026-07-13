abund_silva <- read.csv("results/abundance_table_silva.csv", check.names = FALSE)
abund_rdp   <- read.csv("results/abundance_table_rdp.csv",   check.names = FALSE)

results <- list()

for (mock in c("zymo", "ziel1", "ziel2")) {
    srrs <- readLines(paste0("/tmp/", mock, "_srr.txt"))

    for (db in c("silva", "rdp")) {
        abund <- if (db == "silva") abund_silva else abund_rdp
        cols  <- intersect(srrs, colnames(abund))
        if (length(cols) == 0) next

        counts  <- rowSums(abund[, cols, drop = FALSE])
        present <- abund[counts > 0, c("Genus", "Species")]
        present <- present[!is.na(present$Genus), ]
        present <- unique(present)
        present$mock <- toupper(mock)
        present$db   <- toupper(db)

        results[[paste(mock, db)]] <- present
    }
}

out <- do.call(rbind, results)
out <- out[, c("mock", "db", "Genus", "Species")]
write.table(out, "results/mock_detected_taxa.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
