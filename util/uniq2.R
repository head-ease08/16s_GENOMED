abund <- read.csv("results/abundance_table_silva.csv", check.names = FALSE)

for (mock in c("zymo", "ziel1", "ziel2")) {
    srrs <- readLines(paste0("/tmp/", mock, "_srr.txt"))
    srrs <- intersect(srrs, colnames(abund))
    if (length(srrs) == 0) { cat(mock, ": no matching SRRs\n"); next }

    counts <- rowSums(abund[, srrs, drop = FALSE])
    genera <- unique(abund$Genus[counts > 0 & !is.na(abund$Genus)])

    cat("\n===", toupper(mock), "===\n")
    cat(sort(genera), sep = "\n")
}
