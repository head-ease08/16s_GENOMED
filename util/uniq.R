library(dplyr)

abund <- read.csv("results/abundance_table_silva.csv", check.names = FALSE)
meta  <- read.csv("sra_metadata.csv", sep = ",", check.names = FALSE)

srr_cols <- grep("^SRR", colnames(abund), value = TRUE)
mapping  <- meta[meta[, 1] %in% srr_cols, ]
mapping  <- data.frame(Run = mapping$Run, SampleName = mapping[["Sample Name"]])

mock_srr <- mapping[grepl("Zymo|ZIEL", mapping$SampleName), ]

for (mock_type in c("Zymo", "ZIEL1", "ZIEL2")) {
    srrs <- mock_srr$Run[grepl(mock_type, mock_srr$SampleName)]
    if (length(srrs) == 0) next

    counts   <- rowSums(abund[, srrs, drop = FALSE])
    detected <- unique(abund$Genus[counts > 0 & !is.na(abund$Genus)])

    cat("\n===", mock_type, "===\n")
    cat(sort(detected), sep = "\n")
}
