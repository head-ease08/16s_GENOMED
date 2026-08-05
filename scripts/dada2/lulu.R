log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

if (!requireNamespace("lulu", quietly = TRUE)) {
    remotes::install_github("tobiasgf/lulu")
}
library(lulu)

otu_tab <- read.table(snakemake@input$otu_table,  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, comment.char = "")

# 0 ASVs survived upstream -- otu_table has 0 rows and match_list.txt is a
# 0-byte file in that case, which read.table can't parse ("no lines
# available in input"), and lulu() has nothing to curate anyway.
if (nrow(otu_tab) == 0) {
    saveRDS(NULL, snakemake@output$lulu_result)
    write.table(otu_tab, snakemake@output$curated_table, sep = "\t", quote = FALSE, col.names = NA)
    quit(save = "no", status = 0)
}

match_tab <- read.table(snakemake@input$match_list, header = FALSE, sep = "\t")

curated <- lulu(otu_tab, match_tab)

saveRDS(curated, snakemake@output$lulu_result)
write.table(curated$curated_table, snakemake@output$curated_table, sep = "\t", quote = FALSE, col.names = NA)