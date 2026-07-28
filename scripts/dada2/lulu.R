log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

if (!requireNamespace("lulu", quietly = TRUE)) {
    remotes::install_github("tobiasgf/lulu")
}
library(lulu)

otu_tab   <- read.table(snakemake@input$otu_table,  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
match_tab <- read.table(snakemake@input$match_list, header = FALSE, sep = "\t")

curated <- lulu(otu_tab, match_tab)

saveRDS(curated, snakemake@output$lulu_result)
write.table(curated$curated_table, snakemake@output$curated_table, sep = "\t", quote = FALSE, col.names = NA)