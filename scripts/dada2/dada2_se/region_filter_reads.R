log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

truncLen <- unlist(snakemake@params$truncLen)

saveRDS(
    filterAndTrim(
        snakemake@input$r1,  snakemake@output$r1,
        snakemake@input$r2,  snakemake@output$r2,
        truncLen    = truncLen,
        maxN        = 0,
        maxEE       = c(2, 2),
        truncQ      = 2,
        rm.phix     = TRUE,
        compress    = TRUE,
        multithread = TRUE
    ),
    snakemake@output$stats
)

# filterAndTrim skips writing output entirely when 0 reads survive
# (e.g. sample has no reads for this region in a multiplexed pool) --
# touch valid empty gzips so downstream rules still get their inputs.
for (f in c(snakemake@output$r1, snakemake@output$r2)) {
    if (!file.exists(f)) close(gzfile(f, "wb"))
}
