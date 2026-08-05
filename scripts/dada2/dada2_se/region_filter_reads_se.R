log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dada2)

truncLen <- snakemake@params$truncLen

saveRDS(
    filterAndTrim(
        snakemake@input$r1, snakemake@output$r1,
        truncLen    = truncLen,
        maxN        = 0,
        maxEE       = 2,
        truncQ      = 2,
        rm.phix     = TRUE,
        compress    = TRUE,
        multithread = TRUE
    ),
    snakemake@output$stats
)

# filterAndTrim skips writing output entirely when 0 reads survive
# (e.g. sample has no reads for this region in a multiplexed pool) --
# touch a valid empty gzip so downstream rules still get their input.
if (!file.exists(snakemake@output$r1)) close(gzfile(snakemake@output$r1, "wb"))
