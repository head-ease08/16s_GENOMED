library(dada2)

getN <- function(x) sum(getUniques(x))

sample_names  <- sub("_R1\\.rds$", "", basename(snakemake@input$dada_fwd))
filter_stats  <- do.call(rbind, lapply(snakemake@input$filter_stats, readRDS))
dada_fwd      <- lapply(snakemake@input$dada_fwd, readRDS)
mergers       <- lapply(snakemake@input$mergers, readRDS)
seqtab_nochim <- readRDS(snakemake@input$seqtab_nochim)

track <- data.frame(
    sample   = sample_names,
    input    = filter_stats[, 1],
    filtered = filter_stats[, 2],
    denoised = sapply(dada_fwd, getN),
    merged   = sapply(mergers, getN),
    nonchim  = rowSums(seqtab_nochim)
)

write.csv(track, snakemake@output$summary, row.names = FALSE)
