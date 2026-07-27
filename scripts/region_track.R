log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(dada2)

regions <- snakemake@params$regions
samples <- snakemake@params$samples
rdir    <- snakemake@params$region_dir

# JSON от cutadapt читаем без jsonlite: нужны два числа, регулярка надёжнее
# лишней зависимости в envs/dada2.yaml.
read_cutadapt <- function(path) {
    if (!file.exists(path)) return(c(input = NA, output = NA))
    txt <- paste(readLines(path, warn = FALSE), collapse = " ")
    grab <- function(key) {
        m <- regmatches(txt, regexpr(paste0('"', key, '"\\s*:\\s*[0-9]+'), txt))
        if (length(m) == 0) NA_integer_ else as.integer(sub('.*:\\s*', '', m))
    }
    c(input = grab("input"), output = grab("output"))
}

rows <- list()
for (rg in regions) {
    nochim_path <- file.path(rdir, rg, "seqtab_nochim.rds")
    nochim <- if (file.exists(nochim_path)) readRDS(nochim_path) else NULL

    for (s in samples) {
        ca <- read_cutadapt(file.path("logs/demux_region", rg, paste0(s, ".json")))

        st_path <- file.path(rdir, rg, "filter_stats", paste0(s, ".rds"))
        filt <- NA_integer_
        if (file.exists(st_path)) {
            st <- readRDS(st_path)
            if (!is.null(dim(st)) && "reads.out" %in% colnames(st)) filt <- st[1, "reads.out"]
        }

        final <- NA_integer_
        if (!is.null(nochim) && s %in% rownames(nochim)) final <- sum(nochim[s, ])

        rows[[length(rows) + 1]] <- data.frame(
            region      = rg,
            sample      = s,
            raw_pairs   = ca["input"],
            after_demux = ca["output"],
            after_filter = filt,
            final_reads = final,
            stringsAsFactors = FALSE
        )
    }
}

df <- do.call(rbind, rows)
df$demux_pct <- round(100 * df$after_demux / df$raw_pairs, 2)
df$final_pct <- round(100 * df$final_reads / df$after_demux, 2)

write.table(df, snakemake@output$tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Сводка по регионам в лог — то, ради чего это правило и написано.
cat("\n=== доля ридов, дошедших до ASV, по регионам ===\n")
agg <- aggregate(cbind(raw_pairs, after_demux, final_reads) ~ region, df, sum, na.rm = TRUE)
agg$demux_pct <- round(100 * agg$after_demux / agg$raw_pairs, 1)
agg$final_pct <- round(100 * agg$final_reads / agg$after_demux, 1)
print(agg, row.names = FALSE)

cat("\nЧитать так:\n")
cat("  demux_pct низкий у одного региона  -> праймер не находится, поднимите error_rate\n")
cat("  final_pct низкий у merge-региона   -> не хватает перекрытия, регион надо в SE\n")
cat("  final_pct низкий у SE-региона      -> truncLen выше реальной длины рида\n")

for (i in seq_len(nrow(agg))) {
    if (!is.na(agg$final_pct[i]) && agg$final_pct[i] < 40)
        warning(sprintf("регион %s: до ASV дошло %.1f%% — разберитесь до интерпретации результатов",
                        agg$region[i], agg$final_pct[i]))
}
