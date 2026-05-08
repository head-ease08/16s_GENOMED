forward.qc <- file.path(root.qc, paste0(sample.names, "_R1.fastq.gz"))
reverse.qc <- file.path(root.qc, paste0(sample.names, "_R2.fastq.gz"))
names(forward.qc) <- sample.names
names(reverse.qc) <- sample.names

qc.out <- filterAndTrim(forward.cut, forward.qc,
                        reverse.cut, reverse.qc,
                        truncLen  = c(240, 200),   # tune to your quality plots
                        maxN      = 0,
                        maxEE     = c(2, 2),
                        truncQ    = 2,
                        rm.phix   = TRUE,
                        compress  = TRUE,
                        multithread = TRUE)