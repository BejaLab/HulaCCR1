library(dplyr)
library(tidyr)

tsv_files <- unlist(snakemake@input)
datasets <- unlist(snakemake@params["datasets"])
output_file <- unlist(snakemake@output)
col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data <- lapply(tsv_files, read.table, col.names = col_names) %>%
    setNames(datasets) %>%
    Filter(NROW, .) %>%
    bind_rows(.id = "dataset") %>%
    extract(sseqid, into = "clade", regex = "@([^@]+)$") %>%
    distinct(qseqid, .keep_all = T) %>%
    group_by(dataset, clade) %>%
    summarize(num = n_distinct(qseqid), .groups = "drop")
write.table(data, file = output_file, row.names = F, col.names = F, sep = ",", quote = F)
