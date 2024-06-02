library(dplyr)
library(tidyr)

with(snakemake@input, {
    tsv_files <<- tsv
    clstr_file <<- clstr
})
output_file <- unlist(snakemake@output)

read.cdhit.clstr <- function(fname) 
{
    data.fields <- c("E.Value", "Aln", "Identity")
    read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
        separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
        fill(Cluster) %>%
        filter(!grepl(">", Seq.Num)) %>%
        separate(Col2, into = c("Seq.Len", "Col2"), sep = "aa, >") %>%
        extract(Col2, into = c("Seq.Name", "Is.Representative", "Col2"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>% 
        mutate(Is.Representative = Is.Representative == "*", Col2 = ifelse(Is.Representative, "100%", Col2)) %>% 
        group_by(Cluster) %>%
        mutate(Representative = Seq.Name[which(Is.Representative)]) %>% 
        separate_rows(Col2, sep = ",") %>%
        separate(Col2, into = data.fields, sep = "/", fill = "left", convert = T) %>%
        mutate(Identity = sub("%", "", Identity) %>% as.numeric) %>%
        group_by(Seq.Name) %>% 
        mutate(level.rank = paste0(".", 1:n() - 1), level.rank = ifelse(level.rank == ".0", "", level.rank)) %>%
        pivot_wider(names_from = level.rank, values_from = data.fields, names_sep = "") %>%
        ungroup
}

clstr <- read.cdhit.clstr(clstr_file) %>%
    extract(Seq.Name, into = "qseqid", regex = "\\|(.+):") %>%
    select(qseqid, Representative)

col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data <- lapply(tsv_files, read.table, col.names = col_names) %>%
    lapply(mutate_all, as.character) %>%
    bind_rows %>%
    extract(sseqid, into = "subclade", regex = c(".+\\|.+\\|(.+)\\|")) %>%
    distinct(qseqid, subclade) %>%
    left_join(clstr, by = "qseqid") %>%
    filter(!is.na(Representative)) %>%
    mutate(contig = sub(":.+", "", qseqid)) %>%
    distinct(contig, .keep_all = T) %>%
    group_by(Representative, subclade) %>%
    summarize(n = n()) %>%
    spread(subclade, n, fill = 0)

write.csv(data, file = output_file, row.names = F)
