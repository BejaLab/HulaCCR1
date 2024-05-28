library <- function(...) suppressPackageStartupMessages(base::library(...))
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readxl)

with(snakemake@input, {
    tsv_files <<- tsv
    datasets_file <<- datasets
    clades_file <<- clades
})
with(snakemake@params, {
    min_datasets <<- min_datasets
    min_rnas <<- min_rnas
    cutoff <<- cutoff
    ranks <<- ranks
    include_regex <<- include
    exclude_regex <<- exclude
    ingroup <<- ingroup
})
output_file <- unlist(snakemake@output)

files_sizes <- file.size(tsv_files)
sintax <- tsv_files[files_sizes > 0] %>%
    lapply(read.table, col.names = c("gene", "taxa", "strand", "taxonomy"), sep = "\t") %>%
    bind_rows %>%
    separate(gene, into = c("dataset", "label"), sep = "\\|") %>%
    filter(taxa != "")
sintax_stat <- group_by(sintax, dataset) %>%
    summarize(n_rnas = n())

datasets <- read_excel(datasets_file) %>%
    mutate(dataset = analysis) %>%
    left_join(sintax_stat, by = "dataset") %>%
    filter(n_rnas >= min_rnas) %>%
    arrange(habitat, modified_location) %>%
    mutate(modified_location = factor(modified_location, levels = unique(modified_location)))
subclades <- read.table(clades_file, sep = ",", col.names=c("dataset", "category", "count")) %>%
    filter(grepl(ingroup, category)) %>%
    mutate(data_type = "subclades", rank = "_") %>%
    filter(dataset %in% datasets$analysis) %>%
    arrange(category)

cumpaste = function(x, .sep = "; ") Reduce(function(x1, x2) paste(x1, x2, sep = .sep), x, accumulate = T)

taxa_data <- separate_rows(sintax, taxa, sep = ",") %>%
    extract(taxa, into = c("rank", "taxon", "bootstrap"), regex = "^(.):(.+)\\((.+)\\)", convert = T) %>%
    mutate(rank = factor(rank, levels = unique(rank))) %>%
    group_by(dataset, label) %>%
    mutate(taxonomy = cumpaste(taxon)) %>%
    filter(bootstrap >= cutoff, grepl(include_regex, taxonomy), !grepl(exclude_regex, taxonomy)) %>%
    mutate(taxonomy = sub(".+;([^;]+;[^;]+;[^;]+)$", "\\1", taxonomy)) %>%
    ungroup

data <- filter(taxa_data, dataset %in% datasets$analysis) %>%
    distinct(dataset, taxonomy, rank) %>%
    group_by(taxonomy) %>%
    mutate(n_datasets = n_distinct(dataset)) %>%
    ungroup %>%
    filter(n_datasets >= min_datasets) %>%
    arrange(-n_datasets) %>%
    mutate(data_type = "taxonomy", category = taxonomy) %>%
    filter(rank %in% ranks) %>%
    bind_rows(subclades) %>%
    left_join(datasets, by = "dataset") %>%
    filter(!is.na(modified_location)) %>%
    mutate(dataset = paste0(ifelse(type == "Metagenome", "*", ""), dataset)) %>%
    mutate(category = factor(category, levels = rev(unique(category)))) %>%
    mutate(rank = factor(rank, levels = c(ranks, "_")))

shapes <- c(
    taxonomy = 16,
    subclades = 15
)

p <- ggplot(data) +
    geom_point(aes(x = dataset, y = category, shape = data_type)) +
    facet_grid(rank ~ modified_location, scales = "free", space = "free", switch = "y") +
    scale_shape_manual(values = shapes) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), legend.position = "none")
ggsave(output_file, p, width = 14, height = 7)
