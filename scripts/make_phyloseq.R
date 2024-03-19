


suppressMessages(library("phyloseq",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("ggplot2",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("RColorBrewer",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("patchwork",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("vegan",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tidyverse",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("readr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("dplyr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tibble",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("plyr",quietly = TRUE, warn.conflicts = FALSE))


source('alternative_import_biom.R')

# Import biom file to a phyloseq object
samples <- snakemake@SAMPLES
if (length(samples)>1){
metagenome <- import_biom(snakemake@input[["biom_file"]])
} else {
metagenome <- new_import_biom(snakemake@input[["biom_file"]])
}
metadata <- read.csv(snakemake@input[["metadata"]],sep=';',row.names = 1)
metagenome@sam_data <- sample_data(metadata)

# Change taxon names
colnames(metagenome@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove taxa which is not present in one sample at least
metagenome <- prune_taxa(taxa_sums(metagenome)>0,metagenome)

# Rename taxa
metagenome@tax_table <- substring(metagenome@tax_table,4)

output <- snakemake@output
#Save phyloseq object
saveRDS(metagenome,file=unlist(output),compress=FALSE)
