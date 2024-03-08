


suppressMessages(library("phyloseq",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("ggplot2",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("RColorBrewer",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("patchwork",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("vegan",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tidyverse",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("readr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("dplyr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tibble",quietly = TRUE, warn.conflicts = FALSE))


# Import biom file to a phyloseq object

print(snakemake@input[["biom_file"]])
print(snakemake@input[["metadata"]])
print(snakemake@output[1])
metagenome <- import_biom(file.path(".."+snakemake@input[["biom_file"]]))
metadata <- read.csv(file.path("..",snakemake@input[["metadata"]]),sep=';')
metagenome@sam_data <- sample_data(metadata)

# Change taxon names
colnames(metagenome@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove taxa which is not present in one sample at least
metagenome <- prune_taxa(taxa_sums(metagenome)>0,metagenome)

#Save phyloseq object
save(metagenome,file=file.path("..",snakemake@output))
