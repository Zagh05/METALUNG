


library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("vegan")
library("tidyverse")
library("readr")
library("dyplr")
library("tibble")


# Import biom file to a phyloseq object

metagenome <- import_biom(snakemake@input[["biom_file"]])
metadata <- read.csv(snakemake@input[["metadata"]],sep='\t')
metagenome@sam_data <- sample_data(metadata)

# Change taxon names
colnames(metagenome@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove taxa which is not present in one sample at least
metagenome <- prune_taxa(taxa_sums(metagenome)>0,metagenome)

#Save phyloseq object
save(metagenome, snakemake@output)
