suppressMessages(library("phyloseq",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("ggplot2",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("RColorBrewer",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("patchwork",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("vegan",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tidyverse",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("readr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("dyplr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tibble",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("plyr", quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("ALDEx2"), quietly = TRUE, warn.conflicts = FALSE))

phylo_obj <- snakemake@input
subset <- snakemake@params[["subset"]]

#Load phyloseq object saved in results' directory
metagenome <- load(phylo_obj)


#Use data only on bacteria ?
bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')

#To specify taxonomic levels for analysis (e.g., ‘kingdom’, ‘phylum’, ‘class’, etc.) if subset==0 then taxonomic levels from ‘kingdom’ down to ‘species’ are included

if (subset!=0){
bacteria_meta <- subset_taxa(bacteria_meta, subset)
}


tax_labels <- bacteria_meta@tax_table@.Data #Labels for different taxon

mat_count <- bacteria_meta@otu_table@.Data #Microbial count matrix

group_info <- bacteria_meta@sam_data #Metadata about the samples, such as their experimental conditions or treatment groups


