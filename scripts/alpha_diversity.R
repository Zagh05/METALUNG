
suppressMessages(library("phyloseq",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("ggplot2",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("RColorBrewer",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("patchwork",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("vegan",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tidyverse",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("readr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("dplyr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tibble",quietly = TRUE, warn.conflicts = FALSE))

phylo_obj <- snakemake@input[["phylo_obj"]]
measures <- snakemake@params[["measures"]]
title <- snakemake@params[["title"]]
color <- snakemake@params[["color"]]
shape <- snakemake@params[["shape"]]
output <- unlist(snakemake@output)


metagenome <- readRDS(unlist(phylo_obj))

bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')

bacteria_meta_perc <- transform_sample_counts(bacteria_meta, function(x) x*100 / sum(x))


png(output)
plot_richness(physeq=bacteria_meta, measures=measures, title=title, color=color)

dev.off()