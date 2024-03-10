suppressMessages(library("phyloseq",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("ggplot2",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("RColorBrewer",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("patchwork",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("vegan",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tidyverse",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("readr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("dplyr",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("tibble",quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("plyr", quietly = TRUE, warn.conflicts = FALSE))

phylo_obj <- unlist(snakemake@input)
method <- snakemake@params[["method"]]
distance <- snakemake@params[["distance"]]
color <- snakemake@params[["color"]]
shape <- snakemake@params[["shape"]]
title <- snakemake@params[["title"]]
type <- snakemake@params[["type"]]
wrap <- snakemake@params[["wrap"]]
output <- unlist(snakemake@output)

metagenome <- readRDS(phylo_obj)

bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')

bacteria_meta_perc <- transform_sample_counts(bacteria_meta, function(x) x*100 / sum(x))

bacteria_meta_perc.ord <- ordinate(physeq = bacteria_percentages, method = method, distance=distance)


png(output)

if !(wrap){
    plot_ordination(bacteria_meta_perc, bacteria_meta_perc.ord, type=type, color=color, title=title, shape=shape)

}
else {
    plot_ordination(bacteria_meta_perc, bacteria_meta_perc.ord, type=type, color=color, title=title, shape=shape) +
            facet_wrap(~wrap,scales="free_x")
}


dev.off()

