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


phylo_obj <- snakemake@input
subset <- snakemake@params[["subset"]]
method <- snakemake@params[["method"]]
distance <- snakemake@params[["distance"]]
sample_label <- snakemake@params[["sample_label"]]
taxa_label <- snakemake@params[["taxa_label"]]
wrap <- snakemake@params[["wrap"]]

metagenome <- load(phylo_obj)

bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')

bacteria_meta_perc <- transform_sample_counts(bacteria_meta, function(x) x*100 / sum(x))

bacteria_subset <- subset_taxa(bacteria_meta_perc, subset)

png("")

if !(wrap){
    plot_heatmap(bacteria_subset, method = method, distance = distance, sample.label = sample_label, taxa.label = taxa_label, title = title)

}
else {
    plot_heatmap(bacteria_subset, method = method, distance = distance, sample.label = sample_label, taxa.label = taxa_label, title = title) +
            facet_wrap(~wrap,scales="free_x")
}

dev.off()