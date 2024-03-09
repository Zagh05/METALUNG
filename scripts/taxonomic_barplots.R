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
tax_ranks = snakemake@params[["tax_ranks"]]
abund_thres = snakemake@params[["abundance_threshold"]]
groups = snakemake@params[["groups"]]
abundance = snakemake@params[["abundance"]]
output <- unlist(snakemake@output)

load(phylo_obj)

bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')

bacteria_meta_perc <- transform_sample_counts(bacteria_meta, function(x) x*100 / sum(x))

# Explore samples at specific taxonomic levels

tax_plots <- list()

for (tax in tax_ranks) {

    group_plots <- list()

    for (g in groups) {


# Group all the OTUs that have the same taxonomy at a certain taxonomic rank

    if (abundance=="absolut"){
    bacteria_glom <- tax_glom(bacteria_meta, taxrank=tax)
    }
    else {
    bacteria_glom <- tax_glom(bacteria_meta_perc, taxrank=tax)
    }
# Melt phyloseq object into a dataframe to manipulate them with packages like ggplot2 and vegan

    bacteria_glom_df <- psmelt(bacteria_glom)
    str(bacteria_glom_df)

    bacteria_glom_df <- bacteria_glom_df[bacteria_glom_df$Abundance>abund_thres,]

    bacteria_glom_df$tax <- as.factor(bacteria_glom_df$tax)

# brewer.pal() is a function from the RColorBrewer package that provides color palettes based on ColorBrewer designs.

    colors <- colorRampPalette(brewer.pal(8,'Dark2'))(length(levels(bacteria_glom_df$tax)))

    plot <- ggplot(data=bacteria_glom_df, aes(x=Sample[, y=Abundance, fill=tax))+
        geom_bar(aes(), stat='identity', position='stack')+
        scale_fill_manual(values = colors)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        +labs(title = paste("Taxonomy barplot on",str(tax),"level"))
        facet_wrap(~g, scales="free_x")+



    group_plots[[g]] <- plot

    }

    tax_plots[[tax]] <- group_plots

}

# Save to pdf

for (tax in tax_ranks) {
    tax.pdf <- file.path(output, paste('_',tolower(tax), '.pdf', sep=''))
    pdf(tax.pdf, height=, width=)
    for (p in tax_plots[[tax]]){print(p)}
    dev.off()
    }