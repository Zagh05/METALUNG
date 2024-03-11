#A differential abundance analysis for the comparison of two or more conditions. For example,
#single-organism and meta-RNA-seq high-throughput sequencing assays, or of selected and unselected values from in-vitro sequence selections. Uses a Dirichlet-multinomial model to infer abundance from counts, that has been optimized for three or more experimental replicates. Infers sampling variation and calculates the expected false discovery rate given the biological and sampling
#variation using the Wilcox rank test or Welches t-test (aldex.ttest) or the glm and Kruskal Wallis
#tests (aldex.glm). Reports both P and fdr values calculated by the Benjamini Hochberg correction.





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
suppressMessages(library("ALDEx2"), quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library("zCompositions"), quietly=TRUE, warn.conflicts = FALSE))
suppressMessages(library("compositions"), quietly=TRUE, warn.conflicts = FALSE))

source("CoDaSeq_functions.R")



phylo_obj <- unlist(snakemake@input[["phylo_obj"]])
subset <- snakemake@params[["subset"]]
groups <- snakemake@params[["group"]]
lineage <- snakemake@params[["lineage"]]
lineage_rank <- snakemake@params[["lineage_rank"]]
output_dir <- unlist(snakemake@params[["output_dir"]])

#Load phyloseq object saved in results' directory
metagenome <- readRDS(phylo_obj)


#Use data only on bacteria ?
bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')
sample_data <- bacteria_meta@sam_data

if (lineage!=""){
 bacteria_meta <- subset_taxa(bacteria_meta, lineage_rank == lineage)
}

#To specify taxonomic levels for analysis (e.g., ‘kingdom’, ‘phylum’, ‘class’, etc.) if subset=="" then taxonomic levels from ‘kingdom’ down to ‘species’ are included

if (length(unique(sample_data$group)) != 2){
    warning('Will not do ALDEx2 differential abundance with !=2 groups')}

if (nrow(sample_data) < 3){
    warning('Will not do compositional data analysis with < 3 samples')
} else {
    if(subset!=""){
        do.tax.levels <- subset
    } else {
        do.tax.levels <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    }

    print('Compositional data analysis....')
    tax.pdf <- output_dir
    pdf(tax.pdf, height=8, width=7)
    for (tax.level in do.tax.levels[!(do.tax.levels %in% c('Kingdom'))]){
        print(paste('.....', tax.level))
        bacteria_level <- tax_glom(bacteria_meta, taxrank=tax.level)
        #bacteria_level_norm <- transform_sample_counts(bacteria_level, function(x) x/sum(x))
        #bacteria_level <- prune_taxa(taxa_sums(bacteria_meta_norm)>abundance_threshold, bacteria_meta_norm)
        #bacteria_level <- prune_taxa(taxa_sums(bacteria_meta_norm > 0) / nsamples(bacteria_meta_norm) > otu_cutoff, bacteria_meta_norm)
        mat <- bacteria_level@otu_table@.Data
        # convert to CLR
        x.all <- aldex(mat, conditions=sample_data[[group]], mc.samples=128, test="t", effect=TRUE, denom = 'all', verbose = F)

        #Set up a 1x2 plot layout
        par(mfrow=c(1,2))

        aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference")
        aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",ylab="Difference")

        #Add a main title
        mtext(paste("Differential abundance on",tax.level,"level"), side=3,line=-2,outer=TRUE)
        dev.off()
    }
}