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

source("./scripts/CoDaSeq_functions.R")



phylo_obj <- unlist(snakemake@input[["phylo_obj"]])
subset <- snakemake@params[["subset"]]
groups <- snakemake@params[["group"]]
lineage <- snakemake@params[["lineage"]]
lineage_rank <- snakemake@params[["lineage_rank"]]
output_dir <- unlist(snakemake@params[["output_dir"]])
top_taxa <- snakemake@params[["top_taxa"]]
filter_by_sample <- snakemake@params[["filter_by_sample"]]
bh_fdr_cutoff <- snakemake@params[["bh_fdr_cutoff"]]
found <- snakemake@params[["found"]]
label_significant <- snakemake@params[["label_significant"]]
color_significant <- snakemake@params[["color_significant"]]
#Load phyloseq object saved in results' directory
metagenome <- readRDS(phylo_obj)

#Use data only on bacteria ?
bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')

sample_data <- bacteria_meta@sam_data

taxon <- as.data.frame(bacteria_meta@tax_table@.Data)

if (lineage!=""){
 bacteria_meta <- subset_taxa(bacteria_meta, get(lineage_rank) %in% lineage)
}

#To specify taxonomic levels for analysis (e.g., ‘kingdom’, ‘phylum’, ‘class’, etc.) if subset=="" then taxonomic levels from ‘kingdom’ down to ‘species’ are included

if (length(unique(sample_data$groups)) != 2){
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
        toptax <- topk(toptaxa)
        f1 <- filterfun_sample(toptax)
        wh1 <- genefilter_sample(bacteria_level,f1,A=filter_by_sample)
        bacteria_level <- prune_taxa(wh1,bacteria_level)
        mat <- bacteria_level@otu_table@.Data
        # convert to CLR
        x.all <- aldex(mat, conditions=sample_data[[group]], mc.samples=128, test="t", effect=TRUE, denom = 'all', verbose = F)

        found.by.one<-which(x.all$we.eBH<bh_fdr_cutoff|x.all$wi.eBH<bh_fdr_cutoff)
        found.by.all<-which(x.all$we.eBH<bh_fdr_cutoff&x.all$wi.eBH<bh_fdr_cutoff)

        if (found=="all"){
            found <- found.by.all
        } else if (found=="one"){
            found <- found.by.one
        }
        species.found <- taxon[rownames(taxon)%in%rownames(x.all)[found],]$Species
        genus.found <- taxon[rownames(taxon)%in%rownames(x.all)[found],]$Genus

        if (color_significant=="Genus"){
            colors <- rainbow(length(unique(genus.found)))
            names(colors) <- unique(genus.found)
            colors_labels <- genus.found
        } else if (color_significant=="Species"){
             colors <- rainbow(length(unique(species.found)))
             names(colors) <- unique(species.found)
             colors_labels <- species.found
        }

        if (label_significant=="Genus"){
            labels <- genus.found
        } else if (label_significant=="Species"){
            labels <- species.found
        }

         par(mfrow=c(1,1))

         aldex.plot(x.all, type="MA", test="welch",cutoff.pval=bh_fdr_cutoff,called.cex=1,called.col="red")
         text(x = x.all$rab.all[found],
             y = x.all$diff.btw[found],
             labels = labels,
             pos = 3,
             cex = 0.5,col = colors[colors_labels])
         legend_colors <- colors
         legend_labels <- names(colors)
         title(main = sprintf("Bland Altman plot (q=%0.2f)",bh_fdr_cutoff))
         legend("topright", legend = legend_labels, col = legend_colors, pch = 19, cex = 0.8)

         aldex.plot(x.all, type="MW", test="welch", cutoff.pval=bh_fdr_cutoff,called.cex=1,called.col="red")
         text(x = x.all$diff.win[found],
             y = x.all$diff.btw[found],
             labels = labels,
             pos = 3,
             cex = 0.5,col = colors[colors_labels])
         legend_colors <- colors
         legend_labels <- names(colors)
         title(main = sprintf("Effect plot (q=%0.2f)",bh_fdr_cutoff))
         legend("topright", legend = legend_labels, col = legend_colors, pch = 19, cex = 0.8)

         aldex.plot(x.all, type="volcano", test="welch", cutoff.pval=bh_fdr_cutoff,called.cex=1,called.col="red")
         text(x = x.all$diff.btw[found],
             y = -1*median(log10(x.all$we.eBH[found])),
             labels = labels,
             pos = 3,
             cex = 0.5,col = colors[colors_labels])
         legend_colors <- colors
         legend_labels <- names(colors)
         title(main = sprintf("Volcano plot (q=%0.2f)",bh_fdr_cutoff))
         legend("topright", legend = legend_labels, col = legend_colors, pch = 19, cex = 0.8)

        dev.off()
    }
}