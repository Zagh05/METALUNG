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
source("CoDaSeq_functions.R")



phylo_obj <- unlist(snakemake@input)
subset <- snakemake@params[["subset"]]
group <- snakemake@params[["group"]]
output_dir <- snakemake@params[["output_dir"]]

#Load phyloseq object saved in results' directory
load(phylo_obj)


#Use data only on bacteria ?
bacteria_meta <- subset_taxa(metagenome, Kingdom=='Bacteria')
sample_data <- bacteria_meta@sam_data
#To specify taxonomic levels for analysis (e.g., ‘kingdom’, ‘phylum’, ‘class’, etc.) if subset==0 then taxonomic levels from ‘kingdom’ down to ‘species’ are included
if (length(unique(sample_data$group)) != 2){
    warning('Will not do ALDEx2 differential abundance with !=2 groups')}

if (nrow(sample_data) < 3){
    warning('Will not do compositional data analysis with < 3 samples')
} else {
    pca.plot.list <- list()
    if(segata){
        do.tax.levels <- c('Species')
    } else {
        do.tax.levels <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    }

    print('Compositional data analysis....')
    for (tax.level in do.tax.levels[!(do.tax.levels %in% c('Kingdom'))]){
        print(paste('.....', tax.level))
        bacteria_level <- subset_taxa(bacteria_meta, tax.level)
        use.mat <- bacteria_level@otu_table@.Data
        if(nrow(use.mat) > 2){
            reads.filtered <- tryCatch(codaSeq.filter(use.mat,
                min.reads=codata_min_reads, min.occurrence=codata_otu_cutoff, min.prop=codata_min_prop, samples.by.row=FALSE),
            error=function(e) {
                warning(e)
                matrix(0)
                })
            } else {
                reads.filtered <- matrix(0)
            }
        if(nrow(reads.filtered) > 2){
            # replace 0 values with an estimate of the probability that the zero is not 0
            # only if necessary
            # at this point samples are by row, and variables are by column
            # we can use the GBM or CZM methods. Function from zCompositions.
            if (sum(reads.filtered==0)>0){
                rfz <- t(cmultRepl(t(reads.filtered),  label=0, method="CZM", output="p-counts"))
                # round to integers and round anything low to 1
                rfz <- round(rfz)
                rfz[rfz==0] <- 1
            } else {
                rfz <- reads.filtered
            }
            # convert to CLR
            clr.aldex <- aldex.clr(rfz, conds=rep(NA, ncol(rfz)), mc.samples=128, denom = 'all', verbose = F)
            # get clr from mean of the MC instances
            rfz.clr <- t(sapply(clr.aldex@analysisData , function(x) rowMeans(x)))
            # save CLR matrix
            outf <- file.path(output_dir, paste('clr_values_', tax.level, '.tsv', sep=''))
            write.table(round(rfz.clr, 4), outf, sep='\t', quote=F, row.names = T, col.names = T)

            rfz.clr.pca <- prcomp(rfz.clr)
            rfz.clr.mvar <- mvar(rfz.clr)

            # PCA biplot
            plot.df <- data.frame(rfz.clr.pca$x[,1:2], group=sample_data[sample_data$sample %in% rownames(rfz.clr), 'group'])
            pc1.var <- round(sum(rfz.clr.pca$sdev[1]^2)/rfz.clr.mvar * 100, 1)
            pc2.var <- round(sum(rfz.clr.pca$sdev[2]^2)/rfz.clr.mvar * 100, 1)
            pca.plot <- ggplot(plot.df, aes(x=PC1, y=PC2, col=group)) +
                geom_point(size=3) +
                geom_text(label=rownames(plot.df),
                          nudge_x = sum(abs(range(plot.df$PC1)))/ 50,
                          nudge_y = sum(abs(range(plot.df$PC2)))/ 50) +
                theme_bw() +
                labs(x=paste("PC1: ", pc1.var, "% of variance", sep=""),
                     y=paste("PC2: ", pc2.var, "% of variance", sep=""),
                     title=paste('Compositional PCA,', tax.level)) +
                scale_color_brewer(palette='Set1')
            pca.plot.list[[tax.level]] <- pca.plot


            # need to limit to two groups for aldex
            if (length(unique(sample_data$group)) == 2){
                # get two groups
                group.pos <- unique(sample_data$group)[1]
                group.neg <- unique(sample_data$group)[2]
                s.pos <- sample_data$sample[sample_data$group==group.pos]
                s.neg <- sample_data$sample[sample_data$group==group.neg]

                # aldex differential expression between groups
                aldex.res <- aldex(rfz[,c(s.pos, s.neg)],
                    conditions = c(rep(group.pos, length(s.pos)), rep(group.neg, length(s.neg))),
                    mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE,
                    denom="all", verbose=FALSE)

                # scatterplot
                outf.scatter <- file.path(output_dir, paste('aldex_scatter_', tax.level, '.pdf', sep=''))
                pdf(outf.scatter, width = 7, height = 4)
                mypar(1,2, mar = c(3.5, 3.5, 2.5, 1.1))
                aldex.plot(aldex.res, type="MW", test="wilcox", xlab="Dispersion",ylab="Difference")
                mtext(paste(tax.level, ', wilcox', sep=''), line=0.25, cex=1.25)
                aldex.plot(aldex.res, type="MW", test="welch", xlab="Dispersion",ylab="Difference")
                mtext(paste(tax.level, ', welch', sep=''), line=0.25, cex=1.25)
                dev.off()

                # improve dataframe
                aldex.res <- signif(aldex.res, 4)
                # add sample numbers
                aldex.res$n.pos <- length(s.pos)
                aldex.res$n.neg <- length(s.neg)
                # add in name
                aldex.res$taxa <- rownames(aldex.res)
                # reorder
                aldex.res <- aldex.res[, c('taxa', 'n.pos', 'n.neg', colnames(aldex.res)[1:11])]
                aldex.res <- aldex.res[order(aldex.res$we.eBH),]
                # add abs effect
                aldex.res$abs.effect <- abs(aldex.res$effect)

                # save dataframe of differential results
                outf.res <- file.path(output_dir, paste('aldex_result_', tax.level, '.tsv', sep=''))
                write.table(aldex.res, outf.res, sep='\t', quote=F, row.names = F, col.names = T)

                # make boxplots
                # make plots from everything with this value or less
                plot.thresh <- 0.1
                plot.thresh.col <- 'we.eBH'
                aldex.res.plot <- aldex.res[aldex.res[, plot.thresh.col] <= plot.thresh, ]
                outf.boxplot <- file.path(output_dir, paste('aldex_significant_boxplots_', tax.level, '.pdf', sep=''))
                if (nrow(aldex.res.plot) >0 ){
                    # make clr values
                    m.prop <-  apply(rfz, 2, function(x) x/sum(x))
                    clr.aldex <- aldex.clr(rfz[,c(s.pos, s.neg)],
                        conds = c(rep(group.pos, length(s.pos)), rep(group.neg, length(s.neg))),
                        mc.samples=128, denom = 'all', verbose = F)
                    # get clr from mean of the MC instances
                    rfz.clr <- sapply(clr.aldex@analysisData , function(x) rowMeans(x))

                    # save all significant in one boxplot
                    pdf(file.path(output_dir,'pdf_boxplot.pdf'), height=5, width = 6, onefile=TRUE)
                    toplot <- min(nrow(aldex.res.plot), 50)
                    for ( i in 1:toplot){
                        g <- aldex.res.plot[i, 'taxa']
                        test.df <- data.frame(group = c(rep(group.pos, length(s.pos)), rep(group.neg, length(s.neg))),
                          sample = c(s.pos, s.neg),
                          proportion = c(m.prop[g, s.pos], m.prop[g, s.neg]),
                          clr = c(rfz.clr[g, s.pos], rfz.clr[g, s.neg]))


                        g1 <- ggplot(test.df, aes(x=group, y=proportion, fill=group)) +
                        geom_boxplot() +
                        theme_bw() +
                        labs(subtitle=g) +
                        ylim(c(0,max(test.df$proportion)*1.2)) +
                        stat_compare_means(method = 'wilcox') +
                        stat_compare_means(method = 't.test', label.y.npc = 0.95) + guides(fill=FALSE)

                        g2 <- ggplot(test.df, aes(x=group, y=clr, fill=group)) +
                        geom_boxplot() +
                        theme_bw() +
                        labs(subtitle=g) +
                        ylim(c(min(test.df$clr), max(test.df$clr)*1.2)) +
                        stat_compare_means(method = 'wilcox') +
                        stat_compare_means(method = 't.test', label.y.npc = 0.95)

                        print(ggarrange(g1,g2,widths = c(1,1), common.legend = T))
                    }
                    dev.off()
                }
            }
        }
    }

    # save pca plots from each level
    pdf(file.path(output_dir, 'compositional_PCA_plot.pdf'), height=6.5, width=9)
    for (tn in rev(do.tax.levels)){
        print(pca.plot.list[[tn]])
    }
    trash <- dev.off()
}



message('Done! :D')

