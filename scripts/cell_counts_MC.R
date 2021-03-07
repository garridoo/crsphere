
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# run script silently

sink(type="message")
options(warn=-1)

# cleanup

rm(list = ls())

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")
source("cpcoa.func.R")

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("sva")
library("dunn.test")
library("rcompanion")

# load paths 

source("paths.R")

# files

counts.file <- paste(results.dir, "cell_counts_MC.txt", sep="")

# load data

counts <- read.table(counts.file, header=T, sep="\t")
counts$day <- paste("Day", counts$day, sep="")
counts$day <- factor(gsub("Day", "", counts$day), levels=0:10)

counts$condition <- paste(counts$color, counts$condition)
counts$condition_day <- paste(counts$condition, counts$day)

idx <- counts$color=="blue" &
       !is.na(counts$chlamy_per_mL) & 
       counts$day %in% c("1", "4", "7", "10") 
       T
counts <- counts[idx, ]

# define colors and plot

colors <- data.frame(group=c("blue C", "blue B+C"),
                     color=c(cr_only_color, cr_color))

colors <- colors[colors$group %in% counts$condition, ]
counts$condition <- factor(counts$condition, levels=colors$group)

p <- ggplot(counts, aes(x=day, y=chlamy_per_mL, color=condition)) +
            geom_boxplot(alpha=1, outlier.color="transparent") +
            geom_jitter(position=position_jitterdodge(0.5), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_y_log10(breaks=c(0, 1e+05, 1e+06, 1e+07), limits=c(NA, 1e+07)) +
            annotation_logticks(short=unit(1.5,"mm"), mid=unit(2.5,"mm"), long=unit(3.5,"mm")) +
            labs(x="Day", y="Cell counts") +
            ggtitle("Cell counts mesocosm experiment") +
            main_theme +
            theme(legend.position="none")

ggsave(paste(figures.dir, "cell_counts_MC.pdf", sep=""), p, width=8, height=5.5)

# perform Mann-Whitney test for each timepoint

pvals <- data.frame(day=NULL, p.val=NULL)

for (t in unique(counts$day)) {
    
    a <- counts$chlamy_per_mL[counts$day==t & counts$condition=="blue C"]
    b <- counts$chlamy_per_mL[counts$day==t & counts$condition=="blue B+C"]

    p <- wilcox.test(a, b)$p.value

    pvals <- rbind(pvals, data.frame(day=t, p.val=p))

}

pvals$p.adj <- p.adjust(pvals$p.val, method="hochberg")
pvals <- pvals[pvals$p.adj < 0.05, ]

print(pvals)

