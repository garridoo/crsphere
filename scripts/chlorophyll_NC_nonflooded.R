
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

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

chlorophyll.file <- paste(results.dir, "chlorophyll_NC_nonflooded.txt", sep="")

# load data

chlorophyll <- read.table(chlorophyll.file, header=T, sep="\t")

 
# define colors and plot

colors <- data.frame(group=c("covered", "uncovered"),
                     color=c(covered_color, uncovered_color))

colors <- colors[colors$group %in% chlorophyll$treatment, ]
chlorophyll$treatment <- factor(chlorophyll$treatment, levels=colors$group)

p <- ggplot(chlorophyll, aes(x=treatment, y=fluorescence, color=treatment)) +
            geom_boxplot(alpha=1, outlier.color="transparent") +
            geom_jitter(position=position_jitterdodge(0.5), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Treatment", y="Fluorescence") +
            main_theme +
            theme(legend.position="none")

ggsave(paste(figures.dir, "chlorophyll_NC_nonflooded.pdf", sep=""), p, width=3, height=4)

# perform Mann-Whitney test

a <- chlorophyll$fluorescence[chlorophyll$treatment=="covered"]
b <- chlorophyll$fluorescence[chlorophyll$treatment=="uncovered"]
p <- wilcox.test(a, b)$p.value
print(p)
 
