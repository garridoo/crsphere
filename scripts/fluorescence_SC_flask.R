
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

tecan.file <- paste(results.dir, "fluorescence_SC_flask.txt", sep="")

# load data

tecan <- read.table(tecan.file, header=T, sep="\t")

# cast fluorescence data table into a matrix

dat <- cast(tecan, cond ~ Sample, value="Value")
rownames(dat) <- dat[, 1]
dat <- data.frame(dat[, -1])
dat <- as.matrix(dat)

# batch correction using ComBat

design <- unique(data.frame(Sample=tecan$Sample, split_exp=tecan$split_exp, Vessel=tecan$Vessel, replicate=tecan$replicate, Content=tecan$Content))

mod <- model.matrix(~1, data=design)
batch <- design$split_exp
combat_values <- ComBat(dat=dat, batch=batch, mod=NULL, par.prior=T, mean.only=T, ref.batch=NULL)

# create dataframe for plotting

df <- melt(combat_values)
colnames(df) <- c("condition", "sample", "value")

df <- cbind(df, design[match(df$sample, design$Sample), ])
df$day <- factor(gsub(".*_day", "", df$condition), levels=c(1, 4, 7, 11))

# define colors and plot

colors <- data.frame(group=c("C", "C+SC"),
                     color=c(cr_only_color, cr_color))

colors <- colors[colors$group %in% df$Content, ]
df$Content <- factor(df$Content, levels=colors$group)

p <- ggplot(df, aes(x=day, y=value, color=Content)) +
            geom_boxplot(alpha=1, outlier.color="transparent") +
            geom_jitter(position=position_jitterdodge(0.5), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_y_log10(breaks=c(0, 100, 300, 1000), limits=c(NA, max(df$value))) +
            annotation_logticks(short=unit(1.5,"mm"), mid=unit(2.5,"mm"), long=unit(3.5,"mm")) +
            labs(x="Day", y="Fluorescence") +
            ggtitle("Fluorescence SynCom flask experiment") +
            main_theme +
            theme(legend.position="none")

ggsave(paste(figures.dir, "fluorescence_SC_flask.pdf", sep=""), p, width=8, height=4)

# perform Mann-Whitney test for each timepoint

pvals <- data.frame(day=NULL, p.val=NULL)

for (t in unique(df$day)) {
    
    a <- df$value[df$day==t & df$Content=="C"]
    b <- df$value[df$day==t & df$Content=="C+SC"]

    p <- wilcox.test(a, b)$p.value

    pvals <- rbind(pvals, data.frame(day=t, p.val=p))

}

pvals$p.adj <- p.adjust(pvals$p.val, method="hochberg")
pvals <- pvals[pvals$p.adj < 0.05, ]

print(pvals)

