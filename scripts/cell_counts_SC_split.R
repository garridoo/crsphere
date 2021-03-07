
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

counts.file <- paste(results.dir, "cell_counts_SC_split.txt", sep="")

# load data

counts <- read.table(counts.file, header=T, sep="\t")

# cast fluorescence data table into a matrix

dat <- cast(counts, condition ~ Sample, value="Value", fun.aggregate=sum, fill=NA)
rownames(dat) <- dat[, 1]
dat <- data.frame(dat[, -1])
dat <- as.matrix(dat)


# split conditions not present in every batch (SC only and axenic controls)

idx <-rowSums(is.na(dat))!=0
last <- dat[idx, ]
dat <- dat[!idx, ]

# batch correction using ComBat

design <- unique(data.frame(Sample=counts$Sample, split_exp=counts$split_exp))

mod <- model.matrix(~1, data=design)
batch <- design$split_exp
combat_values <- ComBat(dat=dat, batch=batch, mod=NULL, par.prior=T, mean.only=T, ref.batch=NULL)

combat_values <- rbind(combat_values, last)
dat <- rbind(dat, last)

# create dataframe for plotting

df <- melt(dat)
colnames(df) <- c("condition", "sample", "value")

df <- cbind(df, design[match(df$sample, design$Sample), ])
 
# define colors (by content) and shapes (by neighbor) and plot

colors <- data.frame(group=c("C/TP",    "C/SC",   "C/C+SC",
                             "C+SC/TP", "C+SC/C", "C+SC/SC",
                             "SC/TP",   "SC/C",   "SC/C+SC",
                             "TP/C",    "TP/SC",  "TP/C+SC"),
                     color=c(cr_only_color,  cr_only_color,  cr_only_color,
                             cr_color,       cr_color,       cr_color,
                             bacteria_color, bacteria_color, bacteria_color,
                             c_grey,         c_grey,         c_grey))

shapes <- data.frame(group=c("C/TP",    "C/SC",   "C/C+SC",
                             "C+SC/TP", "C+SC/C", "C+SC/SC",
                             "SC/TP",   "SC/C",   "SC/C+SC",
                             "TP/C",    "TP/SC",  "TP/C+SC"),
                     shape=c(1,         17,       15,
                             1,         16,       18,
                             1,         16,       15,
                             16,        18,       15))     
 
colors <- colors[colors$group %in% df$condition, ]
df$condition <- factor(df$condition, levels=colors$group)

shapes <- shapes[shapes$group %in% df$condition, ]
df$condition <- factor(df$condition, levels=shapes$group)

p <- ggplot(df, aes(x=condition, y=value, color=condition, shape=condition)) +
            geom_boxplot(alpha=1, outlier.color="transparent") +
            geom_jitter(position=position_jitterdodge(2.5), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0, max(df$value))) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="Condition", y="Cell counts") +
            ggtitle("Cell counts SynCom split experiment") +
            main_theme +
            theme(legend.position="none")

ggsave(paste(figures.dir, "cell_counts_SC_split.pdf", sep=""), p, width=12, height=5)

# perform Dunn post hoc test (after rejecting the NULL hypothesis with Kruskall-Wallis)
# only vessels with the same content are compared (e.g. C+SC with C+SC, etc.)
# since fluorescence signal can be dampened due to the presence of bacteria

    df_content <- df

    kw <- kruskal.test(value ~ condition, data=df_content)
    
    dt <- dunn.test(df_content$value, df_content$condition, method="bh")
    l <- cldList(P.adjusted ~ comparisons, data=dt, threshold=0.05)
    
    print(kw)
    print(l)

