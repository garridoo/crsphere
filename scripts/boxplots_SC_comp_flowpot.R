
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# SynCom flask dataset

source("load_data_SC_comp_flowpot.R")

# remove chlamy reads

idx <- rownames(asv_table)!="Chlamydomonas_miseq"
asv_table <- asv_table[idx, ]

# re-normalize otu table

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

# subset samples of interest

idx <- design$system %in% c("FP", "FPL") &
       design$condition %in% c("At+Full", "C+Full", "Full", "Full_MS", "Full_TP") &
       design$fraction %in% c("root", "cells", "soil") &
       T

design$fraction_treatment <- paste(design$fraction, design$treatment)

design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]

# re-order taxonomy table

idx <- match(rownames(asv_table), taxonomy$ID) 
taxonomy <- taxonomy[idx, ]

# colors

colors <- data.frame(group=c("input", "Cr", "At", "soil"),
                     color=c(input_color, cr_color, at_color, soil_color))

shapes <- data.frame(group=c("root Full", "root ICL", "root IRL", "soil Full", "soil ICL", "soil IRL", "cells Full", "input Full", "input IRL", "input ICL"),
                     shape=c(19, 10, 13, 15, 12, 7, 17, 19, 10, 13))

# aggregate RAs according to host of origin

df <- melt(aggregate(asv_table_subset, by=list(taxonomy$collection), FUN=sum))
colnames(df) <- c("SynCom", "SampleID", "RA")

idx <- match(df$SampleID, design$SampleID)
df$host <- design$host[idx]
df$fraction <- design$fraction[idx]
df$fraction_treatment <- paste(design$fraction[idx], design$treatment[idx])

colors <- colors[colors$group %in% df$host, ]
df$host <- factor(df$host, levels=colors$group)

shapes <- shapes[shapes$group %in% df$fraction_treatment, ]
df$fraction_treatment <- factor(df$fraction_treatment, levels=shapes$group)

# plot boxplots

p <- ggplot(df, aes(x=SynCom, y=RA, color=host, shape=fraction_treatment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0.2, 0.8), labels=percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="", y="Relative abundance") +
            ggtitle("Host preference") +
            main_theme +
            theme(legend.position="none")

ggsave(paste(figures.dir, "competition_flowpot_boxplots.pdf", sep=""), p, width=6, height=5)

