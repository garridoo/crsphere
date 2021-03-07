
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# SynCom flask dataset

source("load_data_SC_comp_flask.R")

# remove chlamy reads

idx <- rownames(asv_table)!="Chlamydomonas_miseq"
asv_table <- asv_table[idx, ]

# re-normalize otu table

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

# subset samples of interest

idx <- design$system %in% c("flask") &
       design$treatment %in% c("Full") &
       T

design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]

# re-order taxonomy table

idx <- match(rownames(asv_table), taxonomy$ID) 
taxonomy <- taxonomy[idx, ]
 
# colors

colors <- data.frame(group=c("C+SC", "SC", "input"),
                     color=c(cr_color, bacteria_color, input_color))

# aggregate RAs according to host of origin

df <- melt(aggregate(asv_table_subset, by=list(taxonomy$collection), FUN=sum))
colnames(df) <- c("SynCom", "SampleID", "RA")

idx <- match(df$SampleID, design$SampleID)
df$community <- design$community[idx]

colors <- colors[colors$group %in% df$community, ]
df$community <- factor(df$community, levels=colors$group)
 
# plot boxplots

p <- ggplot(df, aes(x=SynCom, y=RA, color=community)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0.2, 0.8), labels=percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="", y="Relative abundance") +
            ggtitle("Host preference") +
            main_theme +
            theme(legend.position="none")

ggsave(paste(figures.dir, "competition_flask_boxplots.pdf", sep=""), p, width=3, height=5)

