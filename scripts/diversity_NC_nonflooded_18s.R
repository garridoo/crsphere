
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load At Cr natural community dataset

source("load_data_NC_nonflooded_18s.R")

# subset samples of interest

idx <- design$compartment %in% c("Input", "Surface") &
       T
    
design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]
taxonomy_subset <- taxonomy[idx, ]

### plot stacked barplot of superkingdom relative abundances

taxonomy_subset <- taxonomy_subset[idx, ]
asv_table_subset <- asv_table_subset[idx, ]

sp_asv_table_subset <- aggregate(asv_table_subset, by=list(taxonomy_subset$Phylum), FUN=sum)
df <- melt(sp_asv_table_subset)
colnames(df) <- c("clade", "SampleID", "RA")
df$treatment <- design_subset$treatment[match(df$SampleID, design_subset$SampleID)]

df_means <- aggregate(df$RA, by=list(df$treatment, df$clade), FUN=mean)
colnames(df_means) <- c("treatment", "clade", "mean_RA")

df_means <- df_means[!df_means$mean_RA <= 0.01, ]

p <- ggplot(df_means, aes(x=treatment, y=mean_RA, fill=clade)) +
            geom_bar(stat="identity") +
            labs(x="Compartment", y="Relative abundance") +
            theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
            main_theme +
            theme(legend.position="right")

ggsave(paste(figures.dir, "nonflooded_NC_18s_superphylum.pdf", sep=""), p, width=shannon_width+5, height=shannon_height+2)

