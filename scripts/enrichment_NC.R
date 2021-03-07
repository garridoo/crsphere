
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load natural community dataset

source("load_data_NC.R")

# thresholding

threshold <- 0.5
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]
taxonomy_subset <- taxonomy[idx, ]

# parameters

alpha <- 0.05               # significance threshold
p.adj.method <- "none"      # FDR p-value adjustment method

library(edgeR)

# subset samples

idx <- design$system=="soil" &
       design$timepoint %in% c("Day36") &
       design$compartment %in% c("soil", "phycosphere", "root") &
       TRUE

design_subset <- design[idx, ]
asv_table_subset <- asv_table[, idx]
asv_table_subset_norm <- asv_table_norm[, idx]

# remove empty rows from the ASV table

idx <- rowSums(asv_table_subset_norm)!=0
asv_table_subset <- asv_table_subset[idx, ]
asv_table_subset_norm <- asv_table_subset_norm[idx, ]

# create data frame of RAs

df <- melt(asv_table_subset_norm)
colnames(df) <- c("ASV", "SampleID", "RA")
df <- cbind(df, design[match(df$SampleID, design$SampleID), -1])

# determine ASV enrichment status

asvs <- rownames(asv_table_subset_norm)
p.vals_soil <- p.vals_soil_root <- p.vals_root_soil <- p.vals_phycosphere_soil <- data.frame(ASV=asvs, p=NA, p.adj=NA)

for (i in 1:length(asvs)) {

    asv <- asvs[i]
    df_asv <- df[df$ASV==asv, ]

    # get vectors of relative abundances for each subset of samples

    a <- df_asv$RA[df_asv$compartment=="soil"]
    b <- df_asv$RA[df_asv$compartment=="root"]
    c <- df_asv$RA[df_asv$compartment=="phycosphere"]

    if(length(c)==0) c <- 0

    # find species significantly enriched in the input compared to output samples

    p.vals_soil$p[i] <- wilcox.test(a, c(b, c), alternative="greater")$p.value

    # find species enriched in different conditions compared to soil samples

    p.vals_root_soil$p[i] <- wilcox.test(b, a, alternative="greater")$p.value
    p.vals_phycosphere_soil$p[i] <- wilcox.test(c, a, alternative="greater")$p.value

}

p.vals_soil$p.adj <- p.adjust(p.vals_soil$p, method=p.adj.method)
p.vals_soil  <- p.vals_soil[p.vals_soil$p.adj < alpha, ]
soil_ASVs <- as.character(p.vals_soil$ASV)

p.vals_root_soil$p.adj <- p.adjust(p.vals_root_soil$p, method=p.adj.method)
p.vals_root_soil  <- p.vals_root_soil[p.vals_root_soil$p.adj < alpha, ]
root_ASVs <- as.character(p.vals_root_soil$ASV)

p.vals_phycosphere_soil$p.adj <- p.adjust(p.vals_phycosphere_soil$p, method=p.adj.method)
p.vals_phycosphere_soil  <- p.vals_phycosphere_soil[p.vals_phycosphere_soil$p.adj < alpha, ]
phycosphere_ASVs <- as.character(p.vals_phycosphere_soil$ASV)

# create table of enriched ASVs

enriched_asvs <- data.frame(soil_enriched=rownames(asv_table) %in% soil_ASVs,
                            root_enriched=rownames(asv_table) %in% root_ASVs,
                            phycosphere_enriched=rownames(asv_table) %in% phycosphere_ASVs)
row.names(enriched_asvs) <- rownames(asv_table)

