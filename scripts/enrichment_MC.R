
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load mesocosm dataset

source("load_data_MC.R")

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

idx <- design$system=="flask" &
       design$community_AP %in% c("Input_FALSE", "B+C_FALSE", "B_FALSE", "B_TRUE") &
       design$soil %in% c("CAS") &
       design$experiment %in% c("AL1", "AL2", "AL3", "AL4", "AL5", "AL6") &
       design$medium %in% c("TP", "TP+AP", "Input") &
       design$timepoint %in% c("Day0", "Day7") &
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
p.vals_Input_BC <- p.vals_BC_Input <- p.vals_B_Input <- p.vals_BAP_Input <- data.frame(ASV=asvs, p=NA, p.adj=NA)

for (i in 1:length(asvs)) {

    asv <- asvs[i]
    df_asv <- df[df$ASV==asv, ]

    # get vectors of relative abundances for each subset of samples

    a <- df_asv$RA[df_asv$community=="Input"]
    b <- df_asv$RA[df_asv$community=="B+C" & df_asv$medium %in% c("TP", "BD")]
    c <- df_asv$RA[df_asv$community=="B" & df_asv$medium %in% c("TP", "BD")]
    d <- df_asv$RA[df_asv$community=="B" & df_asv$medium %in% c("TP+AP", "BD+AP")]

    if(length(c)==0) c <- 0

    # find species significantly enriched in the input compared to output samples

    p.vals_Input_BC$p[i] <- wilcox.test(a, c(b, d, d), alternative="greater")$p.value

    # find species enriched in different conditions compared to input samples

    p.vals_BC_Input$p[i] <- wilcox.test(b, a, alternative="greater")$p.value
    p.vals_B_Input$p[i] <- wilcox.test(c, a, alternative="greater")$p.value
    p.vals_BAP_Input$p[i] <- wilcox.test(d, a, alternative="greater")$p.value

}

p.vals_Input_BC$p.adj <- p.adjust(p.vals_Input_BC$p, method=p.adj.method)
p.vals_Input_BC  <- p.vals_Input_BC[p.vals_Input_BC$p.adj < alpha, ]
Input_ASVs <- as.character(p.vals_Input_BC$ASV)

p.vals_BC_Input$p.adj <- p.adjust(p.vals_BC_Input$p, method=p.adj.method)
p.vals_BC_Input  <- p.vals_BC_Input[p.vals_BC_Input$p.adj < alpha, ]
BC_ASVs <- as.character(p.vals_BC_Input$ASV)

p.vals_B_Input$p.adj <- p.adjust(p.vals_B_Input$p, method=p.adj.method)
p.vals_B_Input  <- p.vals_B_Input[p.vals_B_Input$p.adj < alpha, ]
B_ASVs <- as.character(p.vals_B_Input$ASV)

p.vals_BAP_Input$p.adj <- p.adjust(p.vals_BAP_Input$p, method=p.adj.method)
p.vals_BAP_Input  <- p.vals_BAP_Input[p.vals_BAP_Input$p.adj < alpha, ]
BAP_ASVs <- as.character(p.vals_BAP_Input$ASV)

# create table of enriched ASVs

enriched_asvs <- data.frame(Input_enriched=rownames(asv_table) %in% Input_ASVs,
                            BC_enriched=rownames(asv_table) %in% BC_ASVs,
                            B_enriched=rownames(asv_table) %in% B_ASVs,
                            BAP_enriched=rownames(asv_table) %in% BAP_ASVs)
row.names(enriched_asvs) <- rownames(asv_table)

