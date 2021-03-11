
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load data

source("IPL_load_data.R")

# normalize OTU tables

otu_table_norm_At_NC <- apply(otu_table_At_NC, 2, function(x) x/sum(x))
otu_table_norm_Cr_NC <- apply(otu_table_Cr_NC, 2, function(x) x/sum(x))
otu_table_norm_Cr_MC <- apply(otu_table_Cr_MC, 2, function(x) x/sum(x))

# remove OTUs with zero counts

idx <- rowSums(otu_table_CrIPL) > 0
otu_table_CrIPL <- otu_table_CrIPL[idx, ]
otu_table_norm_CrIPL <- otu_table_norm_CrIPL[idx, ]

# remove wells without a insufficient number of reads

min_reads <- 100

idx <- colSums(otu_table_CrIPL) > min_reads
otu_table_CrIPL <- otu_table_CrIPL[, idx]
otu_table_norm_CrIPL <- otu_table_norm_CrIPL[, idx]

### calculate recovery rates for culture-independent samples

source("IPL_recovery_rates.R")

