
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")
source("cpcoa.func.R")

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")

# load paths 

source("paths.R")

# files

design.file <- paste(data.dir, "design_NC_nonflooded_18s.txt", sep="")
asv_table.file <- paste(results.dir, "asv_table_NC_nonflooded_18s.txt", sep="")
taxonomy.file <- paste(results.dir, "taxonomy_NC_nonflooded_18s.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
asv_table <- read.table(asv_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T, quote="")

# re-order data matrices

idx <- design$SampleID %in% colnames(asv_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(asv_table))
asv_table <- asv_table[, idx]

idx <- match(rownames(asv_table), taxonomy[, 1])
taxonomy <- taxonomy[idx, ]

# remove land plants and H. sapiens

idx <- !taxonomy$Family %in% c("Embryophyceae_XX", "Mammalia")
taxonomy <- taxonomy[idx, ]
asv_table <- asv_table[idx, ]

# normalize otu tables

design$depth <- colSums(asv_table)
asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

# thresholding

threshold <- 0
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx  <- design$depth >= 1000
design  <- design[idx, ]
asv_table  <- asv_table[, idx]
asv_table_norm  <- asv_table_norm[, idx]

# re-order taxonomy table

idx <- match(rownames(asv_table), taxonomy[, 1]) 
taxonomy <- taxonomy[idx, ]

