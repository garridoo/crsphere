
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
library("sva")

# load paths 

source("paths.R")

# files

design.file <- paste(data.dir, "design_SC_flask.txt", sep="")
asv_table.file <- paste(results.dir, "asv_table_SC_flask.txt", sep="")
taxonomy.file <- paste(data.dir, "crsphere_metadata.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
asv_table <- read.table(asv_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T)

# re-order data matrices

idx <- design$SampleID %in% colnames(asv_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(asv_table))
asv_table <- asv_table[, idx]

# normalize otu tables

design$depth <- colSums(asv_table)
asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

# quantify chlamy ra

idx <- rownames(asv_table)!="Chlamydomonas_miseq"
design$chlamy_reads <- as.numeric(asv_table[!idx, ])

# re-normalize otu table

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))
design$chlamy_ra <- design$chlamy_reads / design$depth

idx  <- design$depth >= 1000
design  <- design[idx, ]
asv_table  <- asv_table[, idx]
asv_table_norm  <- asv_table_norm[, idx]

# re-order taxonomy table

idx <- match(rownames(asv_table), taxonomy$ID) 
taxonomy <- taxonomy[idx, ]
 
